#! /usr/bin/env python3
import os
import sys
import yaml
import logging

# 1. Quickly check if we have a config file in arguments to read OMP threads before heavy imports
threads_per_rank = 1
config_file = "config.yaml"
for i, arg in enumerate(sys.argv):
    if arg in ["--config", "-c"] and i + 1 < len(sys.argv):
        config_file = sys.argv[i + 1]
        break

if os.path.exists(config_file):
    try:
        with open(config_file, "r") as f:
            cfg_data = yaml.safe_load(f)
            threads_per_rank = cfg_data.get("potential", {}).get("threads_per_rank", 1)
    except Exception:
        pass

# Restrict multi-threading to threads_per_rank per rank to prevent thread over-subscription and thrashing
os.environ["OMP_NUM_THREADS"] = str(threads_per_rank)
os.environ["MKL_NUM_THREADS"] = str(threads_per_rank)
os.environ["OPENBLAS_NUM_THREADS"] = str(threads_per_rank)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads_per_rank)
os.environ["NUMEXPR_NUM_THREADS"] = str(threads_per_rank)

# 2. Quiet worker printing to keep terminal stdout pristine
try:
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
except ImportError:
    rank = 0

if rank > 0:
    import builtins
    builtins.print = lambda *args, **kwargs: None
    
    from rich.console import Console
    Console.print = lambda self, *args, **kwargs: None
    Console.log = lambda self, *args, **kwargs: None

    # Redirect sys-level stdout and stderr to prevent C/C++ or logging library pollution
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')

    # Silence python standard logging completely on worker ranks
    logging.getLogger().setLevel(logging.CRITICAL)

import random
import time
import numpy as np

import typer
from rich.console import Console

from pymc_nmr.schema import SimulationConfig
from pymc_nmr.config import Config
from pymc_nmr.job_control import JobControl
from pymc_nmr.species import Species, Element
from pymc_nmr.statistics import Statistics, TypeStatistics
from pymc_nmr.field import Field
from pymc_nmr.basin_hop import BasinHop
from pymc_nmr.airss_style import AirssStyle
from pymc_nmr.monte_carlo import MonteCarlo

app = typer.Typer(help="PyMC-NMR: Modern Monte Carlo & structure search CLI.")
console = Console()

def _aux_log_path(log_file: str, kind: str) -> str:
    # ponytail: rank suffix is already in log_file; just insert _<kind> before extension.
    base, ext = os.path.splitext(log_file)
    return f"{base}_{kind}{ext}"

class LogStream:
    # ponytail: line-buffered wrapper so the existing out_stream.write()/flush()/close() API stays unchanged.
    def __init__(self, path: str, rank: int):
        self._buf = ""
        self._logger = logging.getLogger(f"pymc_nmr.rank{rank}")
        self._logger.setLevel(logging.INFO)
        self._logger.propagate = False
        handler = logging.FileHandler(path, mode="a")
        handler.setFormatter(
            logging.Formatter(
                f"%(asctime)s [rank {rank}] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        self._logger.addHandler(handler)

    def write(self, msg: str):
        self._buf += msg
        *lines, self._buf = self._buf.split("\n")
        for line in lines:
            if line.strip():
                self._logger.info(line)

    def flush(self):
        if self._buf.strip():
            self._logger.info(self._buf)
            self._buf = ""
        for handler in self._logger.handlers:
            handler.flush()

    def close(self):
        self.flush()
        for handler in self._logger.handlers:
            handler.close()
        self._logger.handlers.clear()

def run_simulation(
    config_path: str = "config.yaml",
    structure_path: str = "basin.xyz",
    log_path: str = "mc.log",
    verbose: bool = True
):
    """
    Run PyMC-NMR structure search / MC simulation using YAML config.
    """
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    except ImportError:
        rank = 0
        size = 1
        # ponytail: warn when launched under MPI but mpi4py is missing; otherwise every rank becomes 0.
        if any(os.environ.get(v) for v in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MPIR_CVAR_NUM_VCIS")):
            console.print("[bold yellow]Warning: MPI runtime detected but mpi4py is not installed; running as a single rank.[/bold yellow]")

    # Only print console messages on rank 0
    verbose = verbose and (rank == 0)

    # ponytail: one logger per rank, no central aggregation; combine logs externally if needed.
    logging.getLogger().setLevel(logging.WARNING if rank == 0 else logging.CRITICAL)

    if verbose:
        console.print(f"[bold green]Loading configuration from {config_path}...[/bold green]")
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found at {config_path}")
        
    try:
        with open(config_path, "r") as f:
            raw_yaml = yaml.safe_load(f)
            sim_config = SimulationConfig(**raw_yaml)
    except Exception as e:
        raise ValueError(f"Failed to validate config schema: {e}")

    # Force PyTorch/MACE to use the exact number of threads specified in the YAML config per rank
    try:
        import torch
        torch.set_num_threads(sim_config.potential.threads_per_rank)
        torch.set_num_interop_threads(sim_config.potential.threads_per_rank)
        if verbose:
            console.print(f"[bold blue]Active PyTorch threads per rank: {torch.get_num_threads()}[/bold blue]")
    except ImportError:
        pass

    # Differentiate random seed per rank to explore different trajectories
    seed = 42 + rank
    np.random.seed(seed)
    random.seed(seed)
    # ponytail: also seed torch; MACE/MLIP may draw randomness through torch backends.
    try:
        import torch
        torch.manual_seed(seed)
    except ImportError:
        pass

    # ponytail: round-robin assign ranks to visible GPUs; keep device as "cuda"
    # for MACE compatibility (it only accepts cpu/cuda/mps/xpu, not cuda:N).
    if sim_config.potential.device.lower().startswith("cuda"):
        try:
            import torch
            gpus = torch.cuda.device_count()
            if gpus > 0:
                gpu_id = rank % gpus
                torch.cuda.set_device(gpu_id)
                if verbose:
                    console.print(f"[bold blue]Rank {rank} assigned to GPU {gpu_id} of {gpus}.[/bold blue]")
        except Exception:
            pass

    # Determine initial structure file path
    structure_file = sim_config.job.structure_path
    if structure_path != "basin.xyz":
        structure_file = structure_path

    if verbose:
        console.print(f"[bold green]Loading initial structure from {structure_file}...[/bold green]")
    if not os.path.exists(structure_file):
        raise FileNotFoundError(f"Structure file not found at {structure_file}")

    # Suffix log file names per rank to avoid write conflicts
    if size > 1:
        base, ext = os.path.splitext(log_path)
        log_file = f"{base}_rank{rank}{ext}"
    else:
        log_file = log_path

    out_stream = LogStream(log_file, rank)
    # ponytail: logger already timestamps each line; keep a plain start marker.
    out_stream.write("\n--- Simulation started ---\n")
    out_stream.write(f" rank {rank} of {size}, random seed {seed}\n")

    # Read basin structure first to derive species automatically
    basin = Config()
    try:
        with open(structure_file, "r") as instream:
            restart_iteration, restart_time, restart_energy = basin.read_config(instream)
    except Exception as e:
        raise ValueError(f"Failed to read structure configuration: {e}")

    # Derive unique species dynamically from structure symbols (putting "ghost" at the end if present)
    unique_symbols = []
    for s in basin.symbol:
        if s not in unique_symbols and s != "ghost":
            unique_symbols.append(s)
    if any("ghost" in s for s in basin.symbol):
        unique_symbols.append("ghost")

    # Populate Species object
    spec = Species()
    spec.number_of_elements = len(unique_symbols)
    spec.ele_data = [Element(name=s) for s in unique_symbols]

    # Setup config with derived species
    basin.setup_configuration(spec, out_stream)

    # Populate JobControl from our Pydantic SimulationConfig
    job = JobControl()
    job.load_from_schema(sim_config)
    job.rank_suffix = f"_rank{rank}" if size > 1 else ""

    # Setup Field from schema
    fld = Field()
    fld.load_from_schema(sim_config)
    fld.aux_log = _aux_log_path(log_file, "relax")

    bh = BasinHop()
    mc = MonteCarlo()
    ai = AirssStyle()

    if job.structure_method == "basinhop":
        job.write_bh_control(out_stream)
        bh.initialise(spec, job, out_stream)
    elif job.structure_method == "airss":
        job.write_airss_control(out_stream)
        ai.initialise(spec, job, out_stream)
    elif job.structure_method == "monte":
        job.write_mc_control(out_stream)
        mc.initialise(spec, job, out_stream)
    else:
        out_stream.write("\n*** unrecognised structure search method \n")
        raise ValueError(f"Unrecognised structure search method: {job.structure_method}")

    num_cycles = job.num_cycles
    stats = Statistics()
    type_stats = TypeStatistics()
     
    if verbose:
        console.print(f"[bold blue]Running {num_cycles} cycles of {job.structure_method} structure search...[/bold blue]")

    for cycle in range(num_cycles):
        num_steps = 0
        start_time = time.time()
        initialise = True

        if job.structure_method == "monte":
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Monte Carlo Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")
            mc.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, restart_iteration, out_stream)
        
        elif job.structure_method == "airss":
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" AIRSS Style Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")
            ai.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, initialise, restart_iteration, 
                   restart_energy, out_stream)
        else:
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Basin Hopping Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")
            bh.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, initialise, restart_iteration, restart_energy, 
                   out_stream)

        finish_time = time.time()
        diff = finish_time - start_time

        stats.last_summary(out_stream, num_steps)
        type_stats.last_summary_types(spec, out_stream)

        out_stream.write("\n\n" + " *" * 53 + "\n")
        out_stream.write(f"\n time to Monte Carlo simulation : {diff:.3f} seconds\n")
        out_stream.write("\n *** Monte Carlo Finished. Writing restart\n")
        out_stream.flush()

    out_stream.flush()
    out_stream.close()
    if verbose:
        console.print(f"[bold green]Simulation completed successfully! Log saved to {log_file}.[/bold green]")

@app.command()
def run(
    config_path: str = typer.Option("config.yaml", "--config", "-c", help="Path to YAML configuration file"),
    structure_path: str = typer.Option("basin.xyz", "--structure", "-s", help="Path to initial structure XYZ file"),
    log_path: str = typer.Option("mc.log", "--log", "-l", help="Path to the simulation log file"),
):
    """
    Run PyMC-NMR structure search / MC simulation using YAML config.
    """
    try:
        run_simulation(config_path=config_path, structure_path=structure_path, log_path=log_path)
    except Exception as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        raise typer.Exit(1)

if __name__ == "__main__":
    app()
