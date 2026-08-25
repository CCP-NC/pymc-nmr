"""Validated schema for the ``config.yaml`` input file.

Every key a user may write in ``config.yaml`` is declared here, together with
its default. The models are deliberately strict:

* ``extra="forbid"`` means a mistyped key (``temperatur:``) is reported as an
  error instead of being silently ignored, which would otherwise leave the run
  quietly using the default value.
* Enumerated options use ``Literal``, so an invalid choice is caught when the
  config is read rather than deep inside the run.

Units follow ASE conventions throughout: energies in eV, lengths in Angstrom,
temperatures in Kelvin, and therefore pressure in eV/Angstrom^3.
"""

from typing import List, Literal, Optional, Tuple

from pydantic import BaseModel, ConfigDict, Field


class _Strict(BaseModel):
    """Base model that rejects unknown keys so config typos are not silent."""

    model_config = ConfigDict(extra="forbid")


class PotentialSchema(_Strict):
    """The machine-learned interatomic potential used to evaluate energies."""

    arch: str = "mace_mp"
    model: str
    device: Literal["cpu", "cuda", "mps"] = "cpu"
    precision: Literal["float32", "float64"] = "float64"
    dispersion: bool = False
    threads_per_rank: int = 1


class SwapMove(_Strict):
    """One swap-move type: interchange the positions of two atom species."""

    pairs: List[Tuple[str, str]]
    counterion: Optional[str] = None
    max_distance: float = 3.0
    frequency: int = 100


class MovesSchema(_Strict):
    swaps: List[SwapMove] = Field(default_factory=list)


class MonteCarloSchema(_Strict):
    steps: int = 0
    max_distance: float = 0.001
    distance_update_freq: int = 1000
    distance_ratio: float = 0.37

    vol_move_freq: int = 0
    vol_move_symmetry: Literal["cubic", "tetragonal", "orthorhombic", "vectors"] = "cubic"
    max_vol_displacement: float = 0.1
    vol_update_freq: int = 1000
    vol_ratio: float = 0.37

    atom_move_freq: int = 0
    move_types: List[str] = Field(default_factory=list)

    transmutate_frequency: int = 0
    mutate_types_1: List[str] = Field(default_factory=list)
    mutate_types_2: List[str] = Field(default_factory=list)
    transmute_chem_pot: List[float] = Field(default_factory=list)


class MDSchema(_Strict):
    """Optional Langevin MD "shake" move, used to escape a local basin."""

    move_freq: int = 0
    timestep_fs: float = 2.0
    temperature_k: float = 1000.0
    friction_fs: float = 0.01
    steps: int = 1000


class RelaxSchema(_Strict):
    """Geometry relaxation applied after a move (basin hopping and AIRSS).

    ``style`` selects which cell degrees of freedom are relaxed:
    ``conp`` all cell parameters, ``cona`` cell lengths only, ``conv`` fixed cell.
    """

    method: Literal["lbfgs", "fire", "bfgs"] = "lbfgs"
    steps: int = 1000
    tol: float = 1e-3
    style: Literal["conp", "cona", "conv"] = "conp"


class GridSchema(_Strict):
    """Spatial grid used to find empty sites for counterion placement."""

    x: int = 10
    y: int = 10
    z: int = 10
    cut: float = 2.0


class JobSchema(_Strict):
    structure_path: str = "basin.xyz"
    structure_method: Literal["basinhop", "airss", "monte"] = "basinhop"
    num_cycles: int = 1

    # External pressure for NpT volume moves, in eV/Angstrom^3 (not GPa).
    pressure: float = 0.0
    # Sampling temperature in Kelvin. Required: it sets beta = 1/(kB T) in the
    # Metropolis criterion, so there is no meaningful default and T = 0 would
    # divide by zero.
    temperature: float = Field(gt=0.0)

    wrap: bool = False
    restart: bool = False
    max_force: float = 1.0e-3
    frozen_types: List[str] = Field(default_factory=list)

    sanity_check_freq: int = 1000
    print_freq: int = 1
    dump_archive: bool = False
    archive_frequency: int = 1000
    save_downhill: bool = False
    equil_steps: int = 0
    writestats: bool = False
    writestats_freq: int = 10


class SimulationConfig(_Strict):
    """Top-level model for ``config.yaml``."""

    job: JobSchema
    potential: PotentialSchema
    monte: MonteCarloSchema = Field(default_factory=MonteCarloSchema)
    moves: MovesSchema = Field(default_factory=MovesSchema)
    md: MDSchema = Field(default_factory=MDSchema)
    relax: RelaxSchema = Field(default_factory=RelaxSchema)
    grid: GridSchema = Field(default_factory=GridSchema)
