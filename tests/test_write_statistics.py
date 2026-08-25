# tests for machine-friendly stats output
import io

from pymc_nmr.basin_hop import BasinHop
from pymc_nmr.monte_carlo import MonteCarlo
from pymc_nmr.airss_style import AirssStyle


def _assert_csv_line(f, expected):
    f.seek(0)
    lines = f.read().strip().split("\n")
    assert lines[0] == "step,total_energy,a,b,c,alpha,beta,gamma"
    assert lines[1] == expected


def test_bh_writes_csv_header_and_row():
    bh = BasinHop()
    f = io.StringIO()
    bh.write_statistics(5, -42.0, [1.0, 2.0, 3.0, 90.0, 80.0, 70.0], f)
    _assert_csv_line(f, "5,-42.0,1.0,2.0,3.0,90.0,80.0,70.0")


def test_mc_writes_csv_header_and_row():
    mc = MonteCarlo()
    f = io.StringIO()
    mc.write_statistics(7, -3.14, [4.0, 5.0, 6.0, 60.0, 70.0, 80.0], f)
    _assert_csv_line(f, "7,-3.14,4.0,5.0,6.0,60.0,70.0,80.0")


def test_airss_writes_csv_header_and_row():
    airss = AirssStyle()
    f = io.StringIO()
    airss.write_statistics(9, -1.0, [2.0, 3.0, 4.0, 90.0, 90.0, 90.0], f)
    _assert_csv_line(f, "9,-1.0,2.0,3.0,4.0,90.0,90.0,90.0")


def test_appending_does_not_repeat_header():
    bh = BasinHop()
    f = io.StringIO()
    f.write("existing line\n")
    bh.write_statistics(2, -10.0, [1.0] * 6, f)
    f.seek(0)
    lines = f.read().strip().split("\n")
    assert lines[0] == "existing line"
    assert lines[1] == "2,-10.0,1.0,1.0,1.0,1.0,1.0,1.0"


if __name__ == "__main__":
    test_bh_writes_csv_header_and_row()
    test_mc_writes_csv_header_and_row()
    test_airss_writes_csv_header_and_row()
    test_appending_does_not_repeat_header()
    print("all write_statistics tests passed")
