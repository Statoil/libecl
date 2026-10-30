from __future__ import annotations

from pathlib import Path

import pytest
import resdata.grid._grid as _grid
from cwrap import CFILE
from cwrap import open as copen
from resdata import ResDataType, UnitSystem
from resdata.grid import Grid
from resdata.resfile import ResdataKW
from resdata.util.util import DoubleVector

from ._grid_fixtures import (
    load_egrid_with_single_lgr,
    make_rectangular_grid,
    write_egrid_with_actnum,
    write_egrid_with_single_lgr,
)


@pytest.fixture
def rectangular_grid():
    return make_rectangular_grid(3, 3, 3, 1, 1, 1)


def test_that_grid_can_be_written_as_GRID_and_loaded_back(rectangular_grid, tmp_path):
    filename = tmp_path / "TEST.GRID"
    rectangular_grid.save_GRID(str(filename))

    loaded = Grid(str(filename))
    assert rectangular_grid.equal(loaded, include_lgr=False, include_nnc=False)


def test_that_grid_can_be_written_as_GRDECL(rectangular_grid, tmp_path):
    filename = tmp_path / "TEST.GRDECL"
    with copen(str(filename), "w") as f:
        rectangular_grid.save_grdecl(f, output_unit=UnitSystem.METRIC)
    assert filename.exists()
    assert filename.stat().st_size > 0


def test_that_load_from_file_loads_a_saved_egrid(rectangular_grid, tmp_path):
    filename = tmp_path / "TEST3.EGRID"
    rectangular_grid.save_EGRID(str(filename))
    loaded = Grid.load_from_file(str(filename))
    assert loaded is not None
    assert rectangular_grid.equal(loaded, include_lgr=False, include_nnc=False)


def test_that_grid_keyword_can_be_printed_as_GRDECL(rectangular_grid, tmp_path):
    size = rectangular_grid.get_global_size()
    kw = ResdataKW("PORO", size, ResDataType.RD_FLOAT)
    for i in range(size):
        kw[i] = 0.2 + i * 0.01

    filename = tmp_path / "KW.GRDECL"
    with copen(str(filename), "w") as f:
        rectangular_grid.write_grdecl(kw, f, default_value=-999.0)
    assert filename.exists()
    assert filename.stat().st_size > 0


@pytest.fixture
def grid_with_inactive_cells():
    actnum = [1] * 8
    actnum[0] = 0
    actnum[3] = 0
    return make_rectangular_grid(2, 2, 2, 1, 1, 1, actnum), actnum


def _read_grdecl_kw(path: Path, name: str, size: int, rd_type) -> ResdataKW:
    with copen(str(path), "r") as f:
        return ResdataKW.read_grdecl(f, name, rd_type=rd_type)


def test_that_active_sized_float_keyword_fills_default_in_inactive_indices(
    grid_with_inactive_cells,
    tmp_path,
):
    grid, actnum = grid_with_inactive_cells
    nactive = grid.get_num_active()
    global_size = grid.get_global_size()
    assert nactive < global_size

    kw = ResdataKW("PORO", nactive, ResDataType.RD_FLOAT)
    for a in range(nactive):
        kw[a] = 0.2 + a * 0.01

    filename = tmp_path / "KW_ACTIVE_FLOAT.GRDECL"
    with copen(str(filename), "w") as f:
        grid.write_grdecl(kw, f, default_value=-999.0)

    read_back = _read_grdecl_kw(filename, "PORO", global_size, ResDataType.RD_FLOAT)
    a = 0
    for i in range(global_size):
        if actnum[i]:
            assert read_back[i] == kw[a]
            a += 1
        else:
            assert read_back[i] == pytest.approx(-999.0)


def test_that_active_sized_int_keyword_fills_default_in_inactive_indices(
    grid_with_inactive_cells,
    tmp_path,
):
    grid, actnum = grid_with_inactive_cells
    nactive = grid.get_num_active()
    global_size = grid.get_global_size()

    kw = ResdataKW("SATNUM", nactive, ResDataType.RD_INT)
    for a in range(nactive):
        kw[a] = a + 1

    filename = tmp_path / "KW_ACTIVE_INT.GRDECL"
    with copen(str(filename), "w") as f:
        grid.write_grdecl(kw, f, default_value=-999)

    read_back = _read_grdecl_kw(filename, "SATNUM", global_size, ResDataType.RD_INT)
    a = 0
    for i in range(global_size):
        if actnum[i]:
            assert read_back[i] == kw[a]
            a += 1
        else:
            assert read_back[i] == -999


def test_that_every_grid_can_be_written_to_disk_as_egrid_and_read_back(
    all_grids, tmp_path
):
    for label, grid in all_grids:
        filename = tmp_path / f"WRITE_{label}.EGRID"
        grid.save_EGRID(str(filename))
        assert filename.exists(), f"failed to write {label}"
        reloaded = Grid(str(filename))
        assert reloaded is not None


def test_that_every_grid_can_be_written_as_grdecl(all_grids, tmp_path):
    for label, grid in all_grids:
        filename = tmp_path / f"WRITE_{label}.GRDECL"
        with copen(str(filename), "w") as f:
            grid.save_grdecl(f, output_unit=UnitSystem.METRIC)
        assert filename.exists(), f"failed to write {label}"
        assert filename.stat().st_size > 0


def test_that_get_global_index3_with_out_of_range_ijk_raises_index_error(
    rectangular_grid,
):
    with pytest.raises(IndexError, match="i,j,k"):
        _grid._get_global_index3(rectangular_grid, 99, 0, 0)


def test_that_get_corner_xyz_with_out_of_range_i_raises_index_error(
    rectangular_grid,
):
    with pytest.raises(IndexError, match="invalid i value"):
        _grid._get_corner_xyz(rectangular_grid, -1, 0, 0)


def test_that_get_corner_xyz_with_out_of_range_j_raises_index_error(
    rectangular_grid,
):
    with pytest.raises(IndexError, match="invalid j value"):
        _grid._get_corner_xyz(rectangular_grid, 0, 100, 0)


def test_that_get_corner_xyz_with_out_of_range_k_raises_index_error(
    rectangular_grid,
):
    with pytest.raises(IndexError, match="invalid k value"):
        _grid._get_corner_xyz(rectangular_grid, 0, 0, -1)


def test_that_fetching_lgr_through_an_lgr_subgrid_raises_value_error(tmp_path):
    main_grid = load_egrid_with_single_lgr(
        tmp_path / "MAIN_LGR.EGRID", 3, 3, 3, 2, 2, 2, 1, 1, 1, "LGR1"
    )
    lgr = main_grid.get_lgr("LGR1")
    with pytest.raises(ValueError, match="LGR_grid"):
        _grid._get_named_lgr(lgr, "ANYTHING")


def test_that_grid_value_with_mismatched_kw_size_raises_value_error(
    rectangular_grid,
):
    bad_kw = ResdataKW("BAD", 7, ResDataType.RD_FLOAT)
    with pytest.raises(ValueError, match="incommensurable size"):
        _grid._grid_value(rectangular_grid, bad_kw, 0, 0, 0)


def test_that_grid_value_with_non_numeric_kw_raises_value_error(rectangular_grid):
    size = rectangular_grid.get_global_size()
    char_kw = ResdataKW("CHARKW", size, ResDataType.RD_CHAR)
    with pytest.raises(ValueError, match="can not lookup type"):
        _grid._grid_value(rectangular_grid, char_kw, 0, 0, 0)


def test_that_load_column_with_mismatched_kw_size_raises_value_error(
    rectangular_grid,
):
    bad_kw = ResdataKW("BAD", 7, ResDataType.RD_FLOAT)
    column = DoubleVector()
    with pytest.raises(ValueError, match="incommensurable sizes"):
        _grid._load_column(rectangular_grid, bad_kw, 0, 0, column)


def test_that_load_column_with_non_numeric_kw_raises_value_error(rectangular_grid):
    size = rectangular_grid.get_global_size()
    char_kw = ResdataKW("CHARKW", size, ResDataType.RD_CHAR)
    column = DoubleVector()
    with pytest.raises(ValueError, match="can not lookup type"):
        _grid._load_column(rectangular_grid, char_kw, 0, 0, column)


def test_that_fwrite_grdecl_with_invalid_bool_default_raises_value_error(
    grid_with_inactive_cells, tmp_path
):
    grid, _ = grid_with_inactive_cells
    nactive = grid.get_num_active()
    bool_kw = ResdataKW("BOOLKW", nactive, ResDataType.RD_BOOL)
    out = tmp_path / "BOOL_OUT.GRDECL"
    with copen(str(out), "w") as f:
        cfile = CFILE(f)
        with pytest.raises(ValueError, match="bool interpolation"):
            _grid._fwrite_grdecl(grid, bool_kw, None, cfile, 5.0)


def test_that_fwrite_grdecl_with_mismatched_kw_size_raises_value_error(
    rectangular_grid, tmp_path
):
    bad_kw = ResdataKW("BAD", 7, ResDataType.RD_FLOAT)
    out = tmp_path / "BAD_OUT.GRDECL"
    with copen(str(out), "w") as f:
        cfile = CFILE(f)
        with pytest.raises(ValueError, match="size mismatch"):
            _grid._fwrite_grdecl(rectangular_grid, bad_kw, None, cfile, 0.0)


def test_that_compressed_kw_copy_with_mismatched_sizes_raises_value_error(
    grid_with_inactive_cells,
):
    grid, _ = grid_with_inactive_cells
    target_kw = ResdataKW("T", 3, ResDataType.RD_FLOAT)
    src_kw = ResdataKW("S", 99, ResDataType.RD_FLOAT)
    with pytest.raises(ValueError, match="size mismatch"):
        _grid._compressed_kw_copy(grid, target_kw, src_kw)


def test_that_global_kw_copy_with_mismatched_sizes_raises_value_error(
    grid_with_inactive_cells,
):
    grid, _ = grid_with_inactive_cells
    target_kw = ResdataKW("T", 3, ResDataType.RD_FLOAT)
    src_kw = ResdataKW("S", 99, ResDataType.RD_FLOAT)
    with pytest.raises(ValueError, match="size mismatch"):
        _grid._global_kw_copy(grid, target_kw, src_kw)


def test_that_mapaxes_with_wrong_size_raises_when_loading_egrid(tmp_path):
    bad_size = 7
    filename = tmp_path / "BAD_MAPAXES.EGRID"
    write_egrid_with_single_lgr(
        filename,
        2,
        2,
        2,
        1,
        1,
        1,
        0,
        0,
        0,
        "LGR1",
        mapaxes=[0.0] * bad_size,
    )

    with pytest.raises(ValueError, match="MAPAXES"):
        Grid(str(filename))


@pytest.mark.parametrize("actnum_size", [1, 7, 9])
def test_that_egrid_with_wrongly_sized_actnum_is_rejected(tmp_path, actnum_size):
    filename = tmp_path / "BAD_ACTNUM.EGRID"
    write_egrid_with_actnum(filename, 2, 2, 2, [1] * actnum_size)

    with pytest.raises(ValueError, match="Invalid size of ACTNUM"):
        Grid(str(filename))


@pytest.mark.parametrize("hostnum_size", [1, 7])
def test_that_egrid_with_wrongly_sized_hostnum_is_rejected(tmp_path, hostnum_size):
    filename = tmp_path / "BAD_HOSTNUM.EGRID"
    write_egrid_with_single_lgr(
        filename, 2, 2, 2, 2, 2, 2, 0, 0, 0, "LGR1", hostnum=[1] * hostnum_size
    )

    with pytest.raises(ValueError, match="Invalid size of HOSTNUM"):
        Grid(str(filename))
