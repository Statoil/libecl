from __future__ import annotations

import pytest
from resdata.grid import Grid

from ._grid_fixtures import (
    CELL_ACTIVE_FRACTURE,
    CELL_ACTIVE_MATRIX,
    GLOBAL_STRING,
    load_egrid_dual_porosity,
    load_egrid_with_nested_lgr,
    load_egrid_with_single_lgr,
    load_egrid_with_two_lgrs_and_amalgamated_nnc,
    load_grid_file_main_with_lgr,
    load_grid_file_with_lgr_parent,
    load_grid_file_with_mapaxes,
    load_grid_file_with_nested_lgr,
    write_egrid_dual_porosity,
    write_egrid_with_single_lgr,
)


def test_that_egrid_with_single_lgr_reports_one_lgr(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
    )

    assert grid.get_num_lgr() == 1
    assert grid.has_lgr("LGR1")


def test_that_egrid_with_negative_lgr_number_raises(tmp_path):
    filename = tmp_path / "NEGATIVE_LGR_NUMBER.EGRID"
    write_egrid_with_single_lgr(
        filename,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        "LGR1",
        lgr_nr=-1000,
    )

    with pytest.raises(ValueError, match="Invalid lgr_nr: lgr_nr=-1000"):
        Grid(str(filename))


def test_that_egrid_with_lgr_number_above_upper_bound_raises(tmp_path):
    filename = tmp_path / "TOO_LARGE_LGR_NUMBER.EGRID"
    write_egrid_with_single_lgr(
        filename,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        "LGR1",
        lgr_nr=1_000_001,
    )

    with pytest.raises(
        ValueError, match="LGR number 1000001 exceeds maximum supported value"
    ):
        Grid(str(filename))


def test_that_lgr_grid_has_expected_dimensions(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR_DIMS.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
    )
    lgr = grid.get_lgr("LGR1")
    assert lgr.get_nx() == 2
    assert lgr.get_ny() == 2
    assert lgr.get_nz() == 2


def test_that_main_grid_with_single_lgr_has_no_fracture_cells(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR_FRAC.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
    )
    assert grid.get_active_fracture_index(global_index=0) == -1
    assert grid.get_global_index1F(0) is None


def test_that_actnum_export_matches_per_cell_active_status_with_lgr(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR_ACTNUM.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
    )
    actnum = list(grid.export_actnum())
    for i in range(grid.get_global_size()):
        expected = 1 if grid.active(global_index=i) else 0
        assert actnum[i] == expected


def test_that_grid_with_lgr_can_be_written_as_grid_and_reloaded(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR_RW.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
    )
    grid_filename = tmp_path / "LGR_RW.GRID"
    grid.save_GRID(str(grid_filename))
    assert grid_filename.exists()

    reloaded = Grid(str(grid_filename))
    assert reloaded.get_num_lgr() == 1
    assert reloaded.has_lgr("LGR1")


def test_that_egrid_with_lgr_and_nnc_loads(tmp_path):
    grid = load_egrid_with_single_lgr(
        tmp_path / "LGR_NNC.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
        nncg=[1],
        nncl=[1],
    )
    assert grid.has_lgr("LGR1")


def test_that_egrid_with_nested_lgrs_reports_both_lgrs(tmp_path):
    grid = load_egrid_with_nested_lgr(
        tmp_path / "NESTED.EGRID",
        3,
        3,
        3,
        1,
        1,
        1,
        2,
        2,
        2,
        "OUTER",
        0,
        0,
        0,
        2,
        2,
        2,
        "INNER",
    )
    assert grid.get_num_lgr() == 2
    assert grid.has_lgr("OUTER")
    assert grid.has_lgr("INNER")


def test_that_grid_file_with_empty_parent_lgr_loads(tmp_path):
    grid = load_grid_file_with_lgr_parent(
        tmp_path / "LGR_EMPTY.GRID",
        "LGR1",
        "",
    )
    assert grid.get_num_lgr() == 1
    assert grid.has_lgr("LGR1")


def test_that_grid_file_with_GLOBAL_parent_lgr_loads(tmp_path):
    grid = load_grid_file_with_lgr_parent(
        tmp_path / "LGR_GLOBAL.GRID",
        "LGR1",
        GLOBAL_STRING,
    )
    assert grid.get_num_lgr() == 1
    assert grid.has_lgr("LGR1")


def test_that_grid_file_with_nested_lgr_loads_both(tmp_path):
    grid = load_grid_file_with_nested_lgr(
        tmp_path / "LGR_NESTED.GRID",
        "OUTER",
        "INNER",
    )
    assert grid.get_num_lgr() == 2
    assert grid.has_lgr("OUTER")
    assert grid.has_lgr("INNER")


@pytest.fixture
def egrid_with_mapaxes(tmp_path):
    mapaxes = [10.0, 21.0, 10.0, 20.0, 11.0, 21.0]
    grid = load_egrid_with_single_lgr(
        tmp_path / "MAPAXES.EGRID",
        3,
        3,
        3,
        2,
        2,
        2,
        1,
        1,
        1,
        "LGR1",
        mapaxes=mapaxes,
    )
    return grid, mapaxes


def test_that_export_mapaxes_returns_kw_with_input_values(egrid_with_mapaxes):
    grid, mapaxes = egrid_with_mapaxes
    kw = grid.export_mapaxes()
    assert kw is not None
    assert len(kw) == 6
    for i in range(6):
        assert kw[i] == mapaxes[i]


def test_that_grid_file_with_even_nz_main_and_lgr_is_not_dual(tmp_path):
    grid = load_grid_file_main_with_lgr(
        tmp_path / "EVEN_NZ_WITH_LGR.GRID",
        1,
        1,
        2,
    )
    assert not grid.dual_grid()
    assert grid.has_lgr("LGR1")


def test_that_dual_porosity_egrid_with_nnc_loads_and_drops_invalid_nnc(tmp_path):
    nx, ny, nz = 1, 1, 4
    size = nx * ny * nz
    actnum = [CELL_ACTIVE_MATRIX | CELL_ACTIVE_FRACTURE] * size

    grid = load_egrid_dual_porosity(
        tmp_path / "DUALP_NNC.EGRID",
        nx,
        ny,
        nz,
        actnum,
        nnc1=[1, 3],
        nnc2=[2, 4],
    )
    assert grid is not None
    assert grid.dual_grid()


@pytest.fixture
def dual_porosity_grid(tmp_path):
    nx, ny, nz = 2, 2, 2
    size = nx * ny * nz
    actnum = [CELL_ACTIVE_MATRIX | CELL_ACTIVE_FRACTURE] * size
    actnum[0] = CELL_ACTIVE_MATRIX

    grid = load_egrid_dual_porosity(
        tmp_path / "DUALP_FRAC.EGRID",
        nx,
        ny,
        nz,
        actnum,
    )
    return grid, actnum, size


def test_that_count_of_fracture_active_cells_matches_actnum(dual_porosity_grid):
    grid, _, size = dual_porosity_grid
    assert grid.get_num_active_fracture() == size - 1


def test_that_each_fracture_active_cell_maps_to_a_fracture_index(
    dual_porosity_grid,
):
    grid, actnum, size = dual_porosity_grid
    nactive_fracture = grid.get_num_active_fracture()
    seen = [False] * nactive_fracture
    for g in range(size):
        f = grid.get_active_fracture_index(global_index=g)
        if actnum[g] & CELL_ACTIVE_FRACTURE:
            assert 0 <= f < nactive_fracture
            assert not seen[f]
            seen[f] = True
            assert grid.get_global_index1F(f) == g
        else:
            assert f == -1
    assert all(seen)


def test_that_fracture_index_past_the_last_active_returns_none(
    dual_porosity_grid,
):
    grid, _, _ = dual_porosity_grid
    nactive_fracture = grid.get_num_active_fracture()
    assert grid.get_global_index1F(nactive_fracture) is None


def test_that_egrid_with_amalgamated_lgr_to_lgr_nnc_loads(tmp_path):
    grid = load_egrid_with_two_lgrs_and_amalgamated_nnc(
        tmp_path / "LGR_AMALGAMATED_NNC.EGRID",
        3,
        3,
        3,
        "LGR1",
        0,
        0,
        0,
        "LGR2",
        2,
        2,
        2,
        nna1=[1],
        nna2=[1],
    )
    assert grid.has_lgr("LGR1")
    assert grid.has_lgr("LGR2")


def test_that_nnc1_shorter_than_nnc2_is_rejected(tmp_path):
    actnum = [CELL_ACTIVE_MATRIX | CELL_ACTIVE_FRACTURE] * 4
    filename = tmp_path / "SHORT_NNC1.EGRID"
    write_egrid_dual_porosity(filename, 1, 1, 4, actnum, nnc1=[1], nnc2=[2] * 100000)

    with pytest.raises(ValueError, match=r"NNC1.*NNC2"):
        Grid(str(filename))
