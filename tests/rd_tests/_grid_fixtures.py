from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import pytest
from resdata import ResDataType
from resdata.grid import Grid, GridGenerator
from resdata.resfile import FortIO, ResdataKW, openFortIO

# constants from lib/include/resdata/rd_kw_magic.hpp
FILEHEAD_KW = "FILEHEAD"
GRIDHEAD_KW = "GRIDHEAD"
ENDGRID_KW = "ENDGRID"
ENDLGR_KW = "ENDLGR"
LGR_KW = "LGR"
LGR_PARENT_KW = "LGRPARNT"
HOSTNUM_KW = "HOSTNUM"
COORD_KW = "COORD"
ZCORN_KW = "ZCORN"
ACTNUM_KW = "ACTNUM"
CORSNUM_KW = "CORSNUM"
MAPAXES_KW = "MAPAXES"
NNCHEAD_KW = "NNCHEAD"
NNCG_KW = "NNCG"
NNCL_KW = "NNCL"
NNC1_KW = "NNC1"
NNC2_KW = "NNC2"
NNCHEADA_KW = "NNCHEADA"
NNA1_KW = "NNA1"
NNA2_KW = "NNA2"
DIMENS_KW = "DIMENS"
RADIAL_KW = "RADIAL"
MAPUNITS_KW = "MAPUNITS"
GRIDUNIT_KW = "GRIDUNIT"
CORNERS_KW = "CORNERS"
COORDS_KW = "COORDS"

NNCHEAD_SIZE = 10
NNCHEAD_NUMNNC_INDEX = 0
NNCHEAD_LGR_INDEX = 1
NNCHEADA_ILOC1_INDEX = 0
NNCHEADA_ILOC2_INDEX = 1

GRIDHEAD_SIZE = 100
GRIDHEAD_TYPE_INDEX = 0
GRIDHEAD_NX_INDEX = 1
GRIDHEAD_NY_INDEX = 2
GRIDHEAD_NZ_INDEX = 3
GRIDHEAD_LGR_INDEX = 4
GRIDHEAD_NUMRES_INDEX = 24
GRIDHEAD_GRIDTYPE_CORNERPOINT = 1

FILEHEAD_VERSION_INDEX = 0
FILEHEAD_YEAR_INDEX = 1
FILEHEAD_TYPE_INDEX = 4
FILEHEAD_DUALP_INDEX = 5
FILEHEAD_ORGFORMAT_INDEX = 6
FILEHEAD_GRIDTYPE_CORNERPOINT = 0
FILEHEAD_ORGTYPE_CORNERPOINT = 1
FILEHEAD_SINGLE_POROSITY = 0
FILEHEAD_DUAL_POROSITY = 1

CELL_NOT_ACTIVE = 0
CELL_ACTIVE_MATRIX = 1
CELL_ACTIVE_FRACTURE = 2

GLOBAL_STRING = "GLOBAL"

UNIT_CELL_CORNERS = (
    0.0,
    0.0,
    0.0,
    1.0,
    0.0,
    0.0,
    0.0,
    1.0,
    0.0,
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    1.0,
    1.0,
    0.0,
    1.0,
    0.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
)


def coord_size(nx: int, ny: int) -> int:
    return (nx + 1) * (ny + 1) * 6


def zcorn_size(nx: int, ny: int, nz: int) -> int:
    return nx * ny * nz * 8


def make_int_kw(name: str, values: Sequence[int]) -> ResdataKW:
    kw = ResdataKW(name, len(values), ResDataType.RD_INT)
    for i, v in enumerate(values):
        kw[i] = int(v)
    return kw


def make_float_kw(name: str, values: Sequence[float]) -> ResdataKW:
    kw = ResdataKW(name, len(values), ResDataType.RD_FLOAT)
    for i, v in enumerate(values):
        kw[i] = float(v)
    return kw


def make_char_kw(name: str, values: Sequence[str]) -> ResdataKW:
    kw = ResdataKW(name, len(values), ResDataType.RD_CHAR)
    for i, v in enumerate(values):
        kw[i] = v
    return kw


def make_empty_kw(name: str) -> ResdataKW:
    return ResdataKW(name, 0, ResDataType.RD_INT)


def write_int_kw(fortio: FortIO, name: str, values: Sequence[int]) -> None:
    make_int_kw(name, values).fwrite(fortio)


def write_float_kw(fortio: FortIO, name: str, values: Sequence[float]) -> None:
    make_float_kw(name, values).fwrite(fortio)


def write_char_kw(fortio: FortIO, name: str, values: Sequence[str]) -> None:
    make_char_kw(name, values).fwrite(fortio)


def write_empty_kw(fortio: FortIO, name: str) -> None:
    make_empty_kw(name).fwrite(fortio)


def _filehead_values(dualp_flag: int = FILEHEAD_SINGLE_POROSITY) -> list:
    values = [0] * 100
    values[FILEHEAD_VERSION_INDEX] = 3
    values[FILEHEAD_YEAR_INDEX] = 2007
    values[FILEHEAD_TYPE_INDEX] = FILEHEAD_GRIDTYPE_CORNERPOINT
    values[FILEHEAD_DUALP_INDEX] = dualp_flag
    values[FILEHEAD_ORGFORMAT_INDEX] = FILEHEAD_ORGTYPE_CORNERPOINT
    return values


def _gridhead_values(nx: int, ny: int, nz: int, grid_nr: int) -> list:
    values = [0] * GRIDHEAD_SIZE
    values[GRIDHEAD_TYPE_INDEX] = GRIDHEAD_GRIDTYPE_CORNERPOINT
    values[GRIDHEAD_NX_INDEX] = nx
    values[GRIDHEAD_NY_INDEX] = ny
    values[GRIDHEAD_NZ_INDEX] = nz
    values[GRIDHEAD_NUMRES_INDEX] = 1
    values[GRIDHEAD_LGR_INDEX] = grid_nr
    return values


def _write_filehead(fortio: FortIO, dualp_flag: int = FILEHEAD_SINGLE_POROSITY) -> None:
    write_int_kw(fortio, FILEHEAD_KW, _filehead_values(dualp_flag))


def _write_gridhead(fortio: FortIO, nx: int, ny: int, nz: int, grid_nr: int) -> None:
    write_int_kw(fortio, GRIDHEAD_KW, _gridhead_values(nx, ny, nz, grid_nr))


def _write_grid_body(fortio: FortIO, grid: Grid) -> None:
    grid.export_coord().fwrite(fortio)
    grid.export_zcorn().fwrite(fortio)
    grid.export_actnum().fwrite(fortio)


def _write_nnc_pair(
    fortio: FortIO,
    lgr_nr: int,
    first_kw: str,
    second_kw: str,
    first: Sequence[int],
    second: Sequence[int],
) -> None:
    if not first:
        return
    nnchead = [0] * NNCHEAD_SIZE
    nnchead[NNCHEAD_NUMNNC_INDEX] = len(first)
    nnchead[NNCHEAD_LGR_INDEX] = lgr_nr
    write_int_kw(fortio, NNCHEAD_KW, nnchead)
    write_int_kw(fortio, first_kw, list(first))
    write_int_kw(fortio, second_kw, list(second))


def _write_lgr_egrid_section(
    fortio: FortIO,
    name: str,
    parent: str,
    grid_nr: int,
    lgr_grid: Grid,
    host_global_one_based: int,
    hostnum: Sequence[int] | None = None,
) -> None:
    write_char_kw(fortio, LGR_KW, [name])
    write_char_kw(fortio, LGR_PARENT_KW, [parent])
    _write_gridhead(fortio, lgr_grid.nx, lgr_grid.ny, lgr_grid.nz, grid_nr)
    _write_grid_body(fortio, lgr_grid)

    lgr_size = lgr_grid.get_global_size()
    if hostnum is None:
        hostnum = [host_global_one_based] * lgr_size
    write_int_kw(fortio, HOSTNUM_KW, list(hostnum))
    write_empty_kw(fortio, ENDGRID_KW)
    write_empty_kw(fortio, ENDLGR_KW)


def make_rectangular_grid(
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    dy: float,
    dz: float,
    actnum: Sequence[int] | None = None,
) -> Grid:
    return GridGenerator.create_rectangular((nx, ny, nz), (dx, dy, dz), actnum)


def write_egrid_with_single_lgr(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    lgr_nx: int,
    lgr_ny: int,
    lgr_nz: int,
    host_i: int,
    host_j: int,
    host_k: int,
    lgr_name: str,
    mapaxes: Sequence[float] | None = None,
    nncg: Sequence[int] = (),
    nncl: Sequence[int] = (),
    lgr_nr: int = 1,
    hostnum: Sequence[int] | None = None,
) -> None:
    main_grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0)
    lgr_grid = make_rectangular_grid(
        lgr_nx, lgr_ny, lgr_nz, 1.0 / lgr_nx, 1.0 / lgr_ny, 1.0 / lgr_nz
    )
    host_global = host_i + host_j * nx + host_k * nx * ny

    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f)
        if mapaxes is not None:
            write_float_kw(f, MAPAXES_KW, list(mapaxes))
        _write_gridhead(f, nx, ny, nz, 0)
        _write_grid_body(f, main_grid)
        write_empty_kw(f, ENDGRID_KW)

        _write_lgr_egrid_section(
            f, lgr_name, "", lgr_nr, lgr_grid, host_global + 1, hostnum
        )
        _write_nnc_pair(f, lgr_nr, NNCG_KW, NNCL_KW, nncg, nncl)


def write_egrid_with_nested_lgr(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    host_i: int,
    host_j: int,
    host_k: int,
    outer_nx: int,
    outer_ny: int,
    outer_nz: int,
    outer_name: str,
    inner_host_i: int,
    inner_host_j: int,
    inner_host_k: int,
    inner_nx: int,
    inner_ny: int,
    inner_nz: int,
    inner_name: str,
) -> None:
    main_grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0)
    outer_grid = make_rectangular_grid(
        outer_nx,
        outer_ny,
        outer_nz,
        1.0 / outer_nx,
        1.0 / outer_ny,
        1.0 / outer_nz,
    )
    inner_grid = make_rectangular_grid(
        inner_nx,
        inner_ny,
        inner_nz,
        1.0 / (outer_nx * inner_nx),
        1.0 / (outer_ny * inner_ny),
        1.0 / (outer_nz * inner_nz),
    )

    outer_host_global = host_i + host_j * nx + host_k * nx * ny
    inner_host_global = (
        inner_host_i + inner_host_j * outer_nx + inner_host_k * outer_nx * outer_ny
    )

    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f)
        _write_gridhead(f, nx, ny, nz, 0)
        _write_grid_body(f, main_grid)
        write_empty_kw(f, ENDGRID_KW)

        _write_lgr_egrid_section(
            f,
            outer_name,
            "",
            1,
            outer_grid,
            outer_host_global + 1,
        )
        _write_lgr_egrid_section(
            f,
            inner_name,
            outer_name,
            2,
            inner_grid,
            inner_host_global + 1,
        )


def write_egrid_with_two_lgrs_and_amalgamated_nnc(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    lgr1_name: str,
    host1_i: int,
    host1_j: int,
    host1_k: int,
    lgr2_name: str,
    host2_i: int,
    host2_j: int,
    host2_k: int,
    nna1: Sequence[int],
    nna2: Sequence[int],
) -> None:
    main_grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0)
    lgr1 = make_rectangular_grid(2, 2, 2, 0.5, 0.5, 0.5)
    lgr2 = make_rectangular_grid(2, 2, 2, 0.5, 0.5, 0.5)

    h1 = host1_i + host1_j * nx + host1_k * nx * ny
    h2 = host2_i + host2_j * nx + host2_k * nx * ny

    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f)
        _write_gridhead(f, nx, ny, nz, 0)
        _write_grid_body(f, main_grid)
        write_empty_kw(f, ENDGRID_KW)

        _write_lgr_egrid_section(f, lgr1_name, "", 1, lgr1, h1 + 1)
        _write_lgr_egrid_section(f, lgr2_name, "", 2, lgr2, h2 + 1)

        if nna1:
            nncheada = [0] * NNCHEAD_SIZE
            nncheada[NNCHEADA_ILOC1_INDEX] = 1
            nncheada[NNCHEADA_ILOC2_INDEX] = 2
            write_int_kw(f, NNCHEADA_KW, nncheada)
            write_int_kw(f, NNA1_KW, list(nna1))
            write_int_kw(f, NNA2_KW, list(nna2))


def write_egrid_with_actnum(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    actnum: Sequence[int],
) -> None:
    """Write an EGRID of size nx*ny*nz whose
    ACTNUM has exactly the given values, regardless of the grid size."""
    body_grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0)
    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f)
        _write_gridhead(f, nx, ny, nz, 0)
        body_grid.export_coord().fwrite(f)
        body_grid.export_zcorn().fwrite(f)
        write_int_kw(f, ACTNUM_KW, list(actnum))
        write_empty_kw(f, ENDGRID_KW)


def write_egrid_with_coarse_groups(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    corsnum: Sequence[int],
    actnum: Sequence[int] | None = None,
) -> None:
    grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0, actnum)

    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f)
        _write_gridhead(f, nx, ny, nz, 0)
        _write_grid_body(f, grid)
        write_int_kw(f, CORSNUM_KW, list(corsnum))
        write_empty_kw(f, ENDGRID_KW)


def write_egrid_dual_porosity(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    actnum: Sequence[int],
    nnc1: Sequence[int] = (),
    nnc2: Sequence[int] = (),
    corsnum: Sequence[int] | None = None,
) -> None:
    """Write a dual-porosity EGRID. ACTNUM uses CELL_ACTIVE_MATRIX/FRACTURE bits."""
    grid = make_rectangular_grid(nx, ny, nz, 1.0, 1.0, 1.0)

    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_filehead(f, FILEHEAD_DUAL_POROSITY)
        _write_gridhead(f, nx, ny, nz, 0)
        grid.export_coord().fwrite(f)
        grid.export_zcorn().fwrite(f)
        write_int_kw(f, ACTNUM_KW, list(actnum))
        if corsnum is not None:
            write_int_kw(f, CORSNUM_KW, list(corsnum))
        write_empty_kw(f, ENDGRID_KW)
        _write_nnc_pair(f, 0, NNC1_KW, NNC2_KW, nnc1, nnc2)


def _extruded_unit_corners(i: int, j: int, z0: float, z1: float) -> list:
    return [
        float(i),
        float(j),
        z0,
        float(i + 1),
        float(j),
        z0,
        float(i),
        float(j + 1),
        z0,
        float(i + 1),
        float(j + 1),
        z0,
        float(i),
        float(j),
        z1,
        float(i + 1),
        float(j),
        z1,
        float(i),
        float(j + 1),
        z1,
        float(i + 1),
        float(j + 1),
        z1,
    ]


def _write_grid_dimens(f: FortIO, nx: int, ny: int, nz: int) -> None:
    write_int_kw(f, DIMENS_KW, [nx, ny, nz])


def _write_grid_radial_false(f: FortIO) -> None:
    write_char_kw(f, RADIAL_KW, ["FALSE"])


def _write_grid_mapunits(f: FortIO) -> None:
    write_char_kw(f, MAPUNITS_KW, ["METRES"])


def _write_grid_gridunit(f: FortIO) -> None:
    write_char_kw(f, GRIDUNIT_KW, ["METRES", ""])


def _write_grid_cell(
    f: FortIO,
    i: int,
    j: int,
    k: int,
    global_index: int,
    host_cell: int,
    corners: Sequence[float],
) -> None:
    coords = [i + 1, j + 1, k + 1, global_index + 1, 1, host_cell, 0]
    write_int_kw(f, COORDS_KW, coords)
    write_float_kw(f, CORNERS_KW, list(corners))


def write_grid_file_with_lgrs(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    lgrs: Sequence[dict] = (),
    mapaxes: Sequence[float] | None = None,
) -> None:
    """LGR dict keys: lgr_name, parent_name, emit_parent, nx, ny, nz, host_cell."""
    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_grid_dimens(f, nx, ny, nz)
        _write_grid_mapunits(f)
        if mapaxes is not None:
            write_float_kw(f, MAPAXES_KW, list(mapaxes))
        _write_grid_gridunit(f)
        _write_grid_radial_false(f)

        global_index = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    corners = _extruded_unit_corners(i, j, float(k), float(k + 1))
                    _write_grid_cell(f, i, j, k, global_index, 0, corners)
                    global_index += 1

        for lgr in lgrs:
            if lgr.get("emit_parent", False):
                write_char_kw(
                    f,
                    LGR_KW,
                    [lgr["lgr_name"], lgr.get("parent_name", "")],
                )
            else:
                write_char_kw(f, LGR_KW, [lgr["lgr_name"]])

            _write_grid_dimens(f, lgr["nx"], lgr["ny"], lgr["nz"])
            _write_grid_radial_false(f)
            lgr_global = 0
            for k in range(lgr["nz"]):
                for j in range(lgr["ny"]):
                    for i in range(lgr["nx"]):
                        _write_grid_cell(
                            f,
                            i,
                            j,
                            k,
                            lgr_global,
                            lgr.get("host_cell", 1),
                            UNIT_CELL_CORNERS,
                        )
                        lgr_global += 1


def write_grid_file_with_mapaxes(filename: Path, mapaxes: Sequence[float]) -> None:
    with openFortIO(str(filename), mode=FortIO.WRITE_MODE) as f:
        _write_grid_dimens(f, 1, 1, 1)
        _write_grid_mapunits(f)
        write_float_kw(f, MAPAXES_KW, list(mapaxes))
        _write_grid_gridunit(f)
        _write_grid_radial_false(f)
        write_int_kw(f, COORDS_KW, [1, 1, 1, 1, 1])
        write_float_kw(f, CORNERS_KW, list(UNIT_CELL_CORNERS))


def write_fegrid_minimal(filename: Path) -> None:
    grid = make_rectangular_grid(1, 1, 1, 1.0, 1.0, 1.0)
    with openFortIO(str(filename), mode=FortIO.WRITE_MODE, fmt_file=True) as f:
        _write_filehead(f)
        _write_gridhead(f, 1, 1, 1, 0)
        _write_grid_body(f, grid)
        write_empty_kw(f, ENDGRID_KW)


def _zcorn_index(nx: int, ny: int, i: int, j: int, k: int, c: int) -> int:
    idx = k * 8 * nx * ny + j * 4 * nx + 2 * i
    if c % 2 == 1:
        idx += 1
    if c % 4 == 2:
        idx += 2 * nx
    if c % 4 == 3:
        idx += 2 * nx
    if c >= 4:
        idx += 4 * nx * ny
    return idx


def _set_pillar(
    coord_kw: ResdataKW,
    pillar: int,
    tx: float,
    ty: float,
    tz: float,
    bx: float,
    by: float,
    bz: float,
) -> None:
    off = 6 * pillar
    coord_kw[off + 0] = tx
    coord_kw[off + 1] = ty
    coord_kw[off + 2] = tz
    coord_kw[off + 3] = bx
    coord_kw[off + 4] = by
    coord_kw[off + 5] = bz


def build_single_cell_grid(
    corners: Sequence[Sequence[float]],
    actnum_value: int = 1,
) -> Grid:
    """Build a 1x1x1 corner-point grid from 8 (x,y,z) corner tuples.

    Corner ordering follows EGRID ZCORN order:
       0,1,2,3 = top corners (i_min/i_max x j_min/j_max)
       4,5,6,7 = bottom corners
    """
    nx, ny, nz = 1, 1, 1

    coord_kw = ResdataKW(COORD_KW, coord_size(nx, ny), ResDataType.RD_FLOAT)
    for i in range(coord_size(nx, ny)):
        coord_kw[i] = 0.0
    zcorn_kw = ResdataKW(ZCORN_KW, zcorn_size(nx, ny, nz), ResDataType.RD_FLOAT)
    for i in range(zcorn_size(nx, ny, nz)):
        zcorn_kw[i] = 0.0
    actnum_kw = ResdataKW(ACTNUM_KW, nx * ny * nz, ResDataType.RD_INT)
    actnum_kw[0] = actnum_value

    top_corner = [0, 1, 2, 3]
    bot_corner = [4, 5, 6, 7]
    for j in (0, 1):
        for i in (0, 1):
            pillar = i + j * (nx + 1)
            top = top_corner[i + 2 * j]
            bot = bot_corner[i + 2 * j]
            _set_pillar(
                coord_kw,
                pillar,
                corners[top][0],
                corners[top][1],
                corners[top][2],
                corners[bot][0],
                corners[bot][1],
                corners[bot][2],
            )

    for c in range(8):
        zi = _zcorn_index(nx, ny, 0, 0, 0, c)
        zcorn_kw[zi] = corners[c][2]

    specgrid = make_int_kw("SPECGRID", [nx, ny, nz])
    return Grid.create(specgrid, zcorn_kw, coord_kw, actnum_kw)


def generate_coordkw_grid(
    num_x: int,
    num_y: int,
    num_z: int,
    z_perturbations: Iterable[tuple[int, int, int, int, float]] = (),
    actnum: list[int] | None = None,
) -> Grid:
    """Build a grid via COORD/ZCORN with integer-coordinate corners.

    z_perturbations is a list of (i, j, k, corner, z) overrides applied to
    ZCORN after the regular pattern is laid out.
    """
    coord_kw = ResdataKW(COORD_KW, coord_size(num_x, num_y), ResDataType.RD_FLOAT)
    for i in range(coord_size(num_x, num_y)):
        coord_kw[i] = 0.0
    zcorn_kw = ResdataKW(
        ZCORN_KW,
        zcorn_size(num_x, num_y, num_z),
        ResDataType.RD_FLOAT,
    )

    for j in range(num_y):
        for i in range(num_x):
            _set_pillar(coord_kw, i + j * num_x, i, j, -1, i, j, -1)
            for k in range(num_z):
                for c in range(4):
                    zi1 = _zcorn_index(num_x, num_y, i, j, k, c)
                    zi2 = _zcorn_index(num_x, num_y, i, j, k, c + 4)
                    zcorn_kw[zi1] = float(k)
                    zcorn_kw[zi2] = float(k + 1)

    for i, j, k, c, z in z_perturbations:
        zcorn_kw[_zcorn_index(num_x, num_y, i, j, k, c)] = float(z)

    specgrid = make_int_kw("SPECGRID", [num_x, num_y, num_z])
    actnum_kw = None
    if actnum is not None:
        actnum_kw = make_int_kw("ACTNUM", actnum)
    return Grid.create(specgrid, zcorn_kw, coord_kw, actnum_kw)


def load_egrid_with_single_lgr(filename: Path, *args, **kwargs) -> Grid:
    write_egrid_with_single_lgr(filename, *args, **kwargs)
    return Grid(str(filename))


def load_egrid_with_nested_lgr(filename: Path, *args, **kwargs) -> Grid:
    write_egrid_with_nested_lgr(filename, *args, **kwargs)
    return Grid(str(filename))


def load_egrid_with_two_lgrs_and_amalgamated_nnc(
    filename: Path,
    *args,
    **kwargs,
) -> Grid:
    write_egrid_with_two_lgrs_and_amalgamated_nnc(filename, *args, **kwargs)
    return Grid(str(filename))


def load_egrid_with_coarse_groups(filename: Path, *args, **kwargs) -> Grid:
    write_egrid_with_coarse_groups(filename, *args, **kwargs)
    return Grid(str(filename))


def load_egrid_dual_porosity(filename: Path, *args, **kwargs) -> Grid:
    write_egrid_dual_porosity(filename, *args, **kwargs)
    return Grid(str(filename))


def load_egrid_dual_porosity_with_coarse_groups(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
    actnum: Sequence[int],
    corsnum: Sequence[int],
) -> Grid:
    write_egrid_dual_porosity(filename, nx, ny, nz, actnum, corsnum=corsnum)
    return Grid(str(filename))


def load_grid_file_with_lgr_parent(
    filename: Path,
    lgr_name: str,
    parent_name: str,
) -> Grid:
    write_grid_file_with_lgrs(
        filename,
        1,
        1,
        1,
        [
            {
                "lgr_name": lgr_name,
                "parent_name": parent_name,
                "emit_parent": True,
                "nx": 1,
                "ny": 1,
                "nz": 1,
                "host_cell": 1,
            }
        ],
    )
    return Grid(str(filename))


def load_grid_file_with_mapaxes(filename: Path, mapaxes: Sequence[float]) -> Grid:
    write_grid_file_with_mapaxes(filename, mapaxes)
    return Grid(str(filename))


def load_grid_file_main_with_lgr(
    filename: Path,
    nx: int,
    ny: int,
    nz: int,
) -> Grid:
    write_grid_file_with_lgrs(
        filename,
        nx,
        ny,
        nz,
        [
            {
                "lgr_name": "LGR1",
                "parent_name": "",
                "emit_parent": False,
                "nx": 1,
                "ny": 1,
                "nz": 1,
                "host_cell": 1,
            }
        ],
    )
    return Grid(str(filename))


def load_grid_file_with_nested_lgr(
    filename: Path,
    outer_name: str,
    inner_name: str,
) -> Grid:
    write_grid_file_with_lgrs(
        filename,
        1,
        1,
        1,
        [
            {
                "lgr_name": outer_name,
                "parent_name": "",
                "emit_parent": False,
                "nx": 1,
                "ny": 1,
                "nz": 1,
                "host_cell": 1,
            },
            {
                "lgr_name": inner_name,
                "parent_name": outer_name,
                "emit_parent": True,
                "nx": 1,
                "ny": 1,
                "nz": 1,
                "host_cell": 1,
            },
        ],
    )
    return Grid(str(filename))


def load_fegrid_minimal(filename: Path) -> Grid:
    write_fegrid_minimal(filename)
    return Grid(str(filename))


@pytest.fixture
def all_grids(tmp_path_factory) -> list[tuple[str, Grid]]:
    dirname = tmp_path_factory.mktemp("all_grids")
    grids = []
    grids.append(("rect-2x2x2", make_rectangular_grid(2, 2, 2, 1, 1, 1)))
    grids.append(
        (
            "rect-2x2x2-inactive",
            make_rectangular_grid(2, 2, 2, 1, 1, 1, [1, 1, 1, 1, 0, 1, 1, 1]),
        )
    )
    grids.append(("rect-3x3x3", make_rectangular_grid(3, 3, 3, 1, 1, 1)))
    grids.append(
        (
            "coordkw-perturbed",
            generate_coordkw_grid(2, 2, 2, [(1, 1, 1, 7, 42.0)]),
        )
    )
    grids.append(
        (
            "egrid-single-lgr",
            load_egrid_with_single_lgr(
                dirname / "ALL_GRIDS_LGR.EGRID",
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
            ),
        )
    )
    grids.append(
        (
            "egrid-nested-lgr",
            load_egrid_with_nested_lgr(
                dirname / "ALL_GRIDS_NESTED.EGRID",
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
            ),
        )
    )
    grids.append(
        (
            "egrid-amalgamated-lgrs",
            load_egrid_with_two_lgrs_and_amalgamated_nnc(
                dirname / "ALL_GRIDS_AMALG.EGRID",
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
                [1, 2],
                [1, 2],
            ),
        )
    )
    corsnum = [0] * 8
    corsnum[0] = 1
    corsnum[1] = 1
    grids.append(
        (
            "egrid-coarse",
            load_egrid_with_coarse_groups(
                dirname / "ALL_GRIDS_CORSNUM.EGRID",
                2,
                2,
                2,
                corsnum,
            ),
        )
    )
    actnum = [CELL_ACTIVE_MATRIX | CELL_ACTIVE_FRACTURE] * 8
    grids.append(
        (
            "egrid-dual-porosity",
            load_egrid_dual_porosity(
                dirname / "ALL_GRIDS_DUALP.EGRID",
                2,
                2,
                2,
                actnum,
            ),
        )
    )
    grids.append(
        (
            "egrid-dual-porosity-with-coarse",
            load_egrid_dual_porosity_with_coarse_groups(
                dirname / "ALL_GRIDS_DUALP_CORS.EGRID",
                2,
                2,
                2,
                actnum,
                corsnum,
            ),
        )
    )
    grids.append(
        (
            "fegrid-minimal",
            load_fegrid_minimal(dirname / "ALL_GRIDS.FEGRID"),
        )
    )
    grids.append(
        (
            "grid-mapaxes",
            load_grid_file_with_mapaxes(
                dirname / "ALL_GRIDS_MAPAXES.GRID",
                [1.0, 11.0, 1.0, 1.0, 11.0, 1.0],
            ),
        )
    )
    grids.append(
        (
            "grid-lgr-parent",
            load_grid_file_with_lgr_parent(
                dirname / "ALL_GRIDS_LGRPARENT.GRID",
                "LGR1",
                "GLOBAL",
            ),
        )
    )
    grids.append(
        (
            "grid-main-with-lgr",
            load_grid_file_main_with_lgr(
                dirname / "ALL_GRIDS_MAIN_LGR.GRID",
                1,
                1,
                2,
            ),
        )
    )
    grids.append(
        (
            "grid-nested-lgr",
            load_grid_file_with_nested_lgr(
                dirname / "ALL_GRIDS_GRID_NESTED.GRID",
                "OUTER",
                "INNER",
            ),
        )
    )
    return grids
