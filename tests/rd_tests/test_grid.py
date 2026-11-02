#!/usr/bin/env python
import functools
import itertools
import os.path
import re
from contextlib import suppress
from unittest import skip

import cwrap
import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import given
from resdata import ResDataType, UnitSystem
from resdata.grid import Grid
from resdata.grid import GridGenerator as GridGen
from resdata.resfile import ResdataFile, ResdataKW
from resdata.util.util import DoubleVector, IntVector

from tests import ResdataTest
from tests.util import TestAreaContext

# This dict is used to verify that corners are mapped to the correct
# cell with respect to containment.
CORNER_HOME = {
    (0, 0, 0): 0,
    (0, 0, 1): 9,
    (0, 0, 2): 18,
    (0, 0, 3): 18,
    (0, 1, 0): 3,
    (0, 1, 1): 12,
    (0, 1, 2): 21,
    (0, 1, 3): 21,
    (0, 2, 0): 6,
    (0, 2, 1): 15,
    (0, 2, 2): 24,
    (0, 2, 3): 24,
    (0, 3, 0): 6,
    (0, 3, 1): 15,
    (0, 3, 2): 24,
    (0, 3, 3): 24,
    (1, 0, 0): 1,
    (1, 0, 1): 10,
    (1, 0, 2): 19,
    (1, 0, 3): 19,
    (1, 1, 0): 4,
    (1, 1, 1): 13,
    (1, 1, 2): 22,
    (1, 1, 3): 22,
    (1, 2, 0): 7,
    (1, 2, 1): 16,
    (1, 2, 2): 25,
    (1, 2, 3): 25,
    (1, 3, 0): 7,
    (1, 3, 1): 16,
    (1, 3, 2): 25,
    (1, 3, 3): 25,
    (2, 0, 0): 2,
    (2, 0, 1): 11,
    (2, 0, 2): 20,
    (2, 0, 3): 20,
    (2, 1, 0): 5,
    (2, 1, 1): 14,
    (2, 1, 2): 23,
    (2, 1, 3): 23,
    (2, 2, 0): 8,
    (2, 2, 1): 17,
    (2, 2, 2): 26,
    (2, 2, 3): 26,
    (2, 3, 0): 8,
    (2, 3, 1): 17,
    (2, 3, 2): 26,
    (2, 3, 3): 26,
    (3, 0, 0): 2,
    (3, 0, 1): 11,
    (3, 0, 2): 20,
    (3, 0, 3): 20,
    (3, 1, 0): 5,
    (3, 1, 1): 14,
    (3, 1, 2): 23,
    (3, 1, 3): 23,
    (3, 2, 0): 8,
    (3, 2, 1): 17,
    (3, 2, 2): 26,
    (3, 2, 3): 26,
    (3, 3, 0): 8,
    (3, 3, 1): 17,
    (3, 3, 2): 26,
    (3, 3, 3): 26,
}


def createVolumeTestGridBase(dim, dV, offset=1):
    return [
        GridGen.create_grid(dim, dV, offset=0),
        GridGen.create_grid(dim, dV, offset=offset),
        GridGen.create_grid(dim, dV, offset=offset, irregular_offset=True),
        GridGen.create_grid(dim, dV, offset=offset, concave=True),
        GridGen.create_grid(dim, dV, offset=offset, irregular=True),
        GridGen.create_grid(dim, dV, offset=offset, concave=True, irregular=True),
        GridGen.create_grid(
            dim, dV, offset=offset, irregular_offset=True, concave=True
        ),
        GridGen.create_grid(dim, dV, offset=0, faults=True),
        GridGen.create_grid(dim, dV, offset=offset, faults=True),
        GridGen.create_grid(dim, dV, escape_origo_shift=(100, 100, 0), scale=2),
        GridGen.create_grid(dim, dV, escape_origo_shift=(100, 100, 0), scale=0.5),
        GridGen.create_grid(
            dim, dV, escape_origo_shift=(100, 100, 0), translation=(50, 50, 0)
        ),
        GridGen.create_grid(dim, dV, escape_origo_shift=(100, 100, 0), rotate=True),
        GridGen.create_grid(dim, dV, escape_origo_shift=(100, 100, 0), misalign=True),
        GridGen.create_grid(
            dim,
            dV,
            offset=offset,
            escape_origo_shift=(100, 100, 0),
            irregular_offset=True,
            concave=True,
            irregular=True,
            scale=1.5,
            translation=(5, 5, 0),
            rotate=True,
            misalign=True,
        ),
    ]


def createContainmentTestBase():
    return [
        (3, GridGen.create_grid((6, 6, 6), (1, 1, 1), offset=1)),
        (10, GridGen.create_grid((3, 3, 3), (1, 1, 1), offset=1, concave=True)),
        (4, GridGen.create_grid((10, 10, 1), (1, 1, 1), offset=0.0, misalign=True)),
        (
            3,
            GridGen.create_grid(
                (6, 6, 6),
                (1, 1, 1),
                offset=0.0,
                escape_origo_shift=(100, 100, 0),
                irregular_offset=True,
                concave=True,
                irregular=True,
                scale=1.5,
                translation=(5, 5, 0),
                misalign=True,
            ),
        ),
    ]


def getMinMaxValue(grid):
    corners = [
        grid.get_cell_corner(i, cell)
        for i in range(8)
        for cell in range(grid.get_global_size())
    ]

    return [(min(values), max(values)) for values in zip(*corners)]


def createWrapperGrid(grid):
    """
    Creates a grid that occupies the same space as the given grid,
    but that consists of a single cell.
    """

    x, y, z = grid.get_nx() - 1, grid.get_ny() - 1, grid.get_nz() - 1
    corner_pos = [
        (0, 0, 0),
        (x, 0, 0),
        (0, y, 0),
        (x, y, 0),
        (0, 0, z),
        (x, 0, z),
        (0, y, z),
        (x, y, z),
    ]

    corners = [grid.get_cell_corner(i, ijk=pos) for i, pos in enumerate(corner_pos)]

    return GridGen.create_single_cell_grid(corners)


def average(points):
    p = functools.reduce(lambda a, b: (a[0] + b[0], a[1] + b[1], a[2] + b[2]), points)
    return [elem / float(len(points)) for elem in p]


# This test class should only have test cases which do not require
# external test data. Tests involving Equinor test data are in the
# test_grid_equinor module.
class GridTest(ResdataTest):
    def test_oom_grid(self):
        nx = 2000
        ny = 2000
        nz = 2000

        with self.assertRaises(MemoryError):
            grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))

    def test_posXYEdge(self):
        nx = 10
        ny = 11
        grid = GridGen.create_rectangular((nx, ny, 1), (1, 1, 1))
        self.assertEqual(grid.find_cell_corner_xy(0, 0, 0), 0)
        self.assertEqual(grid.find_cell_corner_xy(nx, 0, 0), nx)
        self.assertEqual(grid.find_cell_corner_xy(0, ny, 0), (nx + 1) * ny)
        self.assertEqual(grid.find_cell_corner_xy(nx, ny, 0), (nx + 1) * (ny + 1) - 1)

        self.assertEqual(grid.find_cell_corner_xy(0.25, 0, 0), 0)
        self.assertEqual(grid.find_cell_corner_xy(0, 0.25, 0), 0)

        self.assertEqual(grid.find_cell_corner_xy(nx - 0.25, 0, 0), nx)
        self.assertEqual(grid.find_cell_corner_xy(nx, 0.25, 0), nx)

        self.assertEqual(grid.find_cell_corner_xy(0, ny - 0.25, 0), (nx + 1) * ny)
        self.assertEqual(grid.find_cell_corner_xy(0.25, ny, 0), (nx + 1) * ny)

        self.assertEqual(
            grid.find_cell_corner_xy(nx - 0.25, ny, 0), (nx + 1) * (ny + 1) - 1
        )
        self.assertEqual(
            grid.find_cell_corner_xy(nx, ny - 0.25, 0), (nx + 1) * (ny + 1) - 1
        )

    def test_dims(self):
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1))
        self.assertEqual(grid.get_global_index1F(0), None)
        self.assertTrue(grid.equal(grid))
        self.assertFalse(grid.dual_grid())
        self.assertEqual(grid.get_num_active_fracture(), 0)
        self.assertEqual(grid.get_active_fracture_index(global_index=0), -1)
        self.assertEqual(
            grid.get_bounding_box_2d(),
            ((0.0, 0.0), (10.0, 0.0), (10.0, 20.0), (0.0, 20.0)),
        )
        self.assertEqual(grid.depth(ijk=(0, 0, 0)), 0.5)
        self.assertEqual(grid.top(0, 0), 0.0)
        self.assertEqual(grid.top_active(0, 0), 0.0)
        self.assertEqual(grid.bottom(0, 0), 30.0)
        self.assertEqual(grid.locate_depth(1.0, 0, 0), 1)
        self.assertEqual(grid.find_cell(1.1, 0.0, 0.0), (1, 0, 0))
        self.assertTrue(grid.cell_regular(ijk=(1, 0, 0)))
        self.assertEqual(grid.get_layer_xyz(1, 0), (1.0, 0.0, 0.0))
        self.assertEqual(grid.distance(0, 1), (-1.0, 0.0, 0.0))
        self.assertEqual(grid.get_num_lgr(), 0)
        self.assertFalse(grid.has_lgr("lgr"))
        self.assertEqual(grid.coarse_groups(), 0)
        self.assertFalse(grid.in_coarse_group(ijk=(0, 0, 0)))

        with self.assertRaises(KeyError):
            grid.get_lgr("lgr")
        with self.assertRaises(KeyError):
            grid.get_lgr(0)
        with self.assertRaises(IndexError):
            grid.get_cell_lgr(ijk=(0, 0, 0))
        self.assertEqual(grid.get_nx(), 10)
        self.assertEqual(grid.get_ny(), 20)
        self.assertEqual(grid.get_nz(), 30)
        self.assertEqual(grid.nx, 10)
        self.assertEqual(grid.ny, 20)
        self.assertEqual(grid.nz, 30)
        self.assertEqual(grid.get_global_size(), 30 * 10 * 20)

        self.assertEqual(grid.get_dims(), (10, 20, 30, 6000))

    def test_global_index(self):
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1))
        with pytest.raises(
            IndexError, match=r"Invalid value global_index:7000  Range: \[0,6000\)"
        ):
            grid.get_active_index(global_index=7000)

    def test_load_column(self):
        column = DoubleVector(2 * 3 * 4)
        grid = GridGen.create_rectangular((2, 3, 4), (1, 1, 1))
        kw = ResdataKW("KW", 2 * 3 * 4, ResDataType.RD_DOUBLE)
        kw[0] = 1.0
        grid.load_column(kw, 0, 0, column)
        assert list(column) == [1.0, 0.0, 0.0, 0.0]

    def test_create(self):
        with self.assertRaises(ValueError):
            grid = GridGen.create_rectangular(
                (10, 20, 30), (1, 1, 1), actnum=[0, 1, 1, 2]
            )

        with self.assertRaises(ValueError):
            grid = GridGen.create_rectangular(
                (10, 20, 30), (1, 1, 1), actnum=IntVector(initial_size=10)
            )
        grid = GridGen.create_rectangular(
            (10, 20, 30), (1, 1, 1)
        )  # actnum=None -> all active
        self.assertEqual(grid.get_num_active(), 30 * 20 * 10)
        actnum = IntVector(default_value=1, initial_size=6000)
        actnum[0] = 0
        actnum[1] = 0
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1), actnum=actnum)
        self.assertEqual(grid.get_num_active(), 30 * 20 * 10 - 2)

    def test_all_iters(self):
        fk = self.createTestPath("local/ECLIPSE/faarikaal/faarikaal1.EGRID")
        grid = Grid(fk)
        cell = grid[3455]
        self.assertEqual(3455, cell.global_index)
        cell = grid[(4, 1, 82)]
        self.assertEqual(3455, cell.global_index)
        self.assertEqual(grid.cell(global_index=3455), grid.cell(active_index=2000))
        self.assertEqual(grid.cell(global_index=3455), grid.cell(i=4, j=1, k=82))

        na = grid.get_num_active()
        self.assertEqual(na, 4160)
        cnt = 0
        for c in grid.cells(active=True):
            cnt += 1
            self.assertTrue(c.active)
        self.assertEqual(cnt, 4160)

        cnt = len([c for c in grid.cells()])
        self.assertEqual(cnt, len(grid))

    def test_repr_and_name(self):
        grid = GridGen.create_rectangular(
            (2, 2, 2), (10, 10, 10), actnum=[0, 0, 0, 0, 1, 1, 1, 1]
        )
        pfx = "Grid("
        rep = repr(grid)
        self.assertEqual(pfx, rep[: len(pfx)])
        self.assertEqual(type(rep), str)
        self.assertEqual(type(grid.get_name()), str)
        with TestAreaContext("python/rd_grid/repr"):
            grid.save_EGRID("CASE.EGRID")
            g2 = Grid("CASE.EGRID")
            r2 = repr(g2)
            self.assertEqual(pfx, r2[: len(pfx)])
            self.assertEqual(type(r2), str)
            self.assertEqual(type(g2.get_name()), str)

    def test_node_pos(self):
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1))
        with self.assertRaises(IndexError):
            grid.get_node_pos(-1, 0, 0)

        with self.assertRaises(IndexError):
            grid.get_node_pos(11, 0, 0)

        p0 = grid.get_node_pos(0, 0, 0)
        self.assertEqual(p0, (0, 0, 0))

        p7 = grid.get_node_pos(10, 20, 30)
        self.assertEqual(p7, (10, 20, 30))

    # The broken file was previously handled by the rd_file_open() call internally
    # in the rd_grid implementation. That will now not fail for a broken file, and then
    # the grid class needs to do more/better checking itself.
    @skip("Needs better error checking inside in the rd_grid")
    def test_truncated_file(self):
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1))
        with TestAreaContext("python/rd_grid/truncated"):
            grid.save_EGRID("TEST.EGRID")

            size = os.path.getsize("TEST.EGRID")
            with open("TEST.EGRID", "r+") as f:
                f.truncate(size // 2)

            with self.assertRaises(IOError):
                Grid("TEST.EGRID")

    def test_posXY1(self):
        nx = 4
        ny = 1
        nz = 1
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        i, j = grid.find_cell_xy(0.5, 0.5, 0)
        self.assertEqual(i, 0)
        self.assertEqual(j, 0)

        i, j = grid.find_cell_xy(3.5, 0.5, 0)
        self.assertEqual(i, 3)
        self.assertEqual(j, 0)

    def test_init_ACTNUM(self):
        nx = 10
        ny = 23
        nz = 7
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        actnum = grid.exportACTNUM()

        self.assertEqual(len(actnum), nx * ny * nz)
        self.assertEqual(actnum[0], 1)
        self.assertEqual(actnum[nx * ny * nz - 1], 1)

        actnum_kw = grid.export_ACTNUM_kw()
        self.assertEqual(len(actnum_kw), len(actnum))
        for a1, a2 in zip(actnum, actnum_kw):
            self.assertEqual(a1, a2)

    def test_posXY(self):
        nx = 10
        ny = 23
        nz = 7
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        with self.assertRaises(IndexError):
            grid.find_cell_xy(1, 1, -1)

        with self.assertRaises(IndexError):
            grid.find_cell_xy(1, 1, nz + 1)

        with pytest.raises(
            ValueError, match=re.escape("Could not find the point:(15,78) in layer:2")
        ):
            grid.find_cell_xy(15, 78, 2)

        i, j = grid.find_cell_xy(1.5, 1.5, 2)
        self.assertEqual(i, 1)
        self.assertEqual(j, 1)

        for i in range(nx):
            for j in range(ny):
                p = grid.find_cell_xy(i + 0.5, j + 0.5, 0)
                self.assertEqual(p[0], i)
                self.assertEqual(p[1], j)

        c = grid.find_cell_corner_xy(0.10, 0.10, 0)
        self.assertEqual(c, 0)

        c = grid.find_cell_corner_xy(0.90, 0.90, 0)
        self.assertEqual(c, (nx + 1) + 1)

        c = grid.find_cell_corner_xy(0.10, 0.90, 0)
        self.assertEqual(c, (nx + 1))

        c = grid.find_cell_corner_xy(0.90, 0.90, 0)
        self.assertEqual(c, (nx + 1) + 1)

        c = grid.find_cell_corner_xy(0.90, 0.10, 0)
        self.assertEqual(c, 1)

    def test_compressed_copy(self):
        nx = 10
        ny = 10
        nz = 10
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        kw1 = ResdataKW("KW", 1001, ResDataType.RD_INT)
        with self.assertRaises(ValueError):
            cp = grid.compressed_kw_copy(kw1)

    def test_bounding_box_2d_invalid_corners(self):
        grid = GridGen.create_rectangular((10, 20, 30), (1, 1, 1))
        with self.assertRaises(ValueError):
            grid.get_bounding_box_2d(lower_left=(0, 1))
        with self.assertRaises(ValueError):
            grid.get_bounding_box_2d(lower_left=(1, 0))
        with self.assertRaises(ValueError):
            grid.get_bounding_box_2d(upper_right=(1, 20))
        with self.assertRaises(ValueError):
            grid.get_bounding_box_2d(upper_right=(10, 1))

    def test_create_3d_with_active_size_kw_and_invalid_size(self):
        nx, ny, nz = 4, 3, 2
        actnum = ResdataKW("ACTNUM", nx * ny * nz, ResDataType.RD_INT)
        for i in range(nx * ny * nz):
            actnum[i] = 1
        actnum[0] = 0
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1), actnum=actnum)
        num_active = grid.get_num_active()
        kw = ResdataKW("PRESSURE", num_active, ResDataType.RD_FLOAT)
        for i in range(num_active):
            kw[i] = float(i + 1)
        arr = grid.create_3d(kw, default=-1.0)
        self.assertEqual(arr.shape, (nx, ny, nz))
        self.assertAlmostEqual(arr[1, 0, 0], 1.0)

        bad = ResdataKW("BAD", num_active + 5, ResDataType.RD_FLOAT)
        with self.assertRaises(ValueError):
            grid.create_3d(bad)

    def test_dxdydz(self):
        nx = 10
        ny = 10
        nz = 10
        grid = GridGen.create_rectangular((nx, ny, nz), (2, 3, 4))

        dx, dy, dz = grid.get_cell_dims(active_index=0)
        self.assertEqual(dx, 2)
        self.assertEqual(dy, 3)
        self.assertEqual(dz, 4)

    def test_create_3d_is_create_kw_inverse(self):
        nx = 10
        ny = 7
        nz = 5
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        kw1 = ResdataKW("SWAT", nx * ny * nz, ResDataType.RD_FLOAT)
        for k, j, i in itertools.product(range(nz), range(ny), range(nx)):
            kw1[i + j * nx + nx * ny * k] = i * j * k
        numpy_3d = grid.create_3d(kw1)
        kw2 = grid.create_kw(numpy_3d, "SWAT", False)
        self.assertEqual(kw2.name, "SWAT")
        assert np.allclose(grid.create_3d(kw2), numpy_3d)

    def test_create_kw(self):
        rng = np.random.default_rng()
        nx = 10
        ny = 7
        nz = 5
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        array1 = rng.integers(0, 100, size=(nx, ny, nz), dtype=np.int32)
        array2 = rng.normal(10, 2, size=(nx, ny, nz))
        _ = grid.create_kw(array1, "SWAT1", False)
        _ = grid.create_kw(array2, "SWAT2", False)

    def test_create_3d_agrees_with_get_value(self):
        nx = 5
        ny = 3
        nz = 2
        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1))
        kw = ResdataKW("SWAT", nx * ny * nz, ResDataType.RD_FLOAT)
        for k, j, i in itertools.product(range(nz), range(ny), range(nx)):
            kw[i + j * nx + nx * ny * k] = i * j * k
        numpy_3d = grid.create_3d(kw)
        for k, j, i in itertools.product(range(nz), range(ny), range(nx)):
            self.assertAlmostEqual(numpy_3d[i, j, k], grid.grid_value(kw, i, j, k))

    def test_len(self):
        nx = 10
        ny = 11
        nz = 12
        actnum = ResdataKW("ACTNUM", nx * ny * nz, ResDataType.RD_INT)
        actnum[0] = 1
        actnum[1] = 1
        actnum[2] = 1
        actnum[3] = 1

        grid = GridGen.create_rectangular((nx, ny, nz), (1, 1, 1), actnum=actnum)
        self.assertEqual(len(grid), nx * ny * nz)
        self.assertEqual(grid.get_num_active(), 4)

    def test_export(self):
        dims = (3, 3, 3)
        coord = GridGen.create_coord(dims, (1, 1, 1))
        zcorn = GridGen.create_zcorn(dims, (1, 1, 1), offset=0)

        grid = Grid.create(dims, zcorn, coord, None)

        self.assertEqual(zcorn, grid.export_zcorn())
        self.assertEqual(coord, grid.export_coord())

    def test_output_units(self):
        n = 10
        a = 1
        grid = GridGen.create_rectangular((n, n, n), (a, a, a))

        with TestAreaContext("python/rd_grid/units"):
            grid.save_EGRID("CASE.EGRID", output_unit=UnitSystem.FIELD)
            f = ResdataFile("CASE.EGRID")
            g = f["GRIDUNIT"][0]
            self.assertEqual(g[0].strip(), "FEET")
            g2 = Grid("CASE.EGRID")
            self.assertFloatEqual(
                g2.cell_volume(global_index=0), 3.28084 * 3.28084 * 3.28084
            )

            grid.save_EGRID("CASE.EGRID")
            f = ResdataFile("CASE.EGRID")
            g = f["GRIDUNIT"][0]
            self.assertEqual(g[0].strip(), "METRES")

            grid.save_EGRID("CASE.EGRID", output_unit=UnitSystem.LAB)
            f = ResdataFile("CASE.EGRID")
            g = f["GRIDUNIT"][0]
            self.assertEqual(g[0].strip(), "CM")
            g2 = Grid("CASE.EGRID")
            self.assertFloatEqual(g2.cell_volume(global_index=0), 100 * 100 * 100)

    def test_volume(self):
        dim = (10, 10, 10)
        dV = (2, 2, 2)

        grids = createVolumeTestGridBase(dim, dV)
        for grid in grids:
            tot_vol = createWrapperGrid(grid).cell_volume(0)
            cell_volumes = [grid.cell_volume(i) for i in range(grid.get_global_size())]
            self.assertTrue(min(cell_volumes) >= 0)
            self.assertFloatEqual(sum(cell_volumes), tot_vol)

    def test_unique_containment(self):
        test_base = createContainmentTestBase()

        for steps_per_unit, grid in test_base:
            wgrid = createWrapperGrid(grid)

            (xmin, xmax), (ymin, ymax), (zmin, zmax) = getMinMaxValue(wgrid)

            x_space = np.linspace(
                xmin - 1, xmax + 1, int(xmax - xmin + 2) * steps_per_unit + 1
            )
            y_space = np.linspace(
                ymin - 1, ymax + 1, int(ymax - ymin + 2) * steps_per_unit + 1
            )
            z_space = np.linspace(
                zmin - 1, zmax + 1, int(zmax - zmin + 2) * steps_per_unit + 1
            )

            # limit amount of points tested by
            # only testing every 3rd point
            x_space = x_space[0:-1:3]
            y_space = y_space[0:-1:3]
            z_space = z_space[0:-1:3]

            for x, y, z in itertools.product(x_space, y_space, z_space):
                hits = [
                    grid.cell_contains(x, y, z, i)
                    for i in range(grid.get_global_size())
                ].count(True)

                self.assertIn(hits, [0, 1])

                expected = 1 if wgrid.cell_contains(x, y, z, 0) else 0
                self.assertEqual(
                    expected,
                    hits,
                    "Expected %d for (%g,%g,%g), got %d" % (expected, x, y, z, hits),
                )

    def test_cell_corner_containment(self):
        n = 4
        d = 10
        grid = GridGen.create_rectangular((n, n, n), (d, d, d))

        for x, y, z in itertools.product(range(0, n * d + 1, d), repeat=3):
            self.assertEqual(
                1, [grid.cell_contains(x, y, z, i) for i in range(n**3)].count(True)
            )

    def test_cell_corner_containment_compatability(self):
        grid = GridGen.create_rectangular((3, 3, 3), (1, 1, 1))

        for x, y, z in itertools.product(range(4), repeat=3):
            for i in range(27):
                if grid.cell_contains(x, y, z, i):
                    self.assertEqual(CORNER_HOME[(x, y, z)], i)

    def test_cell_face_containment(self):
        n = 4
        d = 10
        grid = GridGen.create_rectangular((n, n, n), (d, d, d))

        for x, y, z in itertools.product(range(d // 2, n * d, d), repeat=3):
            for axis, direction in itertools.product(range(3), [-1, 1]):
                p = [x, y, z]
                p[axis] = p[axis] + direction * d / 2
                self.assertEqual(
                    1,
                    [
                        grid.cell_contains(p[0], p[1], p[2], i) for i in range(n**3)
                    ].count(True),
                )

    # This test generates a cell that is concave on ALL 6 sides
    def test_concave_cell_containment(self):
        points = [
            (5, 5, 5),
            (20, 10, 10),
            (10, 20, 10),
            (25, 25, 5),
            (10, 10, 20),
            (25, 5, 25),
            (5, 25, 25),
            (20, 20, 20),
        ]

        grid = GridGen.create_single_cell_grid(points)

        def assertPoint(p):
            return self.assertTrue(grid.cell_contains(p[0], p[1], p[2], 0))

        def assertNotPoint(p):
            return self.assertFalse(grid.cell_contains(p[0], p[1], p[2], 0))

        # Cell center
        assertPoint(average(points))

        # "Side" center
        assertNotPoint(average(points[0:4:]))
        assertNotPoint(average(points[4:8:]))
        assertNotPoint(average(points[1:8:2]))
        assertNotPoint(average(points[0:8:2]))
        assertNotPoint(average(points[0:8:4] + points[1:8:4]))
        assertNotPoint(average(points[2:8:4] + points[3:8:4]))

        # Corners
        for p in points:
            assertPoint(p)

        # Edges
        edges = (
            [(i, i + 1) for i in range(0, 8, 2)]
            + [(i, i + 2) for i in [0, 1, 4, 5]]
            + [(i, i + 4) for i in range(4)]
            + [(1, 2), (2, 7), (1, 7), (4, 7), (2, 4), (4, 1)]
        )
        for a, b in edges:
            assertPoint(average([points[a], points[b]]))

        # Epsilon inside from corners
        middle_point = average(points)
        for p in points:
            assertPoint(average(20 * [p] + [middle_point]))

        # Espilon outside
        middle_point[2] = 0
        for p in points[0:4:]:
            assertNotPoint(average(20 * [p] + [middle_point]))

        middle_point[2] = 30
        for p in points[4:8:]:
            assertNotPoint(average(20 * [p] + [middle_point]))

    # This test generates a cell that is strictly convex on ALL 6 sides
    def test_concvex_cell_containment(self):
        points = [
            (10, 10, 10),
            (25, 5, 5),
            (5, 25, 5),
            (20, 20, 10),
            (5, 5, 25),
            (20, 10, 20),
            (10, 20, 20),
            (25, 25, 25),
        ]

        grid = GridGen.create_single_cell_grid(points)

        def assertPoint(p):
            return self.assertTrue(grid.cell_contains(p[0], p[1], p[2], 0))

        def assertNotPoint(p):
            return self.assertFalse(grid.cell_contains(p[0], p[1], p[2], 0))

        # Cell center
        assertPoint(average(points))

        # "Side" center
        assertPoint(average(points[0:4:]))
        assertPoint(average(points[4:8:]))
        assertPoint(average(points[1:8:2]))
        assertPoint(average(points[0:8:2]))
        assertPoint(average(points[0:8:4] + points[1:8:4]))
        assertPoint(average(points[2:8:4] + points[3:8:4]))

        # Corners
        for p in points:
            assertPoint(p)

        # Edges
        edges = (
            [(i, i + 1) for i in range(0, 8, 2)]
            + [(i, i + 2) for i in [0, 1, 4, 5]]
            + [(i, i + 4) for i in range(4)]
            + [(1, 2), (2, 7), (1, 7), (4, 7), (2, 4), (4, 1)]
        )
        for a, b in edges:
            assertPoint(average([points[a], points[b]]))

        # Epsilon inside from corners
        middle_point = average(points)
        for p in points:
            assertPoint(average(20 * [p] + [middle_point]))

        # Espilon outside
        middle_point[2] = 0
        for p in points[0:4:]:
            assertNotPoint(average(20 * [p] + [middle_point]))

        middle_point[2] = 30
        for p in points[4:8:]:
            assertNotPoint(average(20 * [p] + [middle_point]))


def test_save_grdecl(tmpdir):
    grid = GridGen.create_rectangular((2, 3, 4), (1, 1, 1))
    with tmpdir.as_cwd():
        with cwrap.open("grid.grdecl", "w") as f:
            grid.save_grdecl(f)
        assert grid.equal(Grid.load_from_grdecl("grid.grdecl"))
        assert grid.equal(Grid.load_from_file("grid.grdecl"))


def test_save_grid(tmpdir):
    grid = GridGen.create_rectangular((2, 3, 4), (1, 1, 1))
    with tmpdir.as_cwd():
        grid.save_GRID("grid.GRID")
        assert grid.equal(Grid("grid.GRID"))


def test_write_grdecl(tmpdir):
    kw = ResdataKW("KW", 2 * 3 * 4, ResDataType.RD_DOUBLE)
    grid = GridGen.create_rectangular((2, 3, 4), (1, 1, 1))

    with tmpdir.as_cwd():
        with cwrap.open("kw.grdecl", "w") as f:
            grid.write_grdecl(kw, f)
        with cwrap.open("kw.grdecl", "r") as f:
            kw2 = ResdataKW.read_grdecl(f, "KW")

        assert list(kw) == list(kw2)


def test_create_volume_keyword():
    grid = GridGen.create_rectangular((2, 3, 4), (1, 1, 1))
    assert list(grid.create_volume_keyword()) == [1.0] * (2 * 3 * 4)


invalid_coord_dimensions_grdecl = """
MAPUNITS
 'METRES  '
/
GRIDUNIT
 'METRES  ' '        '
/
SPECGRID
  1  1  1  1  F
/
COORD
   0.00000000E+00
/
ZCORN
   0.00000000E+00
/
"""


def test_that_grdecl_with_invalid_coord_dimensions_raises(tmp_path):
    filename = tmp_path / "invalid.grdecl"
    filename.write_text(invalid_coord_dimensions_grdecl)
    with pytest.raises(ValueError, match="size of COORD keyword = 1"):
        Grid.load_from_file(str(filename))


invalid_zcorn_dimensions_grdecl = """
MAPUNITS
 'METRES  '
/
GRIDUNIT
 'METRES  ' '        '
/
SPECGRID
  1  1  1  1  F
/
COORD
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
/
ZCORN
   0.00000000E+00
/
"""


def test_that_grdecl_with_zcorn_dimensions_raises(tmp_path):
    filename = tmp_path / "invalid.grdecl"
    filename.write_text(invalid_zcorn_dimensions_grdecl)
    with pytest.raises(ValueError, match="size of ZCORN keyword = 1"):
        Grid.load_from_file(str(filename))


invalid_actnum_dimensions_grdecl = """
MAPUNITS
 'METRES  '
/
GRIDUNIT
 'METRES  ' '        '
/
SPECGRID
  1  1  1  1  F
/
COORD
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
/
ZCORN
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
/
ACTNUM
   0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
/
"""


def test_that_grdecl_with_actnum_dimensions_raises(tmp_path):
    filename = tmp_path / "invalid.grdecl"
    filename.write_text(invalid_actnum_dimensions_grdecl)
    with pytest.raises(ValueError, match="size of ACTNUM keyword = 4"):
        Grid.load_from_file(str(filename))


names = st.text(
    min_size=1,
    max_size=7,
    alphabet=st.characters(min_codepoint=40, max_codepoint=126),
).map(lambda x: x.rjust(8))


@st.composite
def kws(
    draw,
    name=names,
    size=st.integers(min_value=1, max_value=32),
    datatype=st.sampled_from(
        [
            ResDataType.RD_INT,
            ResDataType.RD_FLOAT,
            ResDataType.RD_DOUBLE,
            ResDataType.RD_BOOL,
            ResDataType.RD_CHAR,
            ResDataType.RD_STRING,
        ]
    ),
    elements=None,
):
    kw = draw(st.builds(ResdataKW, name, size, datatype))
    for i in range(len(kw)):
        if elements is not None:
            kw[i] = draw(elements)
        elif kw.data_type.is_int():
            kw[i] = draw(st.integers())
        elif kw.data_type.is_float():
            kw[i] = draw(st.floats(width=32, allow_nan=False, allow_infinity=False))
        elif kw.data_type.is_double():
            kw[i] = draw(st.floats(width=64, allow_nan=False, allow_infinity=False))
        elif kw.data_type.is_bool():
            kw[i] = draw(st.booleans())
        else:
            kw[i] = draw(names)
    return kw


@given(
    specgrid=kws(
        name=st.just("SPECGRID"),
        datatype=st.just(ResDataType.RD_INT),
        elements=st.integers(min_value=1, max_value=3),
    ),
    zcorn=kws(name=st.just("ZCORN   "), datatype=st.just(ResDataType.RD_FLOAT)),
    coord=kws(name=st.just("COORD   "), datatype=st.just(ResDataType.RD_FLOAT)),
    actnum=kws(name=st.just("ACTNUM  "), datatype=st.just(ResDataType.RD_INT)),
    mapaxes=kws(name=st.just("MAPAXES "), datatype=st.just(ResDataType.RD_FLOAT)),
)
def test_grid_create(specgrid, zcorn, coord, actnum, mapaxes):
    with suppress(ValueError, IndexError):
        Grid.create(specgrid, zcorn, coord, actnum, mapaxes)
