/**
 * Tests coarse cell groups.
 *
 * A coarse cell group identifies a contiguous block of cells that the
 * simulator treats as a single lumped cell. This is encoded by the
 * CORSNUM keyword in an EGRID file: one int per cell, where 0 means "not in
 * any coarse group" and positive values are group identifiers.
 */

#include <catch2/catch_test_macros.hpp>
#include <resdata/rd_grid.hpp>
#include <resdata/rd_kw.hpp>
#include <resdata/rd_kw_magic.hpp>
#include <limits>
#include <vector>

#include "grid_fixtures.hpp"
#include "tmpdir.hpp"

TEST_CASE_METHOD(Tmpdir, "Load EGRID with coarse cell groups", "[unittest]") {
    GIVEN("An EGRID file on disc with two coarse cell groups") {
        const int nx = 3, ny = 3, nz = 3;
        const int size = nx * ny * nz;

        // Assign two coarse groups. 0 means the cell is not in a
        // coarse group. Group 1 covers the two cells along i at (0,0,0)
        // and (1,0,0); group 2 covers two cells in the last k layer.
        std::vector<int> corsnum(size, 0);
        corsnum[0] = 1;
        corsnum[1] = 1;
        corsnum[size - 1] = 2;
        corsnum[size - 2] = 2;

        auto filename = dirname / "CORSNUM.EGRID";
        auto grid =
            load_egrid_with_coarse_groups(filename, nx, ny, nz, corsnum.data());

        REQUIRE(grid != nullptr);

        THEN("The grid reports that coarse cells are present") {
            REQUIRE(rd_grid_have_coarse_cells(grid.get()));
            REQUIRE(rd_grid_get_num_coarse_groups(grid.get()) == 2);
        }

        THEN("Cells in a coarse group are flagged as such, others are not") {
            REQUIRE(rd_grid_cell_in_coarse_group1(grid.get(), 0));
            REQUIRE(rd_grid_cell_in_coarse_group1(grid.get(), 1));
            REQUIRE_FALSE(rd_grid_cell_in_coarse_group1(grid.get(), 2));
            REQUIRE(rd_grid_cell_in_coarse_group1(grid.get(), size - 1));
            REQUIRE(rd_grid_cell_in_coarse_group1(grid.get(), size - 2));
        }

        THEN("The coarse group objects are accessible by index") {
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 0) != nullptr);
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 1) != nullptr);
        }

        THEN("The grid reports no fracture cells") {
            REQUIRE(rd_grid_get_active_fracture_index1(grid.get(), 0) == -1);
            REQUIRE(rd_grid_get_active_fracture_index1(grid.get(), size / 2) ==
                    -1);
            REQUIRE(rd_grid_get_global_index1F(grid.get(), 0) == -1);
        }

        THEN("The grid reports a consistent (empty) LGR tree") {
            REQUIRE(rd_grid_test_lgr_consistency(grid.get()));
        }

        THEN("rd_grid_init_actnum_data returns the actnum") {
            std::vector<int> actnum(size);
            rd_grid_init_actnum_data(grid.get(), actnum.data());
            for (int i = 0; i < size; i++)
                REQUIRE(actnum[i] ==
                        (rd_grid_cell_active1(grid.get(), i) ? 1 : 0));
        }
    }
}

TEST_CASE_METHOD(Tmpdir, "Load EGRID with non-consecutive coarse group numbers",
                 "[unittest]") {
    // CORSNUM values are allowed to have gap, we return null for skipped slots.
    GIVEN("An EGRID file whose CORSNUM uses groups 1 and 3 but not 2") {
        const int nx = 3, ny = 3, nz = 3;
        const int size = nx * ny * nz;

        std::vector<int> corsnum(size, 0);
        corsnum[0] = 1;
        corsnum[1] = 1;
        corsnum[size - 1] = 3;
        corsnum[size - 2] = 3;

        auto filename = dirname / "CORSNUM_SPARSE.EGRID";
        auto grid =
            load_egrid_with_coarse_groups(filename, nx, ny, nz, corsnum.data());
        REQUIRE(grid != nullptr);

        THEN("rd_grid_get_num_coarse_groups reports the highest group number "
             "even though an intermediate group is unused") {
            REQUIRE(rd_grid_have_coarse_cells(grid.get()));
            REQUIRE(rd_grid_get_num_coarse_groups(grid.get()) == 3);
        }

        THEN("rd_grid_iget_coarse_group returns a non-null object for used "
             "group numbers and null for the skipped one") {
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 0) != nullptr);
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 1) == nullptr);
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 2) != nullptr);
        }
    }
}

TEST_CASE_METHOD(Tmpdir, "Load EGRID with extreme coarse group numbers",
                 "[unittest]") {
    const int nx = 2, ny = 2, nz = 2;
    const size_t size = nx * ny * nz;

    GIVEN("An EGRID file whose CORSNUM contains a very large group number") {
        std::vector<int> corsnum(size, 0);
        corsnum[0] = std::numeric_limits<int>::max();
        corsnum[1] = std::numeric_limits<int>::max();

        auto filename = dirname / "CORSNUM_HUGE.EGRID";
        auto grid =
            load_egrid_with_coarse_groups(filename, nx, ny, nz, corsnum.data());
        REQUIRE(grid != nullptr);

        THEN("only the used group is stored") {
            REQUIRE(rd_grid_get_num_coarse_groups(grid.get()) ==
                    std::numeric_limits<int>::max());
            REQUIRE(rd_grid_iget_coarse_group(
                        grid.get(), std::numeric_limits<int>::max() - 1) !=
                    nullptr);
            REQUIRE(rd_grid_iget_coarse_group(grid.get(), 0) == nullptr);
            REQUIRE(rd_grid_iget_coarse_group(
                        grid.get(), std::numeric_limits<int>::max()) ==
                    nullptr);
        }
    }

    GIVEN("An EGRID file whose CORSNUM contains negative group numbers") {
        std::vector<int> corsnum(size, -4);

        auto filename = dirname / "CORSNUM_NEGATIVE.EGRID";
        auto grid =
            load_egrid_with_coarse_groups(filename, nx, ny, nz, corsnum.data());
        REQUIRE(grid != nullptr);

        THEN("the cells are not part of any coarse group") {
            REQUIRE(rd_grid_get_num_coarse_groups(grid.get()) == 0);
            for (size_t i = 0; i < size; i++)
                REQUIRE(!rd_grid_cell_in_coarse_group1(grid.get(), i));
        }
    }
}
