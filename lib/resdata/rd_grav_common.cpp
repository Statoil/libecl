#include <cstdlib>
#include <cmath>

#include <ert/util/util.hpp>

#include <resdata/rd_kw.hpp>
#include <resdata/rd_file.hpp>
#include <resdata/rd_region.hpp>
#include <resdata/rd_kw_magic.hpp>
#include <resdata/rd_grav_common.hpp>

#include "detail/resdata/rd_grid_cache.hpp"
/*
  This file contains code which is common to both the rd_grav
  implementation for gravity changes, and the rd_subsidence
  implementation for changes in subsidence.
*/

bool *rd_grav_common_alloc_aquifer_cell(const rd::rd_grid_cache &grid_cache,
                                        const rd_file_type *init_file) {
    const rd_kw_type *aquifer_kw = nullptr;
    if (rd_file_has_kw(init_file, AQUIFER_KW)) {
        aquifer_kw = rd_file_iget_named_kw(init_file, AQUIFER_KW, 0);
        if (rd_kw_get_size(aquifer_kw) != grid_cache.size())
            util_abort("%s: %s keyword has %d elements, but the grid has %d "
                       "active cells\n",
                       __func__, AQUIFER_KW, rd_kw_get_size(aquifer_kw),
                       static_cast<int>(grid_cache.size()));
    }

    bool *aquifer_cell =
        (bool *)util_calloc(grid_cache.size(), sizeof *aquifer_cell);

    for (int active_index = 0; active_index < grid_cache.size(); active_index++)
        aquifer_cell[active_index] = false;

    if (aquifer_kw) {
        const int *aquifer_data = rd_kw_get_int_ptr(aquifer_kw);

        for (int active_index = 0; active_index < grid_cache.size();
             active_index++) {
            if (aquifer_data[active_index] < 0)
                aquifer_cell[active_index] = true;
        }
    }

    return aquifer_cell;
}

double rd_grav_common_eval_biot_savart(const rd::rd_grid_cache &grid_cache,
                                       rd_region_type *region,
                                       const bool *aquifer,
                                       const double *weight, double utm_x,
                                       double utm_y, double depth) {
    const auto &xpos = grid_cache.xpos();
    const auto &ypos = grid_cache.ypos();
    const auto &zpos = grid_cache.zpos();
    double sum = 0;
    if (region == NULL) {
        const int size = grid_cache.size();
        int index;
        for (index = 0; index < size; index++) {
            if (!aquifer[index]) {
                double dist_x = (xpos[index] - utm_x);
                double dist_y = (ypos[index] - utm_y);
                double dist_z = (zpos[index] - depth);
                double dist =
                    sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

                /**
            For numerical precision it might be benficial to use the
            util_kahan_sum() function to do a Kahan summation.
        */
                sum += weight[index] * dist_z / (dist * dist * dist);
            }
        }
    } else {
        const int_vector_type *index_vector = rd_region_get_active_list(region);
        const int size = int_vector_size(index_vector);
        const int *index_list = int_vector_get_const_ptr(index_vector);
        int i, index;
        for (i = 0; i < size; i++) {
            index = index_list[i];
            if (!aquifer[index]) {
                double dist_x = (xpos[index] - utm_x);
                double dist_y = (ypos[index] - utm_y);
                double dist_z = (zpos[index] - depth);
                double dist =
                    sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

                sum += weight[index] * dist_z / (dist * dist * dist);
            }
        }
    }
    return sum;
}

static inline double
rd_grav_common_eval_geertsma_kernel(int index, const double *xpos,
                                    const double *ypos, const double *zpos,
                                    double utm_x, double utm_y, double depth,
                                    double poisson_ratio, double seabed) {
    double z = zpos[index];
    z -= seabed;
    double dist_x = xpos[index] - utm_x;
    double dist_y = ypos[index] - utm_y;

    double dist_z1 = z - depth;
    double dist_z2 = dist_z1 - 2 * z;

    double dist1 = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z1 * dist_z1);
    double dist2 = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z2 * dist_z2);

    double cube_dist1 = dist1 * dist1 * dist1;
    double cube_dist2 = dist2 * dist2 * dist2;

    double displacement =
        dist_z1 / cube_dist1 + (3 - 4 * poisson_ratio) * dist_z2 / cube_dist2 -
        6 * depth * (z + depth) * dist_z2 / (dist2 * dist2 * cube_dist2) +
        2 * ((3 - 4 * poisson_ratio) * (z + depth) - depth) / cube_dist2;

    return displacement;
}

double rd_grav_common_eval_geertsma(const rd::rd_grid_cache &grid_cache,
                                    rd_region_type *region, const bool *aquifer,
                                    const double *weight, double utm_x,
                                    double utm_y, double depth,
                                    double poisson_ratio, double seabed) {
    const auto &xpos = grid_cache.xpos();
    const auto &ypos = grid_cache.ypos();
    const auto &zpos = grid_cache.zpos();
    double sum = 0;
    if (region == NULL) {
        const int size = grid_cache.size();
        int index;
        for (index = 0; index < size; index++) {
            if (!aquifer[index]) {

                double displacement = rd_grav_common_eval_geertsma_kernel(
                    index, xpos.data(), ypos.data(), zpos.data(), utm_x, utm_y,
                    depth, poisson_ratio, seabed);

                /**
            For numerical precision it might be benficial to use the
            util_kahan_sum() function to do a Kahan summation.
        */
                sum += weight[index] * displacement;
            }
        }
    } else {
        const int_vector_type *index_vector = rd_region_get_active_list(region);
        const int size = int_vector_size(index_vector);
        const int *index_list = int_vector_get_const_ptr(index_vector);
        int i, index;
        for (i = 0; i < size; i++) {
            index = index_list[i];
            if (!aquifer[index]) {
                double displacement = rd_grav_common_eval_geertsma_kernel(
                    index, xpos.data(), ypos.data(), zpos.data(), utm_x, utm_y,
                    depth, poisson_ratio, seabed);
                sum += weight[index] * displacement;
            }
        }
    }
    return sum;
}
