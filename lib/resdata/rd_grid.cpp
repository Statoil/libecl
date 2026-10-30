#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>

#include <algorithm>
#include <memory>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <string>
#include <optional>
#include <array>
#include <stdexcept>
#include <limits>
#include <filesystem>

#include <fmt/format.h>

#include <ert/util/util.hpp>
#include <ert/util/double_vector.hpp>
#include <ert/util/int_vector.hpp>
#include <ert/util/vector.hpp>
#include <ert/util/stringlist.hpp>

#include <ert/geometry/geo_util.hpp>
#include <ert/geometry/geo_polygon.hpp>

#include <resdata/rd_util.hpp>
#include <resdata/rd_type.hpp>
#include <resdata/rd_kw.hpp>
#include <resdata/rd_kw_grdecl.hpp>
#include <resdata/rd_file.hpp>
#include <resdata/rd_kw_magic.hpp>
#include <resdata/rd_endian_flip.hpp>
#include <resdata/rd_coarse_cell.hpp>
#include <resdata/rd_grid.hpp>
#include <resdata/nnc_info.hpp>

namespace fs = std::filesystem;

/**
  this function implements functionality to load eclispe grid files,
  both .egrid and .grid files - in a transparent fashion.

  observe the following convention:

  global_index:  [0 , nx*ny*nz)
  active_index:  [0 , nactive)

  about indexing
  --------------

  there are three different ways to index/access a cell:

    1. by ijk
    2. by global index, [0 , nx*ny*nz)
    3. by active index, [0 , nactive)

  most of the query functions can take input in several of the
  ways. the expected arguments are indicated as the last part of the
  function name:

    rd_grid_get_pos3()  - 3:  this function expects i,j,k
    rd_grid_get_pos1()  - 1:  this function expects a global index
    rd_grid_get_pos1a() - 1a: this function expects an active index.

*/

/**
   note about lgr
   --------------

   the eclipse local grid refinement (lgr) is organised as follows:

     1. you start with a normal grid.
     2. some of the cells can be subdivided into further internal
        grids, this is the lgr.

   this is illustrated below:


    +--------------------------------------+
    |            |            |            |
    |            |            |            |
    |     x      |      x     |     x      |
    |            |            |            |
    |            |            |            |
    |------------|------------|------------|
    |     |      |      |     |            |
    |     |  x   |   x  |     |            |
    |-----x------|------x-----|     x      |
    |  x  |  x   |   x  |     |            |
    |     |      |      |     |            |
    -------------|------------|------------|
    |            |            |            |
    |            |            |            |
    |     x      |            |            |
    |            |            |            |
    |            |            |            |
    +--------------------------------------+


  the grid above shows the following:

    1. the coarse (i.e. normal) grid has 9 cells, of which 7 marked
       with 'x' are active.

    2. two of the cells have been refined into new 2x2 grids. in the
       refined cells only three and two of the refined cells are
       active.

  in a grid file the keywords for this grid will look like like this:


  .....    __
  coords     \
  corners     |
  coords      |
  corners     |
  coords      |
  corners     |      normal coord / corners kewyords
  coords      |      for the nine coarse cells. observe
  corners     |      that when reading in these cells it is
  coords      |      impossible to know that some of the
  corners     |      cells will be subdivieded in a following
  coords      |      lgr definition.
  corners     |
  coords      |
  corners     |
  coords      |
  corners     |
  coords      |
  corners  __/________________________________________________________
  lgr        \
  lgrilg      |
  dimens      |
  coords      |
  corners     |      first lgr, with some header information,
  coords      |      and then normal coords/corners keywords for
  corners     |      the four refined cells.
  coords      |
  corners     |
  coords      |
  corners  __/
  lgr        \
  lgrilg      |
  dimens      |
  coords      |      second lgr.
  corners     |
  coords      |
  corners     |
  coords      |
  corners     |
  coords      |
  corners  __/


  for egrid files it is essentially the same, except for replacing the
  keywords coords/corners with coord/zcorn/actnum. also the lgr
  headers are somewhat different.

  solution data in restart files comes in a similar way, a restart
  file with lgr can typically look like this:

  .....   __
  .....     \
  startsol   |   all restart data for the ordinary
  pressure   |   grid.
  swat       |
  sgas       |
  ....       |
  endsol  __/____________________________
  lgr       \
  ....       |
  startsol   |   restart data for
  pressure   |   the first lgr.
  sgas       |
  swat       |
  ...        |
  endsol     |
  endlgr  __/
  lgr       \
  ....       |
  startsol   |   restart data for
  pressure   |   the second lgr.
  sgas       |
  swat       |
  ...        |
  endsol     |
  endlgr  __/


  the lgr implementation in is based on the following main principles:

   1. when loading a egrid/grid file one rd_grid_type instance will
      be allocated; this grid will contain the main grid, and all the
      lgr grids.

   2. only one datatype (rd_grid_type) is used both for the main grid
      and the lgr grids.

   3. the main grid will own (memory wise) all the lgr grids, this
      even applies to nested subgrids whose parent is also a lgr.

   4. when it comes to indexing and so on there is no difference
      between lgr grid and the main grid.


      example:
      --------

      {
         rd_file_type * restart_data = rd_file_fread_alloc(restart_filename , true);                      // load some restart info to inspect
         rd_grid_type * grid         = rd_grid_alloc(grid_filename , true);                               // bootstrap rd_grid instance
         stringlist_type * lgr_names  = rd_grid_alloc_name_list( grid );                                   // get a list of all the lgr names.

         printf("grid:%s has %d a total of %d lgr's \n", grid_filename , stringlist_get_size( lgr_names ));
         for (int lgr_nr = 0; lgr_nr < stringlist_get_size( lgr_names); lgr_nr++) {
            rd_grid_type * lgr_grid  = rd_grid_get_lgr( grid , stringlist_iget( lgr_names , lgr_nr ));    // get the rd_grid instance of the lgr - by name.
            rd_kw_type   * pressure_kw;
            int nx,ny,nz,active_size;
            rd_grid_get_dims( lgr_grid , &nx , &ny , &nz , &active_size);                             // get some size info from this lgr.
            printf("lgr:%s has %d x %d x %d elements \n",stringlist_iget(lgr_names , lgr_nr ) , nx , ny , nz);

            // ok - now we want to extract the solution vector (pressure) corresponding to this lgr:
            pressure_kw = rd_file_iget_named_kw( rd_file , "pressure" , rd_grid_get_grid_nr( lgr_grid ));
                                                                                      /|\
                                                                                       |
                                                                                       |
                                                                        we query the lgr_grid instance to find which
                                                                        occurence of the solution data we should look
                                                                        up in the rd_file instance with restart data. puuhh!!

            {
               int center_index = rd_grid_get_global_index3( lgr_grid , nx/2 , ny/2 , nz/2 );          // ask the lgr_grid to get the index at the center of the lgr grid.
               printf("the pressure in the middle of %s is %g \n", stinglist_iget( lgr_names , lgr_nr ) , rd_kw_iget_as_double( pressure_kw , center_index ));
            }
         }
         rd_file_free( restart_data );
         rd_grid_free( grid );
         stringlist_free( lgr_names );
     }

*/

/*
  About coarse groups
  -------------------

  It is possible to get ECLIPSE to join several cells together as one
  coarse cell, to reduce the size of the numerical problem. In this
  implementation this is supported as follows:

    1. Each cell has an integer flag - coarse_group, which points to
       the coarse group that the current cell is part of, or the value
       COARSE_NONE for cells which are not part of a coarsening group.

    2. The details of a coarse cell is implemented in rd_coarse_cell.c.

    3. The rd_grid structure contains a list of rd_coarse_cell
       instances in the coarse_cells vector.

    4. The introduction of coarse groups makes the concept of active
       cells slightly more complex:

         a) All the cells in the coarse group count as one (or zero)
            active cells when asking the grid how many active cells
            there are.

         b) In the case that several of the un-coarsened cells in the
            coarse group are active the mapping between global index
            and active index will no longer be unique - there will be
            several different global indices mapping to the same
            active index.

  The api for coarse related tasks is briefly:

    - int rd_grid_get_num_coarse_groups( const rd_grid_type * main_grid )
    - bool rd_grid_have_coarse_cells( const rd_grid_type * main_grid )
    - rd_coarse_cell_type * rd_grid_iget_coarse_group( const rd_grid_type * rd_grid , int coarse_nr );

  In addition to the API presented by the rd_coarse_cell.c implementation.
*/

/*
  About dual porosity
  -------------------

  Eclipse has support for dual porosity systems where the reservoir is
  made from an interleaved system of matrix blocks and fractures. The
  current implementation has some support for reading such properties
  from the grid files:

    - The active property of the cells is an integer which is a sum of
      the flag values ACTIVE_MATRIX and ACTIVE_FRACTURE.

    - All functions operating on fracture properties have 'fracture'
      as part of the name. The functions operating on the matrix do
      (typically) not have matrix as part of the name. The matrix
      properties are the default properties which apply in the single
      porosity case.

    - In the EGRID files the dual porosity properties are set in the
      ACTNUM field in the file. The numerical values are [0,1,2,3] for
      inactive cells, matrix active, fracture active and
      matrix+fracture active respectively.

    - For the GRID files there is abolutely no metadata to tell that
      this is a dual porosity run (I think ...) - instead the whole
      grid is repeated one more time with cells for the fractures
      following after the matrix cells.

      Naively the GRID file of a dual porosity run will report that it
      contains 2*NZ layers. In the current implementation heuristics
      is used to detect the situation, and the grid will only be
      loaded as consisting of 'NZ' geometric layers.

    - The documentation seems to indicate that the number of active
      fracture cells can in general be different from the number of
      active matrix cells, and the current implementation takes pains
      to support that possibility - but all examples I have come over
      have nactive_fracture == nactive_matrix?

    - Properties and solution data in restart/init/grdecl files are
      not treated here. These properties will just increase (typically
      double) in size - and how to treat them will be a question of
      convention. The following shows a possible solution:

      {
         char fracture_kw[9];
         char matrix_kw[9];
         int  matrix_size   = rd_grid_get_nactive( rd_grid );
         int  fracture_size = rd_grid_get_nactive_fracture( rd_grid );

         swat = rd_file_iget_name_kw( rst_file , "SWAT" , 0);

         snsprintf(fracture_kw , 9 , "F-%6s" , rd_kw_get_header( swat ));
         snsprintf(matrix_kw   , 9 , "M-%6s" , rd_kw_get_header( swat ));

         rd_kw_type * M = rd_kw_alloc_sub_copy( swat , matrix_kw   , 0  , matrix_size );
         rd_kw_type * F = rd_kw_alloc_sub_copy( swat , fracture_kw , matrix_size  , fracture_size );
      }
*/

/*
  About nnc
  ---------

  When loading a grid file the various NNC related keywords are read
  in to assemble information of the NNC. The NNC information is
  organized as follows:

       For cells with NNC's attached the information is kept in a
       nnc_info_type structure. For a particular cell the nnc_info
       structure keeps track of which other cells this particular cell
       is connected to, on a per grid (i.e. LGR) basis.

       In the nnc_info structure the different grids are identified
       through the lgr_nr.



  Example usage:
  --------------

     rd_grid_type * grid = rd_grid_alloc("FILE.EGRID");

     // Get a int_vector instance with all the cells which have nnc info
     // attached.
     const int_vector_type * cells_with_nnc = rd_grid_get_nnc_index_list( grid );

     // Iterate over all the cells with nnc info:
     for (int i=0; i < int_vector_size( cells_with_nnc ); i++) {
         int cell_index =  int_vector_iget( cells_with_nnc , i);
         const nnc_info_type * nnc_info = rd_grid_get_nnc_info1( grid , cell_index);

         // Get all the nnc connections from @cell_index to other cells in the same grid
         {
            const int_vector_type * nnc_list = nnc_info_get_self_index_list( nnc_info );
            for (int j=0; j < int_vector_size( nnc_list ); j++)
               printf("Cell[%d] -> %d  in the same grid \n",cell_index , int_vector_iget(nnc_list , j));
         }


         {
             for (int lgr_index=0; lgr_index < nnc_info_get_size( nnc_info ); lgr_index++) {
                nnc_vector_type * nnc_vector = nnc_info_iget_vector( nnc_info , lgr_index );
                int lgr_nr = nnc_vector_get_lgr_nr( nnc_vector );
                if (lgr_nr != nnc_info_get_lgr_nr( nnc_info )) {
                   const int_vector_type * nnc_list = nnc_vector_get_index_list( nnc_vector );
                   for (int j=0; j < int_vector_size( nnc_list ); j++)
                       printf("Cell[%d] -> %d  in lgr:%d/%s \n",cell_index , int_vector_iget(nnc_list , j) , lgr_nr  , rd_grid_get_lgr_name(rd_grid , lgr_nr));
                }
             }
         }
     }


  Dual porosity and nnc: In ECLIPSE the connection between the matrix
  properties and the fracture properties in a cell is implemented as a
  nnc where the fracture cell has global index in the range [nx*ny*nz,
  2*nz*ny*nz). In ert we we have not implemented this double covering
  in the case of dual porosity models so NNC involving
  fracture cells are not considered.
*/

/*
  about tetraheder decomposition
  ------------------------------

  the table tetraheder_permutations describe how the cells can be
  divided into twelve tetrahedrons. the dimensions in the the table
  are as follows:

   1. the first dimension is the "way" the cell is divided into
      tetrahedrons, there are two different ways. for cells where the
      four point corners on a face are not in the same plane, the two
      methods will not give the same result. which one is "right"??

   2. the second dimension is the tetrahedron number, for each way of
      the two ways there are a total of twelve tetrahedrons.

   3. the third and list dimension is the point number in this
      tetrahedron. when forming a tetrahedron the first input point
      should always be the point corresponding to center of the
      cell. that is not explicit in this table.

   i.e. for instance the third tetrahedron for the first method
   consists of the cells:

        tetraheheder_permutations[0][2] = {0 , 4 , 5}

   in addition to the central point. the value [0..7] correspond the
   the number scheme of the corners in a cell used by eclipse:


       lower layer:   upper layer  (higher value of z - i.e. lower down in resrvoir).

         2---3           6---7
         |   |           |   |
         0---1           4---5


   table entries are ripped from resdatapost code - file: kvpvos.f in
   klib/
*/

/*
About ordering of the corners in the cell
-----------------------------------------

This code reads and builds the grid structure from the
GRID/EGRID files without really considering the question about where
the cells are located in "the real world", the format is quite general
and it should(?) be possible to formulate different conventions
(i.e. handedness and direction of z-axis) with the same format.

The corners in a cell are numbered 0 - 7, where corners 0-3 constitute
one layer and the corners 4-7 consitute the other layer. Observe the
numbering does not follow a consistent rotation around the face:


                                  j
  6---7                        /|\
  |   |                         |
  4---5                         |
                                |
                                o---------->  i
  2---3
  |   |
  0---1

Many grids are left-handed, i.e. the direction of increasing z will
point down towards the center of the earth. Hence in the figure above
the layer 4-7 will be deeper down in the reservoir than layer 0-3, and
also have higher z-value.

Warning: The main author of this code suspects that the coordinate
system can be right-handed as well, giving a z axis which will
increase 'towards the sky'; the safest is probaly to check this
explicitly if it matters for the case at hand.

Method 0 corresponds to a tetrahedron decomposition which will split
the lower layer along the 1-2 diagonal and the upper layer along the
4-7 diagonal, method 1 corresponds to the alternative decomposition
which splits the lower face along the 0-3 diagnoal and the upper face
along the 5-6 diagonal.
*/

static const int tetrahedron_permutations[2][12][3] = {{// K-
                                                        {0, 1, 2},
                                                        {3, 2, 1},
                                                        // J+
                                                        {6, 2, 7},
                                                        {3, 7, 2},
                                                        // I-
                                                        {0, 2, 4},
                                                        {6, 4, 2},
                                                        // I+
                                                        {3, 1, 7},
                                                        {5, 7, 1},
                                                        // J-
                                                        {0, 4, 1},
                                                        {5, 1, 4},
                                                        // K+
                                                        {5, 4, 7},
                                                        {6, 7, 4}},
                                                       {// K-
                                                        {1, 3, 0},
                                                        {2, 0, 3},
                                                        // J+
                                                        {2, 3, 6},
                                                        {7, 6, 3},
                                                        // I-
                                                        {2, 6, 0},
                                                        {4, 0, 6},
                                                        // I+
                                                        {7, 3, 5},
                                                        {1, 5, 3},
                                                        // J-
                                                        {1, 0, 5},
                                                        {4, 5, 0},
                                                        // K+
                                                        {7, 5, 6},
                                                        {4, 6, 5}}};

/*

  the implementation is based on a hierarchy of three datatypes:

   1. rd_grid   - this is the only exported datatype
   2. rd_cell   - internal
   3. point      - internal

*/

typedef struct point_struct point_type;

struct point_struct {
    double x;
    double y;
    double z;
};

static void point_mapaxes_transform(point_type *p, const double origo[2],
                                    const double unit_x[2],
                                    const double unit_y[2]) {
    double new_x = origo[0] + p->x * unit_x[0] + p->y * unit_y[0];
    double new_y = origo[1] + p->x * unit_x[1] + p->y * unit_y[1];

    p->x = new_x;
    p->y = new_y;
}

static void point_mapaxes_invtransform(point_type *p, const double origo[2],
                                       const double unit_x[2],
                                       const double unit_y[2]) {
    double norm = 1.0 / (unit_x[0] * unit_y[1] - unit_x[1] * unit_y[0]);
    double dx = p->x - origo[0];
    double dy = p->y - origo[1];

    double org_x = (dx * unit_y[1] - dy * unit_y[0]) * norm;
    double org_y = (-dx * unit_x[1] + dy * unit_x[0]) * norm;

    p->x = org_x;
    p->y = org_y;
}

static void point_inplace_add(point_type *point, const point_type *add) {
    point->x += add->x;
    point->y += add->y;
    point->z += add->z;
}

static void point_inplace_sub(point_type *point, const point_type *sub) {
    point->x -= sub->x;
    point->y -= sub->y;
    point->z -= sub->z;
}

static void point_inplace_scale(point_type *point, double scale_factor) {
    point->x *= scale_factor;
    point->y *= scale_factor;
    point->z *= scale_factor;
}

/**
   Will compute the vector cross product B x C and store the result in A.
*/

static void point_vector_cross(point_type *A, const point_type *B,
                               const point_type *C) {
    A->x = B->y * C->z - B->z * C->y;
    A->y = -(B->x * C->z - B->z * C->x);
    A->z = B->x * C->y - B->y * C->x;
}

static double point_dot_product(const point_type *v1, const point_type *v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

static void point_compare(const point_type *p1, const point_type *p2,
                          bool *equal) {
    const double tolerance = 0.001;

    double diff_x = (fabs(p1->x - p2->x) / fabs(p1->x + p2->x + 1));
    double diff_y = (fabs(p1->y - p2->y) / fabs(p1->y + p2->y + 1));
    double diff_z = (fabs(p1->z - p2->z) / fabs(p1->z + p2->z + 1));

    if (diff_x + diff_y + diff_z > tolerance)
        *equal = false;
}

static void point_dump_ascii(const point_type *p, FILE *stream,
                             const double *offset) {
    if (offset)
        fprintf(stream, "(%7.2f, %7.2f, %7.2f) ", p->x - offset[0],
                p->y - offset[1], p->z - offset[2]);
    else
        fprintf(stream, "(%7.2f, %7.2f, %7.2f) ", p->x, p->y, p->z);
}

static void point_set(point_type *p, double x, double y, double z) {
    p->x = x;
    p->y = y;
    p->z = z;
}

static void point_shift(point_type *p, double dx, double dy, double dz) {
    p->x += dx;
    p->y += dy;
    p->z += dz;
}

static void point_copy_values(point_type *p, const point_type *src) {
    point_set(p, src->x, src->y, src->z);
}

/*
  Indices used in the cell->active_index[] array.
*/
#define MATRIX_INDEX 0
#define FRACTURE_INDEX 1

#define COARSE_GROUP_NONE -1
#define HOST_CELL_NONE -1

#define CELL_FLAG_VALID                                                        \
    1 /* In the case of GRID files not necessarily all cells geometry values set - in that case this will be left as false. */
#define CELL_FLAG_CENTER                                                       \
    2 /* Has the center value been calculated - this is by default not done to speed up loading a tiny bit. */
#define CELL_FLAG_TAINTED                                                      \
    4 /* Keep invalid cells out of real-world calculations with some heuristics.*/
#define CELL_FLAG_VOLUME 8

typedef struct rd_cell_struct rd_cell_type;

#define GET_CELL_FLAG(cell, flag)                                              \
    (((cell->cell_flags & (flag)) == 0) ? false : true)
#define SET_CELL_FLAG(cell, flag) ((cell->cell_flags |= (flag)))
#define METER_TO_FEET_SCALE_FACTOR 3.28084
#define METER_TO_CM_SCALE_FACTOR 100.0

struct rd_cell_struct {
    point_type center;
    point_type corner_list[8];

    double
        volume; /* Cache volume - whether it is initialized or not is handled by a cell_flags. */
    int active;
    int active_index
        [2]; /* [0]: The active matrix index; [1]: the active fracture index */
    const rd_grid_type *
        lgr; /* if this cell is part of an lgr; this will point to a grid instance for that lgr; NULL if not part of lgr. */
    int host_cell; /* the global index of the host cell for an lgr cell, set to -1 for normal cells. */
    int coarse_group; /* The index of the coarse group holding this cell -1 for non-coarsened cells. */
    int cell_flags;
    nnc_info_ptr nnc_info{nullptr,
                          &nnc_info_free}; /* Non-neighbour connection info*/

    rd_cell_struct() = default;
    ~rd_cell_struct() = default;
    rd_cell_struct(rd_cell_struct &&) = default;
    rd_cell_struct &operator=(rd_cell_struct &&) = default;

    rd_cell_struct(const rd_cell_struct &other) { *this = other; }
    rd_cell_struct &operator=(const rd_cell_struct &other) {
        if (this == &other)
            return *this;
        center = other.center;
        std::copy(std::begin(other.corner_list), std::end(other.corner_list),
                  std::begin(corner_list));
        volume = other.volume;
        active = other.active;
        active_index[0] = other.active_index[0];
        active_index[1] = other.active_index[1];
        lgr = other.lgr;
        host_cell = other.host_cell;
        coarse_group = other.coarse_group;
        cell_flags = other.cell_flags;
        nnc_info.reset(other.nnc_info
                           ? nnc_info_alloc_copy(other.nnc_info.get())
                           : nullptr);
        return *this;
    }
};

static ert_rd_unit_enum
rd_grid_check_unit_system(const rd_kw_type *gridunit_kw);
static const float *rd_grid_get_mapaxes(const rd_grid_type *grid);

#define RD_GRID_ID 991010

struct rd_grid_struct {
    UTIL_TYPE_ID_DECLARATION;
    int lgr_nr; /* EGRID files: corresponds to item 4 in gridhead - 0 for the main grid.
                                          GRID files: 0 for the main grid, then 1 -> number of LGRs in order read from file*/
    std::string
        name; /* the name of the file for the main grid - name of the lgr for lgrs. */
    int ny, nz, nx;
    int size; /* == nx*ny*nz */
    int total_active;
    int total_active_fracture;
    std::vector<bool>
        visited; /* internal helper struct used when searching for index - empty when unused. */
    std::vector<int>
        index_map; /* this a list of nx*ny*nz elements, where value -1 means inactive cell .*/
    std::vector<int>
        inv_index_map; /* this is list of total_active elements - which point back to the index_map. */

    std::vector<int>
        fracture_index_map; /* For fractures: this a list of nx*ny*nz elements, where value -1 means inactive cell .*/
    std::vector<int>
        inv_fracture_index_map; /* For fractures: this is list of total_active elements - which point back to the index_map. */

    std::vector<rd_cell_type> cells;

    std::optional<std::string>
        parent_name; /* the name of the parent for a nested lgr - for the main grid, and also a
                                          lgr descending directly from the main grid this will be nullopt. */
    std::unordered_map<std::string, rd_grid_type *> children;
    const rd_grid_type *
        parent_grid; /* the parent grid for this (lgr) - NULL for the main grid. */
    const rd_grid_type
        *global_grid; /* the global grid - NULL for the main grid. */

    bool coarsening_active;
    std::vector<rd_coarse_cell_ptr> coarse_cells;
    /*
    the two fields below are for *storing* lgr grid instances. observe
    that these fields will be NULL for lgr grids, i.e. grids with
    lgr_nr > 0.
  */
    std::vector<rd_grid_ptr>
        LGR_list; /* a vector of rd_grid instances for LGRs. The index
                     corresponds to the order LGRs are read from file.
                     The vector is empty except for the main grid.*/
    std::vector<int>
        lgr_index_map; /* a vector that maps LGR-nr for EGRID files to index into the LGR_list.*/
    std::unordered_map<std::string, rd_grid_type *>
        LGR_hash; /* a hash of pointers to rd_grid instances - for name based lookup of lgr. */
    int parent_box
        [6]; /* integers i1,i2, j1,j2, k1,k2 of the parent grid region containing this lgr. the indices are inclusive - zero offset */
             /* not used yet .. */

    int dualp_flag;
    bool use_mapaxes;
    double unit_x[2];
    double unit_y[2];
    double origo[2];
    std::optional<std::array<float, 6>> mapaxes;
    rd_kw_ptr coord_kw{nullptr,
                       &rd_kw_free}; /* Retained for writing the grid to file.
                                        In principal it should be possible to
                                        recalculate this from the cell coordinates,
                                        but in cases with skewed cells this has proved
                                        numerically challenging. */

    ert_rd_unit_enum unit_system;
    int eclipse_version;
};

ert_rd_unit_enum rd_grid_get_unit_system(const rd_grid_type *grid) {
    return grid->unit_system;
}

static void rd_cell_compare(const rd_cell_type &c1, const rd_cell_type &c2,
                            bool include_nnc, bool *equal) {
    int i;

    if (c1.active != c2.active)
        *equal = false;

    if (c1.active_index[0] != c2.active_index[0])
        *equal = false;

    if (c1.active_index[1] != c2.active_index[1])
        *equal = false;

    if (c1.coarse_group != c2.coarse_group)
        *equal = false;

    if (c1.host_cell != c2.host_cell)
        *equal = false;

    if (*equal) {
        for (i = 0; i < 8; i++)
            point_compare(&c1.corner_list[i], &c2.corner_list[i], equal);
    }

    if (include_nnc) {
        if (*equal)
            *equal = nnc_info_equal(c1.nnc_info.get(), c2.nnc_info.get());
    }
}

static void rd_cell_assert_center(rd_cell_type &cell);

static void rd_cell_dump_ascii(rd_cell_type &cell, int i, int j, int k,
                               FILE *stream, const double *offset) {
    fprintf(stream,
            "Cell: i:%3d  j:%3d    k:%3d   host_cell:%d  CoarseGroup:%4d "
            "active_nr:%6d  active:%d \nCorners:\n",
            i, j, k, cell.host_cell, cell.coarse_group,
            cell.active_index[MATRIX_INDEX], cell.active);

    rd_cell_assert_center(cell);
    fprintf(stream, "Center   : ");
    point_dump_ascii(&cell.center, stream, offset);
    fprintf(stream, "\n");

    for (int l = 0; l < 8; l++) {
        fprintf(stream, "Corner %d : ", l);
        point_dump_ascii(&cell.corner_list[l], stream, offset);
        fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
}

static void rd_cell_fwrite_GRID(const rd_grid_type *grid,
                                const rd_cell_type &cell, bool fracture_cell,
                                int coords_size, int i, int j, int k,
                                int global_index, const rd_kw_ptr &coords_kw,
                                const rd_kw_ptr &corners_kw,
                                ERT::FortIO &fortio) {
    rd_kw_iset_int(coords_kw.get(), 0, i + 1);
    rd_kw_iset_int(coords_kw.get(), 1, j + 1);
    rd_kw_iset_int(coords_kw.get(), 2, k + 1);
    rd_kw_iset_int(coords_kw.get(), 3, global_index + 1);

    rd_kw_iset_int(coords_kw.get(), 4, 0);
    if (fracture_cell) {
        if (cell.active & CELL_ACTIVE_FRACTURE)
            rd_kw_iset_int(coords_kw.get(), 4, 1);
    } else {
        if (cell.active & CELL_ACTIVE_MATRIX)
            rd_kw_iset_int(coords_kw.get(), 4, 1);
    }

    if (coords_size == 7) {
        rd_kw_iset_int(coords_kw.get(), 5, cell.host_cell + 1);
        rd_kw_iset_int(coords_kw.get(), 6, cell.coarse_group + 1);
    }

    rd_kw_fwrite(coords_kw.get(), fortio.get());
    {
        float *corners =
            static_cast<float *>(rd_kw_get_void_ptr(corners_kw.get()));
        point_type point;

        for (int c = 0; c < 8; c++) {
            point_copy_values(&point, &cell.corner_list[c]);
            if (grid->use_mapaxes)
                point_mapaxes_invtransform(&point, grid->origo, grid->unit_x,
                                           grid->unit_y);

            corners[3 * c] = point.x;
            corners[3 * c + 1] = point.y;
            corners[3 * c + 2] = point.z;
        }
    }
    rd_kw_fwrite(corners_kw.get(), fortio.get());
}

static double max2(double x1, double x2) { return (x1 > x2) ? x1 : x2; }

static double min2(double x1, double x2) { return (x1 < x2) ? x1 : x2; }

static double min4(double x1, double x2, double x3, double x4) {
    return min2(min2(x1, x2), min2(x3, x4));
}

static double max4(double x1, double x2, double x3, double x4) {
    return max2(max2(x1, x2), max2(x3, x4));
}

static double max8(double x1, double x2, double x3, double x4, double x5,
                   double x6, double x7, double x8) {
    return max2(max4(x1, x2, x3, x4), max4(x5, x6, x7, x8));
}

static double min8(double x1, double x2, double x3, double x4, double x5,
                   double x6, double x7, double x8) {
    return min2(min4(x1, x2, x3, x4), min4(x5, x6, x7, x8));
}

static double rd_cell_min_z(const rd_cell_type &cell) {
    return min4(cell.corner_list[0].z, cell.corner_list[1].z,
                cell.corner_list[2].z, cell.corner_list[3].z);
}

static double rd_cell_max_z(const rd_cell_type &cell) {
    return max4(cell.corner_list[4].z, cell.corner_list[5].z,
                cell.corner_list[6].z, cell.corner_list[7].z);
}

/**
   the grid can be rotated so that it is not safe to consider only one
   plane for the x/y min/max.
*/

static double rd_cell_min_x(const rd_cell_type &cell) {
    return min8(cell.corner_list[0].x, cell.corner_list[1].x,
                cell.corner_list[2].x, cell.corner_list[3].x,
                cell.corner_list[4].x, cell.corner_list[5].x,
                cell.corner_list[6].x, cell.corner_list[7].x);
}

static double rd_cell_max_x(const rd_cell_type &cell) {
    return max8(cell.corner_list[0].x, cell.corner_list[1].x,
                cell.corner_list[2].x, cell.corner_list[3].x,
                cell.corner_list[4].x, cell.corner_list[5].x,
                cell.corner_list[6].x, cell.corner_list[7].x);
}

static double rd_cell_min_y(const rd_cell_type &cell) {
    return min8(cell.corner_list[0].y, cell.corner_list[1].y,
                cell.corner_list[2].y, cell.corner_list[3].y,
                cell.corner_list[4].y, cell.corner_list[5].y,
                cell.corner_list[6].y, cell.corner_list[7].y);
}

static double rd_cell_max_y(const rd_cell_type &cell) {
    return max8(cell.corner_list[0].y, cell.corner_list[1].y,
                cell.corner_list[2].y, cell.corner_list[3].y,
                cell.corner_list[4].y, cell.corner_list[5].y,
                cell.corner_list[6].y, cell.corner_list[7].y);
}

/**
   The problem is that some grids have invalid cells. Typically
   the cells accomodating numerical aquifers are located at an utm
   position (0,0).

   Cells which have some pillars located in (0,0) and some cells
   located among the rest of the grid become completely warped - with
   insane volumes, parts of the reservoir volume doubly covered, and
   so on.

   To keep these cells out of the real-world (i.e. involving utm
   coordinates) computations they are marked as 'tainted' in this
   function. The tainting procedure is completely heuristic, and
   probably wrong.

   --------------

   There is second heuristic which marks cells as invalid. In some
   cases (which ??) cells outside the area of interest are just set to
   have all four corner at the same arbitrary depth; these cells are
   inactive and do not affect flow simulations - however the arbitrary
   location of the cells warps visualisation of normal inactive cells
   completely. We therefore try to invalidate such cells here. The
   algorithm used is the same as used for RMS; however RMS will mark
   the cells as inactive - whereas we mark already inactive cells as
   invalid.

   It is very important that this mechanism is NOT used to set cells
   to inactive; as that will break completely when it comes to
   e.g. data like PRESSURE and SWAT which are defined as a number of
   active elements.
 */

static void rd_cell_taint_cell(rd_cell_type *cell) {
    for (auto &p : cell->corner_list) {
        if ((p.x == 0) && (p.y == 0)) {
            SET_CELL_FLAG(cell, CELL_FLAG_TAINTED);
            break;
        }
    }

    /*
    Second heuristic to invalidate cells.
  */
    if (cell->active == CELL_NOT_ACTIVE) {
        if (!GET_CELL_FLAG(cell, CELL_FLAG_TAINTED)) {
            const point_type p0 = cell->corner_list[0];
            int cell_index = 1;
            while (true) {
                const point_type pi = cell->corner_list[cell_index];
                if (pi.z != p0.z)
                    // There is a difference - the cell is certainly valid.
                    break;
                else {
                    cell_index++;
                    if (cell_index == 8) {
                        // They have all been at the same height up until now;
                        // the cell is marked as invalid.
                        SET_CELL_FLAG(cell, CELL_FLAG_TAINTED);
                        break;
                    }
                }
            }
        }
    }
}

static int rd_cell_get_twist(const rd_cell_type &cell) {
    int twist_count = 0;

    for (int c = 0; c < 4; c++) {
        const point_type *p1 = &cell.corner_list[c];
        const point_type *p2 = &cell.corner_list[c + 4];
        if ((p2->z - p1->z) < 0)
            twist_count += 1;
    }
    return twist_count;
}

/**
   Observe that when allocating based on a grid file not all cells are
   necessarily accessed beyond this function. In general not all cells
   will have a coords/corners section in the grid file.
*/

static void rd_cell_init(rd_cell_type *cell, bool init_valid) {
    cell->active = CELL_NOT_ACTIVE;
    cell->lgr = NULL;
    cell->host_cell = HOST_CELL_NONE;
    cell->coarse_group = COARSE_GROUP_NONE;
    cell->cell_flags = 0;
    cell->active_index[MATRIX_INDEX] = -1;
    cell->active_index[FRACTURE_INDEX] = -1;
    if (init_valid)
        cell->cell_flags = CELL_FLAG_VALID;

    cell->nnc_info.reset();
}

/*
#define mod(x,n) ((x) % (n))
static int rd_cell_get_tetrahedron_method( int i , int j , int k) {
  return (1 + (1 - 2*mod(i,2)) * (1 - 2 * mod(j,2)) * (1 - 2*mod(k,2))) / 2;
}
#undef mod
*/

static void rd_cell_set_center(rd_cell_type &cell) {
    point_set(&cell.center, 0, 0, 0);
    for (auto &c : cell.corner_list)
        point_inplace_add(&cell.center, &c);
    point_inplace_scale(&cell.center, 1.0 / 8.0);

    SET_CELL_FLAG((&cell), CELL_FLAG_CENTER);
}

static void rd_cell_assert_center(rd_cell_type &cell) {
    if (!GET_CELL_FLAG((&cell), CELL_FLAG_CENTER))
        rd_cell_set_center(cell);
}

static double C(double *r, int f1, int f2, int f3) {
    if (f1 == 0) {
        if (f2 == 0) {
            if (f3 == 0)
                return r[0]; // 000
            else
                return r[4] - r[0]; // 001
        } else {
            if (f3 == 0)
                return r[2] - r[0]; // 010
            else
                return r[6] + r[0] - r[4] - r[2]; // 011
        }
    } else {
        if (f2 == 0) {
            if (f3 == 0)
                return r[1] - r[0]; // 100
            else
                return r[5] + r[0] - r[4] - r[1]; // 101
        } else {
            if (f3 == 0)
                return r[3] + r[0] - r[2] - r[1]; // 110
            else
                return r[7] + r[4] + r[2] + r[1] - r[6] - r[5] - r[3] -
                       r[0]; // 111
        }
    }
}

static double rd_cell_get_volume(const rd_cell_type &cell) {
    double volume = 0;
    double X[8];
    double Y[8];
    double Z[8];
    for (int c = 0; c < 8; c++) {
        X[c] = cell.corner_list[c].x;
        Y[c] = cell.corner_list[c].y;
        Z[c] = cell.corner_list[c].z;
    }

    for (int pb = 0; pb <= 1; pb++)
        for (int pg = 0; pg <= 1; pg++)
            for (int qa = 0; qa <= 1; qa++)
                for (int qg = 0; qg <= 1; qg++)
                    for (int ra = 0; ra <= 1; ra++)
                        for (int rb = 0; rb <= 1; rb++) {
                            int divisor =
                                (qa + ra + 1) * (pb + rb + 1) * (pg + qg + 1);
                            double dV = C(X, 1, pb, pg) * C(Y, qa, 1, qg) *
                                            C(Z, ra, rb, 1) -
                                        C(X, 1, pb, pg) * C(Z, qa, 1, qg) *
                                            C(Y, ra, rb, 1) -
                                        C(Y, 1, pb, pg) * C(X, qa, 1, qg) *
                                            C(Z, ra, rb, 1) +
                                        C(Y, 1, pb, pg) * C(Z, qa, 1, qg) *
                                            C(X, ra, rb, 1) +
                                        C(Z, 1, pb, pg) * C(X, qa, 1, qg) *
                                            C(Y, ra, rb, 1) -
                                        C(Z, 1, pb, pg) * C(Y, qa, 1, qg) *
                                            C(X, ra, rb, 1);

                            volume += dV / divisor;
                        }

    return fabs(volume);
}

typedef struct tetrahedron_struct tetrahedron_type;

struct tetrahedron_struct {
    point_type p0;
    point_type p1;
    point_type p2;
    point_type p3;
};

/*
  Calculates the volume of a tetrahedron, a normalization of 1/6 is
  omitted.
*/
static inline double tetrahedron_volume6(tetrahedron_type tet) {
    point_type bxc;

    /* vector subtraction */
    tet.p0.x -= tet.p3.x;
    tet.p0.y -= tet.p3.y;
    tet.p0.z -= tet.p3.z;

    tet.p1.x -= tet.p3.x;
    tet.p1.y -= tet.p3.y;
    tet.p1.z -= tet.p3.z;

    tet.p2.x -= tet.p3.x;
    tet.p2.y -= tet.p3.y;
    tet.p2.z -= tet.p3.z;

    /* cross product */
    bxc.x = tet.p1.y * tet.p2.z - tet.p1.z * tet.p2.y;
    bxc.y = -(tet.p1.x * tet.p2.z - tet.p1.z * tet.p2.x);
    bxc.z = tet.p1.x * tet.p2.y - tet.p1.y * tet.p2.x;

    /* dot product */
    return tet.p0.x * bxc.x + tet.p0.y * bxc.y + tet.p0.z * bxc.z;
}

/*
 Returns true if and only if the point p is inside the tetrahedron tet.
*/
static bool tetrahedron_contains(tetrahedron_type tet, const point_type p) {
    const double epsilon = 1e-9;
    double tetra_volume = fabs(tetrahedron_volume6(tet));

    if (tetra_volume < epsilon)
        return false;

    // Decomposes tetrahedron into 4 new tetrahedrons
    point_type tetra_points[4] = {tet.p0, tet.p1, tet.p2, tet.p3};
    double decomposition_volume = 0;
    for (int i = 0; i < 4; ++i) {
        const point_type tmp = tetra_points[i];
        tetra_points[i] = p;

        // Compute volum of decomposition tetrahedron
        tetrahedron_type dec_tet;
        dec_tet.p0 = tetra_points[0];
        dec_tet.p1 = tetra_points[1];
        dec_tet.p2 = tetra_points[2];
        dec_tet.p3 = tetra_points[3];
        decomposition_volume += fabs(tetrahedron_volume6(dec_tet));

        tetra_points[i] = tmp;
    }

    return (fabs(tetra_volume - decomposition_volume) < epsilon);
}

static double parallelogram_area3d(const point_type *p0, const point_type *p1,
                                   const point_type *p2) {
    point_type a = *p1;
    point_type b = *p2;
    point_inplace_sub(&a, p0);
    point_inplace_sub(&b, p0);

    point_type c;
    point_vector_cross(&c, &a, &b);
    return sqrt(point_dot_product(&c, &c));
}

/*
 * Returns true if and only if point p is contained in the triangle denoted by
 * p0, p1 and p2. Note that if the triangle is a line, the function will still
 * return true if the point lies on the line segment.
 */
static bool triangle_contains3d(const point_type *p0, const point_type *p1,
                                const point_type *p2, const point_type *p) {
    double epsilon = 1e-10;
    double vt = parallelogram_area3d(p0, p1, p2);

    double v1 = parallelogram_area3d(p0, p1, p);
    double v2 = parallelogram_area3d(p0, p2, p);
    double v3 = parallelogram_area3d(p1, p2, p);

    // p0, p1, p2 represents a line segment and
    // p lies on the line this segment represents
    if (vt < epsilon && fabs(v1 + v2 + v3) < epsilon) {
        double x_min = util_double_min(p0->x, util_double_min(p1->x, p2->x));
        double x_max = util_double_max(p0->x, util_double_max(p1->x, p2->x));

        double y_min = util_double_min(p0->y, util_double_min(p1->y, p2->y));
        double y_max = util_double_max(p0->y, util_double_max(p1->y, p2->y));

        double z_min = util_double_min(p0->z, util_double_min(p1->z, p2->z));
        double z_max = util_double_max(p0->z, util_double_max(p1->z, p2->z));

        return (x_min - epsilon <= p->x && p->x <= x_max + epsilon &&
                y_min - epsilon <= p->y && p->y <= y_max + epsilon &&
                z_min - epsilon <= p->z && p->z <= z_max + epsilon);
    }

    return (fabs(vt - (v1 + v2 + v3)) < epsilon);
}

/**
       lower layer:   upper layer

         2---3           6---7
         |   |           |   |
         0---1           4---5
*/
static void rd_cell_init_regular(rd_cell_type &cell, const double *offset,
                                 int i, int j, int k, int global_index,
                                 const double *ivec, const double *jvec,
                                 const double *kvec, const int *actnum) {
    point_set(&cell.corner_list[0], offset[0], offset[1],
              offset[2]); // Point 0

    cell.corner_list[1] = cell.corner_list[0]; // Point 1
    point_shift(&cell.corner_list[1], ivec[0], ivec[1], ivec[2]);

    cell.corner_list[2] = cell.corner_list[0]; // Point 2
    point_shift(&cell.corner_list[2], jvec[0], jvec[1], jvec[2]);

    cell.corner_list[3] = cell.corner_list[1]; // Point 3
    point_shift(&cell.corner_list[3], jvec[0], jvec[1], jvec[2]);

    for (int i = 0; i < 4; i++) {
        cell.corner_list[i + 4] = cell.corner_list[i]; // Point 4-7
        point_shift(&cell.corner_list[i + 4], kvec[0], kvec[1], kvec[2]);
    }

    if (actnum != NULL)
        cell.active = actnum[global_index];
    else
        cell.active = CELL_ACTIVE;
}

UTIL_SAFE_CAST_FUNCTION(rd_grid, RD_GRID_ID);
UTIL_IS_INSTANCE_FUNCTION(rd_grid, RD_GRID_ID);

/**
   this function uses heuristics in an attempt to mark cells with
   wrong geometry - see further comments in the function
   rd_cell_taint_cell() which actually does it.
*/

static void rd_grid_taint_cells(rd_grid_type *rd_grid) {
    for (auto &cell : rd_grid->cells) {
        rd_cell_taint_cell(&cell);
    }
}

static bool rd_grid_alloc_cells(rd_grid_type *grid, bool init_valid) {
    try {
        grid->cells.resize(grid->size);
    } catch (const std::bad_alloc &) {
        return false;
    }

    for (auto &cell : grid->cells)
        rd_cell_init(&cell, init_valid);

    return true;
}

/**
   will create a new blank grid instance. if the global_grid argument
   is != NULL the newly created grid instance will copy the mapaxes
   transformations; and set the global_grid pointer of the new grid
   instance. apart from that no further lgr-relationsip initialisation
   is performed.
*/

static rd_grid_type *rd_grid_alloc_empty(rd_grid_type *global_grid,
                                         ert_rd_unit_enum unit_system,
                                         int dualp_flag, int nx, int ny, int nz,
                                         int lgr_nr, bool init_valid) {
    /* Check for int overflow in a roundabout way */
    long size = static_cast<long>(nx) * static_cast<long>(ny);
    if (nx < 0 || ny < 0 || nz < 0 || size > INT_MAX)
        return NULL;
    size *= static_cast<long>(nz);
    if (size > INT_MAX)
        return NULL;

    rd_grid_type *grid = new rd_grid_type();
    UTIL_TYPE_ID_INIT(grid, RD_GRID_ID);
    grid->total_active = 0;
    grid->total_active_fracture = 0;
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->size = nx * ny * nz;
    grid->lgr_nr = lgr_nr;
    grid->global_grid = global_grid;
    grid->coarsening_active = false;

    grid->dualp_flag = dualp_flag;
    grid->unit_system = unit_system;

    if (global_grid != NULL) {
        /*
       this is an lgr instance, and we inherit the global grid
       transformations from the main grid. The underlying mapaxes
       pointer is jest left at NULL for lgr instances.
    */
        grid->unit_x[0] = global_grid->unit_x[0];
        grid->unit_x[1] = global_grid->unit_x[1];
        grid->unit_y[0] = global_grid->unit_y[0];
        grid->unit_y[1] = global_grid->unit_y[1];
        grid->origo[0] = global_grid->origo[0];
        grid->origo[1] = global_grid->origo[1];
        grid->use_mapaxes = global_grid->use_mapaxes;
    } else {
        grid->unit_x[0] = 1;
        grid->unit_x[1] = 0;
        grid->unit_y[0] = 0;
        grid->unit_y[1] = 1;
        grid->origo[0] = 0;
        grid->origo[1] = 0;
        grid->use_mapaxes = false;
    }

    grid->name = "";
    grid->parent_grid = NULL;
    grid->eclipse_version = 0;

    /* This is the large allocation - which can potentially fail. */
    if (!rd_grid_alloc_cells(grid, init_valid)) {
        rd_grid_free(grid);
        grid = NULL;
    }
    return grid;
}

static int rd_grid_get_global_index__(const rd_grid_type *rd_grid, int i, int j,
                                      int k) {
    return i + j * rd_grid->nx + k * rd_grid->nx * rd_grid->ny;
}

static void rd_grid_set_cell_EGRID(rd_grid_type *rd_grid, int i, int j, int k,
                                   double x[4][2], double y[4][2],
                                   double z[4][2], const int *actnum,
                                   const int *corsnum) {

    const int global_index = rd_grid_get_global_index__(rd_grid, i, j, k);
    rd_cell_type &cell = rd_grid->cells.at(global_index);

    for (int iz = 0; iz < 2; iz++) {
        for (int ip = 0; ip < 4; ip++) {
            int c = ip + iz * 4;
            point_set(&cell.corner_list[c], x[ip][iz], y[ip][iz], z[ip][iz]);

            if (rd_grid->use_mapaxes)
                point_mapaxes_transform(&cell.corner_list[c], rd_grid->origo,
                                        rd_grid->unit_x, rd_grid->unit_y);
        }
    }

    /*
    If actnum == NULL that is taken to mean active.

    for normal runs actnum will be 1 for active cells,
    for dual porosity models it can also be 2 and 3.
  */
    if (actnum == NULL)
        cell.active = CELL_ACTIVE;
    else
        cell.active = actnum[global_index];

    if (corsnum != NULL)
        cell.coarse_group = corsnum[global_index] - 1;
}

static void rd_grid_set_cell_GRID(rd_grid_type *rd_grid, int coords_size,
                                  const int *coords, const float *corners) {

    const int i = coords[0] - 1; /* eclipse 1 offset */
    const int j = coords[1] - 1;
    int k = coords[2] - 1;
    int global_index;
    bool matrix_cell = true;
    int active_value = CELL_ACTIVE_MATRIX;

    /*
    This is the rather hysterical treatment of dual porosity in GRID
    files. Cells with k >= nz consitute the fracture part of the
    grid. For these cell the cell properties are not recalculated, but
    the active flag is updated to include the active|inactive
    properties of the fracture.
  */

    if (k >= rd_grid->nz) {
        k -= rd_grid->nz;
        matrix_cell = false;
        active_value = CELL_ACTIVE_FRACTURE;
    }

    global_index = rd_grid_get_global_index__(rd_grid, i, j, k);
    rd_cell_type &cell = rd_grid->cells.at(global_index);

    /* the coords keyword can optionally contain 4,5 or 7 elements:

        coords[0..2] = i,j,k
        coords[3]    = global_cell number (not used here)
        ----
        coords[4]    = 1,0 for active/inactive cells. Does NOT differentiate between matrix and fracture cells.
        coords[5]    = 0 for normal cells, icell of host cell for lgr cell.
        coords[6]    = 0 for normal cells, coarsening group for coarsened cell [not treated yet].

        if coords[4] is not present it is assumed that the cell is active.

        Note about LGRs: The GRID file format has an additional
        keyword called LGRILG which maps out the parent box containing
        an LGR - this only duplicates the information which can be
        learned from the coords[5] element, and is not used in the
        current code.
  */

    {
        int c;

        /*
       The RD_GRID_MAINGRID_LGR_NR != rd_grid->lgr_nr test checks if
       this is a LGR; if this test applies we either have bug - or a
       GRID file with LGRs and only 4/5 elements in the coords keywords.
       In the latter case we must start using the LGRILG keyword.
    */
        if ((RD_GRID_MAINGRID_LGR_NR != rd_grid->lgr_nr) && (coords_size != 7))
            throw std::invalid_argument(
                "Need 7 element coords keywords for LGR - or reimplement to "
                "use LGRILG keyword.");

        switch (coords_size) {
        case 4: /* all cells active */
            cell.active += active_value;
            break;
        case 5: /* only spesific cells active - no lgr */
            cell.active += coords[4] * active_value;
            break;
        case 7:
            cell.active += coords[4] * active_value;
            cell.host_cell = coords[5] - 1;
            cell.coarse_group = coords[6] - 1;
            if (cell.coarse_group >= 0)
                rd_grid->coarsening_active = true;
            break;
        default:
            throw std::invalid_argument(
                fmt::format("coord size:{} unrecognized - should be 4, 5 or 7.",
                            coords_size));
        }

        if (matrix_cell) {
            for (c = 0; c < 8; c++) {
                point_set(&cell.corner_list[c], corners[3 * c],
                          corners[3 * c + 1], corners[3 * c + 2]);

                if (rd_grid->use_mapaxes)
                    point_mapaxes_transform(&cell.corner_list[c],
                                            rd_grid->origo, rd_grid->unit_x,
                                            rd_grid->unit_y);
            }
        }
    }
    SET_CELL_FLAG((&cell), CELL_FLAG_VALID);
}

/**
   The function rd_grid_set_active_index() must be called immediately
   prior to calling this function, to ensure that
   rd_grid->total_active is correct.
*/

static void rd_grid_init_index_map__(rd_grid_type *rd_grid,
                                     std::vector<int> &index_map,
                                     std::vector<int> &inv_index_map,
                                     int active_mask, int type_index) {
    for (int global_index = 0; global_index < rd_grid->size; global_index++) {
        const rd_cell_type &cell = rd_grid->cells.at(global_index);
        if (cell.active & active_mask) {
            index_map[global_index] = cell.active_index[type_index];

            if (cell.coarse_group == COARSE_GROUP_NONE)
                inv_index_map[cell.active_index[type_index]] = global_index;
            //else: In the case of coarse groups the inv_index_map is set below.
        } else
            index_map[global_index] = -1;
    }
}

static void rd_grid_realloc_index_map(rd_grid_type *rd_grid) {
    /* Creating the inverse mapping for the matrix cells. */
    rd_grid->index_map.resize(rd_grid->size);
    rd_grid->inv_index_map.resize(rd_grid->total_active);
    rd_grid_init_index_map__(rd_grid, rd_grid->index_map,
                             rd_grid->inv_index_map, CELL_ACTIVE_MATRIX,
                             MATRIX_INDEX);

    /* Create the inverse mapping for the fractures. */
    if (rd_grid->dualp_flag != FILEHEAD_SINGLE_POROSITY) {
        rd_grid->fracture_index_map.resize(rd_grid->size);
        rd_grid->inv_fracture_index_map.resize(rd_grid->total_active_fracture);
        rd_grid_init_index_map__(rd_grid, rd_grid->fracture_index_map,
                                 rd_grid->inv_fracture_index_map,
                                 CELL_ACTIVE_FRACTURE, FRACTURE_INDEX);
    }

    /* Update the inverse map in the case of coarse cells. Observe that
     in the case of coarse cells with more than one active cell in the
     main grid, the inverse active -> global mapping will map to the
     first active cell in the coarse cell.
  */
    int size = rd_grid_get_num_coarse_groups(rd_grid);
    for (int coarse_group = 0; coarse_group < size; coarse_group++) {
        rd_coarse_cell_type *coarse_cell =
            rd_grid_iget_coarse_group(rd_grid, coarse_group);

        // Coarse cell group numbering is not necessarily in consecutive order
        if (coarse_cell == NULL) {
            continue;
        }

        if (rd_coarse_cell_get_num_active(coarse_cell) > 0) {
            int global_index =
                rd_coarse_cell_iget_active_cell_index(coarse_cell, 0);
            int active_value = rd_coarse_cell_iget_active_value(coarse_cell, 0);
            int active_index = rd_coarse_cell_get_active_index(coarse_cell);
            int active_fracture_index =
                rd_coarse_cell_get_active_fracture_index(coarse_cell);

            if (active_value &
                CELL_ACTIVE_MATRIX) // active->global mapping point to one "random" cell in the coarse group
                rd_grid->inv_index_map[active_index] = global_index;

            if (active_value & CELL_ACTIVE_FRACTURE)
                rd_grid->inv_fracture_index_map[active_fracture_index] =
                    global_index;

            int coarse_size = rd_coarse_cell_get_size(coarse_cell);
            const int_vector_type *global_index_list =
                rd_coarse_cell_get_index_vector(coarse_cell);
            for (int ic = 0; ic < coarse_size; ic++) {
                int gi = int_vector_iget(global_index_list, ic);

                if (active_value &
                    CELL_ACTIVE_MATRIX) // All cells in the coarse group point to the same active index.
                    rd_grid->index_map[gi] = active_index;

                if (active_value & CELL_ACTIVE_FRACTURE)
                    rd_grid->fracture_index_map[gi] = active_fracture_index;
            }
        } // else the coarse cell does not have any active cells.
    }
}

/*
  This function goes through the entire grid and sets the active_index
  of all the cells. The functione rd_grid_realloc_index_map()
  subsequently reads this to create and initialize the index map.
*/

static void rd_grid_set_active_index(rd_grid_type *rd_grid) {
    int active_index = 0;
    int active_fracture_index = 0;

    if (!rd_grid_have_coarse_cells(rd_grid)) {
        /* Keeping a fast path for the 99% most common case of no coarse
       groups and single porosity. */
        for (rd_cell_type &cell : rd_grid->cells) {
            if (cell.active & CELL_ACTIVE_MATRIX) {
                cell.active_index[MATRIX_INDEX] = active_index;
                active_index++;
            }
        }

        if (rd_grid->dualp_flag != FILEHEAD_SINGLE_POROSITY) {
            for (rd_cell_type &cell : rd_grid->cells) {
                if (cell.active & CELL_ACTIVE_FRACTURE) {
                    cell.active_index[FRACTURE_INDEX] = active_fracture_index;
                    active_fracture_index++;
                }
            }
        }
    } else {
        /* --- More involved path in the case of coarsening groups. --- */

        /* 1: Reset the coarse cell active_indices.
          In the involved path rd_coarse_cell_update_index() only updates
          the coarse cells' active_index if it is -1. */
        for (const auto &coarse_cell : rd_grid->coarse_cells) {
            // Coarse cell group numbering is not necessarily in consecutive order
            if (coarse_cell) {
                rd_coarse_cell_reset_active_index(coarse_cell.get());
            }
        }

        /* 2: Go through all the cells and set the active index. In the
          case of coarse cells we only set the common active index of
          the entire coarse cell.
    */
        for (int global_index = 0; global_index < rd_grid->size;
             global_index++) {
            rd_cell_type &cell = rd_grid->cells.at(global_index);
            if (cell.active != CELL_NOT_ACTIVE) {
                if (cell.coarse_group == COARSE_GROUP_NONE) {

                    if (cell.active & CELL_ACTIVE_MATRIX) {
                        cell.active_index[MATRIX_INDEX] = active_index;
                        active_index++;
                    }

                    if (cell.active & CELL_ACTIVE_FRACTURE) {
                        cell.active_index[FRACTURE_INDEX] =
                            active_fracture_index;
                        active_fracture_index++;
                    }

                } else {
                    rd_coarse_cell_type *coarse_cell =
                        rd_grid_iget_coarse_group(rd_grid, cell.coarse_group);
                    rd_coarse_cell_update_index(
                        coarse_cell, global_index, &active_index,
                        &active_fracture_index, cell.active);
                }
            }
        }

        /*
      3: Go through all the coarse cells and set the active index and
         active value of all the cells in the coarse cell to the
         common value for the coarse cell.
    */

        int size = rd_grid_get_num_coarse_groups(rd_grid);
        for (int coarse_group = 0; coarse_group < size; coarse_group++) {
            rd_coarse_cell_type *coarse_cell =
                rd_grid_iget_coarse_group(rd_grid, coarse_group);

            // Coarse cell group numbering is not necessarily in consecutive order
            if (coarse_cell == NULL) {
                continue;
            }

            if (rd_coarse_cell_get_num_active(coarse_cell) > 0) {
                int cell_active_index =
                    rd_coarse_cell_get_active_index(coarse_cell);
                int cell_active_value =
                    rd_coarse_cell_iget_active_value(coarse_cell, 0);
                int group_size = rd_coarse_cell_get_size(coarse_cell);
                const int *coarse_cell_list =
                    rd_coarse_cell_get_index_ptr(coarse_cell);

                for (int i = 0; i < group_size; i++) {
                    int global_index = coarse_cell_list[i];

                    rd_cell_type &cell = rd_grid->cells.at(global_index);

                    if (cell_active_value & CELL_ACTIVE_MATRIX)
                        cell.active_index[MATRIX_INDEX] = cell_active_index;

                    /* Coarse cell and dual porosity - that is probably close to zero measure. */
                    if (cell_active_value & CELL_ACTIVE_FRACTURE) {
                        int cell_active_fracture_index =
                            rd_coarse_cell_get_active_fracture_index(
                                coarse_cell);
                        cell.active_index[FRACTURE_INDEX] =
                            cell_active_fracture_index;
                    }
                }
            }
        }
    }
    rd_grid->total_active = active_index;
    rd_grid->total_active_fracture = active_fracture_index;
}

static void rd_grid_update_index(rd_grid_type *rd_grid) {
    rd_grid_set_active_index(rd_grid);
    rd_grid_realloc_index_map(rd_grid);
}

static rd_coarse_cell_type *
rd_grid_get_or_create_coarse_cell(rd_grid_type *rd_grid, int coarse_nr) {
    while (static_cast<int>(rd_grid->coarse_cells.size()) <= coarse_nr)
        rd_grid->coarse_cells.emplace_back(nullptr, &rd_coarse_cell_free);
    if (!rd_grid->coarse_cells[coarse_nr])
        rd_grid->coarse_cells[coarse_nr].reset(rd_coarse_cell_alloc());

    return rd_grid->coarse_cells[coarse_nr].get();
}

static void rd_grid_init_coarse_cells(rd_grid_type *rd_grid) {
    if (rd_grid->coarsening_active) {
        for (int global_index = 0; global_index < rd_grid->size;
             global_index++) {
            const rd_cell_type &cell = rd_grid->cells.at(global_index);
            if (cell.coarse_group != COARSE_GROUP_NONE) {
                rd_coarse_cell_type *coarse_cell =
                    rd_grid_get_or_create_coarse_cell(rd_grid,
                                                      cell.coarse_group);
                int i, j, k;
                rd_grid_get_ijk1(rd_grid, global_index, &i, &j, &k);
                rd_coarse_cell_update(coarse_cell, i, j, k, global_index);
            }
        }
    }
}

rd_coarse_cell_type *rd_grid_iget_coarse_group(const rd_grid_type *rd_grid,
                                               int coarse_nr) {
    return rd_grid->coarse_cells[coarse_nr].get();
}

bool rd_grid_cell_in_coarse_group1(const rd_grid_type *main_grid,
                                   int global_index) {
    const rd_cell_type &cell = main_grid->cells.at(global_index);
    return cell.coarse_group != COARSE_GROUP_NONE;
}

static void rd_grid_pillar_cross_planes(const point_type *p0, double e_x,
                                        double e_y, double e_z, const double *z,
                                        double *x, double *y) {
    int k;
    if (e_z != 0) {
        for (k = 0; k < 2; k++) {
            double t = (z[k] - p0->z) / e_z;
            x[k] = p0->x + t * e_x;
            y[k] = p0->y + t * e_y;
        }
    } else {
        for (k = 0; k < 2; k++) {
            x[k] = p0->x;
            y[k] = p0->y;
        }
    }
}

/**
   This function must be run before the cell coordinates are
   calculated.  This function is only called for the main grid
   instance, and not for lgrs; the lgrs will inherit the mapaxes
   settings from the global grid.
 */

static void rd_grid_init_mapaxes(rd_grid_type *rd_grid, bool apply_mapaxes,
                                 const float *mapaxes) {
    if (rd_grid->global_grid != NULL)
        throw std::logic_error(
            "GridError: Trying to grid transformation data from mapaxes for a "
            "subgrid(lgr)");
    {
        const double unit_y[2] = {mapaxes[0] - mapaxes[2],
                                  mapaxes[1] - mapaxes[3]};
        const double unit_x[2] = {mapaxes[4] - mapaxes[2],
                                  mapaxes[5] - mapaxes[3]};
        {
            double norm_x =
                1.0 / sqrt(unit_x[0] * unit_x[0] + unit_x[1] * unit_x[1]);
            double norm_y =
                1.0 / sqrt(unit_y[0] * unit_y[0] + unit_y[1] * unit_y[1]);

            rd_grid->unit_x[0] = unit_x[0] * norm_x;
            rd_grid->unit_x[1] = unit_x[1] * norm_x;
            rd_grid->unit_y[0] = unit_y[0] * norm_y;
            rd_grid->unit_y[1] = unit_y[1] * norm_y;
        }
        rd_grid->origo[0] = mapaxes[2];
        rd_grid->origo[1] = mapaxes[3];
        std::array<float, 6> mapaxes_copy;
        memcpy(mapaxes_copy.data(), mapaxes, 6 * sizeof(float));
        rd_grid->mapaxes = mapaxes_copy;

        /*
      If the apply_mapaxes variable is false we will not apply the
      transformation; however we will still internalize the mapaxes
      information, and also output that if the file is saved again.
    */
        rd_grid->use_mapaxes = apply_mapaxes;
    }
}

/**
   this function will add a rd_grid instance as a lgr to the main
   grid. the lgr grid as added to two different structures of the main
   grid:

    1. in the main_grid->LGR_list the lgr instances are inserted in
       order of occurence in the grid file. the following equalities
       should apply:

          occurence number in file == lgr_grid->grid_nr

       This 'mostly' agrees with the GRIDHEAD(4) item, but
       unfortunately not always. Cases have popped up where the series
       of GRIDHEAD(4) values from lgr to lgr have holes :-(

       The LGR_list owns the grid via a unique_ptr with rd_grid_free
       as its deleter, so the grid is destroyed when the vector is
       destroyed.

    2. in the main->lgr_hash the lgr instance is installed with the
       lgrname as key. only a reference is installed in the hash
       table.

    observe that this is in principle somewhat different from the
    install functions below; here the lgr is added to the top level
    grid (i.e. the main grid) which has the storage responsability of
    all the lgr instances. the cell->lgr relationship is established
    in the _install_egrid / install_grid functions further down.
*/

static rd_grid_type *rd_grid_add_lgr(rd_grid_type *main_grid,
                                     rd_grid_ptr lgr_grid) {
    int lgr_nr = lgr_grid->lgr_nr;
    if (lgr_nr < 0)
        throw std::invalid_argument(
            fmt::format("LGR number must be non-negative, got {}", lgr_nr));

    // Bound check for LGR number to avoid excessive memory allocation in
    // lgr_index_map. Value was choses arbitrarily, but should be large enough
    // for any reasonable case.
    constexpr int MAX_LGR_NR = 1'000'000;
    if (lgr_nr > MAX_LGR_NR)
        throw std::invalid_argument(
            fmt::format("LGR number {} exceeds maximum supported value {}",
                        lgr_nr, MAX_LGR_NR));

    auto ptr = lgr_grid.get();
    std::string lgr_name = lgr_grid->name;
    main_grid->LGR_list.push_back(std::move(lgr_grid));
    if (lgr_nr >= static_cast<int>(main_grid->lgr_index_map.size()))
        main_grid->lgr_index_map.resize(lgr_nr + 1, 0);
    main_grid->lgr_index_map[lgr_nr] = main_grid->LGR_list.size() - 1;
    main_grid->LGR_hash[lgr_name] = ptr;
    return ptr;
}

static void rd_grid_install_lgr_common(rd_grid_type *host_grid,
                                       rd_grid_type *lgr_grid) {
    host_grid->children[lgr_grid->name] = lgr_grid;
    lgr_grid->parent_grid = host_grid;
}

/**
   this function will set the lgr pointer of the relevant cells in the
   host grid to point to the lgr_grid. observe that the rd_cell_type
   instances do *not* own the lgr_grid - all lgr_grid instances are
   owned by the main grid.
*/

static void rd_grid_install_lgr_EGRID(rd_grid_type *host_grid,
                                      rd_grid_type *lgr_grid,
                                      const int *hostnum) {
    for (int i = 0; i < lgr_grid->size; i++) {
        int host_index = hostnum[i] - 1;
        rd_cell_type &lgr_cell = lgr_grid->cells.at(i);
        rd_cell_type &host_cell = host_grid->cells.at(host_index);

        host_cell.lgr = lgr_grid;
        lgr_cell.host_cell = host_index;
    }
    rd_grid_install_lgr_common(host_grid, lgr_grid);
}

/**
   similar to rd_grid_install_lgr_egrid for grid based instances.
*/
static void rd_grid_install_lgr_GRID(rd_grid_type *host_grid,
                                     rd_grid_type *lgr_grid) {
    for (const rd_cell_type &lgr_cell : lgr_grid->cells) {
        rd_cell_type &host_cell = host_grid->cells.at(lgr_cell.host_cell);
        host_cell.lgr = lgr_grid;
    }
    rd_grid_install_lgr_common(host_grid, lgr_grid);
}

/**
   sets the name of the lgr and the name of the parent, if this is a
   nested lgr. for normal lgr descending directly from the coarse grid
   the parent_name is set to nullopt.
*/
static void rd_grid_set_lgr_name_EGRID(rd_grid_type *lgr_grid,
                                       const rd_file_type *rd_file,
                                       int grid_nr) {
    rd_kw_type *lgrname_kw =
        rd_file_iget_named_kw(rd_file, LGR_KW, grid_nr - 1);
    lgr_grid->name = rd_kw_iget_stripped_string(lgrname_kw, 0);
    if (rd_file_has_kw(rd_file, LGR_PARENT_KW)) {
        rd_kw_type *parent_kw =
            rd_file_iget_named_kw(rd_file, LGR_PARENT_KW, grid_nr - 1);
        std::string parent = rd_kw_iget_stripped_string(parent_kw, 0);

        if (!parent.empty())
            lgr_grid->parent_name = parent;
        else
            lgr_grid->parent_name = std::nullopt;
    }
}

/**
   sets the name of the lgr and the name of the parent, if this is a
   nested lgr. for lgr descending directly from the parent eclipse
   will supply 'global' (whereas for egrid it will return '' -
   cool?). anyway global -> NULL.
*/

static void rd_grid_set_lgr_name_GRID(rd_grid_type *lgr_grid,
                                      const rd_file_type *rd_file,
                                      int grid_nr) {
    rd_kw_type *lgr_kw = rd_file_iget_named_kw(rd_file, LGR_KW, grid_nr - 1);
    lgr_grid->name = rd_kw_iget_stripped_string(lgr_kw, 0);
    {
        /**
       the lgr keyword can have one or two elements; in the case of two elements
       the second element will be the name of the parent grid - in the case of
       only one element the current lgr is assumed to descend from the main grid
    */
        if (rd_kw_get_size(lgr_kw) == 2) {
            std::string parent = rd_kw_iget_stripped_string(lgr_kw, 1);

            if (!parent.empty() && parent != GLOBAL_STRING)
                lgr_grid->parent_name = parent;
        }
    }
}

/*
  This function will calculate the index of a corner in the the zcorn
  array.  The corner is given with i,j,k indices of the cell, and
  corner number [0...8). The function is mainly a convenience function
  for testing, and should not be used in a tight inner loop unrolling
  the zcorn vector.
*/

int rd_grid_zcorn_index__(int nx, int ny, int i, int j, int k, int c) {
    int zcorn_index = k * 8 * nx * ny + j * 4 * nx + 2 * i;
    if ((c % 2) == 1)
        zcorn_index += 1;

    if ((c % 4) == 2)
        zcorn_index += 2 * nx;

    if ((c % 4) == 3)
        zcorn_index += 2 * nx;

    if (c >= 4)
        zcorn_index += 4 * nx * ny;

    return zcorn_index;
}

static void rd_grid_init_GRDECL_data_jslice(rd_grid_type *rd_grid,
                                            const float *zcorn,
                                            const float *coord,
                                            const int *actnum,
                                            const int *corsnum, int j) {
    const int nx = rd_grid->nx;
    const int ny = rd_grid->ny;
    const int nz = rd_grid->nz;

    for (int i = 0; i < nx; i++) {
        point_type pillars[4][2];
        int pillar_index[4];
        pillar_index[0] = 6 * (j * (nx + 1) + i);
        pillar_index[1] = 6 * (j * (nx + 1) + i + 1);
        pillar_index[2] = 6 * ((j + 1) * (nx + 1) + i);
        pillar_index[3] = 6 * ((j + 1) * (nx + 1) + i + 1);

        for (int ip = 0; ip < 4; ip++) {
            int index = pillar_index[ip];
            point_set(&pillars[ip][0], coord[index], coord[index + 1],
                      coord[index + 2]);

            index += 3;
            point_set(&pillars[ip][1], coord[index], coord[index + 1],
                      coord[index + 2]);
        }

        {
            double ex[4];
            double ey[4];
            double ez[4];

            for (int ip = 0; ip < 4; ip++) {
                ex[ip] = pillars[ip][1].x - pillars[ip][0].x;
                ey[ip] = pillars[ip][1].y - pillars[ip][0].y;
                ez[ip] = pillars[ip][1].z - pillars[ip][0].z;
            }

            for (int k = 0; k < nz; k++) {
                double x[4][2];
                double y[4][2];
                double z[4][2];

                for (int c = 0; c < 2; c++) {
                    z[0][c] = zcorn[k * 8 * nx * ny + j * 4 * nx + 2 * i +
                                    c * 4 * nx * ny];
                    z[1][c] = zcorn[k * 8 * nx * ny + j * 4 * nx + 2 * i + 1 +
                                    c * 4 * nx * ny];
                    z[2][c] = zcorn[k * 8 * nx * ny + j * 4 * nx + 2 * nx +
                                    2 * i + c * 4 * nx * ny];
                    z[3][c] = zcorn[k * 8 * nx * ny + j * 4 * nx + 2 * nx +
                                    2 * i + 1 + c * 4 * nx * ny];
                }

                for (int ip = 0; ip < 4; ip++)
                    rd_grid_pillar_cross_planes(&pillars[ip][0], ex[ip], ey[ip],
                                                ez[ip], z[ip], x[ip], y[ip]);

                rd_grid_set_cell_EGRID(rd_grid, i, j, k, x, y, z, actnum,
                                       corsnum);
            }
        }
    }
}

/*
  2---3
  |   |
  0---1
*/

static void rd_grid_copy_mapaxes(rd_grid_type *target_grid,
                                 const rd_grid_type *src_grid) {
    target_grid->use_mapaxes = src_grid->use_mapaxes;
    target_grid->mapaxes = src_grid->mapaxes;

    for (int i = 0; i < 2; i++) {
        target_grid->unit_x[i] = src_grid->unit_x[i];
        target_grid->unit_y[i] = src_grid->unit_y[i];
        target_grid->origo[i] = src_grid->origo[i];
    }
}

static void rd_grid_copy_content(rd_grid_type *target_grid,
                                 const rd_grid_type *src_grid) {
    for (int i = 0; i < src_grid->size; i++) {
        rd_cell_type *target_cell = &target_grid->cells.at(i);
        const rd_cell_type *src_cell = &src_grid->cells.at(i);

        *target_cell = *src_cell;
    }
    rd_grid_copy_mapaxes(target_grid, src_grid);

    target_grid->parent_name = src_grid->parent_name;
    target_grid->name = src_grid->name;

    target_grid->coarsening_active = src_grid->coarsening_active;
    rd_grid_init_coarse_cells(target_grid);
}

static rd_grid_ptr rd_grid_alloc_copy__(const rd_grid_type *src_grid,
                                        rd_grid_type *main_grid) {
    auto copy_grid =
        rd_grid_ptr(rd_grid_alloc_empty(
                        main_grid, src_grid->unit_system, src_grid->dualp_flag,
                        rd_grid_get_nx(src_grid), rd_grid_get_ny(src_grid),
                        rd_grid_get_nz(src_grid), 0, false),
                    &rd_grid_free);
    if (copy_grid) {
        rd_grid_copy_content(
            copy_grid.get(),
            src_grid); // This will handle everything except LGR relationships which is established in the calling routine
        rd_grid_update_index(copy_grid.get());
    }
    return copy_grid;
}

rd_grid_type *rd_grid_alloc_copy(const rd_grid_type *src_grid) {
    auto copy_grid = rd_grid_alloc_copy__(src_grid, NULL);

    for (const auto &src_lgr : src_grid->LGR_list) {
        // This transfers the storage ownership of the LGR to copy_grid
        auto copy_lgr = rd_grid_alloc_copy__(src_lgr.get(), copy_grid.get());
        auto *lgr = rd_grid_add_lgr(copy_grid.get(), std::move(copy_lgr));

        rd_grid_type *host_grid;
        if (!lgr->parent_name)
            host_grid = copy_grid.get();
        else
            host_grid =
                rd_grid_get_lgr(copy_grid.get(), lgr->parent_name->c_str());

        for (auto &lgr_cell : lgr->cells) {
            rd_cell_type &host_cell = host_grid->cells.at(lgr_cell.host_cell);
            host_cell.lgr = lgr;
        }
        rd_grid_install_lgr_common(host_grid, lgr);
    }

    return copy_grid.release();
}

static const float *
rd_grid_get_mapaxes_from_kw__(const rd_kw_type *mapaxes_kw) {
    if (rd_kw_get_size(mapaxes_kw) != 6)
        throw std::invalid_argument(
            fmt::format("Invalid size of MAPAXES keyword = {}, expected 6",
                        rd_kw_get_size(mapaxes_kw)));
    const float *mapaxes_data = rd_kw_get_float_ptr(mapaxes_kw);

    float x1 = mapaxes_data[2];
    float y1 = mapaxes_data[3];
    float x2 = mapaxes_data[4];
    float y2 = mapaxes_data[5];

    float norm_denominator = x1 * y2 - x2 * y1;

    if (norm_denominator == 0.0) {
        mapaxes_data = NULL;
    }

    return mapaxes_data;
}

static rd_grid_ptr rd_grid_alloc_GRDECL_kw__(
    rd_grid_type *global_grid, int dualp_flag, bool apply_mapaxes,
    const rd_kw_type *gridhead_kw, const rd_kw_type *zcorn_kw,
    const rd_kw_type *coord_kw, const rd_kw_type *gridunit_kw, /* Can be NULL */
    const rd_kw_type *mapaxes_kw,                              /* Can be NULL */
    const rd_kw_type *corsnum_kw,                              /* Can be NULL */
    const int *actnum) {                                       /* Can be NULL */
    int gtype, nx, ny, nz, lgr_nr;
    ert_rd_unit_enum unit_system = RD_METRIC_UNITS;
    gtype = rd_kw_iget_int(gridhead_kw, GRIDHEAD_TYPE_INDEX);
    nx = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NX_INDEX);
    ny = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NY_INDEX);
    nz = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NZ_INDEX);
    lgr_nr = rd_kw_iget_int(gridhead_kw, GRIDHEAD_LGR_INDEX);

    if (nx <= 0 || ny <= 0 || nz <= 0)
        throw std::invalid_argument(fmt::format(
            "Invalid grid dimensions: nx={}, ny={}, nz={}", nx, ny, nz));

    if (lgr_nr < 0)
        throw std::invalid_argument(
            fmt::format("Invalid lgr_nr: lgr_nr={}", lgr_nr));

    const int64_t nx64 = nx;
    const int64_t ny64 = ny;
    const int64_t nz64 = nz;
    const int64_t expected_coord_size = int64_t{6} * (nx64 + 1) * (ny64 + 1);
    const int64_t expected_zcorn_size = int64_t{8} * nx64 * ny64 * nz64;

    if (expected_coord_size > std::numeric_limits<int>::max() ||
        expected_zcorn_size > std::numeric_limits<int>::max())
        throw std::invalid_argument(fmt::format(
            "Grid dimensions too large: nx={}, ny={}, nz={}", nx, ny, nz));

    const int coord_size = static_cast<int>(expected_coord_size);
    const int zcorn_size = static_cast<int>(expected_zcorn_size);

    if (rd_kw_get_size(coord_kw) != coord_size)
        throw std::invalid_argument(
            fmt::format("Invalid size of COORD keyword = {}, expected 6 * "
                        "(nx + 1) * (ny + 1) = {}",
                        rd_kw_get_size(coord_kw), coord_size));

    if (rd_kw_get_size(zcorn_kw) != zcorn_size)
        throw std::invalid_argument(
            fmt::format("Invalid size of ZCORN keyword = {}, expected 8 * "
                        "nx * ny * nz = {}",
                        rd_kw_get_size(zcorn_kw), zcorn_size));

    /*
    The code used to have this test:

       if (grid_nr != rd_kw_iget_int( gridhead_kw , GRIDHEAD_LGR_INDEX))
          util_abort("%s: internal error in grid loader - lgr index mismatch\n",__func__);

    But then suddenly a EGRID file where this did not apply appeared :-(
  */

    if (gtype != GRIDHEAD_GRIDTYPE_CORNERPOINT)
        throw std::invalid_argument(fmt::format(
            "gtype:{} fatal error when loading grid - must have corner point "
            "grid",
            gtype));

    const float *mapaxes = nullptr;
    if (mapaxes_kw)
        mapaxes = rd_grid_get_mapaxes_from_kw__(mapaxes_kw);

    const int *corsnum = nullptr;
    if (corsnum_kw) {
        const int64_t expected_corsnum_size = nx64 * ny64 * nz64;
        if (expected_corsnum_size > std::numeric_limits<int>::max() ||
            rd_kw_get_size(corsnum_kw) != expected_corsnum_size)
            throw std::invalid_argument(
                fmt::format("Invalid size of CORSNUM keyword = {}, expected "
                            "nx * ny * nz = {}",
                            rd_kw_get_size(corsnum_kw), expected_corsnum_size));
        corsnum = rd_kw_get_int_ptr(corsnum_kw);
    }

    if (gridunit_kw)
        unit_system = rd_grid_check_unit_system(gridunit_kw);

    float *zcorn = rd_kw_get_float_ptr(zcorn_kw);
    float *coord = rd_kw_get_float_ptr(coord_kw);
    auto rd_grid =
        rd_grid_ptr(rd_grid_alloc_empty(global_grid, unit_system, dualp_flag,
                                        nx, ny, nz, lgr_nr, true),
                    &rd_grid_free);
    if (rd_grid) {
        if (mapaxes != NULL)
            rd_grid_init_mapaxes(rd_grid.get(), apply_mapaxes, mapaxes);

        if (corsnum != NULL)
            rd_grid->coarsening_active = true;
        rd_grid->coord_kw.reset(
            rd_kw_alloc_new("COORD", coord_size, RD_FLOAT, coord));
        int j;
#pragma omp parallel for
        for (j = 0; j < ny; j++)
            rd_grid_init_GRDECL_data_jslice(rd_grid.get(), zcorn, coord, actnum,
                                            corsnum, j);

        rd_grid_init_coarse_cells(rd_grid.get());
        rd_grid_update_index(rd_grid.get());
        rd_grid_taint_cells(rd_grid.get());
    }
    return rd_grid;
}

static rd_kw_type *rd_grid_alloc_gridhead_kw(int nx, int ny, int nz,
                                             int grid_nr) {
    rd_kw_type *gridhead_kw = rd_kw_alloc(GRIDHEAD_KW, GRIDHEAD_SIZE, RD_INT);
    rd_kw_scalar_set_int(gridhead_kw, 0);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_TYPE_INDEX,
                   GRIDHEAD_GRIDTYPE_CORNERPOINT);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_NX_INDEX, nx);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_NY_INDEX, ny);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_NZ_INDEX, nz);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_NUMRES_INDEX, 1);
    rd_kw_iset_int(gridhead_kw, GRIDHEAD_LGR_INDEX, grid_nr);
    return gridhead_kw;
}

/**
   If you create/load rd_kw instances for the various fields, this
   function can be used to create a GRID instance, without going
   through a GRID/EGRID file. Does not support LGR or coarsening
   hierarchies.
*/

rd_grid_type *
rd_grid_alloc_GRDECL_kw(int nx, int ny, int nz, const rd_kw_type *zcorn_kw,
                        const rd_kw_type *coord_kw,
                        const rd_kw_type *actnum_kw,    /* Can be NULL */
                        const rd_kw_type *mapaxes_kw) { /* Can be NULL */

    const int *actnum_data = NULL;
    if (actnum_kw) {
        const int64_t expected_actnum_size = int64_t{nx} * ny * nz;

        if (expected_actnum_size > std::numeric_limits<int>::max())
            throw std::invalid_argument(fmt::format(
                "Grid dimensions too large: nx={}, ny={}, nz={}", nx, ny, nz));

        const int actnum_size = static_cast<int>(expected_actnum_size);

        if (rd_kw_get_size(actnum_kw) != actnum_size)
            throw std::invalid_argument(
                fmt::format("Invalid size of ACTNUM keyword = {}, expected nz "
                            "* nx * ny = {}",
                            rd_kw_get_size(actnum_kw), actnum_size));
        actnum_data = rd_kw_get_int_ptr(actnum_kw);
    }

    bool apply_mapaxes = true;
    auto gridhead_kw =
        rd_kw_ptr(rd_grid_alloc_gridhead_kw(nx, ny, nz, 0), &rd_kw_free);
    rd_kw_type *gridunit_kw = NULL;
    return rd_grid_alloc_GRDECL_kw__(
               NULL, FILEHEAD_SINGLE_POROSITY, apply_mapaxes, gridhead_kw.get(),
               zcorn_kw, coord_kw, gridunit_kw, mapaxes_kw, NULL, actnum_data)
        .release();
}

static void rd_grid_init_cell_nnc_info(rd_grid_type *rd_grid,
                                       int global_index) {
    rd_cell_type &grid_cell = rd_grid->cells.at(global_index);

    if (!grid_cell.nnc_info)
        grid_cell.nnc_info.reset(nnc_info_alloc(rd_grid->lgr_nr));
}

/*
  The function rd_grid_add_self_nnc() will add a NNC connection
  between two cells in the same grid. Observe that there are two
  peculiarities with this implementation:

   1. In the rd_grid structure the nnc information is distributed
      among the cells. The main purpose of adding the nnc information
      like this is to include the NNC information in the EGRID files
      when writing to disk. Before being written to disk the NNC
      information is serialized into vectors NNC1 and NNC2. It is the
      ordering in the NNC1 and NNC2 vectors which must be correct, and
      that is goverened by the nnc_index argument - i.e. the nnc_index
      serves as an 'ID' for the NNC connections.

      After all NCC information has been entered you should be certain
      that all nnc_index values in the range [0,num_nnc] have been
      set.


   2. To get valid NNC information to load in e.g. Resinsight the
      corresponding TRANNNC keyword must be added to the INIT file,
      i.e. the calling scope must create a rd_kw with
      transmissibility values in parallell with adding NNC information
      to the grid:

         fortio_type * init_file = fortio_open_writer( "CASE.INIT" , ...
         rd_grid_type * grid ...
         rd_kw_type * trannnc_kw = rd_kw_alloc( "TRANNNC" , num_nnc , RD_FLOAT_TYPE );

         for (int i = 0; i < num_nnc; i++) {
             int   g1 = ...
             int   g2 = ..
             float T  = ..

             rd_grid_add_self_nnc( grid , g1 , g2 , i );
             rd_kw_iset( trannnc_kw , i , T );
         }
         ...
         rd_grid_fwrite_EGRID( grid , ... );
         rd_kw_fwrite( trannnc_kw , init_file );

*/

void rd_grid_add_self_nnc(rd_grid_type *grid, int cell_index1, int cell_index2,
                          int nnc_index) {
    const rd_cell_type &grid_cell = grid->cells.at(cell_index1);
    rd_grid_init_cell_nnc_info(grid, cell_index1);
    nnc_info_add_nnc(grid_cell.nnc_info.get(), grid->lgr_nr, cell_index2,
                     nnc_index);
}

/*
  This function populates nnc_info for cells with non neighbour
  connections. For cells C1 and C2 the function will only add the directed link:

      C1 -> C2

  i.e. it is impossible to go from C2 and back to C1. The functions
  will add links as:

    NNC1 -> NNC2   For nnc between cells in the same grid.
    NNCG -> NNCL   For global -> lgr connection
    NNA1 -> NNA2   For links between different LGRs
*/

static void rd_grid_init_nnc_cells(rd_grid_type *grid1, rd_grid_type *grid2,
                                   const rd_kw_type *keyword1,
                                   const rd_kw_type *keyword2) {

    int *grid1_nnc_cells = rd_kw_get_int_ptr(keyword1);
    int *grid2_nnc_cells = rd_kw_get_int_ptr(keyword2);
    int nnc_count = rd_kw_get_size(keyword2);

    for (int nnc_index = 0; nnc_index < nnc_count; nnc_index++) {
        int grid1_cell_index = grid1_nnc_cells[nnc_index] - 1;
        int grid2_cell_index = grid2_nnc_cells[nnc_index] - 1;

        /*
    In the EGRID file, grids with dual porosity are (to some
    extent ...) modeled as two independent grids stacked on top of
    eachother, where the fracture cells have global index in the range
    [nx*ny*nz, 2*nx*ny*nz).

    The physical connection between the matrix and the fractures in
    cell nr c is modelled as an nnc: cell[c] -> cell[c + nx*ny*nz]. In
    the ert ecl library we only have cells in the range [0,nx*ny*nz),
    and fracture is a property of a cell. We therefore do not include
    nnc connections involving fracture cells (i.e. cell_index >=
    nx*ny*nz).
  */
        if ((FILEHEAD_SINGLE_POROSITY != grid1->dualp_flag) &&
            ((grid1_cell_index >= grid1->size) ||
             (grid2_cell_index >= grid2->size)))
            break;

        const rd_cell_type &grid1_cell = grid1->cells.at(grid1_cell_index);
        rd_grid_init_cell_nnc_info(grid1, grid1_cell_index);
        nnc_info_add_nnc(grid1_cell.nnc_info.get(), grid2->lgr_nr,
                         grid2_cell_index, nnc_index);
    }
}

/*
  This function reads the non-neighbour connection data from file and initializes the grid structure with the the nnc data
*/
static void rd_grid_init_nnc(rd_grid_type *main_grid, rd_file_type *rd_file) {
    int num_nnchead_kw = rd_file_get_num_named_kw(rd_file, NNCHEAD_KW);

    /*
    NB: There is a bug in Eclipse version 2015.1, for MPI runs with
        six or more processors (I think ...) the NNC datastructures
        are in an internally inconsistent state; and will lead to a
        hard crash. The issue has been fixed in version 2015.2, but
        unfortunately it is not possible to test for micro version.

        if(num_nnchead_kw > 0 && main_grid->eclipse_version == 2015)
           return;
  */

    for (int i = 0; i < num_nnchead_kw; i++) {
        auto lgr_view = rd_file_view_ptr(
            rd_file_alloc_global_blockview(rd_file, NNCHEAD_KW, i),
            &rd_file_view_free);
        rd_kw_type *nnchead_kw =
            rd_file_view_iget_named_kw(lgr_view.get(), NNCHEAD_KW, 0);
        int lgr_nr = rd_kw_iget_int(nnchead_kw, NNCHEAD_LGR_INDEX);

        if (rd_file_view_has_kw(lgr_view.get(), NNC1_KW)) {
            const rd_kw_type *nnc1 =
                rd_file_view_iget_named_kw(lgr_view.get(), NNC1_KW, 0);
            const rd_kw_type *nnc2 =
                rd_file_view_iget_named_kw(lgr_view.get(), NNC2_KW, 0);

            {
                rd_grid_type *grid =
                    (lgr_nr > 0)
                        ? rd_grid_get_lgr_from_lgr_nr(main_grid, lgr_nr)
                        : main_grid;
                rd_grid_init_nnc_cells(grid, grid, nnc1, nnc2);
            }
        }

        if (rd_file_view_has_kw(lgr_view.get(), NNCL_KW)) {
            const rd_kw_type *nncl =
                rd_file_view_iget_named_kw(lgr_view.get(), NNCL_KW, 0);
            const rd_kw_type *nncg =
                rd_file_view_iget_named_kw(lgr_view.get(), NNCG_KW, 0);
            {
                rd_grid_type *grid =
                    (lgr_nr > 0)
                        ? rd_grid_get_lgr_from_lgr_nr(main_grid, lgr_nr)
                        : main_grid;
                rd_grid_init_nnc_cells(main_grid, grid, nncg, nncl);
            }
        }
    }
}

/*
  This function reads the non-neighbour connection data for
  amalgamated LGRs (that is, non-neighbour connections between two
  LGRs) and initializes the grid structure with the nnc data.
*/
static void rd_grid_init_nnc_amalgamated(rd_grid_type *main_grid,
                                         rd_file_type *rd_file) {
    int num_nncheada_kw = rd_file_get_num_named_kw(rd_file, NNCHEADA_KW);

    for (int i = 0; i < num_nncheada_kw; i++) {
        rd_kw_type *nncheada_kw =
            rd_file_iget_named_kw(rd_file, NNCHEADA_KW, i);
        int lgr_nr1 = rd_kw_iget_int(nncheada_kw, NNCHEADA_ILOC1_INDEX);
        int lgr_nr2 = rd_kw_iget_int(nncheada_kw, NNCHEADA_ILOC2_INDEX);

        rd_grid_type *lgr_grid1 =
            rd_grid_get_lgr_from_lgr_nr(main_grid, lgr_nr1);
        rd_grid_type *lgr_grid2 =
            rd_grid_get_lgr_from_lgr_nr(main_grid, lgr_nr2);

        rd_kw_type *nna1_kw = rd_file_iget_named_kw(rd_file, NNA1_KW, i);
        rd_kw_type *nna2_kw = rd_file_iget_named_kw(rd_file, NNA2_KW, i);

        rd_grid_init_nnc_cells(lgr_grid1, lgr_grid2, nna1_kw, nna2_kw);
    }
}

/**
   Creating a grid based on a EGRID file is a three step process:

   1. Load the file and extracte the keywords.
   2. Call xx_alloc_GRDECL_kw__() to build grid based on keywords.
   3. Call xx_alloc_GRDECL_data__() to build the grid based on keyword data.

   The point is that external scope can create grid based on both a
   list of keywords, and actual data - in addition to the normal
   loading of a full file. Currently the 'shortcuts'
   rd_grid_alloc_grdecl_kw() and rd_grid_alloc_grdecl_data() do not
   support LGRs.
*/

static rd_grid_ptr rd_grid_alloc_EGRID__(rd_grid_type *main_grid,
                                         const rd_file_type *rd_file,
                                         int grid_nr, bool apply_mapaxes,
                                         const int *ext_actnum) {
    rd_kw_type *gridhead_kw =
        rd_file_iget_named_kw(rd_file, GRIDHEAD_KW, grid_nr);
    rd_kw_type *zcorn_kw = rd_file_iget_named_kw(rd_file, ZCORN_KW, grid_nr);
    rd_kw_type *coord_kw = rd_file_iget_named_kw(rd_file, COORD_KW, grid_nr);
    rd_kw_type *corsnum_kw = NULL;
    rd_kw_type *actnum_kw = NULL;
    rd_kw_type *gridunit_kw = NULL;
    rd_kw_type *mapaxes_kw = NULL;
    int dualp_flag;
    int eclipse_version;
    if (grid_nr == 0) {
        rd_kw_type *filehead_kw =
            rd_file_iget_named_kw(rd_file, FILEHEAD_KW, grid_nr);
        dualp_flag = rd_kw_iget_int(filehead_kw, FILEHEAD_DUALP_INDEX);
        eclipse_version = rd_kw_iget_int(filehead_kw, FILEHEAD_YEAR_INDEX);
    } else {
        dualp_flag = main_grid->dualp_flag;
        eclipse_version = main_grid->eclipse_version;
    }

    /*
    If ext_actnum is present that is used as ACTNUM, and the file is not checked
    for an ACTNUM keyword at all.
  */
    const int *actnum_data = NULL;
    if (ext_actnum)
        actnum_data = ext_actnum;
    else {
        if (rd_file_get_num_named_kw(rd_file, ACTNUM_KW) > grid_nr) {
            actnum_kw = rd_file_iget_named_kw(rd_file, ACTNUM_KW, grid_nr);
            const int nx = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NX_INDEX);
            const int ny = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NY_INDEX);
            const int nz = rd_kw_iget_int(gridhead_kw, GRIDHEAD_NZ_INDEX);
            if (nx >= 0 && ny >= 0 && nz >= 0) {
                const int64_t actnum_size = int64_t{nx} * ny * nz;
                if (rd_kw_get_size(actnum_kw) != actnum_size)
                    throw std::invalid_argument(fmt::format(
                        "Invalid size of ACTNUM keyword = {}, expected nz "
                        "* nx * ny = {}",
                        rd_kw_get_size(actnum_kw), actnum_size));
            }
            actnum_data = rd_kw_get_int_ptr(actnum_kw);
        }
    }

    if (grid_nr == 0) {
        /* MAPAXES and COARSENING only apply to the global grid. */
        if (rd_file_has_kw(rd_file, MAPAXES_KW))
            mapaxes_kw = rd_file_iget_named_kw(rd_file, MAPAXES_KW, 0);

        if (rd_file_has_kw(rd_file, CORSNUM_KW))
            corsnum_kw = rd_file_iget_named_kw(rd_file, CORSNUM_KW, 0);

        if (rd_file_has_kw(rd_file, GRIDUNIT_KW))
            gridunit_kw = rd_file_iget_named_kw(rd_file, GRIDUNIT_KW, 0);
    }

    {
        auto rd_grid = rd_grid_alloc_GRDECL_kw__(
            main_grid, dualp_flag, apply_mapaxes, gridhead_kw, zcorn_kw,
            coord_kw, gridunit_kw, mapaxes_kw, corsnum_kw, actnum_data);

        if (RD_GRID_MAINGRID_LGR_NR != grid_nr)
            rd_grid_set_lgr_name_EGRID(rd_grid.get(), rd_file, grid_nr);
        rd_grid->eclipse_version = eclipse_version;
        return rd_grid;
    }
}

static rd_grid_ptr rd_grid_alloc_EGRID_all_grids(const char *grid_file,
                                                 bool apply_mapaxes,
                                                 const int *ext_actnum) {
    rd_file_enum file_type;
    file_type = rd_get_file_type(grid_file, NULL, NULL);
    if (file_type != RD_EGRID_FILE)
        throw std::invalid_argument(fmt::format(
            "{}: wrong file type - expected .EGRID file", grid_file));
    {
        auto rd_file = rd_file_ptr(rd_file_open(grid_file, 0), &rd_file_close);
        if (rd_file) {
            int num_grid = rd_file_get_num_named_kw(rd_file.get(), GRIDHEAD_KW);
            auto main_grid = rd_grid_alloc_EGRID__(nullptr, rd_file.get(), 0,
                                                   apply_mapaxes, ext_actnum);
            int grid_nr;

            for (grid_nr = 1; grid_nr < num_grid; grid_nr++) {
                // The apply_mapaxes argument is ignored for LGR -
                //   it inherits from parent anyway.
                auto lgr_grid = rd_grid_alloc_EGRID__(
                    main_grid.get(), rd_file.get(), grid_nr, false, NULL);
                auto *lgr =
                    rd_grid_add_lgr(main_grid.get(), std::move(lgr_grid));
                {
                    rd_grid_type *host_grid;
                    rd_kw_type *hostnum_kw = rd_file_iget_named_kw(
                        rd_file.get(), HOSTNUM_KW, grid_nr - 1);
                    if (rd_kw_get_size(hostnum_kw) != lgr->size)
                        throw std::invalid_argument(fmt::format(
                            "Invalid size of HOSTNUM keyword = {}, expected "
                            "nx * ny * nz of the LGR = {}",
                            rd_kw_get_size(hostnum_kw), lgr->size));
                    if (!lgr->parent_name)
                        host_grid = main_grid.get();
                    else
                        host_grid = rd_grid_get_lgr(main_grid.get(),
                                                    lgr->parent_name->c_str());

                    rd_grid_install_lgr_EGRID(host_grid, lgr,
                                              rd_kw_get_int_ptr(hostnum_kw));
                }
            }
            main_grid->name = grid_file;
            rd_grid_init_nnc(main_grid.get(), rd_file.get());
            rd_grid_init_nnc_amalgamated(main_grid.get(), rd_file.get());
            return main_grid;
        } else
            throw std::runtime_error(
                fmt::format("Could not open file {}", grid_file));
    }
}

static rd_grid_ptr rd_grid_alloc_EGRID(const char *grid_file,
                                       bool apply_mapaxes) {
    return rd_grid_alloc_EGRID_all_grids(grid_file, apply_mapaxes, NULL);
}

static rd_grid_ptr rd_grid_alloc_GRID_data__(
    rd_grid_type *global_grid, int num_coords, ert_rd_unit_enum unit_system,
    int dualp_flag, bool apply_mapaxes, int nx, int ny, int nz, int grid_nr,
    int coords_size, int **coords, float **corners, const float *mapaxes) {
    if (dualp_flag != FILEHEAD_SINGLE_POROSITY)
        nz = nz / 2;

    auto grid =
        rd_grid_ptr(rd_grid_alloc_empty(global_grid, unit_system, dualp_flag,
                                        nx, ny, nz, grid_nr, false),
                    &rd_grid_free);
    if (grid) {
        if (mapaxes != NULL)
            rd_grid_init_mapaxes(grid.get(), apply_mapaxes, mapaxes);

        for (int i = 0; i < num_coords; i++)
            rd_grid_set_cell_GRID(grid.get(), coords_size, coords[i],
                                  corners[i]);

        rd_grid_init_coarse_cells(grid.get());
        rd_grid_update_index(grid.get());
        rd_grid_taint_cells(grid.get());
    }
    return grid;
}

static int rd_grid_dual_porosity_GRID_check(rd_file_type *rd_file) {
    rd_kw_type *dimens_kw = rd_file_iget_named_kw(rd_file, DIMENS_KW, 0);
    int nx = rd_kw_iget_int(dimens_kw, DIMENS_NX_INDEX);
    int ny = rd_kw_iget_int(dimens_kw, DIMENS_NY_INDEX);
    int nz = rd_kw_iget_int(dimens_kw, DIMENS_NZ_INDEX);

    if ((nz % 2) == 1)
        return FILEHEAD_SINGLE_POROSITY;
    else {
        int dualp_flag = FILEHEAD_DUAL_POROSITY;
        int num_corners = rd_file_get_num_named_kw(rd_file, CORNERS_KW);
        int matrix_index = 0;
        int fracture_index;

        rd_kw_type *matrix_kw;
        rd_kw_type *fracture_kw;

        if (num_corners > nx * ny * nz)
            fracture_index = nx * ny * nz / 2;
        else
            fracture_index = num_corners / 2;

        while (true) {
            matrix_kw =
                rd_file_iget_named_kw(rd_file, CORNERS_KW, matrix_index);
            fracture_kw =
                rd_file_iget_named_kw(rd_file, CORNERS_KW, fracture_index);

            if (!rd_kw_equal(matrix_kw, fracture_kw)) {
                dualp_flag = FILEHEAD_SINGLE_POROSITY;
                break;
            }

            matrix_index++;
            fracture_index++;
            if (fracture_index == nx * ny * nz)
                break;
        }

        return dualp_flag;
    }
}

static rd_grid_ptr rd_grid_alloc_GRID__(rd_grid_type *global_grid,
                                        const rd_file_type *rd_file,
                                        int cell_offset, int grid_nr,
                                        int dualp_flag, bool apply_mapaxes) {
    int nx, ny, nz;
    const float *mapaxes_data = NULL;
    ert_rd_unit_enum unit_system = RD_METRIC_UNITS;

    // 1: Fetching header data from the DIMENS keyword.
    {
        rd_kw_type *dimens_kw =
            rd_file_iget_named_kw(rd_file, DIMENS_KW, grid_nr);
        nx = rd_kw_iget_int(dimens_kw, DIMENS_NX_INDEX);
        ny = rd_kw_iget_int(dimens_kw, DIMENS_NY_INDEX);
        nz = rd_kw_iget_int(dimens_kw, DIMENS_NZ_INDEX);
    }

    // 2: Fetching the mapaxes data from the MAPAXES keyword; the
    //    keyword is optional, and is only applicable to the global grid.
    {
        if ((grid_nr == 0) && (rd_file_has_kw(rd_file, MAPAXES_KW))) {
            const rd_kw_type *mapaxes_kw =
                rd_file_iget_named_kw(rd_file, MAPAXES_KW, 0);
            mapaxes_data = rd_grid_get_mapaxes_from_kw__(mapaxes_kw);
        }
    }

    if ((grid_nr == 0) && (rd_file_has_kw(rd_file, GRIDUNIT_KW)))
        unit_system = rd_grid_check_unit_system(
            rd_file_iget_named_kw(rd_file, GRIDUNIT_KW, 0));

    /*
    The number of COORDS/CORNERS blocks depends on the GRIDFILE option
    used in the ECLIPSE datafile when asking ECLIPSE to save a file in
    GRID format. If the GRIDFILE option is set to 1; only active
    cells are stored. In that case:

      1. We have absolutely no real-world coordinate information about
         the inactive cells.

      2. The GRID file only contains the global information, and can
         not contain LGR information.

    If the GRIDFILE option is set to 2 ECLIPSE will save coordinate
    information for all the cells in the grid, and also lgrs. To
    determine how many coords keywords we should read we use the
    following heuristics:

       1. If global_grid != NULL : this is an lgr, and num_coords =
          nx*ny*nz;

       2. If global_grid == NULL : this is the global grid. We check
          how many COORDS keywords there are in the file:

           a) If there are >= nx*ny*nz keywords we set num_coords ==
              nx*ny*nz.

           b) If there are < nx*ny*nz keywords we set num_coords ==
              #COORDS keywords in the file.

    Possible LGR cells will follow *AFTER* the first nx*ny*nz cells;
    the loop stops at nx*ny*nz. Additionally the LGR cells should be
    discarded (by checking coords[5]) in the rd_grid_set_cell_GRID()
    function.
  */

    int num_coords;

    if (global_grid == NULL) {
        /* This is the main grid - can be both nactive or nx*ny*nz coord elements. */
        int num_coords_kw = rd_file_get_num_named_kw(rd_file, COORDS_KW);
        if (num_coords_kw >= nx * ny * nz)
            num_coords = nx * ny * nz;
        else
            num_coords = num_coords_kw;
    } else
        /* This is an lgr - always nx*ny*nz elements. */
        num_coords = nx * ny * nz;

    // 3: Fetching the main chunk of cell data from the COORDS and
    //    CORNERS keywords.
    int coords_size = -1;

    std::vector<int *> coords(num_coords);
    std::vector<float *> corners(num_coords);

    for (int index = 0; index < num_coords; index++) {
        const rd_kw_type *coords_kw =
            rd_file_iget_named_kw(rd_file, COORDS_KW, index + cell_offset);
        const rd_kw_type *corners_kw =
            rd_file_iget_named_kw(rd_file, CORNERS_KW, index + cell_offset);

        coords[index] = rd_kw_get_int_ptr(coords_kw);
        corners[index] = rd_kw_get_float_ptr(corners_kw);
        coords_size = rd_kw_get_size(coords_kw);
    }
    // Create the grid:
    auto grid = rd_grid_alloc_GRID_data__(
        global_grid, num_coords, unit_system, dualp_flag, apply_mapaxes, nx, ny,
        nz, grid_nr, coords_size, coords.data(), corners.data(), mapaxes_data);

    if (grid_nr > 0)
        rd_grid_set_lgr_name_GRID(grid.get(), rd_file, grid_nr);
    return grid;
}

static rd_grid_ptr rd_grid_alloc_GRID(const char *grid_file,
                                      bool apply_mapaxes) {

    rd_file_enum file_type;
    file_type = rd_get_file_type(grid_file, NULL, NULL);
    if (file_type != RD_GRID_FILE)
        throw std::invalid_argument(fmt::format(
            "{}: wrong file type - expected .GRID file", grid_file));

    int cell_offset = 0;
    rd_file_ptr rd_file{rd_file_open(grid_file, 0), &rd_file_close};
    int num_grid = rd_file_get_num_named_kw(rd_file.get(), DIMENS_KW);
    int grid_nr;
    int dualp_flag;

    dualp_flag = rd_grid_dual_porosity_GRID_check(rd_file.get());
    auto main_grid = rd_grid_alloc_GRID__(NULL, rd_file.get(), cell_offset, 0,
                                          dualp_flag, apply_mapaxes);
    cell_offset += rd_grid_get_global_size(main_grid.get());

    for (grid_nr = 1; grid_nr < num_grid; grid_nr++) {
        auto lgr_grid =
            rd_grid_alloc_GRID__(main_grid.get(), rd_file.get(), cell_offset,
                                 grid_nr, dualp_flag, false);
        cell_offset += rd_grid_get_global_size(lgr_grid.get());
        auto *lgr = rd_grid_add_lgr(main_grid.get(), std::move(lgr_grid));
        {
            rd_grid_type *host_grid;

            if (!lgr->parent_name)
                host_grid = main_grid.get();
            else
                host_grid =
                    rd_grid_get_lgr(main_grid.get(), lgr->parent_name->c_str());

            rd_grid_install_lgr_GRID(host_grid, lgr);
        }
    }
    main_grid->name = grid_file;
    return main_grid;
}

/**
   This function will allocate a new regular grid with dimensions nx x
   ny x nz. The cells in the grid are spanned by the the three unit
   vectors ivec, jvec and kvec.

   The actnum argument should be a pointer to an integer array of
   length nx*ny*nz where actnum[i + j*nx + k*nx*ny] == 1 for active
   cells and 0 for inactive cells. The actnum array can be NULL, in
   which case all cells will be active.
*/
static rd_grid_type *rd_grid_alloc_regular(int nx, int ny, int nz,
                                           const double *ivec,
                                           const double *jvec,
                                           const double *kvec,
                                           const int *actnum) {
    ert_rd_unit_enum unit_system = RD_METRIC_UNITS;
    auto grid = rd_grid_ptr(rd_grid_alloc_empty(NULL, unit_system,
                                                FILEHEAD_SINGLE_POROSITY, nx,
                                                ny, nz, 0, true),
                            &rd_grid_free);
    if (grid) {
        const double grid_offset[3] = {0, 0, 0};

        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int global_index = i + j * nx + k * nx * ny;
                    double offset[3] = {grid_offset[0] + i * ivec[0] +
                                            j * jvec[0] + k * kvec[0],
                                        grid_offset[1] + i * ivec[1] +
                                            j * jvec[1] + k * kvec[1],
                                        grid_offset[2] + i * ivec[2] +
                                            j * jvec[2] + k * kvec[2]};

                    rd_cell_type &cell = grid->cells.at(global_index);
                    rd_cell_init_regular(cell, offset, i, j, k, global_index,
                                         ivec, jvec, kvec, actnum);
                }
            }
        }
        rd_grid_update_index(grid.get());
    }

    return grid.release();
}

/**
   This function will allocate a new rectangular grid with dimensions
   nx x ny x nz. The cells in the grid are rectangular with dimensions
   dx x dy x dz.

   The actnum argument should be a pointer to an integer array of
   length nx*ny*nz where actnum[i + j*nx + k*nx*ny] == 1 for active
   cells and 0 for inactive cells. The actnum array can be NULL, in
   which case all cells will be active.
*/

rd_grid_type *rd_grid_alloc_rectangular(int nx, int ny, int nz, double dx,
                                        double dy, double dz,
                                        const int *actnum) {
    const double ivec[3] = {dx, 0, 0};
    const double jvec[3] = {0, dy, 0};
    const double kvec[3] = {0, 0, dz};

    return rd_grid_alloc_regular(nx, ny, nz, ivec, jvec, kvec, actnum);
}

/**
   This function will allocate a rd_grid instance. As input it takes
   a filename, which can be both a GRID file and an EGRID file (both
   formatted and unformatted).

   When allocating based on an EGRID file the COORDS, ZCORN and ACTNUM
   keywords are extracted, and the rd_grid_alloc_GRDECL() function is
   called with these keywords. This function can be called directly
   with these keywords.
*/

static rd_grid_ptr rd_grid_alloc__(const char *grid_file, bool apply_mapaxes) {
    rd_file_enum file_type;
    rd_grid_ptr rd_grid{nullptr, &rd_grid_free};

    file_type = rd_get_file_type(grid_file, NULL, NULL);
    if (file_type == RD_GRID_FILE)
        rd_grid = rd_grid_alloc_GRID(grid_file, apply_mapaxes);
    else if (file_type == RD_EGRID_FILE)
        rd_grid = rd_grid_alloc_EGRID(grid_file, apply_mapaxes);
    else
        throw std::invalid_argument(fmt::format(
            "Must have .GRID or .EGRID file - {} not recognized", grid_file));

    return rd_grid;
}

rd_grid_type *rd_grid_alloc(const char *grid_file) {
    bool apply_mapaxes = true;
    return rd_grid_alloc__(grid_file, apply_mapaxes).release();
}

/**
   Will load the grid corresponding to the input @input_case;
   depending on the value of @input_case many different paths will be
   tried:

   1 case_input - an existing GRID/EGRID file: Just load the file -
     with no further ado.

   2 case_input - an existing file which is not a grid file;
     if it has definite formatted/unformatted status look only for
     those GRID/EGRID with the same formatted/unformatted status.

   3 case_input is only a basename, look for
     formatted/unformatted files with the correct basename.


   For cases 2 & 3 the function will look for files in the following order:

      BASE.EGRID   BASE.GRID   BASE.FEGRID   BASE.FGRID

   and stop with the first success. Will return NULL if no GRID/EGRID
   files can be found.
*/

static std::optional<std::string>
rd_grid_alloc_case_filename(const char *case_input) {
    rd_file_enum file_type;
    bool fmt_file;
    file_type = rd_get_file_type(case_input, &fmt_file, NULL);
    fs::path case_path(case_input);
    case_path = case_path.parent_path() / case_path.stem();

    if (file_type == RD_GRID_FILE)
        return case_input; /* Case 1 */
    else if (file_type == RD_EGRID_FILE)
        return case_input; /* Case 1 */
    else {
        std::optional<std::string> grid_file = std::nullopt;
        if ((file_type == RD_OTHER_FILE) ||
            (file_type ==
             RD_DATA_FILE)) { /* Case 3 - only basename recognized */
            fs::path EGRID = rd::filename(case_path, RD_EGRID_FILE, false);
            fs::path GRID = rd::filename(case_path, RD_GRID_FILE, false);
            fs::path FEGRID = rd::filename(case_path, RD_EGRID_FILE, true);
            fs::path FGRID = rd::filename(case_path, RD_GRID_FILE, true);

            if (rd::try_exists(EGRID))
                grid_file = EGRID.string();
            else if (rd::try_exists(GRID))
                grid_file = GRID.string();
            else if (rd::try_exists(FEGRID))
                grid_file = FEGRID.string();
            else if (rd::try_exists(FGRID))
                grid_file = FGRID.string();
            /* else: could not find a GRID/EGRID. */
        } else { /* Case 2 - we know the formatted / unformatted status. */
            fs::path EGRID = rd::filename(case_path, RD_EGRID_FILE, fmt_file);
            fs::path GRID = rd::filename(case_path, RD_GRID_FILE, fmt_file);

            if (rd::try_exists(EGRID))
                grid_file = EGRID.string();
            else if (rd::try_exists(GRID))
                grid_file = GRID.string();
        }
        return grid_file;
    }
}

rd_grid_type *rd_grid_load_case__(const char *case_input, bool apply_mapaxes) {
    rd_grid_type *rd_grid = NULL;
    auto grid_file = rd_grid_alloc_case_filename(case_input);
    if (grid_file.has_value() && util_file_exists(grid_file.value().c_str())) {
        rd_grid =
            rd_grid_alloc__(grid_file.value().c_str(), apply_mapaxes).release();
    }
    return rd_grid;
}

rd_grid_type *rd_grid_load_case(const char *case_input) {
    bool apply_mapaxes = true;
    return rd_grid_load_case__(case_input, apply_mapaxes);
}

bool rd_grid_exists(const char *case_input) {
    return rd_grid_alloc_case_filename(case_input).has_value();
}

static bool rd_grid_compare_coarse_cells(const rd_grid_type *g1,
                                         const rd_grid_type *g2, bool verbose) {
    if (g1->coarse_cells.size() == g2->coarse_cells.size()) {
        bool equal = true;

        for (size_t c = 0; c < g1->coarse_cells.size(); c++) {
            const rd_coarse_cell_type *coarse_cell1 = g1->coarse_cells[c].get();
            const rd_coarse_cell_type *coarse_cell2 = g2->coarse_cells[c].get();

            equal = rd_coarse_cell_equal(coarse_cell1, coarse_cell2);
            if (!equal)
                if (verbose)
                    fprintf(stderr, "Difference in coarse cell:%zu \n", c);
        }
        return equal;
    } else
        return false;
}

static bool rd_grid_compare_cells(rd_grid_type *g1, rd_grid_type *g2,
                                  bool include_nnc, bool verbose) {
    int g;
    bool equal = true;
    for (g = 0; g < g1->size; g++) {
        bool this_equal = true;
        rd_cell_type &c1 = g1->cells.at(g);
        rd_cell_type &c2 = g2->cells.at(g);
        rd_cell_compare(c1, c2, include_nnc, &this_equal);

        if (!this_equal) {
            if (verbose) {
                int i, j, k;
                rd_grid_get_ijk1(g1, g, &i, &j, &k);

                printf("Difference in cell: %d : %d,%d,%d  nnc_equal:%d "
                       "Volume:%g \n",
                       g, i, j, k,
                       nnc_info_equal(c1.nnc_info.get(), c2.nnc_info.get()),
                       rd_cell_get_volume(c1));
                printf("-------------------------------------------------------"
                       "----------\n");
                rd_cell_dump_ascii(c1, i, j, k, stdout, NULL);
                printf("-------------------------------------------------------"
                       "----------\n");
                rd_cell_dump_ascii(c2, i, j, k, stdout, NULL);
                printf("-------------------------------------------------------"
                       "----------\n");
            }
            equal = false;
            break;
        }
    }
    return equal;
}

static bool rd_grid_compare_index(const rd_grid_type *g1,
                                  const rd_grid_type *g2, bool verbose) {
    bool equal = true;

    if (g1->total_active != g2->total_active) {
        if (verbose)
            fprintf(stderr, "Difference in total active:%d / %d\n",
                    g1->total_active, g2->total_active);
        equal = false;
    }

    if (equal) {
        if (g1->index_map != g2->index_map) {
            equal = false;
            if (verbose)
                fprintf(stderr, "Difference in index map \n");
        }
    }

    if (equal) {
        if (g1->inv_index_map != g2->inv_index_map) {
            equal = false;
            if (verbose)
                fprintf(stderr, "Difference in inverse index map \n");
        }
    }

    if (equal && (g1->dualp_flag != FILEHEAD_SINGLE_POROSITY)) {
        if (g1->total_active_fracture != g2->total_active_fracture) {
            if (verbose)
                fprintf(stderr, "Difference in toal_active_fracture %d / %d \n",
                        g1->total_active_fracture, g2->total_active_fracture);
            equal = false;
        }

        if (equal) {
            if (g1->fracture_index_map != g2->fracture_index_map) {
                equal = false;
                if (verbose)
                    fprintf(stderr, "Difference in fracture_index_map \n");
            }
        }

        if (equal) {
            if (g1->inv_fracture_index_map != g2->inv_fracture_index_map) {
                equal = false;
                if (verbose)
                    fprintf(stderr, "Difference in inv_fracture_index_map \n");
            }
        }
    }
    return equal;
}

static bool rd_grid_compare_mapaxes(const rd_grid_type *g1,
                                    const rd_grid_type *g2, bool verbose) {
    bool equal = true;
    if (g1->use_mapaxes == g2->use_mapaxes) {
        if (g1->mapaxes != g2->mapaxes)
            equal = false;
    } else
        equal = false;

    if (!equal && verbose)
        fprintf(stderr, "Difference in mapaxes \n");

    return equal;
}

/**
   Return true if grids g1 and g2 are equal, and false otherwise. To
   return true all cells must be identical.
*/

static bool rd_grid_compare__(rd_grid_type *g1, rd_grid_type *g2,
                              bool include_nnc, bool verbose) {

    bool equal = true;
    if (g1->size != g2->size)
        equal = false;

    // The name of the parent grid corresponds to a filename; they can be different.
    if (equal && g1->parent_grid) {
        if (g1->name != g2->name) {
            equal = false;
            if (verbose)
                fprintf(stderr, "Difference in name %s <-> %s \n",
                        g1->name.c_str(), g2->name.c_str());
        }
    }

    /*
    When .GRID files are involved this is hardwired to FILEHEAD_SINGLE_POROSITY.
  */
    if (g1->dualp_flag != g2->dualp_flag) {
        equal = false;
        if (verbose)
            fprintf(stderr, "Dual porosity flags differ: %d / %d \n",
                    g1->dualp_flag, g2->dualp_flag);
    }

    if (equal)
        equal = rd_grid_compare_cells(g1, g2, include_nnc, verbose);

    if (equal)
        equal = rd_grid_compare_index(g1, g2, true /*verbose*/);

    if (equal)
        equal = rd_grid_compare_coarse_cells(g1, g2, verbose);

    if (equal)
        equal = rd_grid_compare_mapaxes(g1, g2, verbose);

    return equal;
}

bool rd_grid_compare(rd_grid_type *g1, rd_grid_type *g2, bool include_lgr,
                     bool include_nnc, bool verbose) {
    bool equal = rd_grid_compare__(g1, g2, include_nnc, verbose);

    if (equal && include_lgr) {
        if (g1->LGR_list.size() == g2->LGR_list.size()) {
            for (size_t grid_nr = 0; grid_nr < g1->LGR_list.size(); grid_nr++) {
                auto &lgr1 = g1->LGR_list.at(grid_nr);
                auto &lgr2 = g2->LGR_list.at(grid_nr);

                printf("Comparing lgr grid:%zu \n", grid_nr);
                equal = rd_grid_compare__(lgr1.get(), lgr2.get(), include_nnc,
                                          verbose);
                if (!equal)
                    break;
            }
        } else
            equal = false;
    }
    return equal;
}

typedef enum {
    NOT_ON_FACE,
    BELONGS_TO_CELL,
    BELONGS_TO_OTHER
} face_status_enum;

/*
    Returns whether the given point is contained within the minimal cube
    encapsulating the cell that has all faces parallel to a coordinate plane.
*/
static bool rd_grid_cube_contains(const rd_cell_type &cell,
                                  const point_type *p) {
    if (p->z < rd_cell_min_z(cell))
        return false;

    if (p->z > rd_cell_max_z(cell))
        return false;

    if (p->x < rd_cell_min_x(cell))
        return false;

    if (p->x > rd_cell_max_x(cell))
        return false;

    if (p->y < rd_cell_min_y(cell))
        return false;

    if (p->y > rd_cell_max_y(cell))
        return false;

    return true;
}

/*
   Returns true if and only if p is on plane "plane" of cell when decomposed by "method".
*/
static bool rd_grid_on_plane(const rd_cell_type &cell, const int method,
                             const int plane, const point_type *p) {
    const point_type *p0 =
        &cell.corner_list[tetrahedron_permutations[method][plane][0]];
    const point_type *p1 =
        &cell.corner_list[tetrahedron_permutations[method][plane][1]];
    const point_type *p2 =
        &cell.corner_list[tetrahedron_permutations[method][plane][2]];
    return triangle_contains3d(p0, p1, p2, p);
}

/*
   Returns true if and only if p is on one of the cells faces and
   "belongs" to this cell. This is done such that every point is contained in at most
   one point.

   Known caveats when using this function:
     - if a point is on the surface of a/many cells, but for all of these cells
       the point is contained on two opposite sides of the cell. Imagine a cake
       being cut as a cake should be cut. To which of the slices does the center
       point of the cake belong? This is a somewhat obscure situation and it is
       not possible to circumvent by only considering the grid cell by cell.
    - if there is a fault and this cell is on the border of the grid.
    - if a cells projection to the xy-plane is concave, this method might give
      false positives.

   Note: The correctness of this function relies *HEAVILY* on the permutation of the
   tetrahedrons in the decompositions.
*/
static face_status_enum rd_grid_on_cell_face(const rd_cell_type &cell,
                                             const int method,
                                             const point_type *p,
                                             const bool max_i, const bool max_j,
                                             const bool max_k) {

    int k_minus = 0, j_pluss = 1, i_minus = 2, i_pluss = 3, j_minus = 4,
        k_pluss = 5;
    bool on[6];
    for (int i = 0; i < 6; ++i) {
        on[i] = (rd_grid_on_plane(cell, method, 2 * i, p) ||
                 rd_grid_on_plane(cell, method, 2 * i + 1, p));
    }

    // Not on any of the cell sides
    if (!on[k_minus] && !on[k_pluss] && !on[j_pluss] && !on[j_minus] &&
        !on[i_minus] && !on[i_pluss])
        return NOT_ON_FACE;

    // Handles side collapses, i.e. the point is contained on opposite sides.
    // Cell passes on the responsibility if not on border of grid.
    bool i_collapse = (on[i_minus] && on[i_pluss]);
    bool j_collapse = (on[j_minus] && on[j_pluss]);
    bool k_collapse = (on[k_minus] && on[k_pluss]);

    for (int i = 0; i < 6; ++i)
        on[i] &= (!i_collapse || max_i) && (!j_collapse || max_j) &&
                 (!k_collapse || max_k);

    on[i_minus] &= !on[i_pluss];
    on[j_minus] &= !on[j_pluss];
    on[k_minus] &= !on[k_pluss];

    // Removed from all sides
    if (!on[k_minus] && !on[k_pluss] && !on[j_pluss] && !on[j_minus] &&
        !on[i_minus] && !on[i_pluss])
        return BELONGS_TO_OTHER;

    // Not on any of the lower priority sides
    if (!on[k_pluss] && !on[j_pluss] && !on[i_pluss])
        return BELONGS_TO_CELL;

    // Contained in cell due to border conditions
    // NOTE: One should read X <= Y as X "implies" Y
    if ((on[i_pluss] <= max_i) && (on[j_pluss] <= max_j) &&
        (on[k_pluss] <= max_k))
        return BELONGS_TO_CELL;

    return BELONGS_TO_OTHER;
}

/*
 Returns true if and only if the tetrahedron defined by p0, p1, p2, p3
 contains p.

 The sole purpose of this functions is to make concave_cell_contains
 more readable.
*/
static bool tetrahedron_by_points_contains(const point_type *p0,
                                           const point_type *p1,
                                           const point_type *p2,
                                           const point_type *p3,
                                           const point_type *p) {

    tetrahedron_type pro_tet;
    pro_tet.p0 = *p0;
    pro_tet.p1 = *p1;
    pro_tet.p2 = *p2;
    pro_tet.p3 = *p3;

    return tetrahedron_contains(pro_tet, *p);
}

/*
 Returns true if and only if the cell "cell" decomposed by "method" contains the point "p".
 This is done by decomposing the cell into 5 tetrahedrons according to the decomposition
 method for the faces.

 Assumes the cell to not be self-intersecting!

 Note: This function relies *HEAVILY* on the permutation of tetrahedron_permutations.
*/
static bool concave_cell_contains(const rd_cell_type &cell, int method,
                                  const point_type *p) {

    const point_type *dia[2][2] = {
        {&cell.corner_list[tetrahedron_permutations[method][0][1]],
         &cell.corner_list[tetrahedron_permutations[method][0][2]]},
        {&cell.corner_list[tetrahedron_permutations[method][10][1]],
         &cell.corner_list[tetrahedron_permutations[method][10][2]]}};

    const point_type *extra[2][2] = {
        {&cell.corner_list[tetrahedron_permutations[method][0][0]],
         &cell.corner_list[tetrahedron_permutations[method][1][0]]},
        {&cell.corner_list[tetrahedron_permutations[method][10][0]],
         &cell.corner_list[tetrahedron_permutations[method][11][0]]}};

    // Test for containment in cell core
    if (tetrahedron_by_points_contains(dia[0][0], dia[1][0], dia[0][1],
                                       dia[1][1], p))
        return true;

    // Test for containment in protrusions
    for (int i = 0; i < 2; ++i) {
        if (tetrahedron_by_points_contains(dia[i][0], dia[i][1],
                                           dia[(i + 1) % 2][0], extra[i][0], p))
            return true;

        if (tetrahedron_by_points_contains(dia[i][0], dia[(i + 1) % 2][1],
                                           dia[i][1], extra[i][1], p))
            return true;
    }

    return false;
}

/*
  Observe the following quirks with this functions:

  - It is quite simple to create a cell where the center point is
    actually *not* inside the cell - that might come as a surprise!

  - Cells with nonzero twist are completely discarded from the search,
    if the point (x,y,z) "should" have been found on the inside of a
    twisted cell the algorithm will incorrectly return false; a
    warning will be printed on stderr if a cell is discarded due to
    twist.

  - See the documentation of rd_grid_on_cell_face for caveats regarding
    containtment of points of cell faces.
*/
static bool rd_grid_cell_contains_xyz3(const rd_grid_type *rd_grid, int i,
                                       int j, int k, double x, double y,
                                       double z) {
    point_type p;
    const rd_cell_type &cell =
        rd_grid->cells.at(rd_grid_get_global_index3(rd_grid, i, j, k));
    point_set(&p, x, y, z);
    int method =
        (i + j + k) %
        2; // Chooses the approperiate decomposition method for the cell

    if (GET_CELL_FLAG((&cell), CELL_FLAG_TAINTED))
        return false;

    // Pruning
    if (!rd_grid_cube_contains(cell, &p))
        return false;

    // Checks if point is on one of the faces of the cell, and if so whether it
    // "belongs" to this cell.
    bool max_i = (i == rd_grid->nx - 1);
    bool max_j = (j == rd_grid->ny - 1);
    bool max_k = (k == rd_grid->nz - 1);
    face_status_enum face_status =
        rd_grid_on_cell_face(cell, method, &p, max_i, max_j, max_k);

    if (face_status != NOT_ON_FACE)
        return face_status == BELONGS_TO_CELL;

    // Twisted cells
    if (rd_cell_get_twist(cell) > 0) {
        fprintf(stderr,
                "** Warning: Point (%g,%g,%g) is in vicinity of twisted cell: "
                "(%d,%d,%d) - function:%s might be mistaken.\n",
                x, y, z, i, j, k, __func__);
        return false;
    }

    // We now check whether the point is strictly inside the cell
    return concave_cell_contains(cell, method, &p);
}

bool rd_grid_cell_contains_xyz1(const rd_grid_type *rd_grid, int global_index,
                                double x, double y, double z) {
    int i, j, k;
    rd_grid_get_ijk1(rd_grid, global_index, &i, &j, &k);
    return rd_grid_cell_contains_xyz3(rd_grid, i, j, k, x, y, z);
}

static void rd_grid_clear_visited(rd_grid_type *grid) {
    grid->visited.assign(grid->size, false);
}

/*
   Box coordinates are not inclusive, i.e. [i1,i2)
*/
static int rd_grid_box_contains_xyz(rd_grid_type *grid, int i1, int i2, int j1,
                                    int j2, int k1, int k2,
                                    const point_type *p) {

    int global_index = -1;
    for (int k = k1; k < k2; k++)
        for (int j = j1; j < j2; j++)
            for (int i = i1; i < i2; i++) {
                global_index = rd_grid_get_global_index3(grid, i, j, k);
                if (!grid->visited[global_index]) {
                    grid->visited[global_index] = true;
                    if (rd_grid_cell_contains_xyz1(grid, global_index, p->x,
                                                   p->y, p->z)) {
                        return global_index;
                    }
                }
            }
    return -1; /* Returning -1; did not find xyz. */
}

/**
 * Search for given xyz coordinate around global start_index in a box of size bx
 */
static int rd_grid_get_global_index_from_xyz_around_box(rd_grid_type *grid,
                                                        double x, double y,
                                                        double z,
                                                        int start_index, int bx,
                                                        point_type *p) {
    /* Try neighbours */
    int i, j, k;
    int i1, i2, j1, j2, k1, k2;
    int nx, ny, nz;
    rd_grid_get_dims(grid, &nx, &ny, &nz, NULL);
    rd_grid_get_ijk1(grid, start_index, &i, &j, &k);

    i1 = util_int_max(0, i - bx);
    j1 = util_int_max(0, j - bx);
    k1 = util_int_max(0, k - bx);

    i2 = util_int_min(nx, i + bx);
    j2 = util_int_min(ny, j + bx);
    k2 = util_int_min(nz, k + bx);

    int global_index =
        rd_grid_box_contains_xyz(grid, i1, i2, j1, j2, k1, k2, p);
    return global_index;
}

/**
   This function will find the global index of the cell containing the
   world coordinates (x,y,z), if no cell can be found the function
   will return -1.

   The function is basically based on scanning through the cells in
   natural (i fastest) order and querying whether the cell[i,j,k]
   contains the (x,y,z) point; not very elegant :-(

   The last argument - 'start_index' - can be used to speed things up
   a bit if you have reasonable guess of where the the (x,y,z) is
   located. The start_index value is used as this:


     start_index == 0: I do not have a clue, start from the beginning
        and scan through the grid linearly.


     start_index != 0:
        1. Check the cell 'start_index'.
        2. Check the neighbours (i +/- 1, j +/- 1, k +/- 1 ).
        3. Give up and do a linear search starting from start_index.

*/
int rd_grid_get_global_index_from_xyz(rd_grid_type *grid, double x, double y,
                                      double z, int start_index) {
    int global_index;
    point_type p;
    point_set(&p, x, y, z);
    rd_grid_clear_visited(grid);

    if (start_index >= 0) {
        /* Try start index */
        if (rd_grid_cell_contains_xyz1(grid, start_index, x, y, z))
            return start_index;
        else {
            /* Try boxes 2, 4, 8, ..., 64  */
            for (int bx = 1; bx <= 6; bx++) {
                global_index = rd_grid_get_global_index_from_xyz_around_box(
                    grid, x, y, z, start_index, 1 << bx, &p);
                if (global_index >= 0)
                    return global_index;
            }
        }
    }

    /*
    OK - the attempted shortcuts did not pay off. Perform full linear search.
  */

    global_index = -1;

    for (int index = 0; index < grid->size; index++) {
        if (rd_grid_cell_contains_xyz1(grid, index, x, y, z))
            return index;
    }
    return -1;
}

static bool rd_grid_sublayer_contanins_xy__(const rd_grid_type *grid, double x,
                                            double y, int k, int i1, int i2,
                                            int j1, int j2,
                                            geo_polygon_type *polygon) {
    geo_polygon_reset(polygon);

    /* Bottom edge */
    for (int i = i1; i < i2; i++) {
        double corner_pos[3];
        rd_grid_get_corner_xyz(grid, i, j1, k, &corner_pos[0], &corner_pos[1],
                               &corner_pos[2]);
        geo_polygon_add_point(polygon, corner_pos[0], corner_pos[1]);
    }

    /* Right edge */
    for (int j = j1; j < j2; j++) {
        double corner_pos[3];
        rd_grid_get_corner_xyz(grid, i2, j, k, &corner_pos[0], &corner_pos[1],
                               &corner_pos[2]);
        geo_polygon_add_point(polygon, corner_pos[0], corner_pos[1]);
    }

    /* Top edge */
    for (int i = i2; i > i1; i--) {
        double corner_pos[3];
        rd_grid_get_corner_xyz(grid, i, j2, k, &corner_pos[0], &corner_pos[1],
                               &corner_pos[2]);
        geo_polygon_add_point(polygon, corner_pos[0], corner_pos[1]);
    }

    /* Left edge */
    for (int j = j2; j > j1; j--) {
        double corner_pos[3];
        rd_grid_get_corner_xyz(grid, i1, j, k, &corner_pos[0], &corner_pos[1],
                               &corner_pos[2]);
        geo_polygon_add_point(polygon, corner_pos[0], corner_pos[1]);
    }
    geo_polygon_close(polygon);
    return geo_polygon_contains_point__(polygon, x, y, true);
}

bool rd_grid_get_ij_from_xy(const rd_grid_type *grid, double x, double y, int k,
                            int *i, int *j) {
    auto polygon = make_geo_polygon("");
    int nx = rd_grid_get_nx(grid);
    int ny = rd_grid_get_ny(grid);
    bool inside = rd_grid_sublayer_contanins_xy__(grid, x, y, k, 0, nx, 0, ny,
                                                  polygon.get());
    if (inside) {
        int i1 = 0;
        int i2 = nx;
        int j1 = 0;
        int j2 = ny;

        while (true) {
            if ((i2 - i1) > 1) {
                int ic = (i1 + i2) / 2;
                if (rd_grid_sublayer_contanins_xy__(grid, x, y, k, i1, ic, j1,
                                                    j2, polygon.get()))
                    i2 = ic;
                else {
                    if (!rd_grid_sublayer_contanins_xy__(grid, x, y, k, ic, i2,
                                                         j1, j2, polygon.get()))
                        throw std::logic_error("point nowhere to be found ...");
                    i1 = ic;
                }
            }

            if ((j2 - j1) > 1) {
                int jc = (j1 + j2) / 2;
                if (rd_grid_sublayer_contanins_xy__(grid, x, y, k, i1, i2, j1,
                                                    jc, polygon.get()))
                    j2 = jc;
                else {
                    if (!rd_grid_sublayer_contanins_xy__(grid, x, y, k, i1, i2,
                                                         jc, j2, polygon.get()))
                        throw std::logic_error("point nowhere to be found ...");
                    j1 = jc;
                }
            }

            if ((i2 - i1) == 1 && (j2 - j1) == 1) {
                *i = i1;
                *j = j1;
                break;
            }
        }
    }
    return inside;
}

void rd_grid_free(rd_grid_type *grid) { delete grid; }

void rd_grid_get_distance(rd_grid_type *grid, int global_index1,
                          int global_index2, double *dx, double *dy,
                          double *dz) {
    rd_cell_type &cell1 = grid->cells.at(global_index1);
    rd_cell_type &cell2 = grid->cells.at(global_index2);

    rd_cell_assert_center(cell1);
    rd_cell_assert_center(cell2);
    {
        *dx = cell1.center.x - cell2.center.x;
        *dy = cell1.center.y - cell2.center.y;
        *dz = cell1.center.z - cell2.center.z;
    }
}

/**
   Only checks that i,j,k are in the required intervals:

      0 <= i < nx
      0 <= j < ny
      0 <= k < nz

*/

bool rd_grid_ijk_valid(const rd_grid_type *grid, int i, int j, int k) {
    bool OK = false;

    if (i >= 0 && i < grid->nx)
        if (j >= 0 && j < grid->ny)
            if (k >= 0 && k < grid->nz)
                OK = true;

    return OK;
}

void rd_grid_get_dims(const rd_grid_type *grid, int *nx, int *ny, int *nz,
                      int *active_size) {
    if (nx != NULL)
        *nx = grid->nx;
    if (ny != NULL)
        *ny = grid->ny;
    if (nz != NULL)
        *nz = grid->nz;
    if (active_size != NULL)
        *active_size = grid->total_active;
}

int rd_grid_get_nz(const rd_grid_type *grid) { return grid->nz; }

int rd_grid_get_nx(const rd_grid_type *grid) { return grid->nx; }

int rd_grid_get_ny(const rd_grid_type *grid) { return grid->ny; }

int rd_grid_get_nactive(const rd_grid_type *grid) { return grid->total_active; }

int rd_grid_get_nactive_fracture(const rd_grid_type *grid) {
    return grid->total_active_fracture;
}

static int rd_grid_get_parent_cell1(const rd_grid_type *grid,
                                    int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    return cell.host_cell;
}

/**
   Converts: (i,j,k) -> global_index. i,j,k are zero offset.
*/

int rd_grid_get_global_index3(const rd_grid_type *rd_grid, int i, int j,
                              int k) {
    if (rd_grid_ijk_valid(rd_grid, i, j, k))
        return rd_grid_get_global_index__(rd_grid, i, j, k);
    else {
        throw std::out_of_range(fmt::format(
            "i,j,k = ({},{},{}) is invalid: nx: [0,{}>  ny: [0,{}>  nz: "
            "[0,{}>",
            i, j, k, rd_grid->nx, rd_grid->ny, rd_grid->nz));
    }
}

/**
   Converts: active_index -> global_index
*/

int rd_grid_get_global_index1A(const rd_grid_type *rd_grid, int active_index) {
    return rd_grid->inv_index_map[active_index];
}

int rd_grid_get_global_index1F(const rd_grid_type *rd_grid,
                               int active_fracture_index) {
    if (rd_grid->inv_fracture_index_map.empty() ||
        active_fracture_index >= rd_grid->total_active_fracture)
        return -1;
    return rd_grid->inv_fracture_index_map[active_fracture_index];
}

/**
   Converts: (i,j,k) -> active_index
   (i,j,k ) are zero offset.

   Will return -1 if the cell is not active.
*/

int rd_grid_get_active_index3(const rd_grid_type *rd_grid, int i, int j,
                              int k) {
    int global_index = rd_grid_get_global_index3(
        rd_grid, i, j, k); /* In range: [0,nx*ny*nz) */
    return rd_grid_get_active_index1(rd_grid, global_index);
}

/**
   Converts: global_index -> active_index.

   Will return -1 if the cell is not active.
*/

int rd_grid_get_active_index1(const rd_grid_type *rd_grid, int global_index) {
    return rd_grid->index_map[global_index];
}

/**
   Converts: global_index -> active_index.

   Will return -1 if the cell is not active.
*/

int rd_grid_get_active_fracture_index1(const rd_grid_type *rd_grid,
                                       int global_index) {
    if (rd_grid->fracture_index_map.empty())
        return -1;
    else
        return rd_grid->fracture_index_map[global_index];
}

/*
  Converts global_index -> (i,j,k)

  This function returns C-based zero offset indices. cell_
*/

void rd_grid_get_ijk1(const rd_grid_type *grid, int global_index, int *i,
                      int *j, int *k) {
    *k = global_index / (grid->nx * grid->ny);
    global_index -= (*k) * (grid->nx * grid->ny);
    *j = global_index / grid->nx;
    global_index -= (*j) * grid->nx;
    *i = global_index;
}

/*
  Converts active_index -> (i,j,k)
*/

void rd_grid_get_ijk1A(const rd_grid_type *rd_grid, int active_index, int *i,
                       int *j, int *k) {
    if (active_index >= 0 && active_index < rd_grid->total_active) {
        int global_index = rd_grid_get_global_index1A(rd_grid, active_index);
        rd_grid_get_ijk1(rd_grid, global_index, i, j, k);
    } else
        throw std::out_of_range(fmt::format(
            "error active_index:{} invalid - grid has only:{} active cells.",
            active_index, rd_grid->total_active));
}

/*
  Functions to get the 'true' (i.e. UTM or whatever) position (x,y,z).

  The cell center is calculated as the plain average of the eight
  corner positions, it is quite simple to construct cells where this
  average position is on the outside of the cell - hence there is no
  guarantee that the (x,y,z) position returned from this function
  actually is on the inside of the cell.
*/

void rd_grid_get_xyz1(rd_grid_type *grid, int global_index, double *xpos,
                      double *ypos, double *zpos) {
    rd_cell_type &cell = grid->cells.at(global_index);
    rd_cell_assert_center(cell);
    {
        *xpos = cell.center.x;
        *ypos = cell.center.y;
        *zpos = cell.center.z;
    }
}

void rd_grid_get_xyz3(rd_grid_type *grid, int i, int j, int k, double *xpos,
                      double *ypos, double *zpos) {
    const int global_index = rd_grid_get_global_index__(grid, i, j, k);
    rd_grid_get_xyz1(grid, global_index, xpos, ypos, zpos);
}

/**
   This function will return (by reference) the x,y,z values of corner
   nr 'corner_nr' in cell 'global_index'. See the documentation of
   tetraheder decomposition for the numbering of the corners.
*/

void rd_grid_get_cell_corner_xyz1(const rd_grid_type *grid, int global_index,
                                  int corner_nr, double *xpos, double *ypos,
                                  double *zpos) {
    if ((corner_nr >= 0) && (corner_nr <= 7)) {
        const rd_cell_type &cell = grid->cells.at(global_index);
        const point_type point = cell.corner_list[corner_nr];
        *xpos = point.x;
        *ypos = point.y;
        *zpos = point.z;
    }
}

void rd_grid_export_cell_corners1(const rd_grid_type *grid, int global_index,
                                  double *x, double *y, double *z) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    for (int i = 0; i < 8; i++) {
        const point_type point = cell.corner_list[i];
        x[i] = point.x;
        y[i] = point.y;
        z[i] = point.z;
    }
}

static void rd_grid_get_cell_corner_xyz3(const rd_grid_type *grid, int i, int j,
                                         int k, int corner_nr, double *xpos,
                                         double *ypos, double *zpos) {
    const int global_index = rd_grid_get_global_index__(grid, i, j, k);
    rd_grid_get_cell_corner_xyz1(grid, global_index, corner_nr, xpos, ypos,
                                 zpos);
}

void rd_grid_get_corner_xyz(const rd_grid_type *grid, int i, int j, int k,
                            double *xpos, double *ypos, double *zpos) {
    if (i < 0 || i > grid->nx)
        throw std::out_of_range(fmt::format(
            "invalid i value:{}  Valid range: [0,{}]", i, grid->nx));

    if (j < 0 || j > grid->ny)
        throw std::out_of_range(fmt::format(
            "invalid j value:{}  Valid range: [0,{}]", j, grid->ny));

    if (k < 0 || k > grid->nz)
        throw std::out_of_range(fmt::format(
            "invalid k value:{}  Valid range: [0,{}]", k, grid->nz));

    {
        int corner_nr = 0;
        if (i == grid->nx) {
            i -= 1;
            corner_nr += 1;
        }

        if (j == grid->ny) {
            j -= 1;
            corner_nr += 2;
        }

        if (k == grid->nz) {
            k -= 1;
            corner_nr += 4;
        }

        rd_grid_get_cell_corner_xyz3(grid, i, j, k, corner_nr, xpos, ypos,
                                     zpos);
    }
}

void rd_grid_get_xyz1A(rd_grid_type *grid, int active_index, double *xpos,
                       double *ypos, double *zpos) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    rd_grid_get_xyz1(grid, global_index, xpos, ypos, zpos);
}

double rd_grid_get_cdepth1(rd_grid_type *grid, int global_index) {
    rd_cell_type &cell = grid->cells.at(global_index);
    rd_cell_assert_center(cell);
    return cell.center.z;
}

double rd_grid_get_cdepth1A(rd_grid_type *grid, int active_index) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    return rd_grid_get_cdepth1(grid, global_index);
}

/**
   Returns the depth of the top surface of the cell.
*/

static double rd_grid_get_top1(const rd_grid_type *grid, int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    double depth = 0;

    for (int i = 0; i < 4; i++)
        depth += cell.corner_list[i].z;

    return depth * 0.25;
}

static double rd_grid_get_top3(const rd_grid_type *grid, int i, int j, int k) {
    const int global_index = rd_grid_get_global_index__(grid, i, j, k);
    return rd_grid_get_top1(grid, global_index);
}

/**
   Returns the depth of the bottom surface of the cell.
*/

static double rd_grid_get_bottom1(const rd_grid_type *grid, int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    double depth = 0;

    for (int i = 0; i < 4; i++)
        depth += cell.corner_list[i + 4].z;

    return depth * 0.25;
}

static double rd_grid_get_bottom3(const rd_grid_type *grid, int i, int j,
                                  int k) {
    const int global_index = rd_grid_get_global_index__(grid, i, j, k);
    return rd_grid_get_bottom1(grid, global_index);
}

int rd_grid_locate_depth(const rd_grid_type *grid, double depth, int i, int j) {
    if (depth < rd_grid_get_top2(grid, i, j))
        return -1;
    else if (depth >= rd_grid_get_bottom2(grid, i, j))
        return -1 * grid->nz;
    else {
        int k = 0;
        double bottom = rd_grid_get_top3(grid, i, j, k);

        while (true) {
            double top = bottom;
            bottom = rd_grid_get_bottom3(grid, i, j, k);

            if ((depth >= top) && (depth < bottom))
                return k;

            k++;
            if (k == grid->nz)
                throw std::logic_error(fmt::format(
                    "internal error when scanning for depth:{:g}", depth));
        }
    }
}

double rd_grid_get_top2(const rd_grid_type *grid, int i, int j) {
    const int global_index = rd_grid_get_global_index__(grid, i, j, 0);
    return rd_grid_get_top1(grid, global_index);
}

double rd_grid_get_bottom2(const rd_grid_type *grid, int i, int j) {
    const int global_index =
        rd_grid_get_global_index__(grid, i, j, grid->nz - 1);
    return rd_grid_get_bottom1(grid, global_index);
}

double rd_grid_get_top1A(const rd_grid_type *grid, int active_index) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    return rd_grid_get_top1(grid, global_index);
}

double rd_grid_get_cell_dz1(const rd_grid_type *grid, int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    double dz = 0;

    for (int i = 0; i < 4; i++)
        dz += (cell.corner_list[i + 4].z - cell.corner_list[i].z);

    return dz * 0.25;
}

double rd_grid_get_cell_dz1A(const rd_grid_type *grid, int active_index) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    return rd_grid_get_cell_dz1(grid, global_index);
}

double rd_grid_get_cell_thickness1(const rd_grid_type *grid, int global_index) {
    return rd_grid_get_cell_dz1(grid, global_index);
}

double rd_grid_get_cell_dx1(const rd_grid_type *grid, int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    double dx = 0;
    double dy = 0;

    for (int c = 1; c < 8; c += 2) {
        dx += cell.corner_list[c].x - cell.corner_list[c - 1].x;
        dy += cell.corner_list[c].y - cell.corner_list[c - 1].y;
    }
    dx *= 0.25;
    dy *= 0.25;

    return sqrt(dx * dx + dy * dy);
}

double rd_grid_get_cell_dx1A(const rd_grid_type *grid, int active_index) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    return rd_grid_get_cell_dx1(grid, global_index);
}

/*
  The current algorithm for calculating the cell dimensions DX,DY and DZ
  reproduces the Eclipse results from the INIT file quite well, relative error
  on the order 1e-4 for DX and DY and 1e-3 for DZ.

  Observe that the DX, DY and DZ values are not tied to the cell volume; i.e.
  the relationship:

       DX * DY * DZ = V

  does generally not hold.
*/

double rd_grid_get_cell_dy1(const rd_grid_type *grid, int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    double dx = 0;
    double dy = 0;

    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 2; i++) {
            int c1 = i + k * 4;
            int c2 = c1 + 2;
            dx += cell.corner_list[c2].x - cell.corner_list[c1].x;
            dy += cell.corner_list[c2].y - cell.corner_list[c1].y;
        }
    }
    dx *= 0.25;
    dy *= 0.25;

    return sqrt(dx * dx + dy * dy);
}

double rd_grid_get_cell_dy1A(const rd_grid_type *grid, int active_index) {
    const int global_index = rd_grid_get_global_index1A(grid, active_index);
    return rd_grid_get_cell_dy1(grid, global_index);
}

const nnc_info_type *rd_grid_get_cell_nnc_info1(const rd_grid_type *grid,
                                                int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    return cell.nnc_info.get();
}

/*
   Global index in [0,...,nx*ny*nz)
*/

bool rd_grid_cell_active1(const rd_grid_type *rd_grid, int global_index) {
    if (rd_grid->index_map[global_index] >= 0)
        return true;
    else
        return false;
}

bool rd_grid_cell_active3(const rd_grid_type *rd_grid, int i, int j, int k) {
    int global_index = rd_grid_get_global_index3(rd_grid, i, j, k);
    return rd_grid_cell_active1(rd_grid, global_index);
}

bool rd_grid_cell_invalid1(const rd_grid_type *rd_grid, int global_index) {
    const rd_cell_type &cell = rd_grid->cells.at(global_index);
    return GET_CELL_FLAG((&cell), CELL_FLAG_TAINTED);
}

bool rd_grid_cell_valid1(const rd_grid_type *rd_grid, int global_index) {
    const rd_cell_type &cell = rd_grid->cells.at(global_index);
    if (GET_CELL_FLAG((&cell), CELL_FLAG_TAINTED))
        return false;
    else
        return (GET_CELL_FLAG((&cell), CELL_FLAG_VALID));
}

static void __assert_main_grid(const rd_grid_type *rd_grid) {
    if (rd_grid->lgr_nr != RD_GRID_MAINGRID_LGR_NR)
        throw std::invalid_argument(
            "tried to get LGR grid from another LGR_grid - only main grid "
            "can be used as first input");
}

/**
   This functon will return a a rd_grid instance corresponding to the
   lgr with name lgr_name. The function will fail HARD if no lgr with
   this name is installed under the present main grid; check first
   with rd_grid_has_lgr() if you are whimp.

   Leading/trailing spaces on lgr_name are stripped prior to the hash lookup.
*/

rd_grid_type *rd_grid_get_lgr(const rd_grid_type *main_grid,
                              const char *__lgr_name) {
    __assert_main_grid(main_grid);
    return main_grid->LGR_hash.at(rd::strip_spaces(__lgr_name));
}

/**
   Returns true/false if the main grid has a a lgr with name
   __lgr_name. Leading/trailing spaces are stripped before checking.
*/

bool rd_grid_has_lgr(const rd_grid_type *main_grid, const char *__lgr_name) {
    if (!__lgr_name)
        return false;

    __assert_main_grid(main_grid);
    return main_grid->LGR_hash.count(rd::strip_spaces(__lgr_name)) > 0;
}

bool rd_grid_has_lgr_nr(const rd_grid_type *main_grid, int lgr_nr) {
    __assert_main_grid(main_grid);
    return static_cast<int>(main_grid->lgr_index_map.size()) > lgr_nr;
}

int rd_grid_get_num_coarse_groups(const rd_grid_type *main_grid) {
    return main_grid->coarse_cells.size();
}

bool rd_grid_have_coarse_cells(const rd_grid_type *main_grid) {
    return !main_grid->coarse_cells.empty();
}

/**
   Return the number of LGR's associated with this main grid
   instance. The main grid is not counted.
*/
int rd_grid_get_num_lgr(const rd_grid_type *main_grid) {
    __assert_main_grid(main_grid);
    return main_grid->LGR_list.size();
}

/**
   The lgr_index has zero offset
   rd_grid_iget_lgr( rd_grid , 0); will return the first lgr
   rd_grid_iget_lgr( rd_grid , 1); will return the second lgr
   The method will fail HARD if lgr_index is out of bounds.
*/

rd_grid_type *rd_grid_iget_lgr(const rd_grid_type *main_grid, int lgr_index) {
    __assert_main_grid(main_grid);
    return main_grid->LGR_list.at(lgr_index).get();
}

/*
   This function returns the lgr with the given lgr_nr. The lgr_nr is
   the fourth element in the GRIDHEAD for EGRID files.  The lgr nr is
   equal to the grid nr if the grid's are consecutive numbered and
   read from file in increasing lgr nr order.

   This method can only be used for EGRID files. For GRID files the
   lgr_nr is 0 for all grids.
*/

rd_grid_type *rd_grid_get_lgr_from_lgr_nr(const rd_grid_type *main_grid,
                                          int lgr_nr) {
    __assert_main_grid(main_grid);
    {
        int lgr_index = main_grid->lgr_index_map.at(lgr_nr);
        return main_grid->LGR_list.at(lgr_index).get();
    }
}

/**
   The following functions will return the LGR subgrid referenced by
   the coordinates given. Observe the following:

   1. The functions will happily return NULL if no LGR is associated
      with the cell indicated - in fact that is (currently) the only
      way to query whether a particular cell has a LGR.

   2. If a certain cell is refined in several levels this function
      will return a pointer to the first level of refinement. The
      return value can can be used for repeated calls to descend
      deeper into the refinement hierarchy.
*/

const rd_grid_type *rd_grid_get_cell_lgr1(const rd_grid_type *grid,
                                          int global_index) {
    const rd_cell_type &cell = grid->cells.at(global_index);
    return cell.lgr;
}

const char *rd_grid_iget_lgr_name(const rd_grid_type *rd_grid, int lgr_index) {
    __assert_main_grid(rd_grid);
    if (lgr_index >= 0 &&
        static_cast<size_t>(lgr_index) < rd_grid->LGR_list.size()) {
        const auto &lgr = rd_grid->LGR_list[lgr_index];
        return rd_grid_get_name(lgr.get());
    } else
        return NULL;
}

const char *rd_grid_get_lgr_name(const rd_grid_type *rd_grid, int lgr_nr) {
    __assert_main_grid(rd_grid);
    if (lgr_nr == 0)
        return rd_grid_get_name(rd_grid);
    {
        int lgr_index = rd_grid->lgr_index_map.at(lgr_nr);
        return rd_grid_iget_lgr_name(rd_grid, lgr_index);
    }
}

int rd_grid_get_lgr_nr_from_name(const rd_grid_type *grid, const char *name) {
    __assert_main_grid(grid);
    if (strcmp(name, grid->name.c_str()) == 0)
        return 0;
    else {
        const rd_grid_type *lgr = rd_grid_get_lgr(grid, name);
        return lgr->lgr_nr;
    }
}

/**
   This function returns the lgr_nr field of the grid; for GRID files, this
   is just the occurence number in the grid file. Starting with 0 at the main
   grid, and then increasing consecutively through the lgr sections.
   For EGRID files, this is the LGR number (fourth element in the
   gridhead).
*/

int rd_grid_get_lgr_nr(const rd_grid_type *rd_grid) { return rd_grid->lgr_nr; }

const char *rd_grid_get_name(const rd_grid_type *rd_grid) {
    return rd_grid->name.empty() ? NULL : rd_grid->name.c_str();
}

int rd_grid_get_global_size(const rd_grid_type *rd_grid) {
    return rd_grid->nx * rd_grid->ny * rd_grid->nz;
}

int rd_grid_get_active_size(const rd_grid_type *rd_grid) {
    return rd_grid_get_nactive(rd_grid);
}

bool rd_grid_cell_regular1(rd_grid_type *rd_grid, int global_index) {
    double x, y, z;
    rd_grid_get_xyz1(rd_grid, global_index, &x, &y, &z);
    return rd_grid_cell_contains_xyz1(rd_grid, global_index, x, y, z);
}

double rd_grid_get_cell_volume1(const rd_grid_type *rd_grid, int global_index) {
    const rd_cell_type &cell = rd_grid->cells.at(global_index);
    int i, j, k;
    rd_grid_get_ijk1(rd_grid, global_index, &i, &j, &k);
    return rd_cell_get_volume(cell);
}

double rd_grid_get_cell_volume1A(const rd_grid_type *rd_grid,
                                 int active_index) {
    int global_index = rd_grid_get_global_index1A(rd_grid, active_index);
    return rd_grid_get_cell_volume1(rd_grid, global_index);
}

/**

   This function is used to translate (with the help of the rd_grid
   functionality) i,j,k to an index which can be used to look up an
   element in the rd_kw instance. It is just a minor convenience
   function.

   * If the rd_kw instance has nx*ny*nz (i,j,k) are translated to a
     global index with rd_grid_get_global_index3(). This is typically
     the case when the rd_kw instance represents a petrophysical
     property which is e.g. loaded from a INIT file.

   * If the rd_kw instance has nactive elements the (i,j,k) indices
     are converted to an active index with
     rd_grid_get_active_index3(). This is typically the case if the
     rd_kw instance is a solution vector which has been loaded from a
     restart file. If you ask for an inactive cell the function will
     return -1.

   * If the rd_kw instance has neither nx*ny*nz nor nactive elements
     the function will fail HARD.

   * The return value is double, irrespective of the type of the
     underlying datatype of the rd_kw instance - the function will
     fail HARD if the underlying type can not be safely converted to
     double, i.e. if it is not in the set [rd_float_type ,
     rd_int_type , rd_double_type].

   * i,j,k: C-based zero offset grid coordinates.

*/

static int rd_grid_get_property_index__(const rd_grid_type *rd_grid,
                                        const rd_kw_type *rd_kw, int i, int j,
                                        int k) {
    int kw_size = rd_kw_get_size(rd_kw);
    int lookup_index = -1;

    if (kw_size == rd_grid->nx * rd_grid->ny * rd_grid->nz)
        lookup_index = rd_grid_get_global_index3(rd_grid, i, j, k);
    else if (kw_size == rd_grid->total_active)
        /* Will be set to -1 if the cell is not active. */
        lookup_index = rd_grid_get_active_index3(rd_grid, i, j, k);
    else
        throw std::invalid_argument(fmt::format(
            "incommensurable size: keyword has {} elements, expected nx*ny*nz "
            "= {} or nactive = {}",
            kw_size, rd_grid->nx * rd_grid->ny * rd_grid->nz,
            rd_grid->total_active));

    return lookup_index;
}

double rd_grid_get_property(const rd_grid_type *rd_grid,
                            const rd_kw_type *rd_kw, int i, int j, int k) {
    rd_data_type data_type = rd_kw_get_data_type(rd_kw);
    if (rd_type_is_numeric(data_type)) {
        int lookup_index =
            rd_grid_get_property_index__(rd_grid, rd_kw, i, j, k);

        if (lookup_index >= 0)
            return rd_kw_iget_as_double(rd_kw, lookup_index);
        else
            return -1; /* Tried to lookup an inactive cell. */

    } else {
        std::unique_ptr<char> type_name(rd_type_alloc_name(data_type));
        throw std::invalid_argument(fmt::format(
            "rd_grid_get_property: can not lookup type:{} - keyword type must "
            "be numeric.",
            type_name.get()));
    }
}

/**
   Will fill the double_vector instance @column with values from
   rd_kw from the column given by (i,j). If @rd_kw has size nactive
   the inactive k values will not be set, i.e. you should make sure
   that the default value of the @column instance has been properly
   set beforehand.

   The column vector will be filled with double values, the content of
   rd_kw will be converted to double in the case INTE,REAL and DOUB
   types, otherwsie it is crash and burn.
*/

void rd_grid_get_column_property(const rd_grid_type *rd_grid,
                                 const rd_kw_type *rd_kw, int i, int j,
                                 double_vector_type *column) {
    rd_data_type data_type = rd_kw_get_data_type(rd_kw);
    if (rd_type_is_numeric(data_type)) {
        int kw_size = rd_kw_get_size(rd_kw);
        bool use_global_index = false;

        if (kw_size == rd_grid->nx * rd_grid->ny * rd_grid->nz)
            use_global_index = true;
        else if (kw_size == rd_grid->total_active)
            use_global_index = false;
        else
            throw std::invalid_argument(fmt::format(
                "incommensurable sizes: nx*ny*nz = {}  nactive = {}  "
                "kw_size = {}",
                rd_grid->size, rd_grid->total_active, rd_kw_get_size(rd_kw)));

        double_vector_reset(column);
        for (int k = 0; k < rd_grid->nz; k++) {
            if (use_global_index) {
                int global_index = rd_grid_get_global_index3(rd_grid, i, j, k);
                double_vector_iset(column, k,
                                   rd_kw_iget_as_double(rd_kw, global_index));
            } else {
                int active_index = rd_grid_get_active_index3(rd_grid, i, j, k);
                if (active_index >= 0)
                    double_vector_iset(
                        column, k, rd_kw_iget_as_double(rd_kw, active_index));
            }
        }
    } else {
        std::unique_ptr<char> type_name(rd_type_alloc_name(data_type));
        throw std::invalid_argument(fmt::format(
            "rd_grid_get_column_property: can not lookup type:{} - keyword "
            "type must be numeric.",
            type_name.get()));
    }
}

void rd_grid_grdecl_fprintf_kw(const rd_grid_type *rd_grid,
                               const rd_kw_type *rd_kw,
                               const char *special_header, FILE *stream,
                               double double_default) {
    int src_size = rd_kw_get_size(rd_kw);
    if (src_size == rd_grid->size) {
        rd_kw_fprintf_grdecl(rd_kw, stream, special_header);
    } else if (src_size == rd_grid->total_active) {
        void *default_ptr = NULL;
        float float_default;
        int int_default;
        int bool_default;
        rd_type_enum rd_type = rd_type_get_type(rd_kw_get_data_type(rd_kw));

        if (rd_type == RD_FLOAT_TYPE) {
            float_default = (float)double_default;
            default_ptr = &float_default;
        } else if (rd_type == RD_INT_TYPE) {
            int_default = (int)double_default;
            default_ptr = &int_default;
        } else if (rd_type == RD_DOUBLE_TYPE) {
            default_ptr = &double_default;
        } else if (rd_type == RD_BOOL_TYPE) {
            int tmp = (int)double_default;
            if (tmp == 1)
                bool_default = RD_BOOL_TRUE_INT;
            else if (tmp == 0)
                bool_default = RD_BOOL_FALSE_INT;
            else
                throw std::domain_error(fmt::format(
                    "Only 0 and 1 are allowed for bool interpolation, got {}",
                    tmp));
            default_ptr = &bool_default;
        }

        if (default_ptr == NULL) {
            std::unique_ptr<char> type_name(
                rd_type_alloc_name(rd_kw_get_data_type(rd_kw)));
            throw std::invalid_argument(fmt::format(
                "rd_grid_grdecl_fprintf_kw: invalid type {}", type_name.get()));
        }

        {
            auto tmp_kw =
                rd_kw_ptr(rd_kw_alloc_scatter_copy(
                              rd_kw, rd_grid->size,
                              rd_grid->inv_index_map.data(), default_ptr),
                          &rd_kw_free);
            rd_kw_fprintf_grdecl(tmp_kw.get(), stream, special_header);
        }
    } else
        throw std::invalid_argument(fmt::format(
            "size mismatch: rd_kw has {} elements, must have either nx*ny*nz "
            "= {} or nactive = {} elements",
            src_size, rd_grid->size, rd_grid->total_active));
}

static bool rd_grid_test_lgr_consistency2(const rd_grid_type *parent,
                                          const rd_grid_type *child) {
    bool consistent = true;
    int child_size = rd_grid_get_global_size(child);
    for (int child_index = 0; child_index < child_size; child_index++) {
        int parent_cell = rd_grid_get_parent_cell1(child, child_index);
        if (parent_cell >= 0) {
            const rd_grid_type *child_test =
                rd_grid_get_cell_lgr1(parent, parent_cell);
            if (child != child_test) {
                fprintf(stderr,
                        "Child parent mapping failure : ptr difference.\n");
                consistent = false;
            }
        } else {
            fprintf(stderr,
                    "Child parent mapping failure : lgr has no parent cell.\n");
            consistent = false;
        }
    }
    return consistent;
}

bool rd_grid_test_lgr_consistency(const rd_grid_type *rd_grid) {
    bool consistent = true;
    for (const auto &lgr_pair : rd_grid->children) {
        const rd_grid_type *lgr = lgr_pair.second;
        consistent &= rd_grid_test_lgr_consistency2(rd_grid, lgr);
        consistent &= rd_grid_test_lgr_consistency(lgr);
    }
    return consistent;
}

/*
 'MAPUNITS'           1 'CHAR'
 'METRES  '

'GRIDUNIT'           2 'CHAR'
 'METRES  ' '        '

*/

/**
   DIMENS
   MAPUNITS
   MAPAXES
   GRIDUNIT
   RADIAL
   COORDS
   CORNERS
   COORDS
   CORNERS
   ....
   DIMENS
   RADIAL
   COORDS
   CORNERS
   COORDS
   CORNERS
   ...
*/

bool rd_grid_use_mapaxes(const rd_grid_type *grid) { return grid->use_mapaxes; }

void rd_grid_init_mapaxes_data_double(const rd_grid_type *grid,
                                      double *mapaxes) {
    if (grid->mapaxes.has_value())
        for (int i = 0; i < 6; i++)
            mapaxes[i] = grid->mapaxes.value()[i];
}

static const float *rd_grid_get_mapaxes(const rd_grid_type *grid) {
    return grid->mapaxes ? grid->mapaxes->data() : NULL;
}

rd_kw_type *rd_grid_alloc_mapaxes_kw(const rd_grid_type *grid) {
    return rd_kw_alloc_new(MAPAXES_KW, 6, RD_FLOAT, rd_grid_get_mapaxes(grid));
}

static rd_kw_type *rd_grid_alloc_mapunits_kw(ert_rd_unit_enum output_unit) {
    rd_kw_type *mapunits_kw = rd_kw_alloc(MAPUNITS_KW, 1, RD_CHAR);

    if (output_unit == RD_FIELD_UNITS)
        rd_kw_iset_string8(mapunits_kw, 0, "FEET");

    if (output_unit == RD_METRIC_UNITS)
        rd_kw_iset_string8(mapunits_kw, 0, "METRES");

    if (output_unit == RD_LAB_UNITS)
        rd_kw_iset_string8(mapunits_kw, 0, "CM");

    return mapunits_kw;
}

static rd_kw_type *rd_grid_alloc_gridunits_kw(ert_rd_unit_enum output_unit) {
    rd_kw_type *gridunits_kw = rd_kw_alloc(GRIDUNIT_KW, 2, RD_CHAR);

    if (output_unit == RD_FIELD_UNITS)
        rd_kw_iset_string8(gridunits_kw, 0, "FEET");

    if (output_unit == RD_METRIC_UNITS)
        rd_kw_iset_string8(gridunits_kw, 0, "METRES");

    if (output_unit == RD_LAB_UNITS)
        rd_kw_iset_string8(gridunits_kw, 0, "CM");

    rd_kw_iset_string8(gridunits_kw, 1, "");
    return gridunits_kw;
}

static ert_rd_unit_enum
rd_grid_check_unit_system(const rd_kw_type *gridunit_kw) {
    const char *length_unit = rd_kw_iget_char_ptr(gridunit_kw, 0);

    if (strncmp(length_unit, "FEET", 4) == 0)
        return RD_FIELD_UNITS;

    if (strncmp(length_unit, "CM", 2) == 0)
        return RD_LAB_UNITS;

    return RD_METRIC_UNITS;
}

float rd_grid_output_scaling(const rd_grid_type *grid,
                             ert_rd_unit_enum output_unit) {
    if (grid->unit_system == output_unit)
        return 1.0;
    else {
        double scale_factor = 1;

        if (grid->unit_system == RD_FIELD_UNITS)
            scale_factor = 1.0 / METER_TO_FEET_SCALE_FACTOR;

        if (grid->unit_system == RD_LAB_UNITS)
            scale_factor = 1.0 / METER_TO_CM_SCALE_FACTOR;

        if (output_unit == RD_FIELD_UNITS)
            scale_factor *= METER_TO_FEET_SCALE_FACTOR;

        if (output_unit == RD_LAB_UNITS)
            scale_factor *= METER_TO_CM_SCALE_FACTOR;

        return scale_factor;
    }
}

static void rd_grid_fwrite_mapaxes(const rd_grid_type *grid,
                                   fortio_type *fortio) {
    rd_kw_type *mapaxes_kw = rd_grid_alloc_mapaxes_kw(grid);
    rd_kw_fwrite(mapaxes_kw, fortio);
    rd_kw_free(mapaxes_kw);
}

static void rd_grid_fwrite_mapunits(fortio_type *fortio,
                                    ert_rd_unit_enum output_unit) {
    rd_kw_type *mapunits_kw = rd_grid_alloc_mapunits_kw(output_unit);
    rd_kw_fwrite(mapunits_kw, fortio);
    rd_kw_free(mapunits_kw);
}

static void rd_grid_fwrite_gridunits(fortio_type *fortio,
                                     ert_rd_unit_enum output_unit) {
    rd_kw_type *gridunits_kw = rd_grid_alloc_gridunits_kw(output_unit);
    rd_kw_fwrite(gridunits_kw, fortio);
    rd_kw_free(gridunits_kw);
}

static void rd_grid_fwrite_main_GRID_headers(const rd_grid_type *rd_grid,
                                             fortio_type *fortio,
                                             ert_rd_unit_enum output_unit) {
    rd_grid_fwrite_mapunits(fortio, output_unit);

    if (rd_grid->use_mapaxes)
        rd_grid_fwrite_mapaxes(rd_grid, fortio);

    rd_grid_fwrite_gridunits(fortio, output_unit);
}

static void rd_grid_fwrite_GRID__(const rd_grid_type *grid, int coords_size,
                                  ERT::FortIO &fortio,
                                  ert_rd_unit_enum output_unit) {
    if (grid->parent_grid != NULL) {
        auto lgr_kw = make_rd_kw(LGR_KW, 1, RD_CHAR);
        rd_kw_iset_string8(lgr_kw.get(), 0, grid->name.c_str());
        rd_kw_fwrite(lgr_kw.get(), fortio.get());
    }

    {
        auto dimens_kw = make_rd_kw(DIMENS_KW, 3, RD_INT);
        rd_kw_iset_int(dimens_kw.get(), 0, grid->nx);
        rd_kw_iset_int(dimens_kw.get(), 1, grid->ny);
        if (grid->dualp_flag == FILEHEAD_SINGLE_POROSITY)
            rd_kw_iset_int(dimens_kw.get(), 2, grid->nz);
        else
            rd_kw_iset_int(dimens_kw.get(), 2, 2 * grid->nz);

        rd_kw_fwrite(dimens_kw.get(), fortio.get());
    }

    if (grid->parent_grid == NULL)
        rd_grid_fwrite_main_GRID_headers(grid, fortio.get(), output_unit);

    {
        auto radial_kw = make_rd_kw(RADIAL_KW, 1, RD_CHAR);
        rd_kw_iset_string8(radial_kw.get(), 0, "FALSE");
        rd_kw_fwrite(radial_kw.get(), fortio.get());
    }

    auto coords_kw = make_rd_kw(COORDS_KW, coords_size, RD_INT);
    auto corners_kw = make_rd_kw(CORNERS_KW, 24, RD_FLOAT);
    for (int k = 0; k < grid->nz; k++) {
        for (int j = 0; j < grid->ny; j++) {
            for (int i = 0; i < grid->nx; i++) {
                int global_index = rd_grid_get_global_index__(grid, i, j, k);
                const rd_cell_type &cell = grid->cells.at(global_index);

                rd_cell_fwrite_GRID(grid, cell, false, coords_size, i, j, k,
                                    global_index, coords_kw, corners_kw,
                                    fortio);
            }
        }
    }

    if (grid->dualp_flag != FILEHEAD_SINGLE_POROSITY) {
        for (int k = grid->nz; k < 2 * grid->nz; k++) {
            for (int j = 0; j < grid->ny; j++) {
                for (int i = 0; i < grid->nx; i++) {
                    int global_index =
                        rd_grid_get_global_index__(grid, i, j, k - grid->nz);
                    const rd_cell_type &cell = grid->cells.at(global_index);

                    rd_cell_fwrite_GRID(grid, cell, true, coords_size, i, j, k,
                                        global_index, coords_kw, corners_kw,
                                        fortio);
                }
            }
        }
    }
}

void rd_grid_fwrite_GRID2(const rd_grid_type *grid, const char *filename,
                          ert_rd_unit_enum output_unit) {
    int coords_size = 5;
    bool fmt_file = false;

    ERT::FortIO fortio(filename, std::ios_base::out, fmt_file, RD_ENDIAN_FLIP);
    if (grid->children.size() > 0)
        coords_size = 7;

    if (grid->coarsening_active)
        coords_size = 7;

    rd_grid_fwrite_GRID__(grid, coords_size, fortio, output_unit);

    for (const auto &igrid : grid->LGR_list) {
        rd_grid_fwrite_GRID__(igrid.get(), coords_size, fortio, output_unit);
    }
}

/*
FILEHEAD          100:INTE
MAPUNITS            1:CHAR
MAPAXES             6:REAL
GRIDUNIT            2:CHAR
GRIDHEAD          100:INTE
COORD           15990:REAL
ZCORN          286720:REAL
ACTNUM          35840:INTE
ENDGRID             0:INTE
*/

/**
   This function does by construction not take a grid instance as
   argument, so that it can be called from an external scope to create
   a standard EGRID header without creating a grid instance first.
*/

static void rd_grid_fwrite_main_EGRID_header(const rd_grid_type *grid,
                                             ERT::FortIO &fortio,
                                             ert_rd_unit_enum output_unit) {
    int EGRID_VERSION = 3;
    int RELEASE_YEAR = 2007;
    int COMPAT_VERSION = 0;
    const float *mapaxes = rd_grid_get_mapaxes(grid);

    {
        auto filehead_kw = make_rd_kw(FILEHEAD_KW, 100, RD_INT);
        rd_kw_scalar_set_int(filehead_kw.get(), 0);

        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_VERSION_INDEX,
                       EGRID_VERSION);
        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_YEAR_INDEX, RELEASE_YEAR);
        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_COMPAT_INDEX,
                       COMPAT_VERSION);
        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_TYPE_INDEX,
                       FILEHEAD_GRIDTYPE_CORNERPOINT);
        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_DUALP_INDEX,
                       grid->dualp_flag);
        rd_kw_iset_int(filehead_kw.get(), FILEHEAD_ORGFORMAT_INDEX,
                       FILEHEAD_ORGTYPE_CORNERPOINT);

        rd_kw_fwrite(filehead_kw.get(), fortio.get());
    }

    rd_grid_fwrite_mapunits(fortio.get(), output_unit);
    if (mapaxes != NULL)
        rd_grid_fwrite_mapaxes(grid, fortio.get());

    rd_grid_fwrite_gridunits(fortio.get(), output_unit);
}

static void rd_grid_fwrite_gridhead_kw(int nx, int ny, int nz, int grid_nr,
                                       fortio_type *fortio) {
    rd_kw_type *gridhead_kw = rd_grid_alloc_gridhead_kw(nx, ny, nz, grid_nr);
    rd_kw_fwrite(gridhead_kw, fortio);
    rd_kw_free(gridhead_kw);
}

/*
  Will scan the halfopen k-interval [k1,k2) to find a cell which has
  valid geometry. If no cell is found the function will return -1.
*/

static int rd_grid_get_valid_index(const rd_grid_type *grid, int i, int j,
                                   int k1, int k2) {
    int global_index;
    int k = k1;
    int delta = (k1 < k2) ? 1 : -1;

    while (true) {
        global_index = rd_grid_get_global_index3(grid, i, j, k);

        const rd_cell_type &cell = grid->cells.at(global_index);
        if (GET_CELL_FLAG((&cell), CELL_FLAG_VALID))
            return global_index;
        else {
            k += delta;
            if (k == k2)
                return -1;
        }
    }
}

static int rd_grid_get_bottom_valid_index(const rd_grid_type *grid, int i,
                                          int j) {
    return rd_grid_get_valid_index(grid, i, j, grid->nz - 1, -1);
}

static int rd_grid_get_top_valid_index(const rd_grid_type *grid, int i, int j) {
    return rd_grid_get_valid_index(grid, i, j, 0, grid->nz);
}

static bool rd_grid_init_coord_section__(const rd_grid_type *grid, int i, int j,
                                         int i_corner, int j_corner,
                                         bool force_set, float *coord_float,
                                         double *coord_double) {

    const int top_index = rd_grid_get_top_valid_index(grid, i, j);
    const int bottom_index = rd_grid_get_bottom_valid_index(grid, i, j);

    if (top_index == -1)
        util_exit("% : no cell with a valid geometry description found in "
                  "(i,j) = %d,%d - then what? \n",
                  __func__, i, j);

    {
        point_type top_point;
        point_type bottom_point;

        const rd_cell_type &bottom_cell = grid->cells.at(bottom_index);
        const rd_cell_type &top_cell = grid->cells.at(top_index);

        /*
      2---3
      |   |
      0---1
    */
        int corner_index = j_corner * 2 + i_corner;
        int coord_offset =
            6 * ((j + j_corner) * (grid->nx + 1) + (i + i_corner));
        {
            point_copy_values(&top_point, &top_cell.corner_list[corner_index]);
            point_copy_values(&bottom_point,
                              &bottom_cell.corner_list[corner_index + 4]);

            if ((top_point.z == bottom_point.z) && (force_set == false)) {
                return false;
            } else {
                if (grid->use_mapaxes) {
                    point_mapaxes_invtransform(&top_point, grid->origo,
                                               grid->unit_x, grid->unit_y);
                    point_mapaxes_invtransform(&bottom_point, grid->origo,
                                               grid->unit_x, grid->unit_y);
                }

                if (coord_float) {
                    coord_float[coord_offset] = top_point.x;

                    coord_float[coord_offset + 1] = top_point.y;
                    coord_float[coord_offset + 2] = top_point.z;
                    coord_float[coord_offset + 3] = bottom_point.x;
                    coord_float[coord_offset + 4] = bottom_point.y;
                    coord_float[coord_offset + 5] = bottom_point.z;
                }

                if (coord_double) {
                    coord_double[coord_offset] = top_point.x;
                    coord_double[coord_offset + 1] = top_point.y;
                    coord_double[coord_offset + 2] = top_point.z;

                    coord_double[coord_offset + 3] = bottom_point.x;
                    coord_double[coord_offset + 4] = bottom_point.y;
                    coord_double[coord_offset + 5] = bottom_point.z;
                }

                return true;
            }
        }
    }
}

static void rd_grid_init_coord_section(const rd_grid_type *grid, int i, int j,
                                       float *coord_float,
                                       double *coord_double) {
    int i_corner = 0;
    int j_corner = 0;

    if (i == grid->nx) {
        i -= 1;
        i_corner = 1;
    }

    if (j == grid->ny) {
        j -= 1;
        j_corner = 1;
    }

    rd_grid_init_coord_section__(grid, i, j, i_corner, j_corner,
                                 /*force_set=*/true, coord_float, coord_double);
}

void rd_grid_init_coord_data(const rd_grid_type *grid, float *coord) {
    /*
    The coord vector contains the points defining the top and bottom
    of the pillars. The vector contains (nx + 1) * (ny + 1) 6 element
    chunks of data, where each chunk contains the coordinates (x,y,z)
    of the top and the bottom of the pillar.
  */
    for (int j = 0; j <= grid->ny; j++) {
        for (int i = 0; i <= grid->nx; i++)
            rd_grid_init_coord_section(grid, i, j, coord, NULL);
    }
}

void rd_grid_init_coord_data_double(const rd_grid_type *grid, double *coord) {
    /*
    The coord vector contains the points defining the top and bottom
    of the pillars. The vector contains (nx + 1) * (ny + 1) 6 element
    chunks of data, where each chunk contains the coordinates (x,y,z)
    f the top and the bottom of the pillar.
  */
    for (int j = 0; j <= grid->ny; j++) {
        for (int i = 0; i <= grid->nx; i++)
            rd_grid_init_coord_section(grid, i, j, NULL, coord);
    }
}

static void rd_grid_assert_coord_kw(rd_grid_type *grid) {
    if (!grid->coord_kw) {
        grid->coord_kw.reset(
            rd_kw_alloc(COORD_KW, rd_grid_get_coord_size(grid), RD_FLOAT));
        rd_grid_init_coord_data(grid, static_cast<float *>(rd_kw_get_void_ptr(
                                          grid->coord_kw.get())));
    }
}

static void rd_grid_init_zcorn_data__(const rd_grid_type *grid,
                                      float *zcorn_float,
                                      double *zcorn_double) {
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                const int cell_index = rd_grid_get_global_index3(grid, i, j, k);
                const rd_cell_type &cell = grid->cells.at(cell_index);

                for (int l = 0; l < 2; l++) {
                    point_type p0 = cell.corner_list[4 * l];
                    point_type p1 = cell.corner_list[4 * l + 1];
                    point_type p2 = cell.corner_list[4 * l + 2];
                    point_type p3 = cell.corner_list[4 * l + 3];

                    int z1 =
                        k * 8 * nx * ny + j * 4 * nx + 2 * i + l * 4 * nx * ny;
                    int z2 = k * 8 * nx * ny + j * 4 * nx + 2 * i + 1 +
                             l * 4 * nx * ny;
                    int z3 = k * 8 * nx * ny + j * 4 * nx + 2 * nx + 2 * i +
                             l * 4 * nx * ny;
                    int z4 = k * 8 * nx * ny + j * 4 * nx + 2 * nx + 2 * i + 1 +
                             l * 4 * nx * ny;

                    if (zcorn_float) {
                        zcorn_float[z1] = p0.z;
                        zcorn_float[z2] = p1.z;
                        zcorn_float[z3] = p2.z;
                        zcorn_float[z4] = p3.z;
                    }

                    if (zcorn_double) {
                        zcorn_double[z1] = p0.z;
                        zcorn_double[z2] = p1.z;
                        zcorn_double[z3] = p2.z;
                        zcorn_double[z4] = p3.z;
                    }
                }
            }
        }
    }
}

void rd_grid_init_zcorn_data(const rd_grid_type *grid, float *zcorn) {
    rd_grid_init_zcorn_data__(grid, zcorn, NULL);
}

void rd_grid_init_zcorn_data_double(const rd_grid_type *grid, double *zcorn) {
    rd_grid_init_zcorn_data__(grid, NULL, zcorn);
}

rd_kw_type *rd_grid_alloc_zcorn_kw(const rd_grid_type *grid) {
    rd_kw_type *zcorn_kw =
        rd_kw_alloc(ZCORN_KW, rd_grid_get_zcorn_size(grid), RD_FLOAT);
    rd_grid_init_zcorn_data(grid, (float *)rd_kw_get_void_ptr(zcorn_kw));
    return zcorn_kw;
}

rd_kw_type *rd_grid_alloc_coord_kw(const rd_grid_type *grid) {
    if (grid->coord_kw)
        return rd_kw_alloc_copy(grid->coord_kw.get());

    rd_kw_type *coord_kw =
        rd_kw_alloc(COORD_KW, RD_GRID_COORD_SIZE(grid->nx, grid->ny), RD_FLOAT);
    rd_grid_init_coord_data(grid, rd_kw_get_float_ptr(coord_kw));

    return coord_kw;
}

int rd_grid_get_coord_size(const rd_grid_type *grid) {
    return RD_GRID_COORD_SIZE(grid->nx, grid->ny);
}

int rd_grid_get_zcorn_size(const rd_grid_type *grid) {
    return RD_GRID_ZCORN_SIZE(grid->nx, grid->ny, grid->nz);
}

void rd_grid_init_actnum_data(const rd_grid_type *grid, int *actnum) {
    for (int i = 0; i < grid->size; i++) {
        const rd_cell_type &cell = grid->cells.at(i);
        if (cell.coarse_group == COARSE_GROUP_NONE)
            actnum[i] = cell.active;
        else {
            /* In the case of coarse cells we must query the coarse cell for
         the original, uncoarsened distribution of actnum values. */
            rd_coarse_cell_type *coarse_cell =
                rd_grid_iget_coarse_group(grid, cell.coarse_group);

            /* 1: Set all the elements in the coarse group to inactive. */
            {
                int group_size = rd_coarse_cell_get_size(coarse_cell);
                const int *index_ptr =
                    rd_coarse_cell_get_index_ptr(coarse_cell);
                int j;
                for (j = 0; j < group_size; j++)
                    actnum[index_ptr[j]] = CELL_NOT_ACTIVE;
            }

            /* 2: Explicitly pick up the active cells from the coarse group
            and mark them correctly. */
            {
                int num_active = rd_coarse_cell_get_num_active(coarse_cell);
                int j;

                for (j = 0; j < num_active; j++) {
                    int global_index =
                        rd_coarse_cell_iget_active_cell_index(coarse_cell, j);
                    int active_value =
                        rd_coarse_cell_iget_active_value(coarse_cell, j);

                    actnum[global_index] = active_value;
                }
            }
        }
    }
}

rd_kw_type *rd_grid_alloc_actnum_kw(const rd_grid_type *grid) {
    rd_kw_type *actnum_kw = rd_kw_alloc(ACTNUM_KW, grid->size, RD_INT);
    rd_grid_init_actnum_data(grid, (int *)rd_kw_get_void_ptr(actnum_kw));
    return actnum_kw;
}

void rd_grid_compressed_kw_copy(const rd_grid_type *grid, rd_kw_type *target_kw,
                                const rd_kw_type *src_kw) {
    if ((rd_kw_get_size(target_kw) == rd_grid_get_nactive(grid)) &&
        (rd_kw_get_size(src_kw) == rd_grid_get_global_size(grid))) {
        int active_index = 0;
        for (int i = 0; i < rd_grid_get_global_size(grid); i++) {
            if (rd_grid_cell_active1(grid, i)) {
                rd_kw_iset(target_kw, active_index, rd_kw_iget_ptr(src_kw, i));
                active_index++;
            }
        }
    } else
        throw std::invalid_argument(fmt::format(
            "rd_grid_compressed_kw_copy: size mismatch target:{}  src:{}  "
            "expected {},{}",
            rd_kw_get_size(target_kw), rd_kw_get_size(src_kw),
            rd_grid_get_nactive(grid), rd_grid_get_global_size(grid)));
}

void rd_grid_global_kw_copy(const rd_grid_type *grid, rd_kw_type *target_kw,
                            const rd_kw_type *src_kw) {
    if ((rd_kw_get_size(src_kw) == rd_grid_get_nactive(grid)) &&
        (rd_kw_get_size(target_kw) == rd_grid_get_global_size(grid))) {
        int active_index = 0;
        for (int i = 0; i < rd_grid_get_global_size(grid); i++) {
            if (rd_grid_cell_active1(grid, i)) {
                rd_kw_iset(target_kw, i, rd_kw_iget_ptr(src_kw, active_index));
                active_index++;
            }
        }
    } else
        throw std::invalid_argument(fmt::format(
            "rd_grid_global_kw_copy: size mismatch target:{}  src:{}  "
            "expected {},{}",
            rd_kw_get_size(target_kw), rd_kw_get_size(src_kw),
            rd_grid_get_global_size(grid), rd_grid_get_nactive(grid)));
}

static void rd_grid_init_hostnum_data(const rd_grid_type *grid, int *hostnum) {
    for (int i = 0; i < grid->size; i++) {
        const rd_cell_type &cell = grid->cells.at(i);
        hostnum[i] = cell.host_cell + 1;
    }
}

static rd_kw_type *rd_grid_alloc_hostnum_kw(const rd_grid_type *grid) {
    rd_kw_type *hostnum_kw = rd_kw_alloc(HOSTNUM_KW, grid->size, RD_INT);
    rd_grid_init_hostnum_data(grid, (int *)rd_kw_get_void_ptr(hostnum_kw));
    return hostnum_kw;
}

static void rd_grid_init_corsnum_data(const rd_grid_type *grid, int *corsnum) {
    for (int i = 0; i < grid->size; i++) {
        const rd_cell_type &cell = grid->cells.at(i);
        corsnum[i] = cell.coarse_group + 1;
    }
}

static rd_kw_type *rd_grid_alloc_corsnum_kw(const rd_grid_type *grid) {
    rd_kw_type *corsnum_kw = rd_kw_alloc(CORSNUM_KW, grid->size, RD_INT);
    rd_grid_init_corsnum_data(grid, (int *)rd_kw_get_void_ptr(corsnum_kw));
    return corsnum_kw;
}

void rd_grid_reset_actnum(rd_grid_type *grid, const int *actnum) {
    const int global_size = rd_grid_get_global_size(grid);
    for (int g = 0; g < global_size; g++) {
        rd_cell_type &cell = grid->cells.at(g);
        if (actnum)
            cell.active = actnum[g];
        else
            cell.active = 1;
    }
    rd_grid_update_index(grid);
}

static void rd_grid_fwrite_self_nnc(const rd_grid_type *grid,
                                    const ERT::FortIO &fortio) {
    const int default_index = 1;
    std::vector<int> g1(0, default_index);
    std::vector<int> g2(0, default_index);
    for (int g = 0; g < rd_grid_get_global_size(grid); g++) {
        const rd_cell_type &cell = grid->cells.at(g);
        auto &nnc_info = cell.nnc_info;
        if (nnc_info) {
            const nnc_vector_type *nnc_vector =
                nnc_info_get_self_vector(nnc_info.get());
            if (!nnc_vector)
                continue;
            for (int i = 0; i < nnc_vector_get_size(nnc_vector); i++) {
                int nnc_index = nnc_vector_iget_nnc_index(nnc_vector, i);
                if (nnc_index < 0)
                    throw std::domain_error(
                        "rd_grid_fwrite_self_nnc: negative nnc_index");

                size_t u_nnc_index = static_cast<size_t>(nnc_index);

                if (u_nnc_index >= g1.size()) {
                    g1.resize(u_nnc_index + 1, default_index);
                    g2.resize(u_nnc_index + 1, default_index);
                }
                g1[nnc_index] = 1 + g;
                g2[nnc_index] = 1 + nnc_vector_iget_grid_index(nnc_vector, i);
            }
        }
    }
    if (g1.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(
            "rd_grid_fwrite_self_nnc: NNC count exceeds INT_MAX");
    int num_nnc = static_cast<int>(g1.size());

    // ensure that g1.data() != nullptr
    if (g1.empty()) {
        g1.resize(1, default_index);
        g2.resize(1, default_index);
    }

    auto nnc1_kw =
        rd_kw_ptr(rd_kw_alloc_new_shared(NNC1_KW, num_nnc, RD_INT, g1.data()),
                  &rd_kw_free);
    auto nnc2_kw =
        rd_kw_ptr(rd_kw_alloc_new_shared(NNC2_KW, num_nnc, RD_INT, g2.data()),
                  &rd_kw_free);
    auto nnchead_kw = make_rd_kw(NNCHEAD_KW, NNCHEAD_SIZE, RD_INT);

    rd_kw_scalar_set_int(nnchead_kw.get(), 0);
    rd_kw_iset_int(nnchead_kw.get(), NNCHEAD_NUMNNC_INDEX, num_nnc);
    rd_kw_iset_int(nnchead_kw.get(), NNCHEAD_LGR_INDEX, grid->lgr_nr);

    rd_kw_fwrite(nnchead_kw.get(), fortio.get());
    rd_kw_fwrite(nnc1_kw.get(), fortio.get());
    rd_kw_fwrite(nnc2_kw.get(), fortio.get());
}

static void rd_grid_fwrite_EGRID__(rd_grid_type *grid, ERT::FortIO &fortio,
                                   ert_rd_unit_enum output_unit) {
    bool is_lgr = true;
    if (grid->parent_grid == NULL)
        is_lgr = false;

    /* Writing header */
    if (!is_lgr) {
        rd_grid_fwrite_main_EGRID_header(grid, fortio, output_unit);
    } else {
        {
            auto lgr_kw = make_rd_kw(LGR_KW, 1, RD_CHAR);
            rd_kw_iset_string8(lgr_kw.get(), 0, grid->name.c_str());
            rd_kw_fwrite(lgr_kw.get(), fortio.get());
        }

        {
            auto lgr_parent_kw = make_rd_kw(LGR_PARENT_KW, 1, RD_CHAR);
            if (grid->parent_name.has_value())
                rd_kw_iset_string8(lgr_parent_kw.get(), 0,
                                   grid->parent_name->c_str());
            else
                rd_kw_iset_string8(lgr_parent_kw.get(), 0, "");

            rd_kw_fwrite(lgr_parent_kw.get(), fortio.get());
        }
    }

    rd_grid_fwrite_gridhead_kw(grid->nx, grid->ny, grid->nz, grid->lgr_nr,
                               fortio.get());
    /* Writing main grid data */
    {
        rd_grid_assert_coord_kw(grid);
        {
            auto coord_kw =
                rd_kw_ptr(rd_kw_alloc_copy(grid->coord_kw.get()), &rd_kw_free);
            auto zcorn_kw =
                rd_kw_ptr(rd_grid_alloc_zcorn_kw(grid), &rd_kw_free);

            if (output_unit != grid->unit_system) {
                double scale_factor = rd_grid_output_scaling(grid, output_unit);
                rd_kw_scale_float(coord_kw.get(), scale_factor);
                rd_kw_scale_float(zcorn_kw.get(), scale_factor);
            }
            rd_kw_fwrite(coord_kw.get(), fortio.get());
            rd_kw_fwrite(zcorn_kw.get(), fortio.get());
        }
        {
            auto actnum_kw =
                rd_kw_ptr(rd_grid_alloc_actnum_kw(grid), &rd_kw_free);
            rd_kw_fwrite(actnum_kw.get(), fortio.get());
        }
        if (is_lgr) {
            auto hostnum_kw =
                rd_kw_ptr(rd_grid_alloc_hostnum_kw(grid), &rd_kw_free);
            rd_kw_fwrite(hostnum_kw.get(), fortio.get());
        }
        if (grid->coarsening_active) {
            auto corsnum_kw =
                rd_kw_ptr(rd_grid_alloc_corsnum_kw(grid), &rd_kw_free);
            rd_kw_fwrite(corsnum_kw.get(), fortio.get());
        }

        {
            auto endgrid_kw = make_rd_kw(ENDGRID_KW, 0, RD_INT);
            rd_kw_fwrite(endgrid_kw.get(), fortio.get());
        }
    }

    if (is_lgr) {
        auto endlgr_kw = make_rd_kw(ENDLGR_KW, 0, RD_INT);
        rd_kw_fwrite(endlgr_kw.get(), fortio.get());
    }
    rd_grid_fwrite_self_nnc(grid, fortio);
}

void rd_grid_fwrite_EGRID2(rd_grid_type *grid, const char *filename,
                           ert_rd_unit_enum output_unit) {
    bool fmt_file = false;
    {
        bool is_fmt;

        if (rd_get_file_type(filename, &is_fmt, NULL) != RD_OTHER_FILE)
            fmt_file = is_fmt;
    }
    ERT::FortIO fortio(filename, std::ios_base::out, fmt_file, RD_ENDIAN_FLIP);

    rd_grid_fwrite_EGRID__(grid, fortio, output_unit);
    for (const auto &igrid : grid->LGR_list) {
        rd_grid_fwrite_EGRID__(igrid.get(), fortio, output_unit);
    }
}

/*
   The construction with rd_grid_fwrite_EGRID() and
   rd_grid_fwrite_EGRID2() is an attempt to create API stability. New
   code should use the rd_grid_fwrite_EGRID2() function.
*/

void rd_grid_fwrite_EGRID(rd_grid_type *grid, const char *filename,
                          bool output_metric) {
    ert_rd_unit_enum output_unit = RD_METRIC_UNITS;

    if (!output_metric)
        output_unit = RD_FIELD_UNITS;

    rd_grid_fwrite_EGRID2(grid, filename, output_unit);
}

/**
   Writes the current grid as grdecl keywords. This function will only write the main grid and not
   possible LGRs which are attached.
*/

void rd_grid_fprintf_grdecl2(rd_grid_type *grid, FILE *stream,
                             ert_rd_unit_enum output_unit) {
    {
        rd_kw_type *mapunits_kw = rd_grid_alloc_mapunits_kw(output_unit);
        rd_kw_fprintf_grdecl(mapunits_kw, stream);
        rd_kw_free(mapunits_kw);
        fprintf(stream, "\n");
    }

    if (grid->mapaxes.has_value()) {
        rd_kw_type *mapaxes_kw = rd_grid_alloc_mapaxes_kw(grid);
        rd_kw_fprintf_grdecl(mapaxes_kw, stream);
        rd_kw_free(mapaxes_kw);
    }

    {
        rd_kw_type *gridunits_kw = rd_grid_alloc_gridunits_kw(output_unit);
        rd_kw_fprintf_grdecl(gridunits_kw, stream);
        rd_kw_free(gridunits_kw);
        fprintf(stream, "\n");
    }

    /*
     The specgrid header is not properly internalized; the fourth and
     fifth elements are just set to hardcoded values.
  */
    {
        int numres = 1;
        char coord_type = 'F';
        fprintf(stream, "%s\n", SPECGRID_KW);
        fprintf(stream, "  %d  %d  %d  %d  %c / \n\n", grid->nx, grid->ny,
                grid->nz, numres, coord_type);
    }

    {
        rd_grid_assert_coord_kw(grid);
        rd_kw_fprintf_grdecl(grid->coord_kw.get(), stream);
        fprintf(stream, "\n");
    }

    {
        rd_kw_type *zcorn_kw = rd_grid_alloc_zcorn_kw(grid);
        rd_kw_fprintf_grdecl(zcorn_kw, stream);
        rd_kw_free(zcorn_kw);
        fprintf(stream, "\n");
    }

    {
        rd_kw_type *actnum_kw = rd_grid_alloc_actnum_kw(grid);
        rd_kw_fprintf_grdecl(actnum_kw, stream);
        rd_kw_free(actnum_kw);
        fprintf(stream, "\n");
    }
}

bool rd_grid_dual_grid(const rd_grid_type *rd_grid) {
    if (rd_grid->dualp_flag == FILEHEAD_SINGLE_POROSITY)
        return false;
    else
        return true;
}

static rd_kw_type *rd_grid_alloc_volume_kw_active(const rd_grid_type *grid) {
    rd_kw_type *volume_kw =
        rd_kw_alloc("VOLUME", rd_grid_get_active_size(grid), RD_DOUBLE);
    {
        double *volume_data = (double *)rd_kw_get_ptr(volume_kw);
        int active_index;
        for (active_index = 0; active_index < rd_grid_get_active_size(grid);
             active_index++) {
            double cell_volume = rd_grid_get_cell_volume1A(grid, active_index);
            volume_data[active_index] = cell_volume;
        }
    }
    return volume_kw;
}

static rd_kw_type *rd_grid_alloc_volume_kw_global(const rd_grid_type *grid) {
    rd_kw_type *volume_kw =
        rd_kw_alloc("VOLUME", rd_grid_get_global_size(grid), RD_DOUBLE);
    {
        double *volume_data = (double *)rd_kw_get_ptr(volume_kw);
        int global_index;
        for (global_index = 0; global_index < rd_grid_get_global_size(grid);
             global_index++) {
            double cell_volume = rd_grid_get_cell_volume1(grid, global_index);
            volume_data[global_index] = cell_volume;
        }
    }
    return volume_kw;
}

rd_kw_type *rd_grid_alloc_volume_kw(const rd_grid_type *grid,
                                    bool active_size) {
    if (active_size)
        return rd_grid_alloc_volume_kw_active(grid);
    else
        return rd_grid_alloc_volume_kw_global(grid);
}

//This function is meant to be used w/ pandas datafram and numpy
//Note: global_index must be allocated w/ rd_grid->total_active or rd_grid->global_size int32 data points
//      index_data must be allocated w/ (4 * rd_grid->total_active or rd_grid->global_size) int32 data points.
void rd_grid_export_index(const rd_grid_type *grid, int *global_index,
                          int *index_data, bool active_only) {
    int pos_indx = 0;
    int pos_data = 0;
    for (int k = 0; k < grid->nz; k++)
        for (int j = 0; j < grid->ny; j++)
            for (int i = 0; i < grid->nx; i++) {
                int g = rd_grid_get_global_index__(grid, i, j, k);
                if (!active_only || grid->cells[g].active_index[0] >= 0) {
                    global_index[pos_indx++] = g;
                    index_data[pos_data++] = i;
                    index_data[pos_data++] = j;
                    index_data[pos_data++] = k;
                    index_data[pos_data++] = grid->cells[g].active_index[0];
                }
            }
}

rd_grid_ptr make_rectangular_grid(int nx, int ny, int nz, double dx, double dy,
                                  double dz, const int *actnum) {
    return {rd_grid_alloc_rectangular(nx, ny, nz, dx, dy, dz, actnum),
            &rd_grid_free};
}

rd_grid_ptr read_grid(const std::filesystem::path &filename) {
    return {rd_grid_alloc(filename.string().c_str()), &rd_grid_free};
}
