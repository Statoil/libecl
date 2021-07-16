#pragma once

#include <map>
#include <string>
#include <vector>
#include <memory>
#include <ert/ecl/enums.h>

#define ECL_SMSPEC_ID 806647

namespace ecl {
class smspec_node;
}

typedef std::map<std::string, const ecl::smspec_node *> node_map;

typedef struct ecl_smspec_struct {
    int __type_id;
    /*
    All the hash tables listed below here are different ways to access
    smspec_node instances. The actual smspec_node instances are
    owned by the smspec_nodes vector;
  */
    node_map field_var_index;
    node_map misc_var_index; /* Variables like 'TCPU' and 'NEWTON'. */
    node_map
        gen_var_index /* This is "everything" - things can either be found as gen_var("WWCT:OP_X") or as well_var("WWCT" , "OP_X") */
        ;

    std::map<std::string, node_map>
        well_var_index; /* Indexes for all well variables:
                                                     {well1: {var1: index1 , var2: index2} , well2: {var1: index1 , var2: index2}} */
    std::map<std::string, node_map>
        group_var_index; /* Indexes for group variables.*/
    std::map<int, node_map>
        region_var_index; /* The stored index is an offset. */
    std::map<int, node_map> block_var_index; /* Block variables like BPR */
    std::map<std::string, std::map<int, node_map>>
        well_completion_var_index; /* Indexes for completion indexes .*/

    std::vector<std::unique_ptr<ecl::smspec_node>> smspec_nodes;
    bool write_mode;
    bool need_nums;
    std::vector<int> index_map;
    std::map<int, int> inv_index_map;
    int params_size;

    int time_seconds;
    int grid_dims[3]; /* Grid dimensions - in DIMENS[1,2,3] */
    int num_regions;
    int Nwells, param_offset;
    std::string
        key_join_string; /* The string used to join keys when building gen_key keys - typically ":" -
                                                      but arbitrary - NOT necessary to be able to invert the joining. */
    std::string
        header_file; /* FULL path to the currenbtly loaded header_file. */

    bool
        formatted; /* Has this summary instance been loaded from a formatted (i.e. FSMSPEC file) or unformatted (i.e. SMSPEC) file. */
    time_t sim_start_time; /* When did the simulation start - worldtime. */

    int time_index; /* The fields time_index, day_index, month_index and year_index */
    int day_index; /* are used by the ecl_sum_data object to locate per. timestep */
    int month_index; /* time information. */
    int year_index;
    bool has_lgr;
    std::vector<float> params_default;

    std::string restart_case;
    ert_ecl_unit_enum unit_system;
    int restart_step;
} ecl_smspec_type;
