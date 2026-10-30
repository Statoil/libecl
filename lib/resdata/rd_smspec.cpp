#include <cstring>
#include <cmath>
#include <ctime>
#include <cerrno>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <filesystem>

#include <ert/util/hash.hpp>
#include <ert/util/util.hpp>
#include <ert/util/float_vector.hpp>
#include <ert/util/stringlist.hpp>
#include "detail/util/path.hpp"

#include <fmt/format.h>

#include <resdata/rd_smspec.hpp>
#include <resdata/rd_file.hpp>
#include <resdata/rd_kw_magic.hpp>
#include <resdata/rd_kw.hpp>
#include <resdata/rd_util.hpp>
#include <resdata/smspec_node.hpp>
#include <resdata/rd_endian_flip.hpp>
#include <resdata/rd_type.hpp>

#ifdef HAVE_FNMATCH
#include <fnmatch.h>
#endif

namespace fs = std::filesystem;

/**
   This file implements the indexing into the summary files.
*/

/*
  Supporting a new variable type:
  -------------------------------

  1. The function smspec_node_alloc() must be updated to return a valid
     rd::smspec_node instance when called with the new var_type.

  2. Update the function rd_smpec_install_gen_key() to install smpec_index
     instances of this particular type. The format of the general key is
     implicitly defined in this function.

  3. The rd_smspec structure supports two different types of lookup:

       a) Lookup based on general key like e.g. WWCT:OP6
       b) Specific lookup based on the variable type, i.e. :

              rd_smspec_get_well_var( xxx , well_name , var).

      Historically everything started with specific lookup as in case b);
      however the general lookup proved to be very convenient, and the
      specfic lookup method has seen less use[*]. The final step in
      supporting a new variable is to update the function
      rd_smspec_fread_header().

      If you want to support specific lookup of the new variable type you
      must add the necessary datastructures to the rd_smspec_struct
      structure and then subsequently fill that structure in the big
      switch() in rd_smspec_fread_header() - if you do not care about
      specific lookup you just have to add an empty case() slot to the
      switch in rd_smspec_fread_header(). The LGR variables, and also
      RD_SMSPEC_SEGMENT_VAR do not support specific lookup.

      [*]: The advantage of the specific lookup is that it is possible
           to supply better error messages (The well 'XX' does not
           exist, instead of just unknown key: 'WWCT:XX'), and it is
           also possible to support queries like: give me all the
           well names.

  4. Mark the variable type as supported with a 'X' in the defintion of
     rd_smspec_var_type in rd_smspec.h.

 */

#define RD_SMSPEC_ID 806647
#define PARAMS_GLOBAL_DEFAULT -99

typedef std::map<std::string, const rd::smspec_node *> node_map;

struct rd_smspec_struct {
    UTIL_TYPE_ID_DECLARATION;
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

    std::vector<std::unique_ptr<rd::smspec_node>> smspec_nodes;
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
    int day_index; /* are used by the rd_sum_data object to locate per. timestep */
    int month_index; /* time information. */
    int year_index;
    bool has_lgr;
    std::vector<float> params_default;

    std::string restart_case;
    ert_rd_unit_enum unit_system;
    int restart_step;
};

/**
About indexing:
---------------

The RDISPE summary files are organised (roughly) like this:

 1. A header-file called xxx.SMPSEC is written, which is common to
    every timestep.

 2. For each timestep the summary vector is written in the form of a
    vector 'PARAMS' of length N with floats. In the PARAMS vector all
    types of data are stacked togeheter, and one must use the header
    info in the SMSPEC file to disentangle the summary data.

Here I will try to describe how the header in SMSPEC is organised, and
how that support is imlemented here. The SMSPEC header is organized
around three character vectors, of length N. To find the position in
the PARAMS vector of a certain quantity, you must consult one, two or
all three of these vectors. The most important vecor - which must
always be consulted is the KEYWORDS vector, then it is the WGNAMES and
NUMS (integer) vectors whcih must be consulted for some variable
types.


Let us a consider a system consisting of:

  * Two wells: P1 and P2 - for each well we have variables WOPR, WWCT
    and WGOR.

  * Three regions: For each region we have variables RPR and RXX(??)

  * We have stored field properties FOPT and FWPT


KEYWORDS = ['TIME','FOPR','FPR','FWCT','WOPR','WOPR,'WWCT','WWCT]
       ....



general_var:
------------
VAR_TYPE:(WELL_NAME|GROUP_NAME|NUMBER):NUMBER

Field var:         VAR_TYPE
Misc var:          VAR_TYPE
Well var:          VAR_TYPE:WELL_NAME
Group var:         VAR_TYPE:GROUP_NAME
Block var:         VAR_TYPE:i,j,k  (where i,j,k is calculated form NUM)
Region var         VAR_TYPE:index  (where index is NOT from the nums vector, it it is just an offset).
Completion var:    VAR_TYPE:WELL_NAME:NUM
....

*/

/*
  The smspec_required_keywords variable contains a list of keywords
  which are *absolutely* required in the SMSPEC file, but observe that
  depending on the content of the "KEYWORDS" array other keywords
  might bre requred as well - this typically includes the NUMS
  keyword. Such 'second-order' dependencies are not accounted for with
  this simple list.
*/

static const size_t num_req_keywords = 5;
static const char *smspec_required_keywords[] = {
    WGNAMES_KW, KEYWORDS_KW, STARTDAT_KW, UNITS_KW, DIMENS_KW};

namespace {

const rd::smspec_node *rd_smspec_get_var_node(const node_map &mp,
                                              const char *var) {
    const auto it = mp.find(var);
    if (it == mp.end())
        return nullptr;

    return it->second;
}

} //end namespace

int rd_smspec_num_nodes(const rd_smspec_type *smspec) {
    return smspec->smspec_nodes.size();
}

/*
  When loading a summary case from file many of the nodes can be ignored, in
  that case the size of PARAMS vector in the data files is larger than the
  number of internalized nodes. Therefor we need to maintain the
  params_size member.
*/

int rd_smspec_get_params_size(const rd_smspec_type *smspec) {
    return smspec->params_size;
}

static rd_smspec_ptr rd_smspec_alloc_empty(bool write_mode,
                                           const std::string &key_join_string) {
    rd_smspec_ptr rd_smspec{new rd_smspec_type(), &rd_smspec_free};
    UTIL_TYPE_ID_INIT(rd_smspec.get(), RD_SMSPEC_ID);

    rd_smspec->sim_start_time = -1;
    rd_smspec->key_join_string = key_join_string;
    rd_smspec->header_file = "";

    rd_smspec->time_index = -1;
    rd_smspec->day_index = -1;
    rd_smspec->year_index = -1;
    rd_smspec->month_index = -1;
    rd_smspec->time_seconds = -1;
    rd_smspec->params_size = -1;

    /*
    The unit system is given as an integer in the INTEHEAD keyword. The INTEHEAD
    keyword is optional, and we have for a long time been completely oblivious
    to the possibility of extracting unit system information from the SMSPEC file.
  */
    rd_smspec->unit_system = RD_METRIC_UNITS;

    rd_smspec->restart_step = -1;
    rd_smspec->write_mode = write_mode;
    rd_smspec->need_nums = false;

    return rd_smspec;
}

std::vector<int> rd_smspec_alloc_mapping(const rd_smspec_type *self,
                                         const rd_smspec_type *other) {
    int params_size = rd_smspec_get_params_size(self);
    std::vector<int> mapping(params_size, -1);

    for (int i = 0; i < rd_smspec_num_nodes(self); i++) {
        const rd::smspec_node &self_node =
            rd_smspec_iget_node_w_node_index(self, i);
        int self_index = self_node.get_params_index();
        const char *key = self_node.get_gen_key1();
        if (rd_smspec_has_general_var(other, key)) {
            const rd::smspec_node &other_node =
                rd_smspec_get_general_var_node(other, key);
            int other_index = other_node.get_params_index();
            mapping[self_index] = other_index;
        }
    }

    return mapping;
}

/**
   Observe that the index here is into the __INTERNAL__ indexing in
   the smspec_nodes vector; and in general widely different from the
   params_index of the returned smspec_node instance.
*/

const rd::smspec_node &
rd_smspec_iget_node_w_node_index(const rd_smspec_type *smspec, int node_index) {
    const auto &node = smspec->smspec_nodes[node_index];
    return *node.get();
}

const rd::smspec_node &
rd_smspec_iget_node_w_params_index(const rd_smspec_type *smspec,
                                   int params_index) {
    int node_index = smspec->inv_index_map.at(params_index);
    return rd_smspec_iget_node_w_node_index(smspec, node_index);
}

/**
 * Returns an ecl data type for which all names will fit. If the maximum name
 * length is at most 8, an RD_CHAR is returned and otherwise a large enough
 * RD_STRING.
 */
static rd_data_type get_wgnames_type(const rd_smspec_type *smspec) {
    size_t max_len = 0;
    for (int i = 0; i < rd_smspec_num_nodes(smspec); ++i) {
        const rd::smspec_node &node =
            rd_smspec_iget_node_w_node_index(smspec, i);
        const char *name = smspec_node_get_wgname(&node);
        if (name)
            max_len = util_size_t_max(max_len, strlen(name));
    }

    return max_len <= RD_STRING8_LENGTH ? RD_CHAR : RD_STRING(max_len);
}

static void rd_smspec_fwrite_INTEHEAD(const rd_smspec_type *smspec,
                                      fortio_type *fortio) {
    rd_kw_ptr intehead = make_rd_kw(INTEHEAD_KW, INTEHEAD_SMSPEC_SIZE, RD_INT);
    rd_kw_iset_int(intehead.get(), INTEHEAD_SMSPEC_UNIT_INDEX,
                   smspec->unit_system);
    /* The simulator type is just hardcoded to ECLIPSE100. */
    rd_kw_iset_int(intehead.get(), INTEHEAD_SMSPEC_IPROG_INDEX,
                   INTEHEAD_ECLIPSE100_VALUE);
    rd_kw_fwrite(intehead.get(), fortio);
}

static void rd_smspec_fwrite_RESTART(const rd_smspec_type *smspec,
                                     fortio_type *fortio) {
    rd_kw_ptr restart_kw =
        make_rd_kw(RESTART_KW, SUMMARY_RESTART_SIZE, RD_CHAR);
    for (int i = 0; i < SUMMARY_RESTART_SIZE; i++)
        rd_kw_iset_string8(restart_kw.get(), i, "");

    if (smspec->restart_case.size() > 0) {
        size_t restart_case_len = smspec->restart_case.size();

        size_t offset = 0;
        for (size_t i = 0; i < SUMMARY_RESTART_SIZE; i++) {
            if (offset < restart_case_len)
                rd_kw_iset_string8(restart_kw.get(), i,
                                   &smspec->restart_case[offset]);
            offset += RD_STRING8_LENGTH;
        }
    }
    rd_kw_fwrite(restart_kw.get(), fortio);
}

static void rd_smspec_fwrite_DIMENS(const rd_smspec_type *smspec,
                                    fortio_type *fortio) {
    rd_kw_ptr dimens_kw = make_rd_kw(DIMENS_KW, DIMENS_SIZE, RD_INT);
    int num_nodes = rd_smspec_num_nodes(smspec);
    rd_kw_iset_int(dimens_kw.get(), DIMENS_SMSPEC_SIZE_INDEX, num_nodes);
    rd_kw_iset_int(dimens_kw.get(), DIMENS_SMSPEC_NX_INDEX,
                   smspec->grid_dims[0]);
    rd_kw_iset_int(dimens_kw.get(), DIMENS_SMSPEC_NY_INDEX,
                   smspec->grid_dims[1]);
    rd_kw_iset_int(dimens_kw.get(), DIMENS_SMSPEC_NZ_INDEX,
                   smspec->grid_dims[2]);
    rd_kw_iset_int(dimens_kw.get(), 4, 0); // Do not know what this is for.
    rd_kw_iset_int(dimens_kw.get(), DIMENS_SMSPEC_RESTART_STEP_INDEX,
                   smspec->restart_step);

    rd_kw_fwrite(dimens_kw.get(), fortio);
}

static void rd_smspec_fwrite_STARTDAT(const rd_smspec_type *smspec,
                                      fortio_type *fortio) {
    auto startdat_kw = make_rd_kw(STARTDAT_KW, STARTDAT_SIZE, RD_INT);
    int second, minute, hour, mday, month, year;
    rd_set_datetime_values(smspec->sim_start_time, &second, &minute, &hour,
                           &mday, &month, &year);

    rd_kw_iset_int(startdat_kw.get(), STARTDAT_DAY_INDEX, mday);
    rd_kw_iset_int(startdat_kw.get(), STARTDAT_MONTH_INDEX, month);
    rd_kw_iset_int(startdat_kw.get(), STARTDAT_YEAR_INDEX, year);
    rd_kw_iset_int(startdat_kw.get(), STARTDAT_HOUR_INDEX, hour);
    rd_kw_iset_int(startdat_kw.get(), STARTDAT_MINUTE_INDEX, minute);
    rd_kw_iset_int(startdat_kw.get(), STARTDAT_MICRO_SECOND_INDEX,
                   second * 1000000);

    rd_kw_fwrite(startdat_kw.get(), fortio);
}

static void rd_smspec_fortio_fwrite(const rd_smspec_type *smspec,
                                    fortio_type *fortio) {
    rd_smspec_fwrite_INTEHEAD(smspec, fortio);
    rd_smspec_fwrite_RESTART(smspec, fortio);
    rd_smspec_fwrite_DIMENS(smspec, fortio);

    int num_nodes = rd_smspec_num_nodes(smspec);
    auto keywords_kw = make_rd_kw(KEYWORDS_KW, num_nodes, RD_CHAR);
    auto units_kw = make_rd_kw(UNITS_KW, num_nodes, RD_CHAR);
    rd_kw_ptr nums_kw{nullptr, &rd_kw_free};

    // If the names_type is an RD_STRING we expect this to be an INTERSECT
    // summary, otherwise an ECLIPSE summary.
    rd_data_type names_type = get_wgnames_type(smspec);
    auto wgnames_kw =
        make_rd_kw(rd_type_is_char(names_type) ? WGNAMES_KW : NAMES_KW,
                   num_nodes, names_type);

    if (smspec->need_nums)
        nums_kw.reset(rd_kw_alloc(NUMS_KW, num_nodes, RD_INT));

    for (int i = 0; i < rd_smspec_num_nodes(smspec); i++) {
        const rd::smspec_node &smspec_node =
            rd_smspec_iget_node_w_node_index(smspec, i);
        /*
      It is possible to add variables with deferred initialisation
      with the rd_sum_add_blank_var() function. Before these
      variables can be actually used for anything interesting they
      must be initialized with the rd_sum_init_var() function.

      If a call to save the smspec file comes before all the
      variable have been initialized things will potentially go
      belly up. This is solved with the following uber-hack:

      o One of the well related keywords is chosen; in
        particular 'WWCT' in this case.

      o The wgname value is set to DUMMY_WELL

      The use of DUMMY_WELL ensures that this field will be
      ignored when/if this smspec file is read in at a later
      stage.
    */
        if (smspec_node.get_var_type() == RD_SMSPEC_INVALID_VAR) {
            rd_kw_iset_string8(keywords_kw.get(), i, "WWCT");
            rd_kw_iset_string8(units_kw.get(), i, "????????");
            rd_kw_iset_string_ptr(wgnames_kw.get(), i, DUMMY_WELL);
        } else {
            rd_kw_iset_string8(keywords_kw.get(), i,
                               smspec_node_get_keyword(&smspec_node));
            rd_kw_iset_string8(units_kw.get(), i,
                               smspec_node_get_unit(&smspec_node));
            {
                const char *wgname = DUMMY_WELL;
                if (smspec_node_get_wgname(&smspec_node))
                    wgname = smspec_node_get_wgname(&smspec_node);
                rd_kw_iset_string_ptr(wgnames_kw.get(), i, wgname);
            }
        }

        if (nums_kw)
            rd_kw_iset_int(nums_kw.get(), i, smspec_node.get_num());
    }
    rd_kw_fwrite(keywords_kw.get(), fortio);
    rd_kw_fwrite(wgnames_kw.get(), fortio);
    if (nums_kw)
        rd_kw_fwrite(nums_kw.get(), fortio);
    rd_kw_fwrite(units_kw.get(), fortio);

    rd_smspec_fwrite_STARTDAT(smspec, fortio);
}

void rd_smspec_fwrite(const rd_smspec_type *smspec, const char *rd_case,
                      bool fmt_file) {
    std::string filename =
        rd::filename(rd_case, RD_SUMMARY_HEADER_FILE, fmt_file, 0).string();
    ERT::FortIO fortio(filename, std::ios_base::out, fmt_file);

    if (!fortio.get())
        throw std::runtime_error(
            fmt::format("Unable to open fortio file {}, error: {}", filename,
                        strerror(errno)));

    rd_smspec_fortio_fwrite(smspec, fortio.get());
}

static rd_smspec_type *
rd_smspec_alloc_writer__(const char *key_join_string, const char *restart_case,
                         int restart_step, time_t sim_start, bool time_in_days,
                         int nx, int ny, int nz) {
    rd_smspec_ptr rd_smspec = rd_smspec_alloc_empty(true, key_join_string);
    /*
    Only a total of 9 * 8 characters is set aside for the restart keyword, if
    the supplied restart case is longer than that we silently ignore it.
  */
    if (restart_case) {
        if (strlen(restart_case) <=
            (SUMMARY_RESTART_SIZE * RD_STRING8_LENGTH)) {
            rd_smspec->restart_case = restart_case;
            rd_smspec->restart_step = restart_step;
        }
    }
    rd_smspec->grid_dims[0] = nx;
    rd_smspec->grid_dims[1] = ny;
    rd_smspec->grid_dims[2] = nz;
    rd_smspec->sim_start_time = sim_start;

    {
        const rd::smspec_node *time_node;

        if (time_in_days) {
            rd_smspec->time_seconds = 3600 * 24;
            time_node = rd_smspec_add_node(rd_smspec.get(), "TIME", "DAYS", 0);
        } else {
            rd_smspec->time_seconds = 3600;
            time_node = rd_smspec_add_node(rd_smspec.get(), "TIME", "HOURS", 0);
        }
        rd_smspec->time_index = time_node->get_params_index();
    }
    return rd_smspec.release();
}

rd_smspec_type *rd_smspec_alloc_restart_writer(
    const char *key_join_string, const char *restart_case, int restart_step,
    time_t sim_start, bool time_in_days, int nx, int ny, int nz) {
    return rd_smspec_alloc_writer__(key_join_string, restart_case, restart_step,
                                    sim_start, time_in_days, nx, ny, nz);
}

rd_smspec_type *rd_smspec_alloc_writer(const char *key_join_string,
                                       time_t sim_start, bool time_in_days,
                                       int nx, int ny, int nz) {
    return rd_smspec_alloc_writer__(key_join_string, NULL, 0, sim_start,
                                    time_in_days, nx, ny, nz);
}

UTIL_SAFE_CAST_FUNCTION(rd_smspec, RD_SMSPEC_ID)

rd_smspec_var_type rd_smspec_identify_var_type(const char *var) {
    return rd::smspec_node::identify_var_type(var);
}

static bool rd_smspec_lgr_var_type(rd_smspec_var_type var_type) {
    if ((var_type == RD_SMSPEC_LOCAL_BLOCK_VAR) ||
        (var_type == RD_SMSPEC_LOCAL_WELL_VAR) ||
        (var_type == RD_SMSPEC_LOCAL_COMPLETION_VAR))

        return true;
    else
        return false;
}

/**
   This function takes a fully initialized smspec_node instance, generates the
   corresponding key and inserts smspec_node instance in the main hash table
   smspec->gen_var_index.

   The format strings used, i.e. VAR:WELL for well based variables is implicitly
   defined through the format strings used in this function.
*/

static void rd_smspec_install_gen_keys(rd_smspec_type *smspec,
                                       const rd::smspec_node &smspec_node) {
    /* Insert the default general mapping. */
    {
        const char *gen_key1 = smspec_node.get_gen_key1();
        if (gen_key1)
            smspec->gen_var_index[gen_key1] = &smspec_node;
    }

    /* Insert the (optional) extra mapping for block related variables and region_2_region variables: */
    {
        const char *gen_key2 = smspec_node.get_gen_key2();
        if (gen_key2)
            smspec->gen_var_index[gen_key2] = &smspec_node;
    }
}

static void rd_smspec_install_special_keys(rd_smspec_type *rd_smspec,
                                           const rd::smspec_node &smspec_node) {
    /**
      This large switch is for installing keys which have custom lookup
      paths, in addition to the lookup based on general keys. Examples
      of this is e.g. well variables which can be looked up through:

      rd_smspec_get_well_var_index( smspec , well_name , var );
  */

    const char *well = smspec_node_get_wgname(&smspec_node);
    const char *group = well;
    const int num = smspec_node_get_num(&smspec_node);
    const char *keyword = smspec_node_get_keyword(&smspec_node);
    rd_smspec_var_type var_type = smspec_node_get_var_type(&smspec_node);

    switch (var_type) {
    case (RD_SMSPEC_COMPLETION_VAR):
        rd_smspec->well_completion_var_index[well][num][keyword] = &smspec_node;
        break;
    case (RD_SMSPEC_FIELD_VAR):
        rd_smspec->field_var_index[keyword] = &smspec_node;
        break;
    case (RD_SMSPEC_GROUP_VAR):
        rd_smspec->group_var_index[group][keyword] = &smspec_node;
        break;
    case (RD_SMSPEC_REGION_VAR):
        rd_smspec->region_var_index[num][keyword] = &smspec_node;
        rd_smspec->num_regions = util_int_max(rd_smspec->num_regions, num);
        break;
    case (RD_SMSPEC_WELL_VAR):
        rd_smspec->well_var_index[well][keyword] = &smspec_node;
        break;
    case (RD_SMSPEC_MISC_VAR):
        /* Misc variable - i.e. date or CPU time ... */
        rd_smspec->misc_var_index[keyword] = &smspec_node;
        break;
    case (RD_SMSPEC_BLOCK_VAR):
        rd_smspec->block_var_index[num][keyword] = &smspec_node;
        break;
        /**
        The variables below are ONLY accesable through the gen_key
        setup; but the must be mentioned in this switch statement,
        otherwise they will induce a hard failure in the default: target
        below.
    */
    case (RD_SMSPEC_LOCAL_BLOCK_VAR):
        break;
    case (RD_SMSPEC_LOCAL_COMPLETION_VAR):
        break;
    case (RD_SMSPEC_LOCAL_WELL_VAR):
        break;
    case (RD_SMSPEC_SEGMENT_VAR):
        break;
    case (RD_SMSPEC_REGION_2_REGION_VAR):
        break;
    case (RD_SMSPEC_AQUIFER_VAR):
        break;
    case (RD_SMSPEC_NETWORK_VAR):
        break;
    default:
        throw std::invalid_argument("Internal error - should not be here \n");
    }
}

bool rd_smspec_equal(const rd_smspec_type *self, const rd_smspec_type *other) {
    if (self->smspec_nodes.size() != other->smspec_nodes.size())
        return false;

    for (size_t i = 0; i < self->smspec_nodes.size(); i++) {
        const rd::smspec_node *node1 = self->smspec_nodes[i].get();
        const rd::smspec_node *node2 = other->smspec_nodes[i].get();

        if (node1->cmp(*node2) != 0)
            return false;
    }

    return true;
}

static void rd_smspec_load_restart(rd_smspec_type *rd_smspec,
                                   const rd_file_type *header) {
    if (!rd_file_has_kw(header, RESTART_KW))
        return;
    const rd_kw_type *restart_kw = rd_file_iget_named_kw(header, RESTART_KW, 0);
    int num_blocks = rd_kw_get_size(restart_kw);
    num_blocks = (num_blocks < 0) ? 0 : num_blocks;
    auto tmp_base = rd::checked_calloc<char>(8 * num_blocks + 1);
    for (int i = 0; i < num_blocks; i++) {
        const char *part = (const char *)rd_kw_iget_ptr(restart_kw, i);
        strncat(tmp_base.get(), part, 8);
    }

    std::string restart_base = rd::strip_spaces(std::string(tmp_base.get()));

    /* We ignore the empty ones. */
    if (!restart_base.size())
        return;

    std::string smspec_header;
    fs::path dir = fs::path(rd_smspec->header_file).parent_path();

    /*
    The conditional block here is to support the following situation:

    1. A simulation with a restart has been performed on Posix with path
       separator '/'.

    2. The simulation is loaded on windows, where the native path
       separator is '\'.

    This code block will translate '/' -> '\' in the restart keyword which
    is read from the summary file.
    */

#ifdef ERT_WINDOWS
    for (int i = 0; i < restart_base.size(); i++) {
        if (restart_base[i] == UTIL_POSIX_PATH_SEP_CHAR)
            restart_base[i] = UTIL_PATH_SEP_CHAR;
    }
#endif

    fs::path restart_path(restart_base);
    if (restart_path.is_relative()) {
        restart_path = dir / restart_path;
    }
    smspec_header = rd::filename(restart_path, RD_SUMMARY_HEADER_FILE,
                                 rd_smspec->formatted, 0)
                        .string();

    std::error_code ec;
    if (!fs::exists(smspec_header, ec) || ec)
        return;

    // Restart from current case is ignored
    if (util_same_file(smspec_header.c_str(), rd_smspec->header_file.c_str()))
        return;

    // restart_path should not exist, but we are checking
    // anyways to keep existing behavior
    if (fs::exists(restart_path, ec) && !ec)
        restart_path = fs::canonical(restart_path);
    rd_smspec->restart_case = restart_path.string();
}

static const rd::smspec_node *
rd_smspec_insert_node(rd_smspec_type *rd_smspec,
                      std::unique_ptr<rd::smspec_node> smspec_node) {
    int params_index = smspec_node->get_params_index();

    /* This indexing must be used when writing. */
    rd_smspec->index_map.push_back(params_index);
    rd_smspec->params_default.resize(params_index + 1, PARAMS_GLOBAL_DEFAULT);
    rd_smspec->params_default[params_index] = smspec_node->get_default();
    rd_smspec->inv_index_map.insert(
        std::make_pair(params_index, rd_smspec->smspec_nodes.size()));

    rd_smspec_install_gen_keys(rd_smspec, *smspec_node.get());
    rd_smspec_install_special_keys(rd_smspec, *smspec_node.get());

    if (smspec_node->need_nums())
        rd_smspec->need_nums = true;

    rd_smspec->smspec_nodes.push_back(std::move(smspec_node));

    if (params_index > rd_smspec->params_size)
        rd_smspec->params_size = params_index + 1;

    if (static_cast<int>(rd_smspec->smspec_nodes.size()) >
        rd_smspec->params_size)
        rd_smspec->params_size = rd_smspec->smspec_nodes.size();

    const auto &node = rd_smspec->smspec_nodes.back();
    return node.get();
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          const char *keyword, int num,
                                          const char *unit,
                                          float default_value) {
    int params_index = rd_smspec->smspec_nodes.size();
    return rd_smspec_insert_node(
        rd_smspec, std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
                       params_index, keyword, num, unit, rd_smspec->grid_dims,
                       default_value, rd_smspec->key_join_string.c_str())));
}

//copy given node with a new index
const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          const rd::smspec_node &node) {
    int params_index = rd_smspec->smspec_nodes.size();
    return rd_smspec_insert_node(rd_smspec,
                                 std::unique_ptr<rd::smspec_node>(
                                     new rd::smspec_node(node, params_index)));
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          const char *keyword, const char *unit,
                                          float default_value) {
    int params_index = rd_smspec->smspec_nodes.size();
    return rd_smspec_insert_node(
        rd_smspec, std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
                       params_index, keyword, unit, default_value)));
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          const char *keyword,
                                          const char *wgname, const char *unit,
                                          float default_value) {
    int params_index = rd_smspec->smspec_nodes.size();
    return rd_smspec_insert_node(
        rd_smspec, std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
                       params_index, keyword, wgname, unit, default_value,
                       rd_smspec->key_join_string.c_str())));
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          const char *keyword,
                                          const char *wgname, int num,
                                          const char *unit,
                                          float default_value) {
    int params_index = rd_smspec->smspec_nodes.size();
    return rd_smspec_insert_node(
        rd_smspec,
        std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
            params_index, keyword, wgname, num, unit, rd_smspec->grid_dims,
            default_value, rd_smspec->key_join_string.c_str())));
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          int params_index, const char *keyword,
                                          const char *wgname, int num,
                                          const char *unit,
                                          float default_value) {
    return rd_smspec_insert_node(
        rd_smspec,
        std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
            params_index, keyword, wgname, num, unit, rd_smspec->grid_dims,
            default_value, rd_smspec->key_join_string.c_str())));
}

const rd::smspec_node *rd_smspec_add_node(rd_smspec_type *rd_smspec,
                                          int params_index, const char *keyword,
                                          const char *wgname, int num,
                                          const char *unit, const char *lgr,
                                          int lgr_i, int lgr_j, int lgr_k,
                                          float default_value) {
    return rd_smspec_insert_node(
        rd_smspec,
        std::unique_ptr<rd::smspec_node>(new rd::smspec_node(
            params_index, keyword, wgname, unit, lgr, lgr_i, lgr_j, lgr_k,
            default_value, rd_smspec->key_join_string.c_str())));
}

const int *rd_smspec_get_index_map(const rd_smspec_type *smspec) {
    return smspec->index_map.data();
}

/**
 * This function is to support the NAMES alias for WGNAMES. If similar
 * situations occur in the future, this is a sane starting point for general
 * support.
 */
static const char *get_active_keyword_alias(rd_file_type *header,
                                            const char *keyword) {
    if (strcmp(keyword, WGNAMES_KW) == 0 || strcmp(keyword, NAMES_KW) == 0)
        return rd_file_has_kw(header, WGNAMES_KW) ? WGNAMES_KW : NAMES_KW;

    return keyword;
}

static bool rd_smspec_check_header(rd_file_type *header) {
    bool OK = true;
    for (size_t i = 0; i < num_req_keywords && OK; i++) {
        OK &= rd_file_has_kw(header, get_active_keyword_alias(
                                         header, smspec_required_keywords[i]));
    }

    return OK;
}

static bool rd_smspec_fread_header(rd_smspec_type *rd_smspec,
                                   const std::string &header_file,
                                   bool include_restart) {
    rd_file_type *header = rd_file_open(header_file.c_str(), 0);
    if (header && rd_smspec_check_header(header)) {
        const char *names_alias = get_active_keyword_alias(header, WGNAMES_KW);
        rd_kw_type *wells = rd_file_iget_named_kw(header, names_alias, 0);
        rd_kw_type *keywords = rd_file_iget_named_kw(header, KEYWORDS_KW, 0);
        rd_kw_type *startdat = rd_file_iget_named_kw(header, STARTDAT_KW, 0);
        rd_kw_type *units = rd_file_iget_named_kw(header, UNITS_KW, 0);
        rd_kw_type *dimens = rd_file_iget_named_kw(header, DIMENS_KW, 0);
        rd_kw_type *nums = NULL;
        rd_kw_type *lgrs = NULL;
        rd_kw_type *numlx = NULL;
        rd_kw_type *numly = NULL;
        rd_kw_type *numlz = NULL;

        int params_index;
        rd_smspec->num_regions = 0;
        rd_smspec->params_size = rd_kw_get_size(keywords);
        if (startdat == NULL)
            throw std::invalid_argument(
                "Could not locate STARTDAT keyword in header");

        if (rd_file_has_kw(header, NUMS_KW))
            nums = rd_file_iget_named_kw(header, NUMS_KW, 0);

        if (rd_file_has_kw(header, INTEHEAD_KW)) {
            const rd_kw_type *intehead =
                rd_file_iget_named_kw(header, INTEHEAD_KW, 0);
            rd_smspec->unit_system = (ert_rd_unit_enum)rd_kw_iget_int(
                intehead, INTEHEAD_SMSPEC_UNIT_INDEX);
            /*
        The second item in the INTEHEAD vector is an integer designating which
        simulator has been used for the current simulation, that is currently
        ignored.
      */
        }

        if (rd_file_has_kw(header,
                           LGRS_KW)) { /* The file has LGR information. */
            lgrs = rd_file_iget_named_kw(header, LGRS_KW, 0);
            numlx = rd_file_iget_named_kw(header, NUMLX_KW, 0);
            numly = rd_file_iget_named_kw(header, NUMLY_KW, 0);
            numlz = rd_file_iget_named_kw(header, NUMLZ_KW, 0);
            rd_smspec->has_lgr = true;
        } else
            rd_smspec->has_lgr = false;

        {
            if (rd_kw_get_size(startdat) < 3) {
                int startdat_size = rd_kw_get_size(startdat);
                rd_file_close(header);
                throw std::invalid_argument(fmt::format(
                    "Invalid size of STARTDAT keyword = {}, expected at "
                    "least 3",
                    startdat_size));
            }
            int *date = rd_kw_get_int_ptr(startdat);
            int year = date[STARTDAT_YEAR_INDEX];
            int month = date[STARTDAT_MONTH_INDEX];
            int day = date[STARTDAT_DAY_INDEX];
            int hour = 0;
            int min = 0;
            int sec = 0;
            if (rd_kw_get_size(startdat) == 6) {
                hour = date[STARTDAT_HOUR_INDEX];
                min = date[STARTDAT_MINUTE_INDEX];
                sec = date[STARTDAT_MICRO_SECOND_INDEX] / 1000000;
            }

            rd_smspec->sim_start_time =
                rd_make_datetime(sec, min, hour, day, month, year);
        }

        rd_smspec->grid_dims[0] =
            rd_kw_iget_int(dimens, DIMENS_SMSPEC_NX_INDEX);
        rd_smspec->grid_dims[1] =
            rd_kw_iget_int(dimens, DIMENS_SMSPEC_NY_INDEX);
        rd_smspec->grid_dims[2] =
            rd_kw_iget_int(dimens, DIMENS_SMSPEC_NZ_INDEX);
        rd_smspec->restart_step =
            rd_kw_iget_int(dimens, DIMENS_SMSPEC_RESTART_STEP_INDEX);

        rd_get_file_type(header_file.c_str(), &rd_smspec->formatted, NULL);

        {
            for (params_index = 0; params_index < rd_kw_get_size(wells);
                 params_index++) {
                float default_value = PARAMS_GLOBAL_DEFAULT;
                int num = SMSPEC_NUMS_INVALID;
                std::string well =
                    rd_kw_iget_stripped_string(wells, params_index);
                std::string kw =
                    rd_kw_iget_stripped_string(keywords, params_index);
                std::string unit =
                    rd_kw_iget_stripped_string(units, params_index);

                rd_smspec_var_type var_type;
                if (nums != NULL)
                    num = rd_kw_iget_int(nums, params_index);
                var_type =
                    rd::smspec_node::valid_type(kw.c_str(), well.c_str(), num);
                if (var_type == RD_SMSPEC_INVALID_VAR) {
                    continue;
                }

                if (rd_smspec_lgr_var_type(var_type)) {
                    int lgr_i = rd_kw_iget_int(numlx, params_index);
                    int lgr_j = rd_kw_iget_int(numly, params_index);
                    int lgr_k = rd_kw_iget_int(numlz, params_index);
                    std::string lgr_name =
                        rd_kw_iget_stripped_string(lgrs, params_index);

                    rd_smspec_insert_node(
                        rd_smspec,
                        std::make_unique<rd::smspec_node>(
                            params_index, kw.c_str(), well.c_str(),
                            unit.c_str(), lgr_name.c_str(), lgr_i, lgr_j, lgr_k,
                            default_value, rd_smspec->key_join_string.c_str()));
                } else
                    rd_smspec_insert_node(
                        rd_smspec,
                        std::make_unique<rd::smspec_node>(
                            params_index, kw.c_str(), well.c_str(), num,
                            unit.c_str(), rd_smspec->grid_dims, default_value,
                            rd_smspec->key_join_string.c_str()));
            }
        }

        rd_smspec->header_file = fs::canonical(header_file).string();
        if (include_restart)
            rd_smspec_load_restart(rd_smspec, header);

        rd_file_close(header);

        return true;
    } else
        return false;
}

rd_smspec_type *rd_smspec_fread_alloc(const std::string &header_file,
                                      const std::string &key_join_string,
                                      bool include_restart) {
    rd_smspec_ptr rd_smspec = rd_smspec_alloc_empty(false, key_join_string);

    if (rd_smspec_fread_header(rd_smspec.get(), header_file, include_restart)) {

        const rd::smspec_node *time_node =
            rd_smspec_get_var_node(rd_smspec->misc_var_index, "TIME");
        if (time_node) {
            const char *time_unit = time_node->get_unit();
            rd_smspec->time_index = time_node->get_params_index();

            if (util_string_equal(time_unit, "DAYS"))
                rd_smspec->time_seconds = 3600 * 24;
            else if (util_string_equal(time_unit, "HOURS"))
                rd_smspec->time_seconds = 3600;
            else
                throw std::invalid_argument(
                    fmt::format("time_unit:{} not recognized", time_unit));
        }

        const rd::smspec_node *day_node =
            rd_smspec_get_var_node(rd_smspec->misc_var_index, "DAY");
        const rd::smspec_node *month_node =
            rd_smspec_get_var_node(rd_smspec->misc_var_index, "MONTH");
        const rd::smspec_node *year_node =
            rd_smspec_get_var_node(rd_smspec->misc_var_index, "YEAR");

        if (day_node != NULL && month_node != NULL && year_node != NULL) {
            rd_smspec->day_index = day_node->get_params_index();
            rd_smspec->month_index = month_node->get_params_index();
            rd_smspec->year_index = year_node->get_params_index();
        }

        if ((rd_smspec->time_index == -1) && (rd_smspec->day_index == -1)) {
            // Unusable configuration.
            // Seems the restart file can also have time specified with
            // 'YEARS' as basic time unit; that mode is not supported.

            throw std::invalid_argument(
                "The SMSPEC file seems to lack all time information, need "
                "either TIME, or DAY/MONTH/YEAR information. Can not "
                "proceed.");
        }
        return rd_smspec.release();
    } else {
        /** Failed to load from disk. */
        return nullptr;
    }
}

/*
   For each type of summary data (according to the types in
   rd_smcspec_var_type there are a set accessor functions:

   xx_get_xx: This function will take the apropriate input, and
   return a double value. The function will fail
   if the rd_smspec object can not recognize the input. THis
   function is not here.

   xxx_has_xx: Ths will return true / false depending on whether the
   rd_smspec object the variable we ask for.

   xxx_get_xxx_index: This function will rerturn an (internal)
   integer index of where the variable in question is stored, this
   index can then be subsequently used for faster lookup. If the
   variable can not be found, the function will return -1.

   In general the index function is the real function, the others are
   only wrappers around this. In addition there are specialized
   functions, like get_well_names() and so on.
*/

namespace {

bool node_exists(const rd::smspec_node *node_ptr) {
    if (node_ptr)
        return true;

    return false;
}

int node_valid_index(const rd::smspec_node *node_ptr) {
    if (node_ptr)
        return node_ptr->get_params_index();

    throw std::out_of_range("Invalid lookup summary object");
}

} // namespace

/* There is a quite wide range of error which are just returned as
   "Not found" (i.e. -1). */
/* Completions not supported yet. */

const rd::smspec_node &
rd_smspec_get_general_var_node(const rd_smspec_type *smspec,
                               const char *lookup_kw) {
    const auto node_ptr =
        rd_smspec_get_var_node(smspec->gen_var_index, lookup_kw);
    if (!node_ptr)
        throw std::out_of_range("No such variable: " + std::string(lookup_kw));

    return *node_ptr;
}

int rd_smspec_get_general_var_params_index(const rd_smspec_type *rd_smspec,
                                           const char *lookup_kw) {
    const auto node_ptr =
        rd_smspec_get_var_node(rd_smspec->gen_var_index, lookup_kw);
    return node_valid_index(node_ptr);
}

bool rd_smspec_has_general_var(const rd_smspec_type *rd_smspec,
                               const char *lookup_kw) {
    const auto node_ptr =
        rd_smspec_get_var_node(rd_smspec->gen_var_index, lookup_kw);
    return node_exists(node_ptr);
}

int rd_smspec_get_time_seconds(const rd_smspec_type *rd_smspec) {
    return rd_smspec->time_seconds;
}

int rd_smspec_get_time_index(const rd_smspec_type *rd_smspec) {
    return rd_smspec->time_index;
}

time_t rd_smspec_get_start_time(const rd_smspec_type *rd_smspec) {
    return rd_smspec->sim_start_time;
}

bool rd_smspec_get_formatted(const rd_smspec_type *rd_smspec) {
    return rd_smspec->formatted;
}

const char *rd_smspec_get_header_file(const rd_smspec_type *rd_smspec) {
    return rd_smspec->header_file.c_str();
}

int rd_smspec_get_restart_step(const rd_smspec_type *rd_smspec) {
    return rd_smspec->restart_step;
}

int rd_smspec_get_first_step(const rd_smspec_type *rd_smspec) {
    if (rd_smspec->restart_step > 0)
        return rd_smspec->restart_step + 1;
    else
        return 1;
}

const char *rd_smspec_get_restart_case(const rd_smspec_type *rd_smspec) {
    if (rd_smspec->restart_case.size() > 0)
        return rd_smspec->restart_case.c_str();
    else
        return NULL;
}

const std::vector<float> &
rd_smspec_get_params_default(const rd_smspec_type *rd_smspec) {
    return rd_smspec->params_default;
}

void rd_smspec_free(rd_smspec_type *rd_smspec) { delete rd_smspec; }

int rd_smspec_get_date_day_index(const rd_smspec_type *smspec) {
    return smspec->day_index;
}

int rd_smspec_get_date_month_index(const rd_smspec_type *smspec) {
    return smspec->month_index;
}

int rd_smspec_get_date_year_index(const rd_smspec_type *smspec) {
    return smspec->year_index;
}

/**
   Fills a stringlist instance with all the gen_key string matching
   the supplied pattern. I.e.

     rd_smspec_alloc_matching_general_var_list( smspec , "WGOR:*");

   will give a list of WGOR for ALL the wells. The function is
   unfortunately not as useful as one might think because ECLIPSE
   will for instance happily give you the WOPR for a
   water injector or WWIR for an oil producer.

   The function can be called several times with different patterns,
   the stringlist is not cleared on startup; the keys in the list are
   unique - keys are not added multiple times. If pattern == NULL all
   keys will match.
*/

void rd_smspec_select_matching_general_var_list(const rd_smspec_type *smspec,
                                                const char *pattern,
                                                stringlist_type *keys) {
    std::set<std::string> ex_keys;
    for (int i = 0; i < stringlist_get_size(keys); i++)
        ex_keys.insert(stringlist_iget(keys, i));

    {
        for (const auto &pair : smspec->gen_var_index) {
            const char *key = pair.first.c_str();

            /*
         The TIME is typically special cased by output and will not
         match the 'all keys' wildcard.
      */
            if (util_string_equal(key, "TIME")) {
                if ((pattern == NULL) || (util_string_equal(pattern, "*")))
                    continue;
            }

            if ((pattern == NULL) || (util_fnmatch(pattern, key) == 0)) {
                if (ex_keys.find(key) == ex_keys.end())
                    stringlist_append_copy(keys, key);
            }
        }
    }

    stringlist_sort(keys, (string_cmp_ftype *)util_strcmp_int);
}

/**
   Allocates a new stringlist and initializes it with the
   rd_smspec_select_matching_general_var_list() function.
*/

stringlist_type *
rd_smspec_alloc_matching_general_var_list(const rd_smspec_type *smspec,
                                          const char *pattern) {
    stringlist_type *keys = stringlist_alloc_new();
    rd_smspec_select_matching_general_var_list(smspec, pattern, keys);
    return keys;
}

/**
    Returns a stringlist instance with all the (valid) well names. It
    is the responsability of the calling scope to free the stringlist
    with stringlist_free();


    If @pattern is different from NULL only wells which 'match' the
    pattern is included; if @pattern == NULL all wells are
    included. The match is done with function fnmatch() -
    i.e. standard shell wildcards.
*/

static stringlist_type *
rd_smspec_alloc_map_list(const std::map<std::string, node_map> &mp,
                         const char *pattern) {
    stringlist_type *map_list = stringlist_alloc_new();

    for (const auto &pair : mp) {
        const char *map_name = pair.first.c_str();

        if (pattern == NULL)
            stringlist_append_copy(map_list, map_name);
        else if (util_fnmatch(pattern, map_name) == 0)
            stringlist_append_copy(map_list, map_name);
    }
    stringlist_sort(map_list, (string_cmp_ftype *)util_strcmp_int);
    return map_list;
}

stringlist_type *rd_smspec_alloc_well_list(const rd_smspec_type *smspec,
                                           const char *pattern) {
    return rd_smspec_alloc_map_list(smspec->well_var_index, pattern);
}

/**
    Returns a stringlist instance with all the (valid) group names. It
    is the responsability of the calling scope to free the stringlist
    with stringlist_free();
*/

stringlist_type *rd_smspec_alloc_group_list(const rd_smspec_type *smspec,
                                            const char *pattern) {
    return rd_smspec_alloc_map_list(smspec->group_var_index, pattern);
}

const int *rd_smspec_get_grid_dims(const rd_smspec_type *smspec) {
    return smspec->grid_dims;
}

ert_rd_unit_enum rd_smspec_get_unit_system(const rd_smspec_type *smspec) {
    return smspec->unit_system;
}
