#pragma

typedef struct ecl_smspec_struct ecl_smspec_type;
typedef struct ecl_sum_data_struct ecl_sum_data_type;
typedef struct ecl_sum_struct ecl_sum_type;

#define ECL_SUM_ID 89067

struct ecl_sum_struct {
    int __type_id;
    ecl_smspec_type *smspec; /* Internalized version of the SMSPEC file. */
    ecl_sum_data_type *data; /* The data - can be NULL. */
    ecl_sum_type *restart_case;

    bool fmt_case;
    bool unified;
    char *key_join_string;
    char *
        path; /* The path - as given for the case input. Can be NULL for cwd. */
    char *abs_path; /* Absolute path. */
    char *base;     /* Only the basename. */
    char *
        ecl_case; /* This is the current case, with optional path component. == path + base*/
    char *
        ext; /* Only to support selective loading of formatted|unformatted and unified|multiple. (can be NULL) */
};
