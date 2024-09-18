void trim_vec(SEXP vec_, int visible_length, int allocated_length);
void set_df_attributes(SEXP df_);
void set_df_attributes_and_trim(SEXP df_, int visible_length, int allocated_length);
void convert_indexing_r_to_c(SEXP ivec_);
void convert_indexing_c_to_r(SEXP ivec_);
void convert_indexing_r_to_c_with_NA(SEXP ivec_);
void convert_indexing_c_to_r_with_NA(SEXP ivec_);
int *create_c_index(SEXP ivec_);
SEXP create_named_list(int n, ...);
bool valid_idx(int x);

#define INVALID_VIDX -999

