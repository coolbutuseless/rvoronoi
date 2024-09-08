void trim_vec(SEXP vec_, int visible_length, int allocated_length);
void set_df_attributes(SEXP df_, int visible_length, int allocated_length);
void convert_indexing_r_to_c(SEXP ivec_);
void convert_indexing_c_to_r(SEXP ivec_);
int *create_c_index(SEXP ivec_);

