

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "utils.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Trim a vector to the given length
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void trim_vec(SEXP vec_, int visible_length, int allocated_length) {
  if (visible_length != allocated_length) {
      SETLENGTH(vec_, visible_length); 
      SET_TRUELENGTH(vec_, allocated_length); 
      SET_GROWABLE_BIT(vec_);
  }
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// @param df_ data.frame object
// @param N number of rows in data.frame
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_df_attributes(SEXP df_, int visible_length, int allocated_length) {
  int nprotect = 0;
  
  if (!isNewList(df_)) {
    error("set_df_attributes(): only accepts 'lists' as input");
  }
  
  if (visible_length > allocated_length) {
    error("set_df_attributes(): visible_length cannot be greated than allocated length");
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Truncate all the vectors to visible length
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (visible_length != allocated_length) {
    for (int col_idx = 0; col_idx < length(df_); col_idx++) {
      SEXP col_ = VECTOR_ELT(df_, col_idx);
      
      SETLENGTH(col_, visible_length); 
      SET_TRUELENGTH(col_, allocated_length); 
      SET_GROWABLE_BIT(col_);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Set row.names
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP rownames = PROTECT(allocVector(INTSXP, 2)); nprotect++;
  SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
  SET_INTEGER_ELT(rownames, 1, -visible_length);
  setAttrib(df_, R_RowNamesSymbol, rownames);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Set as tibble
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP classnames = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(classnames, 0, mkChar("tbl_df"));
  SET_STRING_ELT(classnames, 1, mkChar("tbl"));
  SET_STRING_ELT(classnames, 2, mkChar("data.frame"));
  SET_CLASS(df_, classnames);
  
  UNPROTECT(nprotect);
} 


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Increment integer values by 1 to convert from C 0-indexing to
// R 1-indexing
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void convert_indexing_c_to_r(SEXP ivec_) {
  int *ivec = INTEGER(ivec_);
  
  for (int i = 0; i < length(ivec_); i++) {
    ivec[i]++;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Decrement integer values by 1 to convert from R 1-indexing to
// C 0-indexing
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void convert_indexing_r_to_c(SEXP ivec_) {
  int *ivec = INTEGER(ivec_);
  
  for (int i = 0; i < length(ivec_); i++) {
    ivec[i]--;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Allocate and create a C integer array containing a 0-indexed version
// of the R 1-indexed vector
// The caller is responsible for freeing the returned memory
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int *create_c_index(SEXP ivec_) {
  
  int *ivec = malloc(length(ivec_) * sizeof(int));
  if (ivec == NULL) {
    error("create_c_index(): Couldn't allocate %i members", length(ivec_));
  }
  
  int *ivecp = INTEGER(ivec_);
  for (int i = 0; i < length(ivec_); i++) {
    ivec[i] = ivecp[i] - 1;
  }
  
  return ivec;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Create a named list
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP create_named_list(int n, const char **nms) {
  
  int nprotect = 0;
  SEXP res_ = PROTECT(allocVector(VECSXP, n)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, n)); nprotect++;
  
  for (int i = 0; i < n; i++) {
    SET_STRING_ELT(nms_, i, mkChar(nms[i]));
  }
  
  setAttrib(res_, R_NamesSymbol, nms_);
  
  UNPROTECT(nprotect);
  return res_;
}






