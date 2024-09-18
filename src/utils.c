

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "utils.h"


bool valid_idx(int x) {
  if (x == NA_INTEGER) {
    error("Not expecting NA here!");
  }
  
  if (x == INVALID_VIDX) {
    return false;
  }
  
  if (x < 0) {
    error("Not expecting -ve here");
  }
  
  return true;
}


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
// @param df_ named list object which is to be promoted to data.frame
// @param visible_length trim the data.frame to this number of rows
// @param allocated length how many rows does each list member have?
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_df_attributes_and_trim(SEXP df_, int visible_length, int allocated_length) {
  
  if (!isNewList(df_)) {
    error("set_df_attributes_and_trim(): only accepts 'lists' as input");
  }
  
  if (visible_length > allocated_length) {
    error("set_df_attributes_and_trim(): visible_length (%i) cannot be greater than allocated length (%i)",
          visible_length, allocated_length);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Truncate all the vectors to visible length
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (visible_length != allocated_length) {
    for (int col_idx = 0; col_idx < length(df_); col_idx++) {
      trim_vec(VECTOR_ELT(df_, col_idx), visible_length, allocated_length);
    }
  }
  
  set_df_attributes(df_);
} 



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// @param df_ named list object which is to be promoted to data.frame
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_df_attributes(SEXP df_) {
  int nprotect = 0;
  
  if (!isNewList(df_)) {
    error("set_df_attributes(): only accepts 'lists' as input");
  }
  
  int len = length(VECTOR_ELT(df_, 0));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Set row.names
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP rownames = PROTECT(allocVector(INTSXP, 2)); nprotect++;
  SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
  SET_INTEGER_ELT(rownames, 1, -len);
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
// Increment integer values by 1 to convert from C 0-indexing to
// R 1-indexing
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void convert_indexing_c_to_r_with_NA(SEXP ivec_) {
  int *ivec = INTEGER(ivec_);
  
  for (int i = 0; i < length(ivec_); i++) {
    if (ivec[i] == INVALID_VIDX) {
      ivec[i] = NA_INTEGER;
    } else {
      ivec[i]++;
    }
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Decrement integer values by 1 to convert from R 1-indexing to
// C 0-indexing
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void convert_indexing_r_to_c_with_NA(SEXP ivec_) {
  int *ivec = INTEGER(ivec_);
  
  for (int i = 0; i < length(ivec_); i++) {
    if (ivec[i] == NA_INTEGER) {
      ivec[i] = INVALID_VIDX;
    } else {
      ivec[i]--;
    }
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
SEXP create_named_list(int n, ...) {
  
  va_list args;
  va_start(args, n);
  
  int nprotect = 0;
  SEXP res_ = PROTECT(allocVector(VECSXP, n)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, n)); nprotect++;
  setAttrib(res_, R_NamesSymbol, nms_);
  
  for (int i = 0; i < n; i++) {
    
    const char *nm = va_arg(args, const char *);
    SEXP val_ = va_arg(args, SEXP);
    
    SET_STRING_ELT(nms_, i, mkChar(nm));
    SET_VECTOR_ELT(res_, i, val_);
  }
  
  
  va_end(args);
  UNPROTECT(nprotect);
  return res_;
}






