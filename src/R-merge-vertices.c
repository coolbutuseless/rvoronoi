
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "defs.h"
#include "memory.h"
#include "geometry.h"
#include "voronoi.h"

#include "utils.h"
#include "R-merge-vertices.h"




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Caller is responsible for freeing fv1, fv2
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void merge_vertices_core_(double tol, int nverts, double *x, double *y, int nedges, int *v1, int *v2, int *fnedges,
                          int verbosity) {
  
  
  int *remap1 = malloc((unsigned long)nverts * sizeof(int));
  int *remap2 = malloc((unsigned long)nverts * sizeof(int));
  if (remap1 == NULL || remap2 == NULL) {
    error("merge_vertices_core_(): could not allocate 'remap1' and 'remap2'");
  }
  int rcnt = 0;
  
  
  for (int i = 0; i < nedges; i++) {
    int i1 = v1[i];
    int i2 = v2[i];
    
    if (i1 < 0 || i2 < 0) continue;
    
    if ( (i1 >= nverts) || (i2 >= nverts) || (i1 < 0) || (i2 < 0) ) {
      error("merge_vertices_core_(): [%i] Vertex index out of bounds: %i, %i [0, %i]", 
            i, i1, i2, nverts);
    }
    
    double x1 = x[i1];
    double x2 = x[i2];
    double y1 = y[i1];
    double y2 = y[i2];
    
    double dist = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    
    if (dist < tol) {
      if (i1 > i2) {
        remap1[rcnt] = i2;
        remap2[rcnt] = i1;  // remap2 always holds the higher index
      } else {
        remap1[rcnt] = i1;
        remap2[rcnt] = i2;
      }
      rcnt++;
    } 
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Collapse and remove any edges with duplicated vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%2i, %2i)\n", i, v1[i], v2[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Print remaps
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {
    for (int i = rcnt - 1; i >= 0; i--) {
      Rprintf("Remap: %i: %i == %i\n", i, remap1[i], remap2[i]);
    }
  }
  
  // Track which indices we want to keep
  int ndiscard = 0;
  bool *discard = calloc((unsigned int)nedges, sizeof(bool));
  if (discard == NULL) error("merge_vertices_core_(): could not allocate 'discard'");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Renumber verts
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = rcnt - 1; i >= 0; i--) {
    int vlo = remap1[i];
    int vhi = remap2[i];
    
    // Rprintf("%i: %i -> %i\n", i, vhi, vlo);
    
    for (int j = 0; j < nedges; j++) {
      if ( (v1[j] < 0) || (v2[j] < 0)) continue;
      if (v1[j] == vhi) v1[j] = vlo; 
      if (v2[j] == vhi) v2[j] = vlo; 
      discard[j] = v1[j] == v2[j];
    }
  }
  
  
  for (int i = 0; i < nedges; i++) {
    ndiscard += discard[i];
  }
  
  if (verbosity > 0) {
    Rprintf("Ndistcard: %i\n", ndiscard);
    for (int i = 0; i < nedges; i++) {
      Rprintf("(%i)  %i\n", i, discard[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Collapse and remove any edges with duplicated vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {  
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%2i, %2i)\n", i, v1[i], v2[i]);
    }
  }
  
  
  int src = 0;
  int dst = 0;
  for (dst = 0; dst < nedges - ndiscard; dst++, src++) {
    while(discard[src]) {
      src++;
    }
    // src += discard[src];
    v1[dst] = v1[src];
    v2[dst] = v2[src];
    if (verbosity > 1) {
      Rprintf("src %02i  ->  dst %02i\n", src, dst);
    }
  }
  
  for (; dst < nedges; dst++) {
    v1[dst] = 998;
    v2[dst] = 998;
  }
  
  if (verbosity > 0) {  
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%2i, %2i)\n", i, v1[i], v2[i]);
    }
  }
  
  *fnedges = nedges - ndiscard;
  
  free(remap1);
  free(remap2);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP merge_vertices_(SEXP x_, SEXP y_, SEXP v1_, SEXP v2_, SEXP tol_, SEXP verbosity_) {
  
  int nprotect = 0;
  
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("merge_vertices_() x & y must be equal length");
  }
  if (length(v1_) == 0) {
    warning("merge_vertices_(): No vertices to process\n");
    return R_NilValue;
  }
  if (length(v1_) != length(v2_)) {
    error("merge_vertices_(): v1 & v2 must be equal length");
  }
  
  double *x = REAL(x_);
  double *y = REAL(y_);
  
  int nverts = length(x_);
  int nedges = length(v1_);
  
  SEXP fv1_ = PROTECT(allocVector(INTSXP, nedges)); nprotect++;
  SEXP fv2_ = PROTECT(allocVector(INTSXP, nedges)); nprotect++;
  int *fv1 = INTEGER(fv1_);
  int *fv2 = INTEGER(fv2_);
  
  int *v1 = INTEGER(v1_);
  int *v2 = INTEGER(v2_);
  
  // Convert from R 1-indexing to C 0-indexing
  for (int i = 0; i < nedges; i++) {
    fv1[i] = v1[i] - 1;
    fv2[i] = v2[i] - 1;
  }
  
  
  int fnedges = 0;
  
  merge_vertices_core_(asReal(tol_), nverts, x, y, nedges, fv1, fv2, &fnedges, asInteger(verbosity_));
  
  trim_vec(fv1_, fnedges, nedges);
  trim_vec(fv2_, fnedges, nedges);
  
  // Convert from C 0-indexing to R 1-indexing
  for (int i = 0; i < fnedges; i++) {
    fv1[i] = fv1[i] + 1;
    fv2[i] = fv2[i] + 1;
  }
  
  
  SEXP res_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  
  SET_VECTOR_ELT(res_, 0, fv1_);
  SET_VECTOR_ELT(res_, 1, fv2_);
  
  SET_STRING_ELT(nms_, 0, mkChar("v1"));
  SET_STRING_ELT(nms_, 1, mkChar("v2"));
  
  setAttrib(res_, R_NamesSymbol, nms_);
  
  set_df_attributes(res_, fnedges, fnedges);
  
  
  UNPROTECT(nprotect);
  return res_;
}



