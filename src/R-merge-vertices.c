
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


void remap_insert(int *remap, int i1, int i2) {
  
  // Swap so i1 is always larger
  if (i1 < i2) {
    int tmp = i1;
    i1 = i2;
    i2 = tmp;
  }
  
  // New remapping
  if (remap[i1] < 0) {
    remap[i1] = i2;
    return;
  } else if (remap[i1] == i2) {
    return;
  } else {
    // i1 is already remapped to something!
    if (remap[i1] > i2) {
      remap_insert(remap, remap[i1], i2);
      remap[i1] = i2;
    } else {
      remap_insert(remap, remap[i1], i2);
    }
  }
  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Caller is responsible for freeing fv1, fv2
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void merge_vertices_core_(double tol, 
                          int nverts, double *x, double *y, 
                          int nedges, int *line, int *v1, int *v2, 
                          int *fnedges,
                          int verbosity) {
  
  
  int *remap = malloc((unsigned long)nverts * sizeof(int));
  if (remap == NULL) {
    error("merge_vertices_core_(): could not allocate 'remap'");
  }
  for (int i = 0; i < nverts; i++) {
    remap[i] = -1;
  }
  
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
    
    if (verbosity > 0) {
      Rprintf("Edge: %i (%i -> %i): Dist: %.6f\n", i, i1, i2, dist);
    }
    
    if (dist < tol) {
      if (verbosity > 0) Rprintf("^ Remap\n");
      remap_insert(remap, i1, i2);
    } 
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Collapse and remove any edges with duplicated vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Print remaps
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {
    for (int i = 0; i < nverts; i++) {
      if (remap[i] >= 0) {
        Rprintf("Remap: %i => %i\n", i, remap[i]);
      }
    }
  }
  
  // Track which indices we want to keep
  int ndiscard = 0;
  bool *discard = calloc((unsigned int)nedges, sizeof(bool));
  if (discard == NULL) error("merge_vertices_core_(): could not allocate 'discard'");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Renumber verts
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = nverts - 1; i >= 0; i--) {
    if (remap[i] < 0) continue;
    int vlo = remap[i];
    int vhi = i;
    
    // Rprintf("REMAP: %i -> %i\n", vhi, vlo);
    
    for (int j = 0; j < nedges; j++) {
      if (v1[j] == vhi) v1[j] = vlo; 
      if (v2[j] == vhi) v2[j] = vlo; 
      discard[j] = v1[j] == v2[j];
    }
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Mark duplicate double ended rays for removal
  //  Each double ended ray ends up appearing *twice* in the segment list
  //  so we can delete
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int line_used[nedges];
  // Rprintf("Checking for double-ended ray\n");
  for (int i = 0; i < nedges; i++) line_used[i] = 0;
  for (int i = 0; i < nedges; i++) {
    if (v1[i] >= 0 || v2[i] >= 0) continue;
    if (line_used[line[i]] == 1) {
      // Rprintf("Discard double-ended ray at segment: %i", i);
      discard[i] = 1;
    }
    line_used[line[i]] = 1;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tally the discarded edges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < nedges; i++) {
    ndiscard += discard[i];
  }
  
  if (verbosity > 0) {
    Rprintf("Ndistcard: %i\n", ndiscard);
    for (int i = 0; i < nedges; i++) {
      Rprintf("(%4i)  %4i\n", i, discard[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Collapse and remove any edges with duplicated vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {  
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
    }
  }
  
  
  int src = 0;
  int dst = 0;
  for (dst = 0; dst < nedges - ndiscard; dst++, src++) {
    while(discard[src]) {
      src++;
    }
    // src += discard[src];
    v1[dst]   = v1[src];
    v2[dst]   = v2[src];
    line[dst] = line[src];
    if (verbosity > 1) {
      Rprintf("src %4i  ->  dst %4i\n", src, dst);
    }
  }
  
  // Mark discarded edges with a sentinel value so it's easy to see
  // when something's not right.
  for (; dst < nedges; dst++) {
    v1[dst]   = 998;
    v2[dst]   = 998;
    line[dst] = 998;
  }
  
  if (verbosity > 0) {  
    for (int i = 0; i < nedges; i++) {
      Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
    }
  }
  
  *fnedges = nedges - ndiscard;
  
  free(discard);
  free(remap);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP merge_vertices_(SEXP x_, SEXP y_, 
                     SEXP line_, SEXP v1_, SEXP v2_, 
                     SEXP tol_, SEXP verbosity_) {
  
  int nprotect = 0;
  
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("merge_vertices_() x & y must be equal length");
  }
  if (length(v1_) == 0) {
    warning("merge_vertices_(): No vertices to process\n");
    return R_NilValue;
  }
  if (length(v1_) != length(v2_) || length(v1_) != length(line_)) {
    error("merge_vertices_():line, v1 & v2 must be equal length");
  }
  
  double *x = REAL(x_);
  double *y = REAL(y_);
  
  int nverts = length(x_);
  int nedges = length(v1_);
  
  SEXP fv1_   = PROTECT(allocVector(INTSXP, nedges)); nprotect++;
  SEXP fv2_   = PROTECT(allocVector(INTSXP, nedges)); nprotect++;
  SEXP fline_ = PROTECT(allocVector(INTSXP, nedges)); nprotect++;
  
  // int *fv1   = INTEGER(fv1_);
  // int *fv2   = INTEGER(fv2_);
  // int *fline = INTEGER(fline_);
  
  // These are initially copies of the input vectors 
  memcpy(INTEGER(fv1_)  , INTEGER(v1_)  , nedges * sizeof(int));
  memcpy(INTEGER(fv2_)  , INTEGER(v2_)  , nedges * sizeof(int));
  memcpy(INTEGER(fline_), INTEGER(line_), nedges * sizeof(int));
  
  // Convert them to C 0-indexing for processing
  convert_indexing_r_to_c(fv1_);
  convert_indexing_r_to_c(fv2_);
  convert_indexing_r_to_c(fline_);

  // for (int i = 0; i < length(fv2_); i++) {
  //   Rprintf("(%i) v2[%i] = %i    =   fv2[%i] = %i\n", nedges, i, INTEGER(v2_)[i], i, INTEGER(fv2_)[i]);
  // }
  
  // This will contain the final number of vertices after merging  
  int fnedges = 0;
  
  merge_vertices_core_(asReal(tol_), 
                       nverts, x, y, 
                       nedges, INTEGER(fline_), INTEGER(fv1_), INTEGER(fv2_), 
                       &fnedges, asInteger(verbosity_));
  
  trim_vec(fv1_  , fnedges, nedges);
  trim_vec(fv2_  , fnedges, nedges);
  trim_vec(fline_, fnedges, nedges);
  
  convert_indexing_c_to_r(fv1_);
  convert_indexing_c_to_r(fv2_);
  convert_indexing_c_to_r(fline_);
  
  
  SEXP res_ = PROTECT(
    create_named_list(3, "line", fline_, "v1", fv1_, "v2", fv2_)
  ); nprotect++;
  set_df_attributes(res_, fnedges, fnedges);
  
  
  UNPROTECT(nprotect);
  return res_;
}



