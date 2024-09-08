

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
#include "utils-bbox.h"


void calc_space_for_bound_infinite_edges(int nsegs, int *v1, int *v2, int *nbverts, int *nbsegs) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Count the number of unbound rays and setup storage
  //
  // TODO: add corner vertices and perimeter segments
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  *nbsegs  = 0;
  *nbverts = 0;
  for (int i = 0; i < nsegs; i++) {
    if (v1[i] >= 0 && v2[i] >= 0) continue;
    (*nbsegs)++;
    (*nbverts)++;
    if (v1[i] < 0 && v2[i] < 0) {
      // double unbounded edge adds 2 vertices
      (*nbverts)++;
    }
    
  }
  
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Find infinte edges and intersect them with the boundary
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void bound_infinite_edges(
    bbox_t *bounds,
    int nverts, double *x, double *y, 
    int nsegs , int *li, int *v1, int *v2,
    int nlines, double *a, double *b, double *c,
    int nbverts, double *xb, double *yb,
    int nbsegs, int *rv1, int *rv2) {
  
  
  Rprintf("Found %i rays to bound. Needing %i vertices\n", nbsegs, nbverts);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Hold info about the two intercepts between rays and rectangle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  typedef struct {
    double x;
    double y;
    bool valid;
  } intercept_t;
  
  int ray_idx  = 0;
  int vert_idx = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // For each segment
  //    - find the two intercepts with the bounding rectangle
  //    - add the vertices and segments to the (xb, yb) and (rv1, rv2)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < nsegs; i++) {
    if (v1[i] >= 0 && v2[i] >= 0) continue;
    
    
    int type = 0;
    int v1_anchor = 0;
    double xi = 0;
    double yi = 0;
    int line_idx = 0;
    
    if (v1[i] >= 0) {
      type = 2;
      v1_anchor = v1[i];
    } else if (v2[i] >= 0) {
      type = 1;
      v1_anchor = v2[i];
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // There can only be 2 intercept points which lie on the boundary rectangle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    intercept_t intercept[2] = { 0 };
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intersection with left edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    line_idx = li[i];
    xi = bounds->xmin;
    yi = (c[line_idx]- a[line_idx] * xi) / b[line_idx];
    if (yi >= bounds->ymin && yi <= bounds->ymax) {
      intercept[0].valid = true;
      intercept[0].x     = xi;
      intercept[0].y     = yi;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intercept with right edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xi = bounds->xmax;
    yi = (c[line_idx]- a[line_idx] * xi) / b[line_idx];
    if (yi >= bounds->ymin && yi <= bounds->ymax) {
      intercept[1].valid = true;
      intercept[1].x     = xi;
      intercept[1].y     = yi;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intercept with lower edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double yi1 = bounds->ymin;
    double xi1 = (c[line_idx]- b[line_idx] * yi1) / a[line_idx];
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intercept with upper edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double yi2 = bounds->ymax;
    double xi2 = (c[line_idx]- b[line_idx] * yi2) / a[line_idx];
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // If there's not a left intercept, then pick the left-most out of upper/lower
    // intercepts
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!intercept[0].valid) {
      if (xi1 < xi2) {
        intercept[0].x = xi1;
        intercept[0].y = yi1;
      } else {
        intercept[0].x = xi2;
        intercept[0].y = yi2;
      }
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // If there's not a right intercept, then pick the right-most out of upper/lower
    // intercepts
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!intercept[1].valid) {
      if (xi2 > xi1) {
        intercept[1].x = xi2;
        intercept[1].y = yi2;
      } else {
        intercept[1].x = xi1;
        intercept[1].y = yi1;
      }
    }
    
    Rprintf("Ray: %i   Vert: %i\n", ray_idx, vert_idx);
    if (type == 0) {
      // double ended unbounded ray. 
      // Add vertices
      xb[vert_idx] = intercept[0].x;
      yb[vert_idx] = intercept[0].y;
      vert_idx++;
      xb[vert_idx] = intercept[1].x;
      yb[vert_idx] = intercept[1].y;
      vert_idx++;
      // Add segment
      rv1[ray_idx] = nverts + vert_idx - 2;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else if (type == 1) {
      // ray now bounded on the left
      // Add vertex
      xb[vert_idx] = intercept[0].x;
      yb[vert_idx] = intercept[0].y;
      vert_idx++;
      // Add segment
      rv1[ray_idx] = v1_anchor;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else if (type == 2) {
      // ray not bounded on the right
      // Add vertex
      xb[vert_idx] = intercept[1].x;
      yb[vert_idx] = intercept[1].y;
      vert_idx++;
      // Add segment
      rv1[ray_idx] = v1_anchor;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else {
      error("Impossible 12");
    }
    
    
    
    
  }
  
  
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP bound_infinite_edges_(
    SEXP xmin_, SEXP ymin_, SEXP xmax_, SEXP ymax_,
    SEXP x_, SEXP y_,
    SEXP a_, SEXP b_, SEXP c_,
    SEXP li_, SEXP v1_, SEXP v2_) {
  
  int nprotect = 0;
  
  if (length(x_) != length(y_)) {
    error("bound_infinite_edges_(): bad length for x/y");
  }
  
  if (length(v1_) == 0 || length(v1_) != length(v2_) || 
      length(li_) != length(v1_)) {
    error("bound_infinite_edges_(): bad length for li/v1/v2");
  }
  
  if (length(a_) == 0 || length(a_) != length(b_) || 
      length(a_) != length(c_)) {
    error("bound_infinite_edges_(): bad length for a/b/c");
  }
  
  bbox_t bounds = {
    .xmin = asReal(xmin_),
    .ymin = asReal(ymin_),
    .xmax = asReal(xmax_),
    .ymax = asReal(ymax_)
  };
  
  int nbverts = 0;
  int nbsegs = 0;
  
  calc_space_for_bound_infinite_edges(length(v1_), INTEGER(v1_), INTEGER(v2_), &nbverts, &nbsegs);
  SEXP xb_ = PROTECT(allocVector(REALSXP, nbverts)); nprotect++;
  SEXP yb_ = PROTECT(allocVector(REALSXP, nbverts)); nprotect++;
  
  SEXP rv1_ = PROTECT(allocVector(INTSXP, nbsegs)); nprotect++;
  SEXP rv2_ = PROTECT(allocVector(INTSXP, nbsegs)); nprotect++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert R 1-indexing to C 0-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *v1 = malloc(length(v1_) * sizeof(int));
  int *v2 = malloc(length(v2_) * sizeof(int));
  int *li = malloc(length(li_) * sizeof(int));
  if (v1 == NULL || v2 == NULL || li == NULL) {
    error("bound_infinite_edges_(): li/v1/v2 allocation failed");
  }
  
  int *v1p = INTEGER(v1_);
  int *v2p = INTEGER(v2_);
  int *lip = INTEGER(li_);
  for (int i = 0; i < length(v1_); i++) v1[i] = v1p[i] - 1;
  for (int i = 0; i < length(v1_); i++) v2[i] = v2p[i] - 1;
  for (int i = 0; i < length(v1_); i++) li[i] = lip[i] - 1;

  
  bound_infinite_edges(
    &bounds,
    length(x_), REAL(x_), REAL(y_),
    length(v1_), li, v1, v2,
    length(a_), REAL(a_), REAL(b_), REAL(c_),
    nbverts, REAL(xb_), REAL(yb_),
    nbsegs, INTEGER(rv1_), INTEGER(rv2_)
  );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP verts_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP vnms_  = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(vnms_, 0, mkChar("x"));
  SET_STRING_ELT(vnms_, 1, mkChar("y"));
  setAttrib(verts_, R_NamesSymbol, vnms_);

  SET_VECTOR_ELT(verts_, 0, xb_);
  SET_VECTOR_ELT(verts_, 1, yb_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *rv1 = INTEGER(rv1_);
  int *rv2 = INTEGER(rv2_);
  for (int i = 0; i < nbsegs; i++) rv1[i]++; // convert to R 1-indexing
  for (int i = 0; i < nbsegs; i++) rv2[i]++;
  
  
  SEXP segments_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  SEXP snms_     = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(snms_, 0, mkChar("line"));
  SET_STRING_ELT(snms_, 1, mkChar("v1"));
  SET_STRING_ELT(snms_, 2, mkChar("v2"));
  setAttrib(segments_, R_NamesSymbol, snms_);
  
  SET_VECTOR_ELT(segments_, 0, R_NilValue);
  SET_VECTOR_ELT(segments_, 1, rv1_);
  SET_VECTOR_ELT(segments_, 2, rv2_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // list result
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("vertex"));
  SET_STRING_ELT(nms_, 1, mkChar("segment"));
  setAttrib(res_, R_NamesSymbol, nms_);
  
  SET_VECTOR_ELT(res_, 0, verts_);
  SET_VECTOR_ELT(res_, 1, segments_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(v1);
  free(v2);
  
  UNPROTECT(nprotect);
  return res_;
}



