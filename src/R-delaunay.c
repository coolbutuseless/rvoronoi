
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
#include "R-common.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Delauney Triangulation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP delaunay_(SEXP x_, SEXP y_) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity Check
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("x/y lengths aren't valid");
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialise the calculation context
  // Do delauney? TRUE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int nprotect = 0;
  context_t ctx = { 0 };
  ctx.triangulate = 1;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Track *ALL* the allocations done via 'myalloc()'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctx.alloc_count    = 0;
  ctx.alloc_capacity = 1024;
  ctx.allocs = (void **)calloc((unsigned long)ctx.alloc_capacity, sizeof(void *));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Maximum number of triangles = 2 * n - 2 - b
  // where 'b' is number of vertices on the convex hull
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n = length(x_);
  int max_tris = 2 * n;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertices of each triangles
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP v1_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  SEXP v2_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  SEXP v3_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize context
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctx.ntris = 0;
  ctx.v1 = INTEGER(v1_);
  ctx.v2 = INTEGER(v2_);
  ctx.v3 = INTEGER(v3_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Use Fortune's algo to calculate delaunay triangulation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  freeinit(&ctx.sfl, sizeof *ctx.sites);
  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));
  ctx.siteidx = 0;
  geominit(&ctx);
  voronoi(&ctx, ctx.triangulate);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Indices = data.frame(v1 = integer(), v2 = integer(), v3 = integer())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // const char *idx_nms[3] = {"v1", "v2", "v3"};
  // SEXP idx_ = PROTECT(create_named_list(3, idx_nms)); nprotect++;
  
  SEXP idx_ = PROTECT(create_named_list(3, "v1", v1_, "v2", v2_, "v3", v3_)); nprotect++;
  set_df_attributes(idx_, ctx.ntris, max_tris);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertex coordinates for each triangle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP pidx_ = PROTECT(allocVector(INTSXP , 3 * ctx.ntris)); nprotect++;
  SEXP xs_   = PROTECT(allocVector(REALSXP, 3 * ctx.ntris)); nprotect++;
  SEXP ys_   = PROTECT(allocVector(REALSXP, 3 * ctx.ntris)); nprotect++;
  double *xs = REAL(xs_);
  double *ys = REAL(ys_);
  int *pidx  = INTEGER(pidx_);
  
  double *x = REAL(x_);
  double *y = REAL(y_);
  int *v1 = INTEGER(v1_);
  int *v2 = INTEGER(v2_);
  int *v3 = INTEGER(v3_);
  
  for (int i = 0; i < ctx.ntris; i++) {
    
    *xs++ = x[ v1[i] ];
    *xs++ = x[ v2[i] ];
    *xs++ = x[ v3[i] ];
    
    *ys++ = y[ v1[i] ];
    *ys++ = y[ v2[i] ];
    *ys++ = y[ v3[i] ];
    
    *pidx++ = i + 1;
    *pidx++ = i + 1;
    *pidx++ = i + 1;
  }
  
  SEXP polys_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  SEXP pnms_  = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(pnms_, 0, mkChar("idx"));
  SET_STRING_ELT(pnms_, 1, mkChar("x"));
  SET_STRING_ELT(pnms_, 2, mkChar("y"));
  setAttrib(polys_, R_NamesSymbol, pnms_);
  
  SET_VECTOR_ELT(polys_, 0, pidx_);
  SET_VECTOR_ELT(polys_, 1, xs_);
  SET_VECTOR_ELT(polys_, 2, ys_);
  
  set_df_attributes(polys_, 3 * ctx.ntris, 3 * ctx.ntris);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // List of delaunay results: list(vertex = data.frame(), coords = data.frame())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP lnm_ = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(lnm_, 0, mkChar("segment")); 
  SET_STRING_ELT(lnm_, 1, mkChar("polygon"));
  setAttrib(res_, R_NamesSymbol, lnm_);
  
  SET_VECTOR_ELT(res_, 0, idx_);
  SET_VECTOR_ELT(res_, 1, polys_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' allocations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
  convert_indexing_c_to_r(v1_);
  convert_indexing_c_to_r(v2_);
  convert_indexing_c_to_r(v3_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}





 

