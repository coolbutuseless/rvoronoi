
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


struct Site *next_site(context_t *ctx);
void init_sites(context_t *ctx, double *x, double *y, int n);


SEXP voronoi_(SEXP x_, SEXP y_) {

  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("x/y lengths aren't valid");
  }
  
  int nprotect = 0;
  context_t ctx = { 0 };
  ctx.triangulate = 0;
  
  ctx.alloc_count    = 0;
  ctx.alloc_capacity = 1024;
  ctx.allocs = (void **)calloc(ctx.alloc_capacity, sizeof(void *));
  
  int n = length(x_);
  SEXP vert_x_ = PROTECT(allocVector(REALSXP, 2 * n)); nprotect++;
  SEXP vert_y_ = PROTECT(allocVector(REALSXP, 2 * n)); nprotect++;
  ctx.vert_idx = 0;
  ctx.vert_x = REAL(vert_x_);
  ctx.vert_y = REAL(vert_y_);
  
  SEXP line_a_ = PROTECT(allocVector(REALSXP, 3 * n)); nprotect++;
  SEXP line_b_ = PROTECT(allocVector(REALSXP, 3 * n)); nprotect++;
  SEXP line_c_ = PROTECT(allocVector(REALSXP, 3 * n)); nprotect++;
  ctx.line_idx = 0;
  ctx.line_a = REAL(line_a_);
  ctx.line_b = REAL(line_b_);
  ctx.line_c = REAL(line_c_);
  
  SEXP seg_line_ = PROTECT(allocVector(INTSXP, 3 * n)); nprotect++;
  SEXP seg_v1_   = PROTECT(allocVector(INTSXP, 3 * n)); nprotect++;
  SEXP seg_v2_   = PROTECT(allocVector(INTSXP, 3 * n)); nprotect++;
  ctx.seg_idx = 0;
  ctx.seg_line = INTEGER(seg_line_);
  ctx.seg_v1   = INTEGER(seg_v1_);
  ctx.seg_v2   = INTEGER(seg_v2_);
    
    
  freeinit(&ctx.sfl, sizeof *ctx.sites);

  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));

  ctx.siteidx = 0;
  geominit(&ctx);
  
  voronoi(&ctx, ctx.triangulate, next_site);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertex
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP vert_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP nms_  = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("x"));
  SET_STRING_ELT(nms_, 1, mkChar("y"));
  
  SET_VECTOR_ELT(vert_, 0, vert_x_);
  SET_VECTOR_ELT(vert_, 1, vert_y_);
  
  setAttrib(vert_, R_NamesSymbol, nms_);
  
  //                       visible len, allocated length
  set_df_attributes(vert_, ctx.vert_idx,           2 * n);
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Line
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP line_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  nms_       = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("a"));
  SET_STRING_ELT(nms_, 1, mkChar("b"));
  SET_STRING_ELT(nms_, 2, mkChar("c"));
  
  SET_VECTOR_ELT(line_, 0, line_a_);
  SET_VECTOR_ELT(line_, 1, line_b_);
  SET_VECTOR_ELT(line_, 2, line_c_);
  
  setAttrib(line_, R_NamesSymbol, nms_);
  
  //                       visible len, allocated length
  set_df_attributes(line_, ctx.line_idx,           3 * n);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segments
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP seg_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  nms_      = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("line"));
  SET_STRING_ELT(nms_, 1, mkChar("v1"));
  SET_STRING_ELT(nms_, 2, mkChar("v2"));
  
  SET_VECTOR_ELT(seg_, 0, seg_line_);
  SET_VECTOR_ELT(seg_, 1, seg_v1_);
  SET_VECTOR_ELT(seg_, 2, seg_v2_);
  
  setAttrib(seg_, R_NamesSymbol, nms_);
  
  //                       visible len, allocated length
  set_df_attributes(seg_, ctx.seg_idx,           3 * n);
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Final list of data.frames
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  nms_      = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  
  SET_STRING_ELT(nms_, 0, mkChar("vertex"));
  SET_STRING_ELT(nms_, 1, mkChar("line"));
  SET_STRING_ELT(nms_, 2, mkChar("segment"));
  
  SET_VECTOR_ELT(res_, 0, vert_);
  SET_VECTOR_ELT(res_, 1, line_);
  SET_VECTOR_ELT(res_, 2, seg_);
  
  setAttrib(res_, R_NamesSymbol, nms_);
  
  
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  
  UNPROTECT(nprotect);
  return res_;
}



SEXP delaunay_(SEXP x_, SEXP y_) {
  
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("x/y lengths aren't valid");
  }
  
  int nprotect = 0;
  context_t ctx = { 0 };
  
  ctx.triangulate = 1;
  
  ctx.alloc_count    = 0;
  ctx.alloc_capacity = 1024;
  ctx.allocs = (void **)calloc(ctx.alloc_capacity, sizeof(void *));
  
  int n = length(x_);
  SEXP v1_ = PROTECT(allocVector(INTSXP, 2 * n)); nprotect++; 
  SEXP v2_ = PROTECT(allocVector(INTSXP, 2 * n)); nprotect++; 
  SEXP v3_ = PROTECT(allocVector(INTSXP, 2 * n)); nprotect++; 
  
  ctx.idel = 0;
  
  ctx.v1 = INTEGER(v1_);
  ctx.v2 = INTEGER(v2_);
  ctx.v3 = INTEGER(v3_);
  
  
  
  freeinit(&ctx.sfl, sizeof *ctx.sites);
  
  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));
  
  ctx.siteidx = 0;
  geominit(&ctx);
  
  voronoi(&ctx, ctx.triangulate, next_site);
  
  
  SEXP res_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("v1"));
  SET_STRING_ELT(nms_, 1, mkChar("v2"));
  SET_STRING_ELT(nms_, 2, mkChar("v3"));
  
  SET_VECTOR_ELT(res_, 0, v1_);
  SET_VECTOR_ELT(res_, 1, v2_);
  SET_VECTOR_ELT(res_, 2, v3_);
  
  setAttrib(res_, R_NamesSymbol, nms_);
  
  //                      visible len, allocated length
  set_df_attributes(res_,    ctx.idel,           2 * n);
  
  
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  UNPROTECT(nprotect);
  return res_;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* sort ctx->sites on y, then x, coord */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int site_comparison(const void *v1, const void *v2) {
  
  struct Point *s1 = (struct Point *)v1;
  struct Point *s2 = (struct Point *)v2;
  
  if (s1->y < s2->y)
    return (-1);
  if (s1->y > s2->y)
    return (1);
  if (s1->x < s2->x)
    return (-1);
  if (s1->x > s2->x)
    return (1);
  return (0);
}

 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* return a single in-storage site */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Site *next_site(context_t *ctx) {
  struct Site *s;
  if (ctx->siteidx < ctx->nsites) {
    s = &ctx->sites[ctx->siteidx];
    ctx->siteidx += 1;
    return (s);
  } else
    return ((struct Site *)NULL);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* read all ctx->sites, sort, and compute xmin, xmax, ymin, ymax */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void init_sites(context_t *ctx, double *x, double *y, int n) {
  int i;
  
  ctx->nsites = n;
  ctx->sites = (struct Site *)myalloc(ctx, ctx->nsites * sizeof *ctx->sites);
  for (i = 0; i < n; i++) {
    ctx->sites[i].coord.x = x[i];
    ctx->sites[i].coord.y = y[i];
    ctx->sites[i].sitenbr = i;
    ctx->sites[i].refcnt  = 0;
  }
  qsort(ctx->sites, ctx->nsites, sizeof *ctx->sites, site_comparison);
  ctx->xmin = ctx->sites[0].coord.x;
  ctx->xmax = ctx->sites[0].coord.x;
  for (i = 1; i < ctx->nsites; i += 1) {
    if (ctx->sites[i].coord.x < ctx->xmin)
      ctx->xmin = ctx->sites[i].coord.x;
    if (ctx->sites[i].coord.x > ctx->xmax)
      ctx->xmax = ctx->sites[i].coord.x;
  };
  ctx->ymin = ctx->sites[0].coord.y;
  ctx->ymax = ctx->sites[ctx->nsites - 1].coord.y;
}
