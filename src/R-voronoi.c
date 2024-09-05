
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
#include "R-extract-polygons.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Helper function used with sort()
// sort ctx->sites on y, then x, coord 
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
// Convert coordinates to ctx->sites, sort, and compute xmin, xmax, ymin, ymax 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void init_sites(context_t *ctx, double *x, double *y, int n) {
  
  ctx->nsites = n;
  ctx->sites = (struct Site *)myalloc(ctx, (unsigned int)ctx->nsites * sizeof *ctx->sites);
  for (int i = 0; i < n; i++) {
    ctx->sites[i].coord.x = x[i];
    ctx->sites[i].coord.y = y[i];
    ctx->sites[i].sitenbr = i;
    ctx->sites[i].refcnt  = 0;
  }
  
  qsort(ctx->sites, (unsigned long)ctx->nsites, sizeof *ctx->sites, site_comparison);
  ctx->xmin = ctx->sites[0].coord.x;
  ctx->xmax = ctx->sites[0].coord.x;
  
  for (int i = 1; i < ctx->nsites; i++) {
    if (ctx->sites[i].coord.x < ctx->xmin)
      ctx->xmin = ctx->sites[i].coord.x;
    if (ctx->sites[i].coord.x > ctx->xmax)
      ctx->xmax = ctx->sites[i].coord.x;
  };
  
  ctx->ymin = ctx->sites[0].coord.y;
  ctx->ymax = ctx->sites[ctx->nsites - 1].coord.y;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Voronoi Tesselation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP voronoi_(SEXP x_, SEXP y_) {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity Check
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("x/y lengths aren't valid");
  }
  
  int nprotect = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialise the calculation context
  // Do delauney? FALSE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  context_t ctx = { 0 };
  ctx.triangulate = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Track *ALL* the allocations done via 'myalloc()'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctx.alloc_count    = 0;
  ctx.alloc_capacity = 1024;
  ctx.allocs = (void **)calloc((unsigned long)ctx.alloc_capacity, sizeof(void *));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Space for vertices
  // If there are 'n' points then 
  // maximum number of vertices in the voronoi tesselation is 2n - 5
  // maximum number of edges is 3n - 6
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n = length(x_);
  int max_verts = 2 * n;
  int max_edges = 3 * n;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Voronoi Vertices (x, y)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP vert_x_ = PROTECT(allocVector(REALSXP, max_verts)); nprotect++;
  SEXP vert_y_ = PROTECT(allocVector(REALSXP, max_verts)); nprotect++;
  ctx.nverts = 0;
  ctx.vert_x = REAL(vert_x_);
  ctx.vert_y = REAL(vert_y_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Voronoi lines (a, b, c) => ax + by = c
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP line_a_ = PROTECT(allocVector(REALSXP, max_edges)); nprotect++;
  SEXP line_b_ = PROTECT(allocVector(REALSXP, max_edges)); nprotect++;
  SEXP line_c_ = PROTECT(allocVector(REALSXP, max_edges)); nprotect++;
  ctx.nlines = 0;
  ctx.line_a = REAL(line_a_);
  ctx.line_b = REAL(line_b_);
  ctx.line_c = REAL(line_c_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // The segments within each line
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP seg_line_ = PROTECT(allocVector(INTSXP, max_edges)); nprotect++;
  SEXP seg_v1_   = PROTECT(allocVector(INTSXP, max_edges)); nprotect++;
  SEXP seg_v2_   = PROTECT(allocVector(INTSXP, max_edges)); nprotect++;
  ctx.nsegs = 0;
  ctx.seg_line = INTEGER(seg_line_);
  ctx.seg_v1   = INTEGER(seg_v1_);
  ctx.seg_v2   = INTEGER(seg_v2_);
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate voronoi tesselation using Fortune's Sweep algorithm 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  freeinit(&ctx.sfl, sizeof *ctx.sites);
  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));
  ctx.siteidx = 0;
  geominit(&ctx);
  voronoi(&ctx, ctx.triangulate);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Merge vertices in voronoi segments to prep for polygon building
  // void merge_vertices_core_(double tol, 
  //                           int nverts, double *x, double *y, 
  //                           int nedges, int *v1, int *v2, 
  //                           int *fnedges,
  //                           int verbosity)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int fnedges = 0; // Final number of edges after merging

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Set up the "working area" for the merging of vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *v1 = malloc((unsigned long)ctx.nsegs * sizeof(int));
  int *v2 = malloc((unsigned long)ctx.nsegs * sizeof(int));
  if (v1 == NULL || v2 == NULL) error("voronoi_(): Couldn't allocate for 'v1' and 'v2'");

  for (int i = 0; i < ctx.nsegs; i++) {
    v1[i] = ctx.seg_v1[i]; 
    v2[i] = ctx.seg_v2[i];
  }

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Merge close vertices which are an artefact of the tessellation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  merge_vertices_core_(1e-5, ctx.nverts, ctx.vert_x, ctx.vert_y,
                       ctx.nsegs, v1, v2, &fnedges, 0);

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extract the polygons using the (temporary) merged vertices
  // 
  // SEXP extract_polygons_core_(int vert_n, double *vert_x, double *vert_y,
  //                                      int seg_n, int *seg_v1, int *seg_v2)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP polys_ = PROTECT(
    extract_polygons_core_(ctx.nverts, ctx.vert_x, ctx.vert_y, fnedges, v1, v2)
  ); nprotect++;

  free(v1);
  free(v2);
  
  
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertices:   data.frame(x = ..., y = ...)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP vert_ = PROTECT(allocVector(VECSXP, 2)); nprotect++;
  SEXP nms_  = PROTECT(allocVector(STRSXP, 2)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("x"));
  SET_STRING_ELT(nms_, 1, mkChar("y"));
  
  SET_VECTOR_ELT(vert_, 0, vert_x_);
  SET_VECTOR_ELT(vert_, 1, vert_y_);
  
  setAttrib(vert_, R_NamesSymbol, nms_);
  set_df_attributes(vert_, ctx.nverts, max_verts);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Lines: data.frame(a = ..., b = ..., c = ...)
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
  set_df_attributes(line_, ctx.nlines, max_edges);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segments: data.frame(line = integer(), v1 = integer(), v2 = integer())
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
  set_df_attributes(seg_, ctx.nsegs, max_edges);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extents: list(xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric())
  // Covers all seed points and voronoi vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double xmin =  INFINITY;
  double xmax = -INFINITY;
  double ymin =  INFINITY;
  double ymax = -INFINITY;
  for (int i = 0; i < ctx.nverts; i++) {
    if (ctx.vert_x[i] > xmax) xmax = ctx.vert_x[i];
    if (ctx.vert_x[i] < xmin) xmin = ctx.vert_x[i];
    if (ctx.vert_y[i] > ymax) ymax = ctx.vert_y[i];
    if (ctx.vert_y[i] < ymin) ymin = ctx.vert_y[i];
  }
  
  double *x = REAL(x_);
  double *y = REAL(y_);
  for (int i = 0; i < length(x_); i++) {
    if (x[i] > xmax) xmax = x[i];
    if (x[i] < xmin) xmin = x[i];
    if (y[i] > ymax) ymax = y[i];
    if (y[i] < ymin) ymin = y[i];
  }
  
  SEXP xmin_ = PROTECT(ScalarReal(xmin)); nprotect++;
  SEXP xmax_ = PROTECT(ScalarReal(xmax)); nprotect++;
  SEXP ymin_ = PROTECT(ScalarReal(ymin)); nprotect++;
  SEXP ymax_ = PROTECT(ScalarReal(ymax)); nprotect++;
  
  SEXP ext_ = PROTECT(allocVector(VECSXP, 4)); nprotect++;
  nms_      = PROTECT(allocVector(STRSXP, 4)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("xmin"));
  SET_STRING_ELT(nms_, 1, mkChar("xmax"));
  SET_STRING_ELT(nms_, 2, mkChar("ymin"));
  SET_STRING_ELT(nms_, 3, mkChar("ymax"));
  
  SET_VECTOR_ELT(ext_, 0, xmin_);
  SET_VECTOR_ELT(ext_, 1, xmax_);
  SET_VECTOR_ELT(ext_, 2, ymin_);
  SET_VECTOR_ELT(ext_, 3, ymax_);
  
  setAttrib(ext_, R_NamesSymbol, nms_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Final result: named list of data.frames
  //  list(
  //     vertex  = data.frame()
  //     line    = data.frame()
  //     segment = data.frame()
  //     extents = list()
  // )
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, 5)); nprotect++;
  nms_      = PROTECT(allocVector(STRSXP, 5)); nprotect++;
  
  SET_STRING_ELT(nms_, 0, mkChar("vertex"));
  SET_STRING_ELT(nms_, 1, mkChar("line"));
  SET_STRING_ELT(nms_, 2, mkChar("segment"));
  SET_STRING_ELT(nms_, 3, mkChar("extents"));
  SET_STRING_ELT(nms_, 4, mkChar("polygons"));
  
  SET_VECTOR_ELT(res_, 0, vert_);
  SET_VECTOR_ELT(res_, 1, line_);
  SET_VECTOR_ELT(res_, 2, seg_);
  SET_VECTOR_ELT(res_, 3, ext_);
  SET_VECTOR_ELT(res_, 4, polys_);
  
  setAttrib(res_, R_NamesSymbol, nms_);
  setAttrib(res_, R_ClassSymbol, mkString("vor"));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' memory
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *seg_line = INTEGER(seg_line_);
  int *seg_v1   = INTEGER(seg_v1_);
  int *seg_v2   = INTEGER(seg_v2_);
  for (int i = 0; i < length(seg_line_); i++) seg_line[i]++;
  for (int i = 0; i < length(seg_v1_  ); i++) seg_v1  [i]++;
  for (int i = 0; i < length(seg_v2_  ); i++) seg_v2  [i]++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}



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
  // Final result = data.frame(v1 = integer(), v2 = integer(), v3 = integer())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, 3)); nprotect++;
  SEXP nms_ = PROTECT(allocVector(STRSXP, 3)); nprotect++;
  SET_STRING_ELT(nms_, 0, mkChar("v1"));
  SET_STRING_ELT(nms_, 1, mkChar("v2"));
  SET_STRING_ELT(nms_, 2, mkChar("v3"));
  
  SET_VECTOR_ELT(res_, 0, v1_);
  SET_VECTOR_ELT(res_, 1, v2_);
  SET_VECTOR_ELT(res_, 2, v3_);
  
  setAttrib(res_, R_NamesSymbol, nms_);
  set_df_attributes(res_, ctx.ntris, max_tris);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' allocations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *v1 = INTEGER(v1_);
  int *v2 = INTEGER(v2_);
  int *v3 = INTEGER(v3_);
  for (int i = 0; i < ctx.ntris; i++) v1[i]++;
  for (int i = 0; i < ctx.ntris; i++) v2[i]++;
  for (int i = 0; i < ctx.ntris; i++) v3[i]++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}





 

