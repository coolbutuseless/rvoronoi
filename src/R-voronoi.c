
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
#include "R-common.h"
#include "R-merge-vertices.h"
#include "R-extract-polygons.h"
#include "R-infinite-edges.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Voronoi Tesselation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP voronoi_(SEXP x_, SEXP y_, SEXP match_polygons_) {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity Check
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) <= 1) {
    error("Must be at least 2 seed points");
  }
  
  if (length(x_) != length(y_)) {
    error("x & y are not the same length");
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
  // bbox for bounding the unbounded polygos
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bbox_t bounds = bbox_new();
  bbox_add(&bounds, length(x_), REAL(x_), REAL(y_));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Space for vertices
  // If there are 'n' points then 
  //   maximum number of vertices in the voronoi tesselation is 2n - 5
  //   maximum number of edges is 3n - 6
  //
  // But the above result does not account for infinite rays.
  //   If there are 'n' points, then the maximum number of infinite rays = 'n'
  //   If these infinite rays are being clipped by intersection with a rectangular
  //     boundary, then
  //     - there are 4 extra vertices for the corners of the boundary
  //     - there is 1 extra vertex where the ray intersects the boundary.
  //     - for the special case of a double-ended infinite ray, then there
  //       are two new vertices created.
  //     - for a sequence of collinear points, there will be multiple
  //       double ended rays.  For 'n' collinear points there is a maximum
  //       of 'n-1' double ended rays, which would mean '2n-2' boundary 
  //       intersection vertices
  //
  //  The maximum number of boundary vertices = 2n - 2
  //
  // Boundary intersections mean that there are now segments along the 
  // perimeter with which to make new bounding polygons for these exterior points
  //   - There are 4 segments needed to represent the sides of the boundary
  //   - Each boundary intersection adds another segment.
  //
  // The maximum number of boundary segments 
  //   = 4 + (num boundary intersections) 
  //   = 4 + (2n - 2)
  //   = 2n + 2
  //
  // I'm going to allocate maximum space both vertices and edges 
  // Note: it is impossible to have a voronoi where the maximum interior 
  //        vertices *and* maximum exterior vertices simultaneously.
  //        A better analysis here would give me a better upper bound, 
  //        but I don't need a tight bound at the moment as size of
  //        voronois is not expected to be must about 1000 seed points.
  //
  // Add '10' just to avoid off-by-one errors in my thinking.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n = length(x_); // number of seed points
  
  int max_interior_verts = 2 * n - 5;
  int max_interior_edges = 3 * n - 6;
  
  int max_exterior_verts = 2 * n + 4;
  int max_exterior_edges = 2 * n + 4;
  
  int max_verts = max_interior_verts + 10; 
  int max_edges = max_interior_edges + 10; 
  
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
  // Add voronoi vertices to bounding box. Expand by 10%
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bbox_add(&bounds, ctx.nverts, ctx.vert_x, ctx.vert_y);
  bbox_expand(&bounds, 0.10);
  
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
  SEXP v1m_   = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
  SEXP v2m_   = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
  SEXP linem_ = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
  
  int *v1m   = INTEGER(v1m_);
  int *v2m   = INTEGER(v2m_);
  int *linem = INTEGER(linem_);
  
  memcpy(v1m  , ctx.seg_v1  , ctx.nsegs * sizeof(int));
  memcpy(v2m  , ctx.seg_v2  , ctx.nsegs * sizeof(int));
  memcpy(linem, ctx.seg_line, ctx.nsegs * sizeof(int));

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Merge close vertices which are an artefact of the tessellation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  merge_vertices_core_(1e-10, 
                       ctx.nverts, ctx.vert_x, ctx.vert_y,
                       ctx.nsegs, linem, v1m, v2m, 
                       &fnedges, 0);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate the extra vertices and segments produced if we 
  // bound the infinite rays
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define DOBOUNDING
  
#ifdef DOBOUNDING
  int nbverts = 0;
  int nbsegs  = 0;
  
  // Use the merged segments to find infinite rays (single and double ended)
  // and then calculate the actual number of vertices and segments needed
  // to store the bounded external polygons
  calc_space_for_bound_infinite_edges(fnedges, v1m, v2m, &nbverts, &nbsegs);
  
  // Rprintf("DoBounding: fnedges = %i, nbverts = %i,  nbsegs = %i\n", 
          // fnedges, nbverts, nbsegs);
  
  // void bound_infinite_edges(
  //     bbox_t *bounds,
  //     int nverts, double *x, double *y,
  //     int nsegs , int *li, int *v1, int *v2,
  //     int nlines, double *a, double *b, double *c,
  //     int nbverts, double *xb, double *yb,
  //     int nbsegs, int *rv1, int *rv2);
  
  double *xb = calloc(nbverts , sizeof(double));
  double *yb = calloc(nbverts , sizeof(double));
  int *rv1   = calloc(nbsegs  , sizeof(int));
  int *rv2   = calloc(nbsegs  , sizeof(int));
  
  bound_infinite_edges(
    &bounds,
    ctx.nverts, ctx.vert_x, ctx.vert_y,
    fnedges, linem, v1m, v2m,
    ctx.nlines, ctx.line_a, ctx.line_b, ctx.line_c,
    &nbverts, xb, yb,
    &nbsegs, rv1, rv2
  );
#endif
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Create a new MEGA vertex list by concatenting the 
  //   voronoi vertices and the exterior vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef DOBOUNDING
  SEXP xf_ = PROTECT(allocVector(REALSXP, ctx.nverts + nbverts)); nprotect++;
  SEXP yf_ = PROTECT(allocVector(REALSXP, ctx.nverts + nbverts)); nprotect++;
  double *xf = REAL(xf_);
  double *yf = REAL(yf_);
  
  memcpy(xf + 0         , ctx.vert_x, ctx.nverts * sizeof(double));
  memcpy(yf + 0         , ctx.vert_y, ctx.nverts * sizeof(double));
  memcpy(xf + ctx.nverts,         xb,    nbverts * sizeof(double));
  memcpy(yf + ctx.nverts,         yb,    nbverts * sizeof(double));
  
  SEXP exterior_ = PROTECT(create_named_list(2, "x", xf_, "y", yf_)); nprotect++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Append the new exterior segments to the merged vertex list
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memcpy(v1m   + fnedges, rv1, nbsegs * sizeof(int));
  memcpy(v2m   + fnedges, rv2, nbsegs * sizeof(int));
  for (int i = 0; i < nbsegs; i++) {
    (linem + fnedges)[i] = -99;
  }
  // memset(linem + fnedges, -99, nbsegs * sizeof(int));

  free(xb);
  free(yb);
  free(rv1);
  free(rv2);
  fnedges += nbsegs;
#else 
  SEXP exterior_ = R_NilValue;
#endif
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extract the polygons using the (temporary) merged vertices
  // 
  // SEXP extract_polygons_core(
  //    int vert_n, double *vert_x, double *vert_y,  // The voronoi vertices
  //    int seg_n, int *seg_v1, int *seg_v2,         // Voronoi edges
  //    int seed_n, double *seed_x, double *seed_y   // voronoi seed points
  //  )
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP polys_ = R_NilValue;
  
#ifdef DOBOUNDING
  if (asLogical(match_polygons_)) {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts + nbverts, xf, yf,       // Voronoi vertices + perimeter intersections
        fnedges, v1m, v2m,                  // Voronoi edges + perimeter edges
        length(x_), REAL(x_), REAL(y_)      // Seed points
      )
    ); nprotect++;
  } else {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts + nbverts, xf, yf,    // Voronoi vertices + perimeter intersections
        fnedges, v1m, v2m,                  // Voronoi edges
        0, NULL, NULL                       // Seed points
      )
    ); nprotect++;
  }
#else 
  if (asLogical(match_polygons_)) {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts, ctx.vert_x, ctx.vert_y, // Voronoi vertices only
        fnedges, v1m, v2m,                  // Voronoi edges + perimeter edges
        length(x_), REAL(x_), REAL(y_)      // Seed points
      )
    ); nprotect++;
  } else {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts, ctx.vert_x, ctx.vert_y, // Voronoi vertices only
        fnedges, v1m, v2m,                  // Voronoi edges
        0, NULL, NULL                       // Seed points
      )
    ); nprotect++;
  } 
#endif
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Trim the merged indices to size 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP msegs_ = PROTECT(
    create_named_list(3, "line", linem_, "v1", v1m_, "v2", v2m_)
  ); nprotect++;
  set_df_attributes(msegs_, fnedges, length(v1m_));
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertices:   data.frame(x = ..., y = ...)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP vert_ = PROTECT(
    create_named_list(2, "x", vert_x_, "y", vert_y_)
  ); nprotect++;
  set_df_attributes(vert_, ctx.nverts, max_verts);
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Lines: data.frame(a = ..., b = ..., c = ...)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP line_ = PROTECT(
    create_named_list(3, "a", line_a_, "b", line_b_, "c", line_c_)
  ); nprotect++;
  set_df_attributes(line_, ctx.nlines, max_edges);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segments: data.frame(line = integer(), v1 = integer(), v2 = integer())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP seg_ = PROTECT(
    create_named_list(3, "line", seg_line_, "v1", seg_v1_, "v2", seg_v2_)
  ); nprotect++;
  set_df_attributes(seg_, ctx.nsegs, max_edges);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extents: list(xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric())
  // Covers all seed points and voronoi vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP xmin_ = PROTECT(ScalarReal(bounds.xmin)); nprotect++;
  SEXP xmax_ = PROTECT(ScalarReal(bounds.xmax)); nprotect++;
  SEXP ymin_ = PROTECT(ScalarReal(bounds.ymin)); nprotect++;
  SEXP ymax_ = PROTECT(ScalarReal(bounds.ymax)); nprotect++;
  
  SEXP ext_ = PROTECT(
    create_named_list(4, "xmin", xmin_, "xmax", xmax_, "ymin", ymin_, "ymax", ymax_)
  ); nprotect++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Final result: named list of data.frames
  //  list(
  //     vertex  = data.frame()
  //     line    = data.frame()
  //     segment = data.frame()
  //     extents = list()
  // )
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(
    create_named_list(7, "vertex", vert_, "line", line_, "segment", seg_, 
                      "extents", ext_, "polygons", polys_, "msegments", msegs_,
                      "exterior", exterior_)
  ); nprotect++;
  setAttrib(res_, R_ClassSymbol, mkString("vor"));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' memory
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free_all_myalloc(&ctx);

  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  convert_indexing_c_to_r(seg_line_);
  convert_indexing_c_to_r(seg_v1_);
  convert_indexing_c_to_r(seg_v2_);
  
  convert_indexing_c_to_r(linem_);
  convert_indexing_c_to_r(v1m_);
  convert_indexing_c_to_r(v2m_);
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}


