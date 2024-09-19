
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
#include "R-infinite-segments.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Voronoi Tessellation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP voronoi_(SEXP x_, SEXP y_, SEXP calc_polygons_, SEXP match_sites_, SEXP bound_segments_,
              SEXP merge_tolerance_) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity Check
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) != length(y_)) {
    error("x & y are not the same length");
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // setup
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int nprotect = 0;
  bool calc_polygons  = asLogical(calc_polygons_);
  bool match_sites    = asLogical(match_sites_);
  bool bound_segments = asLogical(bound_segments_);
  double merge_tolerance = asReal(merge_tolerance_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Can't form any polygons if there are no sites!
  // We still go through the motions of this function in order
  // to return data to R in the right form
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) == 0) {
    calc_polygons  = false;
    match_sites    = false;
    bound_segments = false;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialise the calculation context
  // Do delauney triangulation? FALSE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  context_t ctx = { 0 };
  ctx.triangulate = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Track *ALL* the allocations done via 'myalloc()'
  // These will be freed in bulk by calling 'free_all_myalloc()'
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
  // Memory allocation for vertices in voronoi diagram.
  //
  // Note: the upper bound on number of vertices includes 
  //   (a) vertices and segments in the diagram itself
  //   (b) extra vertices and segments created when we bound the voronoi within 
  //       a rectangular boundary.  
  //       NOTE: The boundary is always encompasses all sites and all finite
  //             voronoi vertices
  //
  // If there are 'n' sites then 
  //   maximum number of vertices in the voronoi tessellation is 2n - 5
  //   maximum number of segments is 3n - 6
  //
  // However, we also need to account for infinite and semi-infinite segments
  //
  // If there are 'n' sites, an upper bound on the number of 
  // semi-infinite/infinite segments = 'n'
  // 
  // If these infinite rays segents being clipped by intersection with a rectangular
  //     boundary, then
  //     - there are 4 extra vertices for the corners of the boundary.
  //     - there is 1 extra vertex where the semi-infinite segment intersects the boundary.
  //     - for the special case of an infinite segment, then there
  //       are two intersection vertices created.
  //     - for a sequence of co-linear sites, there will be multiple
  //       infinite segments.  For 'n' co-linear sites there is a maximum
  //       of 'n-1' infinite segments, which would mean '2n-2' boundary 
  //       intersection vertices
  //
  // The maximum number of vertices which lie on a rectangular boundary vertices = 2n - 2
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
  //        A better analysis here would give me a lower upper bound, 
  //        but I don't need a tight bound at the moment as size of
  //        voronois is not expected to be must about 1000 seed points.
  //
  // Add a margin of '10' just to widely avoid off-by-one errors in my thinking.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n = length(x_); // number of seed points
  
  int max_interior_verts = 2 * n - 5;
  int max_interior_edges = 3 * n - 6;
  
  // int max_exterior_verts = 2 * n + 4;
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
  
  ctx.nsegs    = 0;
  ctx.seg_line = INTEGER(seg_line_);
  ctx.seg_v1   = INTEGER(seg_v1_);
  ctx.seg_v2   = INTEGER(seg_v2_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate voronoi tessellation using Fortune's Sweep algorithm 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  freeinit(&ctx.sfl, sizeof *ctx.sites);
  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));
  ctx.siteidx = 0;
  geominit(&ctx);
  voronoi(&ctx, ctx.triangulate);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Add voronoi vertices to bounding box. Expand size by 10%
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bbox_add(&bounds, ctx.nverts, ctx.vert_x, ctx.vert_y);
  bbox_expand(&bounds, 0.10);
  
  
  SEXP polys_  = R_NilValue; // Extracted polygons
  SEXP msegs_  = R_NilValue; // Segments after merging and polygon extraction
  SEXP mverts_ = R_NilValue; // vertices after merging and polygon extraction
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate the extra vertices and segments produced if we 
  // bound the infinite rays
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int nbverts = 0;
  int nbsegs  = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Merge vertices in voronoi segments to prep for polygon building
  // void merge_vertices_core_(double tol, 
  //                           int nverts, double *x, double *y, 
  //                           int nedges, int *v1, int *v2, 
  //                           int *fnedges,
  //                           int verbosity)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int fnedges = 0; // Final number of edges after merging
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Final vertices after merging and polygon extraction
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP xf_ = R_NilValue;
  SEXP yf_ = R_NilValue;
  
  double *xf = NULL;
  double *yf = NULL;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Final segments after merging and polygon extraction
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP v1m_   = R_NilValue;
  SEXP v2m_   = R_NilValue;
  SEXP linem_ = R_NilValue;
  
  int *v1m   = NULL;
  int *v2m   = NULL;
  int *linem = NULL;

  
  if (calc_polygons || bound_segments) {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set up the "working area" for the merging of vertices
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    v1m_   = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
    v2m_   = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
    linem_ = PROTECT(allocVector(INTSXP, ctx.nsegs + max_exterior_edges)); nprotect++;
    
    msegs_ = PROTECT(
      create_named_list(3, "line", linem_, "v1", v1m_, "v2", v2m_)
    ); nprotect++;
    
    v1m   = INTEGER(v1m_);
    v2m   = INTEGER(v2m_);
    linem = INTEGER(linem_);
    
    memcpy(v1m  , ctx.seg_v1  , ctx.nsegs * sizeof(int));
    memcpy(v2m  , ctx.seg_v2  , ctx.nsegs * sizeof(int));
    memcpy(linem, ctx.seg_line, ctx.nsegs * sizeof(int));
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Merge close vertices which are an artefact of the tessellation
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    merge_vertices_core_(merge_tolerance, 
                         ctx.nverts, ctx.vert_x, ctx.vert_y,
                         ctx.nsegs, linem, v1m, v2m, 
                         &fnedges, 0);
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Use the merged segments to find infinite segments (infinite and semi-infinite)
    // and then calculate the actual number of vertices and segments needed
    // to store the bounded external polygons
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calc_space_for_bound_infinite_segments(fnedges, v1m, v2m, &nbverts, &nbsegs);
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Allocate the bounded vertices and segment vertex indices
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double *xb = calloc(nbverts , sizeof(double));
    double *yb = calloc(nbverts , sizeof(double));
    int *rv1   = calloc(nbsegs  , sizeof(int));
    int *rv2   = calloc(nbsegs  , sizeof(int));
    
    bound_infinite_segments(
      &bounds,
      ctx.nverts, ctx.vert_x, ctx.vert_y,
      fnedges, linem, v1m, v2m,
      ctx.nlines, ctx.line_a, ctx.line_b, ctx.line_c,
      &nbverts, xb, yb,
      &nbsegs, rv1, rv2
    );
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Create a new MEGA vertex list by concatenting
    //   * the voronoi vertices
    //   * the exterior vertices
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xf_ = PROTECT(allocVector(REALSXP, ctx.nverts + nbverts)); nprotect++;
    yf_ = PROTECT(allocVector(REALSXP, ctx.nverts + nbverts)); nprotect++;
    xf  = REAL(xf_);
    yf  = REAL(yf_);
    
    memcpy(xf + 0         , ctx.vert_x, ctx.nverts * sizeof(double));
    memcpy(yf + 0         , ctx.vert_y, ctx.nverts * sizeof(double));
    memcpy(xf + ctx.nverts,         xb,    nbverts * sizeof(double));
    memcpy(yf + ctx.nverts,         yb,    nbverts * sizeof(double));
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Append the new exterior segments to the merged vertex list
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    memcpy(v1m + fnedges, rv1, nbsegs * sizeof(int));
    memcpy(v2m + fnedges, rv2, nbsegs * sizeof(int));
    for (int i = 0; i < nbsegs; i++) {
      (linem + fnedges)[i] = INVALID_VIDX;
    }
    
    free(xb);
    free(yb);
    free(rv1);
    free(rv2);
    fnedges += nbsegs;
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Trim the merged indices to size (and make into a data.frame)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    set_df_attributes_and_trim(msegs_, fnedges, length(v1m_));
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // All the vertices including those from bounding infinite segments 
    // and perimeter
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mverts_ = PROTECT(
      create_named_list(
        2, 
        "x", xf_,
        "y", yf_
      )
    ); nprotect++;
    set_df_attributes(mverts_);
  } // End of calculation of merged/bound vertices
  
  
  if (calc_polygons) {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the polygons using the (temporary) merged/bound vertices/segments
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (match_sites) {
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
          fnedges, v1m, v2m,               // Voronoi edges
          0, NULL, NULL                    // Seed points
        )
      ); nprotect++;
    }
    
    convert_indexing_c_to_r_with_NA(linem_);
    convert_indexing_c_to_r_with_NA(v1m_);
    convert_indexing_c_to_r_with_NA(v2m_);
    
  }  // End: if (calc_polygons)
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertices:   data.frame(x = ..., y = ...)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP vert_ = PROTECT(
    create_named_list(2, "x", vert_x_, "y", vert_y_)
  ); nprotect++;
  set_df_attributes_and_trim(vert_, ctx.nverts, max_verts);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Lines: data.frame(a = ..., b = ..., c = ...)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP line_ = PROTECT(
    create_named_list(3, "a", line_a_, "b", line_b_, "c", line_c_)
  ); nprotect++;
  set_df_attributes_and_trim(line_, ctx.nlines, max_edges);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segments: data.frame(line = integer(), v1 = integer(), v2 = integer())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP seg_ = PROTECT(
    create_named_list(3, "line", seg_line_, "v1", seg_v1_, "v2", seg_v2_)
  ); nprotect++;
  set_df_attributes_and_trim(seg_, ctx.nsegs, max_edges);
  
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
  //     sites = ..., vertices = ..., ...
  // )
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP sites_ = PROTECT(
    create_named_list(2, "x", x_, "y", y_)  // original site/seed points  
  ); nprotect++;
  set_df_attributes(sites_);
  
  SEXP res_ = PROTECT(
    create_named_list(
      8, 
      "sites"    , sites_,
      "vertices" , vert_, 
      "segments" , seg_,
      "polygons" , polys_, 
      "lines"    , line_, 
      "extents"  , ext_, 
      "mvertices", mverts_,
      "msegments", msegs_
    )
  ); nprotect++;
  setAttrib(res_, R_ClassSymbol, mkString("vor"));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' memory
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free_all_myalloc(&ctx);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  convert_indexing_c_to_r_with_NA(seg_line_);
  convert_indexing_c_to_r_with_NA(seg_v1_);
  convert_indexing_c_to_r_with_NA(seg_v2_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}


