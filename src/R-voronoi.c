
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


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Voronoi Tesselation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP voronoi_(SEXP x_, SEXP y_, SEXP match_polygons_) {

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
  // bbox for bounding the unbounded polygos
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bbox_t bounds = bbox_new();
  bbox_add(&bounds, length(x_), REAL(x_), REAL(y_));
  
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
  SEXP v1m_   = PROTECT(allocVector(INTSXP, ctx.nsegs)); nprotect++;
  SEXP v2m_   = PROTECT(allocVector(INTSXP, ctx.nsegs)); nprotect++;
  SEXP linem_ = PROTECT(allocVector(INTSXP, ctx.nsegs)); nprotect++;
  
  int *v1m   = INTEGER(v1m_);
  int *v2m   = INTEGER(v2m_);
  int *linem = INTEGER(linem_);
  
  memcpy(v1m  , ctx.seg_v1  , ctx.nsegs * sizeof(int));
  memcpy(v2m  , ctx.seg_v2  , ctx.nsegs * sizeof(int));
  memcpy(linem, ctx.seg_line, ctx.nsegs * sizeof(int));

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Merge close vertices which are an artefact of the tessellation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  merge_vertices_core_(1e-5, 
                       ctx.nverts, ctx.vert_x, ctx.vert_y,
                       ctx.nsegs, linem, v1m, v2m, 
                       &fnedges, 0);
  
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
  
  if (asLogical(match_polygons_)) {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts, ctx.vert_x, ctx.vert_y, // Voronoi vertices
        fnedges, v1m, v2m,                  // Voronoi edges
        length(x_), REAL(x_), REAL(y_)      // Seed points
      )
    ); nprotect++;
  } else {
    polys_ = PROTECT(
      extract_polygons_internal(
        ctx.nverts, ctx.vert_x, ctx.vert_y, // Voronoi vertices
        fnedges, v1m, v2m,                  // Voronoi edges
        0, NULL, NULL                       // Seed points
      )
    ); nprotect++;
  }
  
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
  bbox_add(&bounds, ctx.nverts, ctx.vert_x, ctx.vert_y);
  bbox_expand(&bounds, 0.10);
  
  
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
    create_named_list(6, "vertex", vert_, "line", line_, "segment", seg_, 
                      "extents", ext_, "polygons", polys_, "msegments", msegs_)
  ); nprotect++;
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


