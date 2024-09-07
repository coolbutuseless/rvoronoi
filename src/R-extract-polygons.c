
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
#include "R-extract-polygons.h"
#include "R-polygon-matching.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// An 'edge' is a directed edge with the tesselation
// For every bounded segment (ctx.seg_v1, ctx.seg_v2) add the two 
// directed edges  (v1, v2) and (v2, v1)
// 
// For every edge, include the angle in radians it makes with the horizontal.
// values should be in range [0, 2pi)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  int v1;
  int v2;
  double theta;
} edge_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// A 'wedge' is a sequence of 3 vertices defining an angle at 'v2'
// v1,v2,v3 are indexes into (ctx.vert_x, ctx.vert_y)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  int v1;
  int v2;
  int v3;
  bool used;
} wedge_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// qsort() helper for edges
//  Sort by 'v1' then 'theta'
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int edge_comparison(const void *p1, const void *p2) {
  
  edge_t *e1 = (edge_t *)p1;
  edge_t *e2 = (edge_t *)p2;
  
  if (e1->v1 < e2->v1) {
    return -1;
  } else if (e1->v1 > e2->v1) {
    return  1;
  } else {
    if (e1->theta < e2->theta) {
      return -1;
    } else {
      return 1;
    }
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// qsort() helper for wedges
//  Sort by 'v1' then 'v2'
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int wedge_comparison(const void *p1, const void *p2) {
  
  wedge_t *w1 = (wedge_t *)p1;
  wedge_t *w2 = (wedge_t *)p2;
  
  if (w1->v1 < w2->v1) {
    return -1;
  } else if (w1->v1 > w2->v1) {
    return  1;
  } else {
    if (w1->v2 < w2->v2) {
      return -1;
    } else {
      return 1;
    }
  }
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Free the C array-of-structs holding the poly_t polygon definition
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void free_polys(int npolys, poly_t *polys) {
  for (int i = 0; i < npolys; i++) {
    free(polys[i].v);
    free(polys[i].x);
    free(polys[i].y);
  }
  free(polys);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Given the vertices and edge information (after running voronoi())
// build the list of regions/polygons
//
// The method used in this package to convert a list of uncorrelated edges to 
// a list of polygons is from Jiang & Bunke's "An optimal 
// algorithm for extracting the regions of a plane graph" (Pattern Recognition
// Letters 14 (1993), p553-558).
//
//  * Phase 1: Find all the wedges
//     1. Split undirected edge into two directed edges
//     2. Add angle 
//     3. Sort by v1, theta
//     4. Form wedges
//  * Phase 2:  Group wedges into regions
//     1. Sort the wedge list by v0, v1
//     2. Mark all wedges as unused
//     3. Find the next unused wedge - mark as 'used', add to region. If no unsed wedge, then algorithm complete
//     4. Search for matching continuing wedge. E.g. if initial wedge is c(1, 4, 7), look for wedge c(4, 7, x). Add to region
//     5. If new wedge matches original wedge, then: region is extracted. Go to 3 else: go to 4
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
poly_t *extract_polygons_core(int vert_n, double *vert_x, double *vert_y, int seg_n, int *seg_v1, int *seg_v2, int *npolysp) {
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // How many bounded edges are there - with neither vertex index being NA
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n_undir_edges = 0;
  for (int i = 0; i < seg_n; i++) {
    n_undir_edges += seg_v1[i] >= 0 && seg_v2[i] >= 0;
  }
  // Rprintf("poly: %i/%i edges are bounded\n", n_undir_edges, seg_n);
  int n_dir_edges = 2 * n_undir_edges;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Allocate space for TWICE the number of bounded edges (undirected)
  // 'edge' is going to hold **directed** edges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  edge_t *edge = calloc((unsigned long)n_dir_edges, sizeof(edge_t));
  if (edge == NULL) error("Could not allocate edge memory");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Copy edges in primary direction and calculate angle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int idx = 0;
  for (int i = 0; i < seg_n; i++) {
    if (seg_v1[i] < 0 || seg_v2[i] < 0)
      continue;
    
    if (seg_v1[i] == seg_v2[i]) {
      error("dupe vertex at C index %i. Use `merge_vertices()`", i);
    }
    
    int v1 = seg_v1[i]; // Convert from Rs 1-index to C's 0-indexing
    int v2 = seg_v2[i];
    
    double x1 = vert_x[v1];
    double x2 = vert_x[v2];
    double y1 = vert_y[v1];
    double y2 = vert_y[v2];
    
    // Rprintf("(%2i, %2i) (%.3f, %.3f) (%.3f, %.3f)\n", v1, v2, x1, y1, x2, y2);
    
    double theta = atan2(y2 - y1, x2 - x1);
    theta = theta < 0 ? theta + 2 * M_PI : theta;
    
    edge[idx].v1    = v1;
    edge[idx].v2    = v2;
    edge[idx].theta = theta;
    idx++;
  }
  
  if (idx != n_undir_edges) {
    error("extract_polygons() poly: sanity idx != nedges.  %i != %i", idx, n_undir_edges);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Copy edges in the alternate direction. Flip the angle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < n_undir_edges; i++) {
    edge[n_undir_edges + i].v1    = edge[i].v2;
    edge[n_undir_edges + i].v2    = edge[i].v1;
    double theta = edge[i].theta + M_PI;
    theta = theta >= 2 * M_PI ? theta - 2 * M_PI : theta;
    edge[n_undir_edges + i].theta = theta;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sort edges by (v1, theta)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  qsort(edge, (size_t)n_dir_edges, sizeof(edge_t), edge_comparison);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Create wedges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  wedge_t *wedge = calloc((unsigned long)n_dir_edges, sizeof(wedge_t));
  if (wedge == NULL) error("Could not allocated wedge memory");
  idx = 0;
  int v1_group = edge[0].v1;
  int v2_first = edge[0].v2;
  
  // Create a wedge from each edge
  int i;
  for (i = 0; i < n_dir_edges - 1; i++) {
    
    // If the next 'v1' matches the current 'v1_group' just create the wedge
    // from  (v2-current, v1_group, v2-next)
    if (edge[i+1].v1 == v1_group) {
      wedge[i].v1 = edge[i].v2;
      wedge[i].v2 = v1_group;
      wedge[i].v3 = edge[i+1].v2;
    } else {
      // Create the final wedge in a 'v1-group' by looping back to the first 
      // edge in this group
      wedge[i].v1 = edge[i].v2;
      wedge[i].v2 = v1_group;
      wedge[i].v3 = v2_first;
      
      v1_group = edge[i+1].v1;
      v2_first = edge[i+1].v2;
    }
  }
  
  // Close last group
  wedge[i].v1 = edge[i].v2;
  wedge[i].v2 = v1_group;
  wedge[i].v3 = v2_first;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sort wedges by (v1, v2)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  qsort(wedge, (size_t)n_dir_edges, sizeof(wedge_t), wedge_comparison);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Figure out search bounds for each wedge starting with a particular
  // vertex
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  typedef struct {
    int first; // At which index is this vertex first seen in 'wedge'?
    int last;  // At which index is this vertex last  seen in 'wedge'?
  } bounds_t;
  bounds_t *bounds = NULL;

  // Find maximum vertex which appears - which tells me the memory 
  // allocation for the number of 'bounds'
  int max_vert_num = wedge[n_dir_edges - 1].v1;
  int nbounds = max_vert_num + 1;

  // Allocate the 'bounds'
  bounds = calloc(nbounds, sizeof(bounds_t));
  if (bounds == NULL) error("Failed to allocate 'bounds'");

  // Initialise the cvert (current vertex) to be the first seen vertex
  int cvert = wedge[0].v1;
  
  // Note that the first seen vertex starts at '0'
  bounds[cvert].first = 0;

  // Find the start and end of all indices
  for (int i = 1; i < n_dir_edges; i++) {
    if (wedge[i].v1 != cvert) {
      bounds[cvert].last = i - 1;
      cvert = wedge[i].v1;
      bounds[cvert].first = i;
    }
  }

  // last time the last vertex index is seen is the last edge (obviously :)
  bounds[cvert].last = n_dir_edges - 1;

  // Verbosity/debugging
  // for (int i = 0; i < n_dir_edges; i++) {
  //   Rprintf("W %3i %i: (%3i, %3i, %3i)\n", i, wedge[i].used, wedge[i].v1, wedge[i].v2, wedge[i].v3);
  // }
  //
  // for (int i = 0; i < nbounds; i++) {
  //   Rprintf("B %i  [%i, %i]\n", i, bounds[i].first, bounds[i].last);
  // }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  * Phase 2:  Group wedges into regions
  //     3. Find the next unused wedge - mark as 'used', add to region. 
  //        If no unused wedge, then algorithm complete
  //     4. Search for matching continuing wedge. E.g. if initial wedge is 
  //        c(1, 4, 7), look for wedge c(4, 7, x). Add to region
  //     5. If new wedge matches original wedge, 
  //        Then: Complete region has been extracted. Go to 3 
  //        ELse: go to 4 so find next matching wedge
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Current polygons in C form: pointer to arrays of integer vertex indices
  poly_t *polys;
  int poly_capacity = 32;
  int npolys = 0;
  polys = calloc((unsigned long)poly_capacity, sizeof(poly_t));
  if (polys == NULL) error("Couldn't allocate polys");
  
  // Temp storage for the current set of vertex indices
  int vidx[1024];
  int nvert = 0;
  
  while(1) {
    
    nvert = 0;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Find first unused wedge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int first = 0;
    bool found_unused = false;
    for (first = 0; first < n_dir_edges; first++) {
      if (!wedge[first].used) {
        wedge[first].used = true;
        found_unused = true;
        break;
      }
    }
    
    if (!found_unused) {
      // No more unused wedges. Job complete!
      break;
    }
    
    vidx[nvert++] = wedge[first].v1;
    vidx[nvert++] = wedge[first].v2;
    
    int i = first;  // current wedge
    bool capturing = true;
    
    
    // Find all subsequent wedge matches for the initial wedge
    while (true) {
      
      if (capturing) {
        vidx[nvert++] = wedge[i].v3;
      }
      
      if (nvert > 1000) {
        error("'nvert' > 1000");
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Linear search for matching wedge
      // Original paper uses binary search to achieve O(nlogn)
      // My method: index the start/end index of each v1 within wedge 
      // and just search within there.
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int j = 0; // index of matching edge
      bool match_found = false;
      
      
      // Old linear search for next wedge
      // for (j = 0; j < n_dir_edges; j++) {
      //   if (!wedge[j].used &&
      //       wedge[j].v1 == wedge[i].v2 &&
      //       wedge[j].v2 == wedge[i].v3) {
      //     wedge[j].used = true;
      //     match_found = true;
      //     break;
      //   }
      // }
      
      // New sub-limear search for next matching wedge.
      // just searches between the known start/end index in 'wedge' where
      int this_v1 = wedge[i].v2;
      for (j = bounds[this_v1].first; j <= bounds[this_v1].last; j++) {
        if (!wedge[j].used &&
            wedge[j].v1 == wedge[i].v2 &&
            wedge[j].v2 == wedge[i].v3) {
          wedge[j].used = true;
          match_found = true;
          break;
        }
      }
      
      if (!match_found) {
        error("Internal error: no 'match' wedge found");
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // If the end of this matching wedge is identical to the 
      // start of the first wedge, then we have completed a polygon!
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (wedge[j].v3 == wedge[first].v1) {
        capturing = false;
      }
      
      if (wedge[j].v2 == wedge[first].v1 &&
          wedge[j].v3 == wedge[first].v2) {
        break;
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Continue search to find the next matching wedge
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      i = j;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // record this set of temporary verts as a polygon
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (npolys == poly_capacity) {
      // Exaand storage if we've reached capacity
      poly_capacity *= 2;
      polys  = realloc(polys , (unsigned long)poly_capacity * sizeof(poly_t));
    } 
    
    polys[npolys].v = malloc((unsigned long)nvert * sizeof(int));
    if (polys[npolys].v == NULL) error("polys[npolys] failed allocation");
    memcpy(polys[npolys].v, &vidx, (unsigned long)nvert * sizeof(int));
    polys[npolys].nvert = nvert;
    
    npolys++;
  
  } // while(1): Find next unused wedge
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Dump info about all polygons
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // for (int i = 0; i < npolys; i++) {
  //   int nvert = nverts[i];
  //   int *vidx = polys[i];
  //   
  //   Rprintf("[%2i]  ", i);
  //   for (int j = 0; j < nvert; j++) {
  //     Rprintf("%2i ", vidx[j]);
  //   }
  //   Rprintf("\n");
  //   
  // }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate polygon coordinates, bbox etc
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < npolys; i++) {
    int nvert = polys[i].nvert;
    int *vidx = polys[i].v;
    
    polys[i].polygon_idx = i;
    polys[i].taken = 0;
    polys[i].deleted = 0;
    polys[i].point_idx = -1;
    polys[i].cx = 0;
    polys[i].cy = 0;
    
    polys[i].bbox = bbox_new();
    
    polys[i].x = malloc(nvert * sizeof(double));
    polys[i].y = malloc(nvert * sizeof(double));
    if (polys[i].x == NULL || polys[i].y == NULL) {
      error("Could not allocate mem for polys[i].x and .y");
    }
    
    for (int j = 0; j < nvert; j++) {
      polys[i].x[j] = vert_x[ vidx[j] ];
      polys[i].y[j] = vert_y[ vidx[j] ];
      
      polys[i].cx += polys[i].x[j];
      polys[i].cy += polys[i].y[j];
      
      // polys[i].bbox.xmin = polys[i].x[j] < polys[i].bbox.xmin ? polys[i].x[j] : polys[i].bbox.xmin;
      // polys[i].bbox.xmax = polys[i].x[j] > polys[i].bbox.xmax ? polys[i].x[j] : polys[i].bbox.xmax;
      // polys[i].bbox.ymin = polys[i].y[j] < polys[i].bbox.ymin ? polys[i].y[j] : polys[i].bbox.ymin;
      // polys[i].bbox.ymax = polys[i].y[j] > polys[i].bbox.ymax ? polys[i].y[j] : polys[i].bbox.ymax;
    }
    
    bbox_add(&polys[i].bbox, nvert, polys[i].x, polys[i].y);
    
    
    
    // compute centroid
    polys[i].cx /= nvert;
    polys[i].cy /= nvert;
    
    // compute bbox area
    polys[i].bbox_area = (polys[i].bbox.xmax - polys[i].bbox.xmin) * (polys[i].bbox.ymax - polys[i].bbox.ymin);
  }
  
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(bounds);
  free(wedge);
  free(edge);
  
  *npolysp = npolys;
  return polys;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Create 'poly_t *polys' object and then translate into final 
// R SEXP nested list representation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP extract_polygons_internal(int vert_n, double *vert_x, double *vert_y, 
                               int seg_n, int *seg_v1, int *seg_v2,
                               int seed_n, double *seed_x, double *seed_y) {
  
  int nprotect = 0;
  int npolys = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extract polygons
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  poly_t *polys = extract_polygons_core(vert_n, vert_x, vert_y, seg_n, seg_v1, seg_v2, &npolys);
  
  if (npolys == 0) {
    return R_NilValue;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find largest polygon by bbox
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int max_bbox_area_idx = -1;
  double max_bbox_area  = -1;
  for (int i = 0; i < npolys; i++) {
    if (polys[i].bbox_area > max_bbox_area) {
      max_bbox_area = polys[i].bbox_area;
      max_bbox_area_idx = i;
    }
  }
  polys[max_bbox_area_idx].deleted = true;
  int n_valid_polys = npolys - 1;  
  
  SEXP polys_ = R_NilValue;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // If seed points given then reorganise polygons in the same order as 
  // the seed points.  Otherwise a generic compact representation is used
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (seed_n > 0) {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Seed / Polygon matching
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Rprintf("extract_polygons_internal: n seeds = %i\n", seed_n);
    for (int i = 0; i < seed_n; i++) {
      match_polygons_to_seed_points(i, seed_x[i], seed_y[i], npolys, polys);
    }
    
    polys_ = PROTECT(allocVector(VECSXP, seed_n)); nprotect++;
    
    for (int i = 0; i < npolys; i++) {
      if (polys[i].deleted) continue;
      
      SEXP ll_  = PROTECT(allocVector(VECSXP, 2)); 
      SEXP nms_ = PROTECT(allocVector(STRSXP, 2)); 
      
      SEXP x_ = PROTECT(allocVector(REALSXP, polys[i].nvert)); 
      SEXP y_ = PROTECT(allocVector(REALSXP, polys[i].nvert)); 
      
      SET_STRING_ELT(nms_, 0, mkChar("x"));
      SET_STRING_ELT(nms_, 1, mkChar("y"));
      
      SET_VECTOR_ELT(ll_, 0, x_);
      SET_VECTOR_ELT(ll_, 1, y_);
      
      setAttrib(ll_, R_NamesSymbol, nms_);
      
      if (polys[i].point_idx < 0) {
        warning("Poly [%i] has a point index of %i\n", i, polys[i].point_idx);  
      } else {
        SET_VECTOR_ELT(polys_, polys[i].point_idx, ll_);
      }
      
      // Copy from the poly_t struct
      memcpy(REAL(x_), polys[i].x, polys[i].nvert * sizeof(double));
      memcpy(REAL(y_), polys[i].y, polys[i].nvert * sizeof(double));
      
      
      UNPROTECT(4); // everything is protected as they're now members of list 'res_'
    }
    
  } else {
    // No seed points. Just a compact list of polygons
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Convert C polygons to an R list of lists of coordinates:
    //  list(list(x = ..., y = ...), list(x = ..., y = ...))
    // 
    // Remove any deleted polygons while doing this.
    //  e.g. may have been deleted for being the bounding region rather than a
    //       single polygon
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    polys_ = PROTECT(allocVector(VECSXP, n_valid_polys)); nprotect++;
    
    int llidx = 0;  
    for (int i = 0; i < npolys; i++) {
      if (polys[i].deleted) continue;
      
      SEXP ll_  = PROTECT(allocVector(VECSXP, 2)); 
      SEXP nms_ = PROTECT(allocVector(STRSXP, 2)); 
      
      SEXP x_ = PROTECT(allocVector(REALSXP, polys[i].nvert)); 
      SEXP y_ = PROTECT(allocVector(REALSXP, polys[i].nvert)); 
      
      SET_STRING_ELT(nms_, 0, mkChar("x"));
      SET_STRING_ELT(nms_, 1, mkChar("y"));
      
      SET_VECTOR_ELT(ll_, 0, x_);
      SET_VECTOR_ELT(ll_, 1, y_);
      
      setAttrib(ll_, R_NamesSymbol, nms_);
      
      SET_VECTOR_ELT(polys_, llidx, ll_);
      
      // Copy from the poly_t struct
      memcpy(REAL(x_), polys[i].x, polys[i].nvert * sizeof(double));
      memcpy(REAL(y_), polys[i].y, polys[i].nvert * sizeof(double));
      
      
      UNPROTECT(4); // everything is protected as they're now members of list 'res_'
      llidx++;
    }
  }

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free_polys(npolys, polys);
  UNPROTECT(nprotect);
  return polys_;
}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shimx`
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP extract_polygons_(SEXP x_, SEXP y_, SEXP v1_, SEXP v2_) {

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
  
  int *v1 = malloc((unsigned long)length(v1_) * sizeof(int));
  int *v2 = malloc((unsigned long)length(v2_) * sizeof(int));
  if (v1 == NULL || v2 == NULL) error("extract_polygons_() v1/v2 failed allocation");
  
  int *v1_p = INTEGER(v1_);
  int *v2_p = INTEGER(v2_);
  
  // Convert from R 1-indexing to C 0-indexing
  for (int i = 0; i < length(v1_); i++) {
    v1[i] = v1_p[i] - 1;
    v2[i] = v2_p[i] - 1;
  }
  
  
  // SEXP extract_polygons_internal(int vert_n, double *vert_x, double *vert_y, int seg_n, int *seg_v1, int *seg_v2)
  return extract_polygons_internal(
    length(x_), REAL(x_), REAL(y_), // voronoi vertices
    length(v1_), v1, v2,            // voronoi edges
    0, NULL, NULL                   // seed points 
  );
}



