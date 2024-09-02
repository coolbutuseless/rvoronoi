
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
#include "R-polygon.h"


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
SEXP extract_polygons(context_t *ctx) {
  
  int nprotect = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // How many bounded edges are there - with neither vertex index being NA
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n_undir_edges = 0;
  for (int i = 0; i < ctx->seg_idx; i++) {
    n_undir_edges += ctx->seg_v1[i] != NA_INTEGER && ctx->seg_v2[i] != NA_INTEGER;
  }
  // Rprintf("poly: %i/%i edges are bounded\n", n_undir_edges, ctx->seg_idx);
  int n_dir_edges = 2 * n_undir_edges;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Allocate space for TWICE the number of bounded edges (undirected)
  // 'edge' is going to hold **directed** edges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  edge_t *edge = calloc(n_dir_edges, sizeof(edge_t));
  if (edge == NULL) error("Could not allocate edge memory");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Copy edges in primary direction and calculate angle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int idx = 0;
  for (int i = 0; i < ctx->seg_idx; i++) {
    if (ctx->seg_v1[i] == NA_INTEGER || ctx->seg_v2[i] == NA_INTEGER)
      continue;
    
    int v1 = ctx->seg_v1[i] - 1; // Convert from Rs 1-index to C's 0-indexing
    int v2 = ctx->seg_v2[i] - 1;
    
    double x1 = ctx->vert_x[v1];
    double x2 = ctx->vert_x[v2];
    double y1 = ctx->vert_y[v1];
    double y2 = ctx->vert_y[v2];
    
    // Rprintf("(%2i, %2i) (%.3f, %.3f) (%.3f, %.3f)\n", v1, v2, x1, y1, x2, y2);
    
    double theta = atan2(y2 - y1, x2 - x1);
    theta = theta < 0 ? theta + 2 * M_PI : theta;
    
    edge[idx].v1    = v1;
    edge[idx].v2    = v2;
    edge[idx].theta = theta;
    idx++;
  }
  
  if (idx != n_undir_edges) {
    error("poly: sanity idx != nedges.  %i != %i", idx, n_undir_edges);
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
  qsort(edge, n_dir_edges, sizeof(edge_t), edge_comparison);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Create wedges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  wedge_t *wedge = calloc(n_dir_edges, sizeof(wedge_t));
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
  qsort(wedge, n_dir_edges, sizeof(wedge_t), wedge_comparison);
  
  
  // for (int i = 0; i < n_dir_edges; i++) {
  //   Rprintf("%3i %i: (%3i, %3i, %3i)\n", i, wedge[i].used, wedge[i].v1, wedge[i].v2, wedge[i].v3);
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
  int **polys;
  int poly_capacity = 32;
  int poly_idx = 0;
  polys = calloc(poly_capacity, sizeof(int *));
  if (polys == NULL) error("Couldn't allocate polys");
  
  // Keep track of the number of vertices within each polygon
  int *nverts;
  nverts = calloc(poly_capacity, sizeof(int));
  if (nverts == NULL) error("Couldn't allocate nverts");
  
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
      // My method: index the start of each v1_group and just 
      // jump to there to start the linear search.
      // TODO: Implement jump-started linear search
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int j = 0; // index of matching edge
      bool match_found = false;
      for (j = 0; j < n_dir_edges; j++) {
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
    if (poly_idx == poly_capacity) {
      // Exaand storage if we've reached capacity
      poly_capacity *= 2;
      polys = realloc(polys, poly_capacity * sizeof(int *));
      nverts = realloc(nverts, poly_capacity * sizeof(int));
    } 
    polys[poly_idx] = malloc(nvert * sizeof(int));
    if (polys[poly_idx] == NULL) error("polys[poly_idx] failed allocation");
    memcpy(polys[poly_idx], &vidx, nvert * sizeof(int));
    nverts[poly_idx] = nvert;
    
    poly_idx++;
  
  } // while(1): Find next unused wedge
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Dump info about all polygons
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // for (int i = 0; i < poly_idx; i++) {
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
  // Convert C polygons to an R list of lists of coordinates:
  //  list(list(x = ..., y = ...), list(x = ..., y = ...))
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(allocVector(VECSXP, poly_idx)); nprotect++;
  double *bbox;
  bbox = malloc(poly_idx * sizeof(double));
  if (bbox == NULL) error("Couldn't allocate bbox");
  
  for (int i = 0; i < poly_idx; i++) {
    int nvert = nverts[i];
    int *vidx = polys[i];
    
    SEXP ll_  = PROTECT(allocVector(VECSXP, 2)); 
    SEXP nms_ = PROTECT(allocVector(STRSXP, 2)); 
    
    SEXP x_ = PROTECT(allocVector(REALSXP, nvert)); 
    SEXP y_ = PROTECT(allocVector(REALSXP, nvert)); 
    
    SET_STRING_ELT(nms_, 0, mkChar("x"));
    SET_STRING_ELT(nms_, 1, mkChar("y"));
    
    SET_VECTOR_ELT(ll_, 0, x_);
    SET_VECTOR_ELT(ll_, 1, y_);
    
    setAttrib(ll_, R_NamesSymbol, nms_);
    
    SET_VECTOR_ELT(res_, i, ll_);
    
    double *x = REAL(x_);
    double *y = REAL(y_);
    
    double xmin =  INFINITY;
    double xmax = -INFINITY;
    double ymin =  INFINITY;
    double ymax = -INFINITY;
    
    for (int j = 0; j < nvert; j++) {
      x[j] = ctx->vert_x[ vidx[j] ];
      y[j] = ctx->vert_y[ vidx[j] ];
      
      xmin = x[j] < xmin ? x[j] : xmin;
      xmax = x[j] > xmax ? x[j] : xmax;
      ymin = y[j] < ymin ? y[j] : ymin;
      ymax = y[j] > ymax ? y[j] : ymax;
      
    }
    
    bbox[i] = (xmax - xmin) * (ymax - ymin);
    
    UNPROTECT(4);
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find max bounding box areas so we can remove the largest bounded region
  // The largest bounded region is the boundary around the outside, and
  // does represent an interior polygon
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_filtered_ = res_;
  if (poly_idx > 1) {
    int largest_bbox_idx = -1;
    double largest_bbox  = -1;
    
    for (int i = 0; i < poly_idx; i++) {
      // Rprintf("%i -> %.5f\n", i, bbox[i]);
      if (bbox[i] > largest_bbox) {
        largest_bbox     = bbox[i];
        largest_bbox_idx = i;
      }
    }
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Create a filtered version of the result with the polygon with largest
    // bounding box removed
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    res_filtered_ = PROTECT(allocVector(VECSXP, poly_idx - 1)); nprotect++;
    
    int out_idx = 0;
    for (int i = 0; i < poly_idx; i++) {
      if (i == largest_bbox_idx) {
        // do nothing
      } else {
        SET_VECTOR_ELT(res_filtered_, out_idx, VECTOR_ELT(res_, i));
        out_idx++;
      }
    }
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(bbox);
  for (int i = 0; i < poly_idx; i++) {
    free(polys[i]);
  }
  free(polys);
  free(nverts);
  
  free(wedge);
  free(edge);
  UNPROTECT(nprotect);
  return res_filtered_;
}

