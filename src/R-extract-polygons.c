
#define R_NO_REMAP

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
// An 'edge' is a directed edge
// 
// For every segment (line, v1, v2)
//    add the two directed edges  (v1, v2) and (v2, v1)
// 
// For every edge, include the angle in radians it makes with the horizontal.
// Angles should be in range [0, 2pi)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  int v1;
  int v2;
  double theta; // Angle in radians. In range [0, 2pi)
} edge_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// A 'wedge' is a sequence of 3 vertices defining the angle at 'v2'
// v1,v2,v3 are indexes into (ctx.xvor, ctx.yvor)
//
// 'used' keeps track of whether the wedge has already been used as part
// of a polygon
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
// Remember that there are actually "npolys + 1" polygons allocated.
// One polygon is always ignored as it is the "boudning region" and 
// not a voronoi cell
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void free_polys(int npolys, poly_t *polys) {
  for (int i = 0; i < npolys + 1; i++) {
    free(polys[i].v);
    free(polys[i].x);
    free(polys[i].y);
  }
  free(polys);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Fortune's Voronoi returns vertices and segments (line, v1, v2)
//
// This function builds those elements into a set of polygons.
//
// The method used in this package to convert a list of uncorrelated segments to 
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
//
// 'poly_t' is defined in 'R-polygon-matching.h'
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
poly_t *extract_polygons_core(int n_vor_verts, double *xvor, double *yvor, 
                              int n_vor_segs, int *v1vor, int *v2vor, 
                              int *npolys) {
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // How many finite segments are there? 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n_finite_segments = 0;
  for (int i = 0; i < n_vor_segs; i++) {
    n_finite_segments += valid_idx(v1vor[i]) && valid_idx(v2vor[i]);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Check there are some finite segments to process. 
  // If there are no segments then no polygons are possible!
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (n_finite_segments == 0) {
    *npolys = 0;
    return NULL;
  }
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // For each finite segment in the voronoi, two directed edges need
  // to be created.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n_dir_edges = 2 * n_finite_segments;

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Allocate space for directed edges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  edge_t *edge = calloc((unsigned long)n_dir_edges, sizeof(edge_t));
  if (edge == NULL) Rf_error("extract_polygons_core(): Could not allocate edge memory");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Copy edges in primary direction and calculate angle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int idx = 0;
  for (int i = 0; i < n_vor_segs; i++) {
    
    // The vertex indicies at each end of this segment
    int v1 = v1vor[i]; 
    int v2 = v2vor[i];
    
    // Don't process semi-infinite or infinite segments
    if (!valid_idx(v1) || !valid_idx(v2))
      continue;
    
    // Sanity check. There should not be any zero-length segments 
    if (v1 == v2) {
      Rf_error("extract_polygons_core(): dupe vertex at C index %i. Use `merge_vertices()`?", i);
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate the angle of this segment (in radians) and
    // ensure that it lies within the range [0, 2pi)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double x1 = xvor[v1];
    double x2 = xvor[v2];
    double y1 = yvor[v1];
    double y2 = yvor[v2];
    double theta = atan2(y2 - y1, x2 - x1);
    theta = theta < 0 ? theta + 2 * M_PI : theta;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Create this directed edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    edge[idx].v1    = v1;
    edge[idx].v2    = v2;
    edge[idx].theta = theta;
    idx++;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity check 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (idx != n_finite_segments) {
    Rf_error("extract_polygons() poly: sanity idx != nedges.  %i != %i", idx, n_finite_segments);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Copy edges in the alternate direction and flip the angle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < n_finite_segments; i++) {
    edge[n_finite_segments + i].v1    = edge[i].v2;
    edge[n_finite_segments + i].v2    = edge[i].v1;
    double theta = edge[i].theta + M_PI;
    theta = theta >= 2 * M_PI ? theta - 2 * M_PI : theta;
    edge[n_finite_segments + i].theta = theta;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sort edges by (v1, theta)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  qsort(edge, (size_t)n_dir_edges, sizeof(edge_t), edge_comparison);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Create wedges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  wedge_t *wedge = calloc((unsigned long)n_dir_edges, sizeof(wedge_t));
  if (wedge == NULL) Rf_error("Could not allocated wedge memory");
  
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
  // For each wedge (i1, i2, i3), we must find the matching wedge 
  // (i2, i3, *) which has not already been claimed.
  //
  // Since the wedges are sorted by (i1, i2), the original paper does a 
  // binary search to find the matching wedge.
  // 
  // However! For a generic random (non-pathological) voronoi each vertex 
  // will only take part in a handful of wedges.
  //
  // For this implementation of the algorithm, I will pre-calculate 
  // the start/stop index of each 'i1' value in the sorted wedge list.  
  // 
  // This means that a search for (i1, i2, *) can lookup the wedge index
  // for where (i1, *, *) starts and ends, and onlly check for matches within
  // those bounds.
  //
  // My totally untested assumption is that this index lookup method 
  // is cache friendlier than a binary search.
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
  bounds = calloc((size_t)nbounds, sizeof(bounds_t));
  if (bounds == NULL) Rf_error("Failed to allocate %i 'bounds' for %i undirected edges", nbounds, n_dir_edges);

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

  // The last time the final vertex index is seen is the last edge 
  bounds[cvert].last = n_dir_edges - 1;
  
  
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
  poly_t *polys;
  int poly_capacity = 32;
  polys = calloc((unsigned long)poly_capacity, sizeof(poly_t));
  if (polys == NULL) Rf_error("Couldn't allocate polys");
  
  // Temp storage for the current set of vertex indices
#define MAX_VERTS_PER_POLY 4096
  int vidx[MAX_VERTS_PER_POLY];
  int nvert = 0; // Number of vertices in this particular polygon
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Keep searching for polygons until we run out of wedges
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  while(1) {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reset the vertex count for this new polygon
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // If no unsused wedges were found, it means we've extracted all 
    // the polygons from this data.  Job complete!
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!found_unused) {
      break;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Record the v1 and v2 for the wedge found (v1, v2, *)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vidx[nvert++] = wedge[first].v1;
    vidx[nvert++] = wedge[first].v2;
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // i = current wedge index
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int i = first;  
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // When we loop around to the first wedge, we stop capturing vertices
    // for this polygon, but we need to keep marking wedges as used.
    // For any particular polygon extraction, capturing = false for the 
    // last two wedges in the sequence
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool capturing = true;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Find subsequent wedge match for the current wedge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while (true) {
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Store the 'v3' vertex from the current wedge
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (capturing) {
        vidx[nvert++] = wedge[i].v3;
        if (nvert == MAX_VERTS_PER_POLY) {
          Rf_error("Number of vertices for a single polygon exceeds limit of 16384. Please file an issue");
        }
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Linear search for matching wedge
      // Original paper uses binary search to achieve O(nlogn)
      // My method: index the start/end index of each v1 within wedge 
      // and just search within there.
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int j = 0; // index of matching edge
      bool match_found = false;
      
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Old linear search for next wedge
      // I.e. exhaustive search through *all* wedges to find an unused 
      // wedge which starts with (v2, v3, *) to match current (v1, v2, v3)
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // for (j = 0; j < n_dir_edges; j++) {
      //   if (!wedge[j].used &&
      //       wedge[j].v1 == wedge[i].v2 &&
      //       wedge[j].v2 == wedge[i].v3) {
      //     wedge[j].used = true;
      //     match_found = true;
      //     break;
      //   }
      // }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New lineear search for next matching wedge but with pre-calculated 
      // offset to to where (v2, *, *) starts in the wedge list 
      // I.e. only search linearly through (v2, *, *) wedges.
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int this_v1 = wedge[i].v2;
      for (j = bounds[this_v1].first; j <= bounds[this_v1].last; j++) {
        if (!wedge[j].used &&
            wedge[j].v1 == wedge[i].v2 &&   // this check should always return true. leaving for sanity reasons
            wedge[j].v2 == wedge[i].v3) {
          wedge[j].used = true;
          match_found = true;
          break;
        }
      }
      
      if (!match_found) {
        Rf_error("Internal error: no 'match' wedge found");
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // If the end of this matching wedge is identical to the 
      // start of the first wedge, then we have completed a polygon!
      // 
      //  - don't capture any more vertices (capturing = FALSE)
      //  - but keep processing so we mark this wedge and the following wedge
      //    as used.
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (wedge[j].v3 == wedge[first].v1) {
        capturing = false;
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // We have not found wedge (*, v1i, v2i) which matches the 
      // starting wedge (v1i, v2i, v3i).  
      // Our polygon is complete!
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (wedge[j].v2 == wedge[first].v1 &&
          wedge[j].v3 == wedge[first].v2) {
        break;
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Continue search to find the next matching wedge
      // Set the current wedge index to be this latest matching wedge index and
      // continue
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      i = j;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // record this set of temporary verts as a polygon
    // Exaand storage if we've reached capacity for storing polygons
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (*npolys == poly_capacity) {
      poly_capacity *= 2;
      polys  = realloc(polys , (unsigned long)poly_capacity * sizeof(poly_t));
    } 
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Allocate room for the vertices in this polygon
    // Copy over the "working memory" where the polygon vertices have been
    // accumulated (i.e. 'vidx') and set the vertex count for this poly_t struct
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    polys[*npolys].v = malloc((unsigned long)nvert * sizeof(int));
    if (polys[*npolys].v == NULL) Rf_error("polys[npolys] failed allocation");
    memcpy(polys[*npolys].v, &vidx, (unsigned long)nvert * sizeof(int));
    polys[*npolys].nvert = nvert;
    
    // Increase the polygon count then go found the next one
    (*npolys)++;
  
  } // while(1): Find next unused wedge
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find largest polygon by bbox and mark for deletion.
  // The polygon extraction algorithm  extracts the exterior polyon of 
  // the outside boundary.  It is the largest polygon in the extracted polygons
  // and does not represent a voronoi cell.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int max_bbox_area_idx = -1;
  double max_bbox_area  = -1;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate polygon coordinates, bbox etc
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < *npolys; i++) {
    int nvert = polys[i].nvert;
    int *vidx = polys[i].v;
    
    polys[i].polygon_idx =  i; // Polygon number in order it was found
    polys[i].taken       =  0; // Has polygon been matched to site
    polys[i].site_idx    = -1; // Which of the original sites does this polygon match
    polys[i].cx          =  0; // Centroid x
    polys[i].cy          =  0; // Centroid y
    
    polys[i].bbox = bbox_new(); // Bounding box - used to accelerate matching to initial sites
    
    polys[i].x = malloc((size_t)nvert * sizeof(double));
    polys[i].y = malloc((size_t)nvert * sizeof(double));
    if (polys[i].x == NULL || polys[i].y == NULL) {
      Rf_error("Could not allocate mem for polys[i].x and .y");
    }
    
    // Calculate actual coordinates of each polygon vertex
    for (int j = 0; j < nvert; j++) {
      polys[i].x[j] = xvor[ vidx[j] ];
      polys[i].y[j] = yvor[ vidx[j] ];
      
      polys[i].cx += polys[i].x[j];
      polys[i].cy += polys[i].y[j];
    }
    
    // Expand bounding box to encompas all polygon vertices
    bbox_add(&polys[i].bbox, nvert, polys[i].x, polys[i].y);
    
    // compute polygon centroid
    polys[i].cx /= nvert;
    polys[i].cy /= nvert;
    
    // compute bbox area
    polys[i].bbox_area = (polys[i].bbox.xmax - polys[i].bbox.xmin) * (polys[i].bbox.ymax - polys[i].bbox.ymin);
    
    // Keep track of largest
    if (polys[i].bbox_area > max_bbox_area) {
      max_bbox_area = polys[i].bbox_area;
      max_bbox_area_idx = i;
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Overwrite the largest polygon with the last polygon
  // Adjust the 'npoly' count so the final poly in 'polys' is ignored.
  // Still need to keep it around when we free it in 'free_polys'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  *npolys -= 1;
  poly_t tmp = polys[max_bbox_area_idx];
  polys[max_bbox_area_idx] = polys[*npolys];
  polys[*npolys] = tmp;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(bounds);
  free(wedge);
  free(edge);
  
  return polys;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Create 'poly_t *polys' object and then translate into final 
// R SEXP nested list representation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP extract_polygons_internal(int n_vor_verts, double *xvor, double *yvor, 
                               int n_vor_segs, int *v1vor, int *v2vor,
                               int n_sites, double *xsite, double *ysite) {
  
  int nprotect = 0;
  int npolys = 0;
  
  bool matching_polygons_to_sites = n_sites > 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Extract polygons
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  poly_t *polys = extract_polygons_core(n_vor_verts, xvor, yvor, 
                                        n_vor_segs, v1vor, v2vor, 
                                        &npolys);
  
  if (npolys == 0) {
    return R_NilValue;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Ensure n_sites (if given is the same as the number of valid polys)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (matching_polygons_to_sites && n_sites != npolys) {
    Rf_error("extract-polygons: internal error: n_sites != n_valid_polys: %i != %i\n", n_sites, npolys);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // If seed points given then we will re-organise polygons in the same order as
  // the seed points.  To do so, for each polygon's "poly_t" struct
  // annotate "poly.site_idx" with the matching site index for each polygon
  //
  // If no sites given, just return polygons in the order
  // they were found
  //
  // Remove any deleted polygons while doing this.
  //  e.g. may have been deleted for being the bounding region rather than a
  //       single polygon
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (matching_polygons_to_sites) {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // For each site, find the polygon it is contained within.
    // The polygon 'site_idx' is set to this site index.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (int i = 0; i < n_sites; i++) {
      int res = find_matching_polygon_for_site(i, xsite[i], ysite[i], npolys, polys);
      if (res < 0) {
        Rf_warning("extract_polygons_internal(): Seed [%i] unmatched\n", i);
      }
    }
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // if    matching_polygons_to_sites
  // then  Calculate cumulative vertices prior to this one
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cumulative_verts[npolys];
  
  for (int i = 0; i < npolys; i++) {
    cumulative_verts[i] = 0;
  }

  int tt = 0;
  
  if (matching_polygons_to_sites) {
    for (int i = 0; i < npolys; i++) {
      if (polys[i].site_idx < 0) continue;
      cumulative_verts[polys[i].site_idx] = polys[i].nvert;
      tt += polys[i].nvert;
    }
  }  else {
    // not matching polygons to sites
    for (int i = 0; i < npolys; i++) {
      cumulative_verts[i] = polys[i].nvert;
      tt += polys[i].nvert;
    }
  }
  
  // calculate all vertices up to each polygon
  // this will be used to access the start location to fill
  int total_verts = 0;
  for (int i = 0; i < npolys; i++) {
    int tmp = cumulative_verts[i];
    cumulative_verts[i] = total_verts;
    total_verts += tmp;
  }
  
  if (tt != total_verts) {
    Rprintf("total_verts != tt.  %i != %i\n", total_verts, tt);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Data.frame of coordinates of all polygons
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP idx_ = PROTECT(Rf_allocVector(INTSXP , total_verts)); nprotect++;
  SEXP xv_  = PROTECT(Rf_allocVector(REALSXP, total_verts)); nprotect++;
  SEXP yv_  = PROTECT(Rf_allocVector(REALSXP, total_verts)); nprotect++;
  SEXP v_   = PROTECT(Rf_allocVector(INTSXP , total_verts)); nprotect++;
  
  SEXP coords_ = PROTECT(
    create_named_list(
      4,
      "idx", idx_,
      "x"  , xv_,
      "y"  , yv_,
      "v"  , v_
    )
  ); nprotect++;
  set_df_attributes(coords_);
  
  int *idx   = INTEGER(idx_);
  double *xv = REAL(xv_);
  double *yv = REAL(yv_);
  int *v     = INTEGER(v_);


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // For each polygon, slot its vertices into the right location in
  // (xv, yv).  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < npolys; i++) {
    
    int site_idx;
    if (matching_polygons_to_sites) {
      site_idx = polys[i].site_idx;
    } else {
      site_idx = i;
    }
    
    // how many vertices come before this polygon?
    int loc = cumulative_verts[site_idx];
    
    // copy poly_t data into R vectors at the location starting *after*
    // all the vertices which come before it.
    memcpy(xv + loc, polys[i].x, (size_t)polys[i].nvert * sizeof(double));
    memcpy(yv + loc, polys[i].y, (size_t)polys[i].nvert * sizeof(double));
    memcpy(v  + loc, polys[i].v, (size_t)polys[i].nvert * sizeof(int));
    
    for (int j = 0; j < polys[i].nvert; j++) {
      idx[loc + j] = site_idx ;
    }
  }
  
  convert_indexing_c_to_r(idx_);
  convert_indexing_c_to_r(v_);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Bounding boxes
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP bbox_idx_ = PROTECT(Rf_allocVector(INTSXP , npolys)); nprotect++;
  SEXP xmin_     = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  SEXP ymin_     = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  SEXP xmax_     = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  SEXP ymax_     = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  SEXP area_     = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  
  SEXP bbox_ = PROTECT(
    create_named_list(
      6,
      "idx", bbox_idx_,
      "xmin", xmin_,
      "ymin", ymin_,
      "xmax", xmax_,
      "ymax", ymax_,
      "area", area_
    )
  ); nprotect++;
  set_df_attributes(bbox_);
  
  int *bbox_idx = INTEGER(bbox_idx_);
  double *xmin = REAL(xmin_);
  double *ymin = REAL(ymin_);
  double *xmax = REAL(xmax_);
  double *ymax = REAL(ymax_);
  double *area = REAL(area_);
  
  for (int i = 0; i < npolys; i++) {
    int site_idx;
    if (matching_polygons_to_sites) {
      site_idx = polys[i].site_idx;
    } else {
      site_idx = i;
    }
    
    bbox_idx[i] = i + 1;
    xmin[site_idx] = polys[i].bbox.xmin;
    ymin[site_idx] = polys[i].bbox.ymin;
    xmax[site_idx] = polys[i].bbox.xmax;
    ymax[site_idx] = polys[i].bbox.ymax;
    area[site_idx] = polys[i].bbox_area;
    
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Centroids
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP centroid_idx_ = PROTECT(Rf_allocVector(INTSXP, npolys)); nprotect++;
  SEXP cx_ = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  SEXP cy_ = PROTECT(Rf_allocVector(REALSXP, npolys)); nprotect++;
  
  SEXP centroids_ = PROTECT(
    create_named_list(
      3,
      "idx", centroid_idx_,
      "x"  , cx_,
      "y"  , cy_
    )
  ); nprotect++;
  set_df_attributes(centroids_);
  
  int *centroid_idx = INTEGER(centroid_idx_);
  double *cx = REAL(cx_);
  double *cy = REAL(cy_);
  
  for (int i = 0; i < npolys; i++) {
    int site_idx;
    if (matching_polygons_to_sites) {
      site_idx = polys[i].site_idx;
    } else {
      site_idx = i;
    }
    
    centroid_idx[i] = i + 1;
    cx[site_idx] = polys[i].cx;
    cy[site_idx] = polys[i].cy;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Combine all polygon information in a list
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP npolys_ = PROTECT(Rf_ScalarInteger(npolys)); nprotect++;
  
  SEXP polys_ = PROTECT(
    create_named_list(
      4,
      "npolygons" , npolys_,
      "coords"    , coords_,
      "bbox"      , bbox_,
      "centroids" , centroids_
    )
  ); nprotect++;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free_polys(npolys, polys);
  UNPROTECT(nprotect);
  return polys_;
}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP extract_polygons_(SEXP x_, SEXP y_, SEXP v1_, SEXP v2_, SEXP xseed_, SEXP yseed_, SEXP verbosity_) {

  if (Rf_length(x_) == 0 || Rf_length(x_) != Rf_length(y_)) {
    Rf_error("merge_vertices_() x & y must be equal length");
  }
  if (Rf_length(v1_) == 0) {
    Rf_warning("merge_vertices_(): No vertices to process\n");
    return R_NilValue;
  }
  if (Rf_length(v1_) != Rf_length(v2_)) {
    Rf_error("merge_vertices_(): v1 & v2 must be equal length");
  }
  
  int *v1 = create_c_index(v1_);
  int *v2 = create_c_index(v2_);

  SEXP polys_ = R_NilValue;
  
  if (Rf_isNull(xseed_) || Rf_isNull(yseed_)) {  
    // Calculate compact unordered polygons
    polys_ = extract_polygons_internal(
      Rf_length(x_), REAL(x_), REAL(y_), // voronoi vertices
      Rf_length(v1_), v1, v2,            // voronoi edges
      0, NULL, NULL                   // seed points 
    );
  } else {
    // Calculate polygons to match the seed points
    polys_ = extract_polygons_internal(
      Rf_length(x_), REAL(x_), REAL(y_),            // voronoi vertices
      Rf_length(v1_), v1, v2,                       // voronoi edges
      Rf_length(xseed_), REAL(xseed_), REAL(yseed_) // seed points 
    );
  }
  
  
  free(v1);
  free(v2);
  return polys_;
}



