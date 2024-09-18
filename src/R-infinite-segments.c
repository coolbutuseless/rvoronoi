

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
#include "R-infinite-segments.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Struct for organizing points along the perimeter 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  double x; 
  double y;
  int v;
} point_t;



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Helper function used with sort()
// sort point_t along y-direction only
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int horizontal_comparison(const void *v1, const void *v2) {
  
  point_t *s1 = (point_t *)v1;
  point_t *s2 = (point_t *)v2;

  if (s1->x < s2->x)
    return (-1);
  if (s1->x > s2->x)
    return (1);
  return (0);
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Helper function used with sort()
// sort point_t along y-direction only
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int vertical_comparison(const void *v1, const void *v2) {
  
  point_t *s1 = (point_t *)v1;
  point_t *s2 = (point_t *)v2;
  
  if (s1->y < s2->y)
    return (-1);
  if (s1->y > s2->y)
    return (1);
  return (0);
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate the upper limit of bounded vertices and bounded segments
// required if the semi-infinite and infinite segments from Fortune's 
// algorithm are to be bound within rectangular limits
//
// @param nsegs number of voronoi segements (in total)
// @param v1,v2 vertex indices for the ends of each segment
// @param nbverts (return value) number of bounded vertices required for
//        intersection of infinite rays with rectangular boundary
// @param nbsegs (return value) number of segments which will be generated by
//        intersection infinite segments with rectangular boundary
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void calc_space_for_bound_infinite_segments(int nsegs, int *v1, int *v2, int *nbverts, int *nbsegs) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate the upper limit on the number of bounded segments 
  // and vertices required if infinite segments are intersected with
  // rectangular boundary
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  *nbsegs  = 0;
  *nbverts = 0;
  
  for (int i = 0; i < nsegs; i++) {
    // Ignore finite segments
    if (valid_idx(v1[i]) && valid_idx(v2[i])) continue; 
    
    // Each infinite segment adds one segment and one vertex
    (*nbsegs)++;
    (*nbverts)++;
    
    // If this is a truly infinite segment (infinite in both directions)
    // then it adds a second intersection vertex
    if (!valid_idx(v1[i]) && !valid_idx(v2[i])) {
      // Infinite segment adds 2 vertices
      (*nbverts)++;
    }
  }
  
  
  // Rprintf("Found %2i rays to bound. Needing %2i vertices\n", *nbsegs, *nbverts);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Around the perimeter, need to add 
  //   -  4 vertices for each of the corners of the boundary
  //   -  (4 + nbverts) new segments around the boundary
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  (*nbsegs)  += 4 + *nbverts;
  (*nbverts) += 4;
  
  
  // Rprintf("Total external segements = %2i. External Vertices = %2i\n", *nbsegs, *nbverts);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Find the closest distance between a single point (x, y) and 
// all the 'n' points in (xs, ys)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double closest_distance(double x, double y, int n, double *xs, double *ys) {
  
  double res = INFINITY;
  
  for (int i = 0; i < n; i++) {
    
    double d = (xs[i] - x) * (xs[i] - x) + (ys[i] - y) * (ys[i] - y);
    
    if (d < res) res = d;
  }
  
  return res;
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Find infinte segments and intersect them with the boundary
//
//  User must pre-allocat xb,yb,rv1,rv2
//
// @param bounds bounding box which encompases all sites and voronoi vertices
// @param nverts number of voronoi vertices
// @param x,y vertex coordinates
// @param nsegs total number of voronoi segments
// @param li line index for each of the segements
// @param v1,v2 vertex indicies for each end of the segment
// @param nlines number of line definitions
// @param a,b,c line definition of the form ax + by = c
// @param nbverts (input) total allocated size of xb,yb.
//                (output) total actual number of boundary points
// @param xb,yb coordinates of intersection points with the rectangular boundary
// @param nbsegs (input) total allocated size of rv1, rv2
// @param rv1,rv2 the finite segments created by the truncaction of the 
//                semi-infinite and infinite segments where they intersect
//                the rectangular boundary
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void bound_infinite_segments(
    bbox_t *bounds,
    int nverts, double *x, double *y,            // Voronoi vertices
    int nsegs , int *li, int *v1, int *v2,       // Voronoi edges
    int nlines, double *a, double *b, double *c, // Voronoi lines
    int *nbverts, double *xb, double *yb,         // boundary intersection points (and boundary corners)
    int *nbsegs, int *rv1, int *rv2) {            // bounded rays, and perimeter segments
  
  int verbosity = 0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Hold info about the two possible intercepts between rays and rectangle
  // 
  // Infinite segments always intersect the boundary at two points
  // Semi-infinite rays intersect at a single point
  // 
  // Method:
  //    Intersect the line with the 4 segments making up the boundary
  //    Keep the ones whih are actually internal to the region
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
  //    - add the vertices to (xb, yb) 
  //    - add the segments to (rv1, rv2)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < nsegs; i++) {
    
    // Ignore finite segments
    if (valid_idx(v1[i]) && valid_idx(v2[i])) continue;
    
    int type      = 0;
    int v1_anchor = 0;
    double xi     = 0;
    double yi     = 0;
    int line_idx  = 0;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set the intersection type depending on whether the known vertex
    // is the first or second vertex on the segment.
    // if this is a true infinite segment with neither vertex known then type = 0
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (valid_idx(v1[i])) {
      type = 2;
      v1_anchor = v1[i];
    } else if (valid_idx(v2[i])) {
      type = 1;
      v1_anchor = v2[i];
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // There can only be 2 intercept points which lie on the boundary rectangle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    intercept_t intercept[2] = { 0 };
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intersection with left boundary edge
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
    // Intercept with right boundary edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xi = bounds->xmax;
    yi = (c[line_idx]- a[line_idx] * xi) / b[line_idx];
    if (yi >= bounds->ymin && yi <= bounds->ymax) {
      intercept[1].valid = true;
      intercept[1].x     = xi;
      intercept[1].y     = yi;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intercept with lower boundary edge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double yi1 = bounds->ymin;
    double xi1 = (c[line_idx]- b[line_idx] * yi1) / a[line_idx];
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Intercept with upper boundary edge
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Depending on how the ray is anchored (i.e. left, right or both)
    // add vertex/vertices for boundary intersection and add segment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (verbosity > 0) {
      Rprintf("Ray: %i   Vert: %i\n", ray_idx, vert_idx);
    }
    
    if (type == 0) {
      // double ended unbounded ray. 
      // Add vertices for two interesections with boundary
      xb[vert_idx] = intercept[0].x;
      yb[vert_idx] = intercept[0].y;
      vert_idx++;
      xb[vert_idx] = intercept[1].x;
      yb[vert_idx] = intercept[1].y;
      vert_idx++;
      // Add segment between these two intersections
      rv1[ray_idx] = nverts + vert_idx - 2;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else if (type == 1) {
      // ray now bounded on the left
      // Add vertex for boundary intersection
      xb[vert_idx] = intercept[0].x;
      yb[vert_idx] = intercept[0].y;
      vert_idx++;
      // Add segment between anchor point and intersection
      rv1[ray_idx] = v1_anchor;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else if (type == 2) {
      // ray not bounded on the right
      // Add vertex for boundary intersection
      xb[vert_idx] = intercept[1].x;
      yb[vert_idx] = intercept[1].y;
      vert_idx++;
      // Add segment between anchor point and intersection
      rv1[ray_idx] = v1_anchor;
      rv2[ray_idx] = nverts + vert_idx - 1;
      ray_idx++;
    } else {
      error("Impossible 12");
    }
  } // next seg
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Add vertices for boundary.
  // Note: for pathological cases a point can intersect the corner
  //       of the rectangular boundary.  
  //       The initial code naively added all points on each boundary edge in 
  //       sequence, which lead to corner points being added *TWICE*
  // THe new method checks whether each of the boundary points are 
  // already in the list of new vertices before adding them.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double d = 0;
  d = closest_distance(bounds->xmin, bounds->ymin, vert_idx, xb, yb);
  if (d > 0) {
    xb[vert_idx] = bounds->xmin; yb[vert_idx] = bounds->ymin; vert_idx++;
  }
  
  d = closest_distance(bounds->xmax, bounds->ymin, vert_idx, xb, yb);
  if (d > 0) {
    xb[vert_idx] = bounds->xmax; yb[vert_idx] = bounds->ymin; vert_idx++;
  }
  
  d = closest_distance(bounds->xmax, bounds->ymax, vert_idx, xb, yb);
  if (d > 0) {
    xb[vert_idx] = bounds->xmax; yb[vert_idx] = bounds->ymax; vert_idx++;
  }
  
  d = closest_distance(bounds->xmin, bounds->ymax, vert_idx, xb, yb);
  if (d > 0) {
    xb[vert_idx] = bounds->xmin; yb[vert_idx] = bounds->ymax; vert_idx++;
  }
  
  
  if (vert_idx > *nbverts) {
    error("Expecting vert_idx <= nbverts :: %i == %i", vert_idx, *nbverts);
  }
  
  // How many verts did we actaully add
  *nbverts = vert_idx;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find all points along top boundary edge
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int npoints = 0;
  point_t points[*nbverts];
  
  for (int i = 0; i < *nbverts; i++) {
    if (yb[i] == bounds->ymax) {
      points[npoints].x = xb[i];
      points[npoints].y = yb[i];
      points[npoints].v = nverts + i;
      npoints++;
    }
  }
  if (verbosity > 0) {
    Rprintf("Points on top perimeter: %i\n", npoints);
  }
  qsort(points, npoints, sizeof(point_t), horizontal_comparison);
  
  // Add segments
  for (int i = 0; i < npoints - 1; i++) {
    rv1[ray_idx] = points[i    ].v;
    rv2[ray_idx] = points[i + 1].v;
    ray_idx++;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find all points along bottom boundary edge
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npoints = 0;
  
  for (int i = 0; i < *nbverts; i++) {
    if (yb[i] == bounds->ymin) {
      points[npoints].x = xb[i];
      points[npoints].y = yb[i];
      points[npoints].v = nverts + i;
      npoints++;
    }
  }
  if (verbosity > 0) {
    Rprintf("Points on bottom perimeter: %i\n", npoints);
  }
  qsort(points, npoints, sizeof(point_t), horizontal_comparison);
  
  // Add segments
  for (int i = 0; i < npoints - 1; i++) {
    rv1[ray_idx] = points[i    ].v;
    rv2[ray_idx] = points[i + 1].v;
    ray_idx++;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find all points along left boundary edge
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npoints = 0;
  
  for (int i = 0; i < *nbverts; i++) {
    if (xb[i] == bounds->xmin) {
      points[npoints].x = xb[i];
      points[npoints].y = yb[i];
      points[npoints].v = nverts + i;
      npoints++;
    }
  }
  if (verbosity > 0) {
    Rprintf("Points on left perimeter: %i\n", npoints);
  }
  qsort(points, npoints, sizeof(point_t), vertical_comparison);
  
  // Add segments
  for (int i = 0; i < npoints - 1; i++) {
    rv1[ray_idx] = points[i    ].v;
    rv2[ray_idx] = points[i + 1].v;
    ray_idx++;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find all points along right boundary edge
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npoints = 0;
  
  for (int i = 0; i < *nbverts; i++) {
    if (xb[i] == bounds->xmax) {
      points[npoints].x = xb[i];
      points[npoints].y = yb[i];
      points[npoints].v = nverts + i;
      npoints++;
    }
  }
  if (verbosity > 0) {
    Rprintf("Points on right perimeter: %i\n", npoints);
  }
  qsort(points, npoints, sizeof(point_t), vertical_comparison);
  
  // Add segments
  for (int i = 0; i < npoints - 1; i++) {
    rv1[ray_idx] = points[i    ].v;
    rv2[ray_idx] = points[i + 1].v;
    ray_idx++;
  }
  
  
  if (ray_idx > *nbsegs) {
    error("Expecting ray_idx <= nbsegs :: %i == %i", ray_idx, *nbsegs);
  }
  
  *nbsegs  = ray_idx;
  *nbverts = vert_idx;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP bound_infinite_segments_(
    SEXP xmin_, SEXP ymin_, SEXP xmax_, SEXP ymax_,
    SEXP x_, SEXP y_,
    SEXP a_, SEXP b_, SEXP c_,
    SEXP li_, SEXP v1_, SEXP v2_) {
  
  int nprotect = 0;
  
  if (length(x_) != length(y_)) {
    error("bound_infinite_segments_(): bad length for x/y");
  }
  
  if (length(v1_) == 0 || length(v1_) != length(v2_) || 
      length(li_) != length(v1_)) {
    error("bound_infinite_segments_(): bad length for li/v1/v2");
  }
  
  if (length(a_) == 0 || length(a_) != length(b_) || 
      length(a_) != length(c_)) {
    error("bound_infinite_segments_(): bad length for a/b/c");
  }
  
  bbox_t bounds = {
    .xmin = asReal(xmin_),
    .ymin = asReal(ymin_),
    .xmax = asReal(xmax_),
    .ymax = asReal(ymax_)
  };
  
  int nbverts = 0;
  int nbsegs = 0;
  
  calc_space_for_bound_infinite_segments(length(v1_), INTEGER(v1_), INTEGER(v2_), &nbverts, &nbsegs);
  SEXP xb_ = PROTECT(allocVector(REALSXP, nbverts)); nprotect++;
  SEXP yb_ = PROTECT(allocVector(REALSXP, nbverts)); nprotect++;
  
  SEXP rv1_ = PROTECT(allocVector(INTSXP, nbsegs)); nprotect++;
  SEXP rv2_ = PROTECT(allocVector(INTSXP, nbsegs)); nprotect++;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert R 1-indexing to C 0-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *v1 = create_c_index(v1_);
  int *v2 = create_c_index(v2_);
  int *li = create_c_index(li_);

  
  bound_infinite_segments(
    &bounds,
    length(x_), REAL(x_), REAL(y_),
    length(v1_), li, v1, v2,
    length(a_), REAL(a_), REAL(b_), REAL(c_),
    &nbverts, REAL(xb_), REAL(yb_),
    &nbsegs, INTEGER(rv1_), INTEGER(rv2_)
  );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New vertices
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP verts_ = PROTECT(
    create_named_list(2, "x", xb_, "y", yb_)
  ); nprotect++;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Segment
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  convert_indexing_c_to_r(rv1_);
  convert_indexing_c_to_r(rv2_);

  
  SEXP segments_ = PROTECT(
    create_named_list(3, "line", R_NilValue, "v1", rv1_, "v2", rv2_)
  ); nprotect++;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // list result
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(
    create_named_list(2, "vertex", verts_, "segment", segments_)
  ); nprotect++;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(v1);
  free(v2);
  free(li);
  
  UNPROTECT(nprotect);
  return res_;
}



