
#define R_NO_REMAP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "utils-bbox.h"
#include "R-polygon-matching.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// is the given point (x, y) to the left of the line  (x1,y1) - (x2,y2) ?
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static bool point_left_of_line(double x, double y, 
                               double x1, double y1, double x2, double y2) {
  return (x2 - x1) * (y  - y1) - (y2 - y1) * (x - x1) > 0;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Is the point within the polygon?
//
// (For this test, the handedness of the polygon doesn't matter)
//
// Method:
// Check each edge around the polygon returns the same value for 
// 'point_left_of_line' in regards to the test point.
// IF all the same, then the point is inside the polygon!
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool point_in_polygon_core_(double x, double y, int N, double *xp, double *yp) {
  
  bool status = point_left_of_line(x, y, xp[0], yp[0], xp[1], yp[1]);
  bool new_status = false;
  
  for (int i = 1; i < N - 1; i++) {
    new_status = point_left_of_line(x, y, xp[i], yp[i], xp[i+1], yp[i+1]);
    if (new_status != status) return false;
    status = new_status;
  }
  
  // Wrap around to first point
  new_status = point_left_of_line(x, y, xp[N-1], yp[N-1], xp[0], yp[0]);
  if (new_status != status) return false;
  
  return true;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP points_in_convex_polygon_(SEXP x_, SEXP y_, SEXP xp_, SEXP yp_) {
  
  if (Rf_length(x_) != Rf_length(y_)) {
    Rf_error("points_in_convex_polygon_(): x and y must be same length");
  }
  if (Rf_length(x_) == 0) {
    Rf_error("points_in_convex_polygon_(): need at least one point");
  }
  
  if (Rf_length(xp_) != Rf_length(yp_)) {
    Rf_error("points_in_convex_polygon_(): xp and yp must be same length");
  }
  if (Rf_length(xp_) < 3) {
    Rf_error("points_in_convex_polygon_(): not a polygon");
  }
  
  SEXP res_ = PROTECT(Rf_allocVector(LGLSXP, Rf_length(x_)));
  int *resp = LOGICAL(res_);
  double *x = REAL(x_);
  double *y = REAL(y_);
  double *xp = REAL(xp_);
  double *yp = REAL(yp_);
  
  for (int i = 0; i < Rf_length(x_); i++) {
    resp[i] = point_in_polygon_core_(x[i], y[i], Rf_length(xp_), xp, yp);
  }
  
  UNPROTECT(1);
  return res_;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Match a point to a polygon
// @param site_idx site index
// @param x,y coords of site point
// @param npolys number of polygons
// @param polys vector of poly_t structs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int find_matching_polygon_for_site(int site_idx, double x, double y, int npolys, poly_t *polys) {
  
  // loop over all polygons
  //     if point inside bounding box
  //         if point inside polygon
  //              mark polygon with the 'point_idx'
  //              return polygon_idx
  // otherwise 
  //    return -1
  for (int i = 0; i < npolys; i++) {
    
    // Rprintf("PIP: %i in %i [%i]\n", site_idx, poly.polygon_idx, poly.deleted);
    
    // No need to search this polygon if it has been deleted, or if it
    // has been claimed by another point.  
    // Remember: For a voronoi tessellation a seed point can only match one
    // polygon (and vice versa)
    if (polys[i].taken) continue;
    
    bbox_t bbox = polys[i].bbox;
    
    if (x > bbox.xmin && x < bbox.xmax && y > bbox.ymin && y < bbox.ymax) {
      // The point intersects bounding box for this polygon!
      // Now check if it lies inside the actual polygon
      bool pip = point_in_polygon_core_(x, y, polys[i].nvert, polys[i].x, polys[i].y);
      if (pip) {
        // point is in this polygon!
        polys[i].taken = true;
        polys[i].site_idx = site_idx;
        // Rprintf("Point %3i => Polygon %3i\n", site_idx, polys[i].polygon_idx);
        return polys[i].polygon_idx;
      }
    }
    
  } // next polygon
  
  return -1;
}













