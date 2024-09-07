
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
// is the given point to the left of the line?
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static bool point_left_of_line(double x, double y, 
                               double x1, double y1, double x2, double y2) {
  return (x2 - x1) * (y  - y1) - (y2 - y1) * (x - x1) > 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Is the point within the polygon?
//
// The handedness of the polygon doesn't matter
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool point_in_polygon_core_(double x, double y, int N, double *xp, double *yp) {
  
  bool status = point_left_of_line(x, y, xp[0], yp[0], xp[1], yp[1]);
  bool new_status = false;
  
  for (int i = 1; i < N - 1; i++) {
    new_status = point_left_of_line(x, y, xp[i], yp[i], xp[i+1], yp[i+1]);
    if (new_status != status) return false;
    status = new_status;
  }
  
  // Wrap around to start
  new_status = point_left_of_line(x, y, xp[N-1], yp[N-1], xp[0], yp[0]);
  if (new_status != status) return false;
  
  return true;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP points_in_convex_polygon_(SEXP x_, SEXP y_, SEXP xp_, SEXP yp_) {
  
  if (length(x_) != length(y_)) {
    error("points_in_convex_polygon_(): x and y must be same length");
  }
  if (length(x_) == 0) {
    error("points_in_convex_polygon_(): need at least one point");
  }
  
  if (length(xp_) != length(yp_)) {
    error("points_in_convex_polygon_(): xp and yp must be same length");
  }
  if (length(xp_) < 3) {
    error("points_in_convex_polygon_(): not a polygon");
  }
  
  SEXP res_ = PROTECT(allocVector(LGLSXP, length(x_)));
  int *resp = LOGICAL(res_);
  double *x = REAL(x_);
  double *y = REAL(y_);
  double *xp = REAL(xp_);
  double *yp = REAL(yp_);
  
  for (int i = 0; i < length(x_); i++) {
    resp[i] = point_in_polygon_core_(x[i], y[i], length(xp_), xp, yp);
  }
  
  UNPROTECT(1);
  return res_;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Match a point to a polygon
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int match_polygons_to_seed_points(int point_idx, double x, double y, int npolys, poly_t *polys) {
  
  // Rprintf("PIP ---------------------------------\n");
  
  // loop over all polygons
  for (int i = 0; i < npolys; i++) {
    
    // Rprintf("PIP: %i in %i [%i]\n", point_idx, poly.polygon_idx, poly.deleted);
    
    if (polys[i].deleted) continue;
    
    bbox_t bbox = polys[i].bbox;
    
    if (x > bbox.xmin && x < bbox.xmax && y > bbox.ymin && y < bbox.ymax) {
      // point intersects bounding box for this polygon
      // bool point_in_polygon_core_(double x, double y, int N, double *xp, double *yp)
      bool pip = point_in_polygon_core_(x, y, polys[i].nvert, polys[i].x, polys[i].y);
      if (pip) {
        // point is in this polygon!
        if (polys[i].taken) {
          error("match_polygons_to_seed_points(): Point (%f, %f) matched polygon [%i] already taken by point [%i]", 
                x, y, polys[i].polygon_idx, polys[i].point_idx);
        }
        polys[i].taken = true;
        polys[i].point_idx = point_idx;
        // Rprintf("Point %3i => Polygon %3i\n", point_idx, polys[i].polygon_idx);
        return polys[i].polygon_idx;
      }
    }
    
  } // next polygon
  
  // Rprintf("match_polygons_to_seed_points(): Point %i (%f, %f) did not match any polygons!\n", point_idx, x, y);
  return -1;
}













