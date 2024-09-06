
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


static bool point_left_of_line(double x, double y, 
                               double x1, double y1, double x2, double y2) {
  return (x2 - x1) * (y  - y1) - (y2 - y1) * (x - x1) > 0;
}


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


SEXP point_in_polygon_(SEXP x_, SEXP y_, SEXP xp_, SEXP yp_) {
  
  return ScalarLogical(
    point_in_polygon_core_(asReal(x_), asReal(y_), 
                           length(xp_), REAL(xp_), REAL(yp_))
  );
  
}
