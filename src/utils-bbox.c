
#define R_NO_REMAP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "utils.h"
#include "utils-bbox.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialise an empty boudning box
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bbox_t bbox_new(void) {
  bbox_t bbox = {
    .xmin =  INFINITY,
    .ymin =  INFINITY,
    .xmax = -INFINITY,
    .ymax = -INFINITY
  };
  
  return bbox;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Add new points to a bounding box
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void bbox_add(bbox_t *bbox, int n, double *x, double *y) {
  
  for (int i = 0; i < n; i++) {
    
    if (x[i] < bbox->xmin) bbox->xmin = x[i];
    if (x[i] > bbox->xmax) bbox->xmax = x[i];
    if (y[i] < bbox->ymin) bbox->ymin = y[i];
    if (y[i] > bbox->ymax) bbox->ymax = y[i];
  }
  
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Expand a bbox
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void bbox_expand(bbox_t *bbox, double factor) {
  
  if (factor <= 0) {
    Rf_error("box_expand(): factor must be > 0");
  }
  
  double w = bbox->xmax - bbox->xmin;
  double h = bbox->ymax - bbox->ymin;
  
  double maxd = w > h ? w : h;
  if (maxd == 0) {
    // bbox has no area!
    maxd = 1;
  }
  
  if (w == 0) {
    bbox->xmin = -maxd / 2;
    bbox->xmax =  maxd / 2;
  } else {
    bbox->xmin -= w * factor / 2;
    bbox->xmax += w * factor / 2;
  }
  
  if (h == 0) {
    bbox->ymin = -maxd / 2;
    bbox->ymax =  maxd / 2;
  } else {
    bbox->ymin -= h * factor / 2;
    bbox->ymax += h * factor / 2;
  }
  
}


