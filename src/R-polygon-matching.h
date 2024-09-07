
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Bounding box type
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  double xmin;
  double ymin;
  double xmax;
  double ymax;
} bbox_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Store information about the voronoi polygons 
// This will be used during the matching process where we line up 
// polygons with the original seed points
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  bbox_t bbox;      // bounding box
  double bbox_area; // bounding box area
  int *v;           // indices into voronoi vertices
  int nvert;        // Number of vertices
  double *x;        // x coordinates
  double *y;        // y coordinates
  int polygon_idx;  // index of this polygon
  int point_idx;    // index of this point
  bool taken;       // has this polygon been claimed?
  double cx;        // x centroid
  double cy;        // y centroid
  bool deleted;     // has this polygon been deleted
} poly_t;

int match_polygons_to_seed_points(int point_idx, double x, double y, int npolys, poly_t *polys);

