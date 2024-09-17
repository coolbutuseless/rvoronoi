
// #define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern SEXP voronoi_ (SEXP x_, SEXP y_, SEXP calc_polygons_, SEXP match_sites_);
extern SEXP delaunay_(SEXP x_, SEXP y_, SEXP calc_polygons_, SEXP calc_areas_, SEXP calc_segments_);

extern SEXP merge_vertices_(SEXP x_, SEXP y_, 
                            SEXP line_, SEXP v1_, SEXP v2_, 
                            SEXP tol_, SEXP verbosity_);

extern SEXP extract_polygons_(SEXP x_, SEXP y_, SEXP v1_, SEXP v2_, SEXP xseed_, SEXP yseed_, SEXP verbosity_);

extern SEXP points_in_convex_polygon_ (SEXP x_, SEXP y_, SEXP xp_, SEXP yp_);

extern SEXP bound_infinite_segments_(
    SEXP xmin_, SEXP ymin_, SEXP xmax_, SEXP ymax_,
    SEXP x_, SEXP y_,
    SEXP a_, SEXP b_, SEXP c_,
    SEXP li_, SEXP v1_, SEXP v2_);

static const R_CallMethodDef CEntries[] = {
  
  {"voronoi_" , (DL_FUNC) &voronoi_ , 4},
  {"delaunay_", (DL_FUNC) &delaunay_, 5},
  
  {"merge_vertices_"  , (DL_FUNC) &merge_vertices_  , 7},
  {"extract_polygons_", (DL_FUNC) &extract_polygons_, 7},
  
  {"points_in_convex_polygon_", (DL_FUNC) &points_in_convex_polygon_, 4},
  {"bound_infinite_segments_", (DL_FUNC) &bound_infinite_segments_, 12},
  
  {NULL , NULL, 0}
};


void R_init_rvoronoi(DllInfo *info) {
  R_registerRoutines(
    info,      // DllInfo
    NULL,      // .C
    CEntries,  // .Call
    NULL,      // Fortran
    NULL       // External
  );
  R_useDynamicSymbols(info, FALSE);
}



