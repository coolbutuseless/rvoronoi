
// #define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern SEXP voronoi_ (SEXP x_, SEXP y_);
extern SEXP delaunay_(SEXP x_, SEXP y_);

extern SEXP merge_vertices_  (SEXP x_, SEXP y_, SEXP v1_, SEXP v2_, SEXP tol_, SEXP verbosity_);
extern SEXP extract_polygons_(SEXP x_, SEXP y_, SEXP v1_, SEXP v2_,            SEXP verbosity_);

extern SEXP points_in_convex_polygon_ (SEXP x_, SEXP y_, SEXP xp_, SEXP yp_);

static const R_CallMethodDef CEntries[] = {
  
  {"voronoi_" , (DL_FUNC) &voronoi_ , 2},
  {"delaunay_", (DL_FUNC) &delaunay_, 2},
  
  {"merge_vertices_"  , (DL_FUNC) &merge_vertices_  , 6},
  {"extract_polygons_", (DL_FUNC) &extract_polygons_, 5},
  
  {"points_in_convex_polygon_", (DL_FUNC) &points_in_convex_polygon_, 4},
  
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



