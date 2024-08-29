
// #define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern SEXP voronoi_ (SEXP x_, SEXP y_);
extern SEXP delaunay_(SEXP x_, SEXP y_);

static const R_CallMethodDef CEntries[] = {
  
  {"voronoi_" , (DL_FUNC) &voronoi_ , 2},
  {"delaunay_", (DL_FUNC) &delaunay_, 2},
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



