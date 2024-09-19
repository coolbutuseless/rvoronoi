
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
#include "R-common.h"


typedef struct {
  int v1;
  int v2;
} vpair_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int vpair_comparison(const void *v1, const void *v2) {
  
  vpair_t *s1 = (vpair_t *)v1;
  vpair_t *s2 = (vpair_t *)v2;
  
  if (s1->v1 < s2->v1) {
    return (-1);
  }
  if (s1->v1 > s2->v1) {
    return (1);
  }
  if (s1->v2 < s2->v2) {
    return (-1);
  }
  if (s1->v2 > s2->v2) {
    return (1);
  }
  
  return(0);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Delauney Triangulation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP delaunay_(SEXP x_, SEXP y_, SEXP calc_polygons_, SEXP calc_areas_, SEXP calc_segments_) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sanity Check
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(x_) != length(y_)) {
    error("x/y lengths must be the same");
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialise the calculation context
  // Do delauney? TRUE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int nprotect = 0;
  context_t ctx = { 0 };
  ctx.triangulate = 1;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Track *ALL* the allocations done via 'myalloc()'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctx.alloc_count    = 0;
  ctx.alloc_capacity = 1024;
  ctx.allocs = (void **)calloc((unsigned long)ctx.alloc_capacity, sizeof(void *));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Maximum number of triangles = 2 * n - 2 - b
  // where 'b' is number of vertices on the convex hull
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int n = length(x_);
  int max_tris = 2 * n;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertices of each triangles
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP v1_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  SEXP v2_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  SEXP v3_ = PROTECT(allocVector(INTSXP, max_tris)); nprotect++; 
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize context
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ctx.ntris = 0;
  ctx.v1 = INTEGER(v1_);
  ctx.v2 = INTEGER(v2_);
  ctx.v3 = INTEGER(v3_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Use Fortune's algo to calculate delaunay triangulation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  freeinit(&ctx.sfl, sizeof *ctx.sites);
  init_sites(&ctx, REAL(x_), REAL(y_), length(x_));
  ctx.siteidx = 0;
  geominit(&ctx);
  voronoi(&ctx, ctx.triangulate);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Indices = data.frame(v1 = integer(), v2 = integer(), v3 = integer())
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP idx_  = R_NilValue;
  SEXP area_ = R_NilValue;
  if (asLogical(calc_areas_)) {
    area_ = PROTECT(allocVector(REALSXP, ctx.ntris)); nprotect++;
    idx_  = PROTECT(create_named_list(4, "v1", v1_, "v2", v2_, "v3", v3_, "area", area_)); nprotect++;
  } else {
    idx_ = PROTECT(create_named_list(3, "v1", v1_, "v2", v2_, "v3", v3_)); nprotect++;
  }
  set_df_attributes_and_trim(idx_, ctx.ntris, max_tris);
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Area
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (asLogical(calc_areas_)) {
    double *area  = REAL(area_);
    
    double *x = REAL(x_);
    double *y = REAL(y_);
    
    for (int i = 0; i < ctx.ntris; i++) {
      area[i] = 0.5 * (
        x[ ctx.v1[i] ] * (y[ ctx.v2[i] ] - y[ ctx.v3[i] ]) +
          x[ ctx.v2[i] ] * (y[ ctx.v3[i] ] - y[ ctx.v1[i] ]) +
          x[ ctx.v3[i] ] * (y[ ctx.v1[i] ] - y[ ctx.v2[i] ]) 
      );
      area[i] = fabs(area[i]);
    }
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Vertex coordinates for each triangle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP polys_ = R_NilValue;
  
  if (asLogical(calc_polygons_)) {
    SEXP pidx_ = PROTECT(allocVector(INTSXP , 3 * ctx.ntris)); nprotect++;
    SEXP xs_   = PROTECT(allocVector(REALSXP, 3 * ctx.ntris)); nprotect++;
    SEXP ys_   = PROTECT(allocVector(REALSXP, 3 * ctx.ntris)); nprotect++;
    SEXP vs_   = PROTECT(allocVector(INTSXP , 3 * ctx.ntris)); nprotect++;
    
    double *xs = REAL(xs_);
    double *ys = REAL(ys_);
    int *pidx  = INTEGER(pidx_);
    int *vidx  = INTEGER(vs_);
    
    double *x = REAL(x_);
    double *y = REAL(y_);
    int *v1 = INTEGER(v1_);
    int *v2 = INTEGER(v2_);
    int *v3 = INTEGER(v3_);
    
    for (int i = 0; i < ctx.ntris; i++) {
      
      *xs++ = x[ v1[i] ];
      *xs++ = x[ v2[i] ];
      *xs++ = x[ v3[i] ];
      
      *ys++ = y[ v1[i] ];
      *ys++ = y[ v2[i] ];
      *ys++ = y[ v3[i] ];
      
      *pidx++ = i + 1;
      *pidx++ = i + 1;
      *pidx++ = i + 1;
      
      *vidx++ = v1[i] + 1;
      *vidx++ = v2[i] + 1;
      *vidx++ = v3[i] + 1;
    }
    
    
    polys_ = PROTECT(
      create_named_list(4, "idx", pidx_, "x", xs_, "y", ys_, "vidx", vs_)
    ); nprotect++;
    set_df_attributes(polys_);
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Each segment (each segment should only appear once even though 
  // it might take part in 2 triangles).
  // Only keep vertex pairs where va < vb
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP segs_ = R_NilValue;
  
  if (asLogical(calc_segments_)) {
    vpair_t *vpairs = malloc((long unsigned int)ctx.ntris * 3 * sizeof(vpair_t));
    if (vpairs == NULL) {
      error("couldn't allocate vpairs");
    }
    
    
    int nvpairs = 0;
    for (int i = 0; i < ctx.ntris; i++) {
      if (ctx.v1[i] < ctx.v2[i]) {
        vpairs[nvpairs].v1 = ctx.v1[i];
        vpairs[nvpairs].v2 = ctx.v2[i];
      } else {
        vpairs[nvpairs].v1 = ctx.v2[i];
        vpairs[nvpairs].v2 = ctx.v1[i];
      }
      nvpairs++;
      
      if (ctx.v2[i] < ctx.v3[i]) {
        vpairs[nvpairs].v1 = ctx.v2[i];
        vpairs[nvpairs].v2 = ctx.v3[i];
      } else {
        vpairs[nvpairs].v1 = ctx.v3[i];
        vpairs[nvpairs].v2 = ctx.v2[i];
      }
      nvpairs++;
      
      if (ctx.v3[i] < ctx.v1[i]) {
        vpairs[nvpairs].v1 = ctx.v3[i];
        vpairs[nvpairs].v2 = ctx.v1[i];
      } else {
        vpairs[nvpairs].v1 = ctx.v1[i];
        vpairs[nvpairs].v2 = ctx.v3[i];
      }
      nvpairs++;
    }
    
    qsort(vpairs, (size_t)nvpairs, sizeof(vpair_t), vpair_comparison);
    
    int nsegs = 1;
    for (int i=1; i < nvpairs; i++) {
      // Rprintf("vpair[%2i] %2i -> %2i\n", i, vpairs[i].v1, vpairs[i].v2);
      if (vpairs[i].v1 == vpairs[i-1].v1 && vpairs[i].v2 == vpairs[i-1].v2) continue;
      nsegs++;
    }
    
    
    SEXP sv1_  = PROTECT(allocVector(INTSXP , nsegs)); nprotect++;
    SEXP sv2_  = PROTECT(allocVector(INTSXP , nsegs)); nprotect++;
    SEXP x1_   = PROTECT(allocVector(REALSXP, nsegs)); nprotect++;
    SEXP y1_   = PROTECT(allocVector(REALSXP, nsegs)); nprotect++;
    SEXP x2_   = PROTECT(allocVector(REALSXP, nsegs)); nprotect++;
    SEXP y2_   = PROTECT(allocVector(REALSXP, nsegs)); nprotect++;
    SEXP dist_ = PROTECT(allocVector(REALSXP, nsegs)); nprotect++;
    segs_ = PROTECT(
      create_named_list(
        7,
        "v1", sv1_,
        "v2", sv2_,
        "x1", x1_,
        "y1", y1_,
        "x2", x2_,
        "y2", y2_,
        "len", dist_
      )
    ); nprotect++;
    set_df_attributes(segs_);
    
    int *sv1 = INTEGER(sv1_);
    int *sv2 = INTEGER(sv2_);
    double *x1 = REAL(x1_);
    double *y1 = REAL(y1_);
    double *x2 = REAL(x2_);
    double *y2 = REAL(y2_);
    double *dist = REAL(dist_);
    
    double *x = REAL(x_);
    double *y = REAL(y_);
    
    sv1[0]  = vpairs[0].v1;
    sv2[0]  = vpairs[0].v2;
    x1[0]   = x[ sv1[0] ];
    y1[0]   = y[ sv1[0] ];
    x2[0]   = x[ sv2[0] ];
    y2[0]   = y[ sv2[0] ];
    dist[0] = sqrt(
      (x2[0] - x1[0]) * (x2[0] - x1[0]) + 
        (y2[0] - y1[0]) * (y2[0] - y1[0])
    );
    nsegs = 1;
    
    for (int i=1; i < nvpairs; i++) {
      if (vpairs[i].v1 == vpairs[i-1].v1 && vpairs[i].v2 == vpairs[i-1].v2) continue;
      
      sv1[nsegs] = vpairs[i].v1;
      sv2[nsegs] = vpairs[i].v2;
      
      x1[nsegs] = x[ sv1[nsegs] ];
      y1[nsegs] = y[ sv1[nsegs] ];
      x2[nsegs] = x[ sv2[nsegs] ];
      y2[nsegs] = y[ sv2[nsegs] ];
      
      dist[nsegs] = sqrt(
        (x2[nsegs] - x1[nsegs]) * (x2[nsegs] - x1[nsegs]) +
          (y2[nsegs] - y1[nsegs]) * (y2[nsegs] - y1[nsegs])
      );
      
      nsegs++;
    }
    
    
    free(vpairs);
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // data.frame of initial sites
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP sites_ = PROTECT(
    create_named_list(
      2,
      "x", x_,
      "y", y_
    )
  ); nprotect++;
  set_df_attributes(sites_);
  
  SEXP ntris_ = PROTECT(ScalarInteger(ctx.ntris)); nprotect++;
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // List of delaunay result
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SEXP res_ = PROTECT(
    create_named_list(
      5, 
      "sites"   , sites_,
      "ntris"   , ntris_,
      "tris"    , idx_, 
      "polygons", polys_,
      "segments", segs_
    )
  ); nprotect++;
  
  setAttrib(res_, R_ClassSymbol, mkString("del"));
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Free all the 'myalloc()' allocations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < ctx.alloc_count; i++) {
    free(ctx.allocs[i]);
  }
  free(ctx.allocs);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Convert C 0-indexing to R 1-indexing
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
  convert_indexing_c_to_r(v1_);
  convert_indexing_c_to_r(v2_);
  convert_indexing_c_to_r(v3_);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tidy and return
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  UNPROTECT(nprotect);
  return res_;
}







