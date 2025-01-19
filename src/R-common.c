
#define R_NO_REMAP

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


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Helper function used with sort()
// sort ctx->sites on y, then x, coord 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int site_comparison(const void *v1, const void *v2) {
  
  struct Point *s1 = (struct Point *)v1;
  struct Point *s2 = (struct Point *)v2;
  
  if (s1->y < s2->y)
    return (-1);
  if (s1->y > s2->y)
    return (1);
  if (s1->x < s2->x)
    return (-1);
  if (s1->x > s2->x)
    return (1);
  return (0);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Convert coordinates to ctx->sites, sort, and compute xmin, xmax, ymin, ymax 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void init_sites(context_t *ctx, double *x, double *y, int n) {
  
  ctx->nsites = n;
  ctx->sites = (struct Site *)myalloc(ctx, (unsigned int)ctx->nsites * sizeof *ctx->sites);
  for (int i = 0; i < n; i++) {
    ctx->sites[i].coord.x = x[i];
    ctx->sites[i].coord.y = y[i];
    ctx->sites[i].sitenbr = i;
    ctx->sites[i].refcnt  = 0;
  }
  
  qsort(ctx->sites, (unsigned long)ctx->nsites, sizeof *ctx->sites, site_comparison);
  ctx->xmin = ctx->sites[0].coord.x;
  ctx->xmax = ctx->sites[0].coord.x;
  
  for (int i = 1; i < ctx->nsites; i++) {
    if (ctx->sites[i].coord.x < ctx->xmin)
      ctx->xmin = ctx->sites[i].coord.x;
    if (ctx->sites[i].coord.x > ctx->xmax)
      ctx->xmax = ctx->sites[i].coord.x;
    
    if (ctx->sites[i].coord.x == ctx->sites[i - 1].coord.x &&
        ctx->sites[i].coord.y == ctx->sites[i - 1].coord.y) {
      free_all_myalloc(ctx);
      Rf_error("Input points contain duplicates. Not allowed [mem: %i]", ctx->total_alloc);
    }
    
  };
  
  ctx->ymin = ctx->sites[0].coord.y;
  ctx->ymax = ctx->sites[ctx->nsites - 1].coord.y;
}

