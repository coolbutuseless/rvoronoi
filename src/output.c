
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "defs.h"
#include "output.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_bisector(context_t *ctx, struct Edge *e) {
  if (!ctx->triangulate) {
    // Rprintf("l %f %f %f\n", e->a, e->b, e->c);
    ctx->line_a[ctx->nlines] = e->a;
    ctx->line_b[ctx->nlines] = e->b;
    ctx->line_c[ctx->nlines] = e->c;
    ctx->nlines++;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_ep(context_t *ctx, struct Edge *e) {
  if (!ctx->triangulate) {
    // Rprintf("e %d", e->edgenbr);
    // Rprintf(" %d ", e->ep[le] != NULL ? e->ep[le]->sitenbr : -1);
    // Rprintf("%d\n", e->ep[re] != NULL ? e->ep[re]->sitenbr : -1);
    
    ctx->seg_line[ctx->nsegs] = e->edgenbr; // convert to R 1-indexing
    ctx->seg_v1  [ctx->nsegs] = e->ep[le] != NULL ? e->ep[le]->sitenbr : -999;
    ctx->seg_v2  [ctx->nsegs] = e->ep[re] != NULL ? e->ep[re]->sitenbr : -999;
    
    ctx->nsegs++;
    
  };
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_vertex(context_t *ctx, struct Site *v) {
  if (!ctx->triangulate) {
    ctx->vert_x[ctx->nverts] = v->coord.x;
    ctx->vert_y[ctx->nverts] = v->coord.y;
    ctx->nverts++;
    // Rprintf("v %f %f\n", v->coord.x, v->coord.y);
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_site(context_t *ctx, struct Site *s) {
  if (!ctx->triangulate) {
    // Rprintf("s %f %f\n", s->coord.x, s->coord.y);
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_triple(context_t *ctx, struct Site *s1, struct Site *s2, struct Site *s3) {
  if (ctx->triangulate) {
    ctx->v1[ctx->ntris] = s1->sitenbr;
    ctx->v2[ctx->ntris] = s2->sitenbr;
    ctx->v3[ctx->ntris] = s3->sitenbr;
    ctx->ntris++;
    
    // Rprintf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
  }
}


