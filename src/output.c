
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
    ctx->line_a[ctx->line_idx] = e->a;
    ctx->line_b[ctx->line_idx] = e->b;
    ctx->line_c[ctx->line_idx] = e->c;
    ctx->line_idx++;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_ep(context_t *ctx, struct Edge *e) {
  if (!ctx->triangulate) {
    // Rprintf("e %d", e->edgenbr);
    // Rprintf(" %d ", e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1);
    // Rprintf("%d\n", e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1);
    
    ctx->seg_line[ctx->seg_idx] = e->edgenbr + 1; // convert to R 1-indexing
    ctx->seg_v1  [ctx->seg_idx] = e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr + 1 : NA_INTEGER;
    ctx->seg_v2  [ctx->seg_idx] = e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr + 1 : NA_INTEGER;
    
    ctx->seg_idx++;
    
  };
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void out_vertex(context_t *ctx, struct Site *v) {
  if (!ctx->triangulate) {
    ctx->vert_x[ctx->vert_idx] = v->coord.x;
    ctx->vert_y[ctx->vert_idx] = v->coord.y;
    ctx->vert_idx++;
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
    ctx->v1[ctx->idel] = s1->sitenbr + 1;
    ctx->v2[ctx->idel] = s2->sitenbr + 1;
    ctx->v3[ctx->idel] = s3->sitenbr + 1;
    ctx->idel++;
    
    // Rprintf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
  }
}


