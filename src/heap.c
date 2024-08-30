
#include "defs.h"
#include "heap.h"
#include "geometry.h"
#include "memory.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PQinsert(context_t *ctx, struct Halfedge *he, struct Site *v, double offset) {
  struct Halfedge *last, *next;

  he->vertex = v;
  ref(v);
  he->ystar = v->coord.y + offset;
  last = &ctx->PQhash[PQbucket(ctx, he)];
  while ((next = last->PQnext) != NULL &&
         (he->ystar > next->ystar ||
          (he->ystar == next->ystar && v->coord.x > next->vertex->coord.x))) {
    last = next;
  };
  he->PQnext = last->PQnext;
  last->PQnext = he;
  ctx->PQcount++;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PQdelete(context_t *ctx, struct Halfedge *he) {
  struct Halfedge *last;

  if (he->vertex != NULL) {
    last = &ctx->PQhash[PQbucket(ctx, he)];
    while (last->PQnext != he)
      last = last->PQnext;
    last->PQnext = he->PQnext;
    ctx->PQcount--;
    deref(ctx, he->vertex);
    he->vertex = NULL;
  };
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int PQbucket(context_t *ctx, struct Halfedge *he) {
  int bucket;

  bucket = (he->ystar - ctx->ymin) / ctx->deltay * ctx->PQhashsize;
  if (bucket < 0)
    bucket = 0;
  if (bucket >= ctx->PQhashsize)
    bucket = ctx->PQhashsize - 1;
  if (bucket < ctx->PQmin)
    ctx->PQmin = bucket;
  return (bucket);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int PQempty(context_t *ctx) { 
  return (ctx->PQcount == 0); 
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Point PQ_min(context_t *ctx) {
  struct Point answer;

  while (ctx->PQhash[ctx->PQmin].PQnext == NULL) {
    ctx->PQmin++;
  };
  answer.x = ctx->PQhash[ctx->PQmin].PQnext->vertex->coord.x;
  answer.y = ctx->PQhash[ctx->PQmin].PQnext->ystar;
  return (answer);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Halfedge *PQextractmin(context_t *ctx) {
  struct Halfedge *curr;

  curr = ctx->PQhash[ctx->PQmin].PQnext;
  ctx->PQhash[ctx->PQmin].PQnext = curr->PQnext;
  ctx->PQcount--;
  return (curr);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PQinitialize(context_t *ctx) {

  ctx->PQcount = 0;
  ctx->PQmin = 0;
  ctx->PQhashsize = 4 * ctx->sqrt_nsites;
  ctx->PQhash = (struct Halfedge *)myalloc(ctx, ctx->PQhashsize * sizeof *ctx->PQhash);
  for (int i = 0; i < ctx->PQhashsize; i++)
    ctx->PQhash[i].PQnext = NULL;
}


