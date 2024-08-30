
#include "defs.h"
#include "edgelist.h"
#include "memory.h"
#include "geometry.h"
#include "output.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ELinitialize(context_t *ctx) {
  freeinit(&ctx->hfl, sizeof **ctx->ELhash);
  ctx->ELhashsize = 2 * ctx->sqrt_nsites;
  ctx->ELhash = (struct Halfedge **)myalloc(ctx, sizeof *ctx->ELhash * ctx->ELhashsize);
  for (int i = 0; i < ctx->ELhashsize; i++)
    ctx->ELhash[i] = NULL;
  ctx->ELleftend = HEcreate(ctx, NULL, 0);
  ctx->ELrightend = HEcreate(ctx, NULL, 0);
  ctx->ELleftend->ELleft = NULL;
  ctx->ELleftend->ELright = ctx->ELrightend;
  ctx->ELrightend->ELleft = ctx->ELleftend;
  ctx->ELrightend->ELright = NULL;
  ctx->ELhash[0] = ctx->ELleftend;
  ctx->ELhash[ctx->ELhashsize - 1] = ctx->ELrightend;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Halfedge *HEcreate(context_t *ctx, struct Edge *e, int pm) {
  struct Halfedge *answer;
  answer = (struct Halfedge *)getfree(ctx, &ctx->hfl);
  answer->ELedge = e;
  answer->ELpm = pm;
  answer->PQnext = NULL;
  answer->vertex = NULL;
  answer->ELrefcnt = 0;
  return (answer);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ELinsert(struct Halfedge *lb, struct Halfedge *new) {
  new->ELleft = lb;
  new->ELright = lb->ELright;
  (lb->ELright)->ELleft = new;
  lb->ELright = new;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Get entry from hash table, pruning any deleted nodes 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Halfedge *ELgethash(context_t *ctx, int b) {
  struct Halfedge *he;

  if (b < 0 || b >= ctx->ELhashsize)
    return (NULL);
  he = ctx->ELhash[b];
  if (he == NULL || he->ELedge != (struct Edge *)DELETED)
    return (he);

  /* Hash table points to deleted half edge.  Patch as necessary. */
  ctx->ELhash[b] = NULL;
  if ((he->ELrefcnt--) == 0)
    makefree((struct Freenode *)he, (struct Freelist *)&ctx->hfl);
  return (NULL);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Halfedge *ELleftbnd(context_t *ctx, struct Point *p) {
  struct Halfedge *he;

  /* Use hash table to get close to desired halfedge */
  int bucket = (p->x - ctx->xmin) / ctx->deltax * ctx->ELhashsize;
  if (bucket < 0)
    bucket = 0;
  if (bucket >= ctx->ELhashsize)
    bucket = ctx->ELhashsize - 1;
  he = ELgethash(ctx, bucket);
  if (he == NULL) {
    for (int i = 1; 1; i++) {
      if ( ((he = ELgethash(ctx, bucket - i)) != NULL) || ((he = ELgethash(ctx, bucket + i)) != NULL) )
        break;
    };
  };
  /* Now search linear list of halfedges for the corect one */
  if (he == ctx->ELleftend || (he != ctx->ELrightend && right_of(he, p))) {
    do {
      he = he->ELright;
    } while (he != ctx->ELrightend && right_of(he, p));
    he = he->ELleft;
  } else
    do {
      he = he->ELleft;
    } while (he != ctx->ELleftend && !right_of(he, p));

  /* Update hash table and reference counts */
  if (bucket > 0 && bucket < ctx->ELhashsize - 1) {
    if (ctx->ELhash[bucket] != NULL)
      ctx->ELhash[bucket]->ELrefcnt--;
    ctx->ELhash[bucket] = he;
    ctx->ELhash[bucket]->ELrefcnt++;
  };
  return (he);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// This delete routine can't reclaim node, since pointers from hash
// table may be present.   
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ELdelete(struct Halfedge *he) {
  (he->ELleft)->ELright = he->ELright;
  (he->ELright)->ELleft = he->ELleft;
  he->ELedge = (struct Edge *)DELETED;
}


struct Halfedge *ELright(struct Halfedge *he) {
  return (he->ELright);
}


struct Halfedge *ELleft(struct Halfedge *he) {
  return (he->ELleft);
}


struct Site *leftreg(context_t *ctx, struct Halfedge *he) {
  if (he->ELedge == NULL)
    return (ctx->bottomsite);
  return (he->ELpm == le ? he->ELedge->reg[le] : he->ELedge->reg[re]);
}


struct Site *rightreg(context_t *ctx, struct Halfedge *he) {
  if (he->ELedge == NULL)
    return (ctx->bottomsite);
  return (he->ELpm == le ? he->ELedge->reg[re] : he->ELedge->reg[le]);
}
