
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include "defs.h"
#include "voronoi.h"
#include "heap.h"
#include "output.h"
#include "edgelist.h"
#include "geometry.h"

/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* return a single in-storage site */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Site *nextsite(context_t *ctx) {
  struct Site *s;
  if (ctx->siteidx < ctx->nsites) {
    s = &ctx->sites[ctx->siteidx];
    ctx->siteidx++;
    return (s);
  } else
    return NULL;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
// deltax, deltay (can all be estimates).
// Performance suffers if they are wrong; better to make nsites,
// deltax, and deltay too big than too small.  (?) */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void voronoi(context_t *ctx, int triangulate) {
  struct Site *newsite, *bot, *top, *temp, *p;
  struct Site *v;
  struct Point newintstar = { 0 };
  int pm;
  struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
  struct Edge *e;

  PQinitialize(ctx);
  ctx->bottomsite = nextsite(ctx);
  out_site(ctx, ctx->bottomsite);
  ELinitialize(ctx);

  newsite = nextsite(ctx);
  while (1) {
    if (!PQempty(ctx))
      newintstar = PQ_min(ctx);

    if (newsite != NULL &&
        (PQempty(ctx) || newsite->coord.y < newintstar.y ||
         (newsite->coord.y == newintstar.y &&
          newsite->coord.x < newintstar.x))) { /* new site is smallest */
      out_site(ctx, newsite);
      lbnd = ELleftbnd(ctx, &(newsite->coord));
      rbnd = ELright(lbnd);
      bot = rightreg(ctx, lbnd);
      e = bisect(ctx, bot, newsite);
      bisector = HEcreate(ctx, e, le);
      ELinsert(lbnd, bisector);
      if ((p = intersect(ctx, lbnd, bisector)) != NULL) {
        PQdelete(ctx, lbnd);
        PQinsert(ctx, lbnd, p, dist(p, newsite));
      };
      lbnd = bisector;
      bisector = HEcreate(ctx, e, re);
      ELinsert(lbnd, bisector);
      if ((p = intersect(ctx, bisector, rbnd)) != NULL) {
        PQinsert(ctx, bisector, p, dist(p, newsite));
      };
      newsite = nextsite(ctx);
    } else if (!PQempty(ctx))
    /* intersection is smallest */
    {
      lbnd = PQextractmin(ctx);
      llbnd = ELleft(lbnd);
      rbnd = ELright(lbnd);
      rrbnd = ELright(rbnd);
      bot = leftreg(ctx, lbnd);
      top = rightreg(ctx, rbnd);
      out_triple(ctx, bot, top, rightreg(ctx, lbnd));
      v = lbnd->vertex;
      makevertex(ctx, v);
      endpoint(ctx, lbnd->ELedge, lbnd->ELpm, v);
      endpoint(ctx, rbnd->ELedge, rbnd->ELpm, v);
      ELdelete(lbnd);
      PQdelete(ctx, rbnd);
      ELdelete(rbnd);
      pm = le;
      if (bot->coord.y > top->coord.y) {
        temp = bot;
        bot = top;
        top = temp;
        pm = re;
      }
      e = bisect(ctx, bot, top);
      bisector = HEcreate(ctx, e, pm);
      ELinsert(llbnd, bisector);
      endpoint(ctx, e, re - pm, v);
      deref(ctx, v);
      if ((p = intersect(ctx, llbnd, bisector)) != NULL) {
        PQdelete(ctx, llbnd);
        PQinsert(ctx, llbnd, p, dist(p, bot));
      };
      if ((p = intersect(ctx, bisector, rrbnd)) != NULL) {
        PQinsert(ctx, bisector, p, dist(p, bot));
      };
    } else
      break;
  };

  for (lbnd = ELright(ctx->ELleftend); lbnd != ctx->ELrightend; lbnd = ELright(lbnd)) {
    e = lbnd->ELedge;
    out_ep(ctx, e);
  };
}
