
#include <math.h>

#include "defs.h"
#include "geometry.h"
#include "memory.h"
#include "output.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void geominit(context_t *ctx) {
  struct Edge e;
  double sn;

  freeinit(&ctx->efl, sizeof e);
  ctx->nvertices = 0;
  ctx->nedges = 0;
  sn = ctx->nsites + 4;
  ctx->sqrt_nsites = (int)ceil(sqrt((double)sn));
  ctx->deltay = ctx->ymax - ctx->ymin;
  ctx->deltax = ctx->xmax - ctx->xmin;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Edge *bisect(context_t *ctx, struct Site *s1, struct Site *s2) {
  double dx, dy, adx, ady;
  struct Edge *newedge;

  newedge = (struct Edge *)getfree(ctx, &ctx->efl);

  newedge->reg[0] = s1;
  newedge->reg[1] = s2;
  ref(s1);
  ref(s2);
  newedge->ep[0] = NULL;
  newedge->ep[1] = NULL;

  dx = s2->coord.x - s1->coord.x;
  dy = s2->coord.y - s1->coord.y;
  adx = dx > 0 ? dx : -dx;
  ady = dy > 0 ? dy : -dy;
  newedge->c = s1->coord.x * dx + s1->coord.y * dy + (dx * dx + dy * dy) * 0.5;
  if (adx > ady) {
    newedge->a = 1.0;
    newedge->b = dy / dx;
    newedge->c /= dx;
  } else {
    newedge->b = 1.0;
    newedge->a = dx / dy;
    newedge->c /= dy;
  };

  newedge->edgenbr = ctx->nedges;
  out_bisector(ctx, newedge);
  ctx->nedges++;
  return (newedge);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Site *intersect(context_t *ctx, struct Halfedge *el1, struct Halfedge *el2) {
  struct Edge *e1, *e2, *e;
  struct Halfedge *el;
  double d, xint, yint;
  int right_of_site;
  struct Site *v;

  e1 = el1->ELedge;
  e2 = el2->ELedge;
  if (e1 == NULL || e2 == NULL)
    return (NULL);
  if (e1->reg[1] == e2->reg[1])
    return (NULL);

  d = e1->a * e2->b - e1->b * e2->a;
  if (-1.0e-10 < d && d < 1.0e-10)
    return (NULL);

  xint = (e1->c * e2->b - e2->c * e1->b) / d;
  yint = (e2->c * e1->a - e1->c * e2->a) / d;

  if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
      (e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
       e1->reg[1]->coord.x < e2->reg[1]->coord.x)) {
    el = el1;
    e = e1;
  } else {
    el = el2;
    e = e2;
  };
  right_of_site = xint >= e->reg[1]->coord.x;
  if ((right_of_site && el->ELpm == le) || (!right_of_site && el->ELpm == re))
    return (NULL);

  v = (struct Site *)getfree(ctx, &ctx->sfl);
  v->refcnt = 0;
  v->coord.x = xint;
  v->coord.y = yint;
  return (v);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// returns 1 if p is to right of halfedge e 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int right_of(struct Halfedge *el, struct Point *p) {
  struct Edge *e;
  struct Site *topsite;
  int right_of_site, above, fast;
  double dxp, dyp, dxs, t1, t2, t3, yl;

  e = el->ELedge;
  topsite = e->reg[1];
  right_of_site = p->x > topsite->coord.x;
  if (right_of_site && el->ELpm == le)
    return (1);
  if (!right_of_site && el->ELpm == re)
    return (0);

  if (e->a == 1.0) {
    dyp = p->y - topsite->coord.y;
    dxp = p->x - topsite->coord.x;
    fast = 0;
    if ((!right_of_site & (e->b < 0.0)) || (right_of_site & (e->b >= 0.0))) {
      above = dyp >= e->b * dxp;
      fast = above;
    } else {
      above = p->x + p->y * e->b > e->c;
      if (e->b < 0.0)
        above = !above;
      if (!above)
        fast = 1;
    };
    if (!fast) {
      dxs = topsite->coord.x - (e->reg[0])->coord.x;
      above = e->b * (dxp * dxp - dyp * dyp) <
              dxs * dyp * (1.0 + 2.0 * dxp / dxs + e->b * e->b);
      if (e->b < 0.0)
        above = !above;
    };
  } else /*e->b==1.0 */
  {
    yl = e->c - e->a * p->x;
    t1 = p->y - yl;
    t2 = p->x - topsite->coord.x;
    t3 = yl - topsite->coord.y;
    above = t1 * t1 > t2 * t2 + t3 * t3;
  };
  return (el->ELpm == le ? above : !above);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void endpoint(context_t *ctx, struct Edge *e, int lr, struct Site *s) {
  e->ep[lr] = s;
  ref(s);
  if (e->ep[re - lr] == NULL)
    return;
  out_ep(ctx, e);
  deref(ctx, e->reg[le]);
  deref(ctx, e->reg[re]);
  makefree((struct Freenode *)e, (struct Freelist *)&ctx->efl);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double dist(struct Site *s, struct Site *t) {
  double dx, dy;
  dx = s->coord.x - t->coord.x;
  dy = s->coord.y - t->coord.y;
  return (sqrt(dx * dx + dy * dy));
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void makevertex(context_t *ctx, struct Site *v) {
  v->sitenbr = ctx->nvertices;
  ctx->nvertices++;
  out_vertex(ctx, v);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void deref(context_t *ctx, struct Site *v) {
  v->refcnt--;
  if (v->refcnt == 0)
    makefree((struct Freenode *)v, (struct Freelist *)&ctx->sfl);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ref(struct Site *v) { 
  v->refcnt++; 
}


