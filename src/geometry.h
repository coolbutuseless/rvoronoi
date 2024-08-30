
void geominit(context_t *ctx);
struct Edge *bisect(context_t *ctx, struct Site *s1, struct Site *s2);
struct Site *intersect(context_t *ctx, struct Halfedge *el1, struct Halfedge *el2);
int right_of(struct Halfedge *el, struct Point *p);
void endpoint(context_t *ctx, struct Edge *e, int lr, struct Site *s);
double dist(struct Site *s, struct Site *t);
void makevertex(context_t *ctx, struct Site *v);
void deref(context_t *ctx, struct Site *v);
void ref(struct Site *v);
