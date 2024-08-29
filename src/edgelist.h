
void ELinitialize(context_t *ctx);
struct Halfedge *HEcreate(context_t *ctx, struct Edge *e, int pm);
void ELinsert(struct Halfedge *lb, struct Halfedge *new);
struct Halfedge *ELgethash(context_t *ctx, int b);
struct Halfedge *ELleftbnd(context_t *ctx, struct Point *p);
void ELdelete(struct Halfedge *he);
struct Halfedge *ELright(struct Halfedge *he);
struct Halfedge *ELleft(struct Halfedge *he);
struct Site *leftreg(context_t *ctx, struct Halfedge *he);
struct Site *rightreg(context_t *ctx, struct Halfedge *he);
