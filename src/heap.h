
void PQinsert(context_t *ctx, struct Halfedge *he, struct Site *v, float offset);
void PQdelete(context_t *ctx, struct Halfedge *he) ;
int PQbucket(context_t *ctx, struct Halfedge *he);
int PQempty(context_t *ctx);
struct Point PQ_min(context_t *ctx);
struct Halfedge *PQextractmin(context_t *ctx);
void PQinitialize(context_t *ctx);
