void freeinit(struct Freelist *fl, int size);
char *getfree(context_t *ctx, struct Freelist *fl);
void makefree(struct Freenode *curr, struct Freelist *fl);
char *myalloc(context_t *ctx, unsigned n);
