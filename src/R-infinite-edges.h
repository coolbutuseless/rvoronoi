void calc_space_for_bound_infinite_edges(int nsegs, int *v1, int *v2, int *nbverts, int *nbsegs);
void bound_infinite_edges(
    bbox_t *bounds,
    int nverts, double *x, double *y, 
    int nsegs , int *li, int *v1, int *v2,
    int nlines, double *a, double *b, double *c,
    int *nbverts, double *xb, double *yb,
    int *nbsegs, int *rv1, int *rv2);

