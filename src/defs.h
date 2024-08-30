#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

#define le 0
#define re 1


struct Freenode {
  struct Freenode *nextfree;
};


struct Freelist {
  struct Freenode *head;
  int nodesize;
};


struct Point {
  double x, y;
};


/* structure used both for sites and for vertices */
struct Site {
  struct Point coord;
  int sitenbr;
  int refcnt;
};


struct Edge {
  double a, b, c;
  struct Site *ep[2];
  struct Site *reg[2];
  int edgenbr;
};


struct Halfedge {
  struct Halfedge *ELleft, *ELright;
  struct Edge *ELedge;
  int ELrefcnt;
  char ELpm;
  struct Site *vertex;
  double ystar;
  struct Halfedge *PQnext;
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// All global context is now here
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  
  int triangulate;
  
  int total_alloc;
  
  double xmin, xmax, ymin, ymax, deltax, deltay;
  
  struct Site *sites;
  
  struct Freelist sfl;
  struct Site *bottomsite;
  
  int nsites;
  int siteidx;
  int sqrt_nsites;
  int nvertices;
  
  int nedges;
  struct Freelist efl;
  
  struct Freelist hfl;
  struct Halfedge *ELleftend, *ELrightend;
  int ELhashsize;
  struct Halfedge **ELhash;
  
  
  struct Halfedge *PQhash;
  int PQhashsize;
  int PQcount;
  int PQmin;
  
  
  // Delauncy indices to return to R
  int *v1;
  int *v2;
  int *v3;
  int idel;
  
  // Voronoi Information to return to R
  int line_idx;
  int vert_idx;
  int seg_idx;
  double *line_a, *line_b, *line_c;
  double *vert_x, *vert_y;
  int *seg_line, *seg_v1, *seg_v2;
  
  // Half edge allocations
  void **allocs;
  int alloc_count;
  int alloc_capacity;
  
} context_t;




