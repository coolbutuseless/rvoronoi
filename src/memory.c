#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "defs.h"
#include "memory.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void freeinit(struct Freelist *fl, int size) {
  fl->head = NULL;
  fl->nodesize = size;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
char *getfree(context_t *ctx, struct Freelist *fl) {
  struct Freenode *t;
  if (fl->head == NULL) {
    t = (struct Freenode *)myalloc(ctx, ctx->sqrt_nsites * fl->nodesize);
    for (int i = 0; i < ctx->sqrt_nsites; i++)
      makefree((struct Freenode *)((char *)t + i * fl->nodesize), fl);
  };
  t = fl->head;
  fl->head = (fl->head)->nextfree;
  return ((char *)t);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void makefree(struct Freenode *curr, struct Freelist *fl) {
  curr->nextfree = fl->head;
  fl->head = curr;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
char *myalloc(context_t *ctx, unsigned n) {
  char *t;
  if ((t = malloc(n)) == (char *)0) {
    error("myalloc(): Insufficient memory processing site %d (%d bytes in use)\n",
            ctx->siteidx, ctx->total_alloc);
  };
  
  ctx->total_alloc += n;
  
  // keep track of memory allocations to free at end
  if (ctx->alloc_count == ctx->alloc_capacity) {
    ctx->alloc_capacity *= 2;
    ctx->allocs = realloc(ctx->allocs, ctx->alloc_capacity * sizeof(void *));
  }
  
  ctx->allocs[ctx->alloc_count] = (void *)t;
  ctx->alloc_count++;
  
  
  return (t);
}
