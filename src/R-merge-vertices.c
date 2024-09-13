
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
#include "geometry.h"
#include "voronoi.h"

#include "utils.h"
#include "R-merge-vertices.h"



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Insert a new index remapping into the 'remap' data structure
//
// Given 'n' vertices, some of these may be close enough to merge into a 
// single vertex.  The 'segment' data returned vrom Fortune's voronoi references
// two vertices - thus making a line segment.
// If one of these vertices is deleted (as it is deemed equivalant to another
// vertex), then any segment which references this vertex must be renumbered
// to instead reference the other vertex.
// 
// This remapping has some corner cases.
//   * e.g. if vertex '4' remaps to vertex '2', and also vertex '4' should
//     remap to vertex '3', then a third remapping is implied i.e. 
//     vertex '3' should also remap to vertex '2'
//   * e.g. if vertex '4' remaps to vertex '3', but vertex '3' also remaps to
//     vertex '2', then this implies that vertex '4' should ultimately remap
//     to vertex '2'
//
// The 'remap' structure is just a vector of 'n' target vertices - one for
// each index.  It indicates the target index this vertex should be remapped 
// to (and no remapping is done if the target index is -1)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define NOT_REMAPPED -1
void remap_insert(int *remap, int i1, int i2) {
  // Swap so i1 is always larger
  if (i1 < i2) {
    int tmp = i1;
    i1 = i2;
    i2 = tmp;
  }
  
  // New remapping
  if (remap[i1] == NOT_REMAPPED) {
    // vertex is not currently remapped to anything, so we can remap it
    // unconditionally to the given value
    remap[i1] = i2;
    return;
  } else if (remap[i1] == i2) {
    // Vertex is already remapped to this other vertex. Do nothing.
    return;
  } else {
    // Vertex index 'i1' has already been mapped to something else.
    // Figure out which of the new vertex remapping or the current remapping
    // is smaller.  Remap 'i1' to the smaller vertex and then 
    // cascade the other remapping
    if (remap[i1] > i2) {
      remap_insert(remap, remap[i1], i2);
      remap[i1] = i2;
    } else {
      remap_insert(remap, remap[i1], i2);
    }
  }
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Merge voronoi vertices where:  (dist between vertices)^2 < tol
//
// @param tol tolerance
// @param nverts number of voronoi vertices
// @param x,y voronoi vertices
// @param nsegs number of voronoi segments
// @param line the index of the voronoi line on which this segment lies
// @param v1,v2 the indices of the vertices defining the extents of the segment
//        on the given line
// @param nsegs_final the final number of edges remaining after the merge 
//
// NOTE: 'nsegs_final', 'line', 'v1' and 'v2' are MODIFIED IN PLACE
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void merge_vertices_core_(double tol, 
                          int nverts, double *x, double *y, 
                          int nsegs, int *line, int *v1, int *v2, 
                          int *nsegs_final,
                          int verbosity) {
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialise the data structure for remapping a vertex onto another 
  // vertex.  As vertices which are very close are found, the higher numbered
  // vertex is remapped onto the lower number vertex.
  // See 'remap_insert()' for more information
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int *remap = malloc((unsigned long)nverts * sizeof(int));
  if (remap == NULL) {
    error("merge_vertices_core_(): could not allocate 'remap'");
  }
  for (int i = 0; i < nverts; i++) {
    remap[i] = NOT_REMAPPED;
  }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Track edges to be discarded.
  // disard[j] = true when:
  //    1. segment has zero-length i.e. same vertex on each end
  //    2. segment is a repeated version of an infinite segment
  //       Fortune's voronoi outputs segments e.g. in the
  //       case of just 2 points on the plane.  The algorithm however outputs
  //       *TWO* copies of the infinite segments and one of them
  //       MUST be deleted for the polygon extraction to work.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int ndiscard = 0;
  bool *discard = calloc((unsigned int)nsegs, sizeof(bool));
  if (discard == NULL) error("merge_vertices_core_(): could not allocate 'discard'");
  
  
  
  if (tol > 0) {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // For the data returned by Fortune's Voronoi, the only vertices which need
    // to be compared are the two ends of each segment.  I.e. it is not necessary
    // to do a distance check between all pairs of points
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (int i = 0; i < nsegs; i++) {
      
      // The vertex indices of the two ends of this segment
      int i1 = v1[i];
      int i2 = v2[i];
      
      // If this is an infinite 'ray', then the distance between endpoints is INF
      if (i1 < 0 || i2 < 0) continue;
      
      // Sanity check
      if ( (i1 >= nverts) || (i2 >= nverts) || (i1 < 0) || (i2 < 0) ) {
        error("merge_vertices_core_(): [%i] Vertex index out of bounds: %i, %i [0, %i]", 
              i, i1, i2, nverts);
      }
      
      // Distance between two ends of the segment
      double x1 = x[i1];
      double x2 = x[i2];
      double y1 = y[i1];
      double y2 = y[i2];
      
      double dist = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
      
      if (verbosity > 0) {
        Rprintf("Edge: %i (%i -> %i): Dist: %.6f\n", i, i1, i2, dist);
      }
      
      if (dist < tol) {
        if (verbosity > 0) Rprintf("^ Remap\n");
        remap_insert(remap, i1, i2);
      } 
    }
    
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Collapse and remove any edges with duplicated vertices
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (verbosity > 0) {
      for (int i = 0; i < nsegs; i++) {
        Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
      }
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Print remaps
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (verbosity > 0) {
      for (int i = 0; i < nverts; i++) {
        if (remap[i] >= 0) {
          Rprintf("Remap: %i => %i\n", i, remap[i]);
        }
      }
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Renumber verts from highest-to-lowest
    //
    // For each vertex:
    //    has it been remapped?
    //    If yes:
    //       for each segment, 
    //          if either end of this segment matches this vertex
    //              remap it to the target vertex
    //          if both ends of this segment are the same, then delete this segment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (int i = nverts - 1; i >= 0; i--) {
      if (remap[i] == NOT_REMAPPED) continue;
      int vlo = remap[i];
      int vhi = i;
      
      for (int j = 0; j < nsegs; j++) {
        if (v1[j] == vhi) v1[j] = vlo; 
        if (v2[j] == vhi) v2[j] = vlo; 
        
        discard[j] = v1[j] == v2[j];
      }
    }
    
  } // End: if (tol > 0)
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Mark duplicate 'infinite segments' for removal
  //  Each 'infinite segment' ends up appearing *twice* in the segment list
  //  so we can delete.
  //  An example of repeated infinite segments
  //     line     v1      v2
  //  --------------------------
  //        1      1       2
  //        2    -999    -999
  //        3      4       5
  //        2    -999    -999      REPEATED
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bool line_used[nsegs];
  for (int i = 0; i < nsegs; i++) line_used[i] = false;
  for (int i = 0; i < nsegs; i++) {
    if (v1[i] >= 0 || v2[i] >= 0) continue;
    if (line_used[line[i]]) {
      discard[i] = true;
    }
    line_used[line[i]] = true;
  }
  
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tally the discarded segments
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int i = 0; i < nsegs; i++) {
    ndiscard += discard[i];
  }
  
  if (verbosity > 0) {
    Rprintf("Ndistcard: %i\n", ndiscard);
    for (int i = 0; i < nsegs; i++) {
      Rprintf("(%4i)  %4i\n", i, discard[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Remove segments marked for 'discard' by compacting the 
  // vectors of v1, v2 and line
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (verbosity > 0) {  
    for (int i = 0; i < nsegs; i++) {
      Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
    }
  }
  
  
  int src = 0;
  int dst = 0;
  for (dst = 0; dst < nsegs - ndiscard; dst++, src++) {
    while(discard[src]) {
      src++;
    }
    // src += discard[src];
    v1[dst]   = v1[src];
    v2[dst]   = v2[src];
    line[dst] = line[src];
    if (verbosity > 1) {
      Rprintf("src %4i  ->  dst %4i\n", src, dst);
    }
  }
  
  // Mark out-of-bounds segments with a sentinel value so it's easier to see
  // when something's not right.
  for (; dst < nsegs; dst++) {
    v1[dst]   = 998;
    v2[dst]   = 998;
    line[dst] = 998;
  }
  
  if (verbosity > 0) {  
    for (int i = 0; i < nsegs; i++) {
      Rprintf("%2i: (%4i, %4i)\n", i, v1[i], v2[i]);
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Set the number of segments which are left are merging
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  *nsegs_final = nsegs - ndiscard;
  
  free(discard);
  free(remap);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R shim used for testing vertex merging
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEXP merge_vertices_(SEXP x_, SEXP y_, 
                     SEXP line_, SEXP v1_, SEXP v2_, 
                     SEXP tol_, SEXP verbosity_) {
  
  int nprotect = 0;
  
  if (length(x_) == 0 || length(x_) != length(y_)) {
    error("merge_vertices_() x & y must be equal length");
  }
  if (length(v1_) == 0) {
    warning("merge_vertices_(): No vertices to process\n");
    return R_NilValue;
  }
  if (length(v1_) != length(v2_) || length(v1_) != length(line_)) {
    error("merge_vertices_():line, v1 & v2 must be equal length");
  }
  
  double *x = REAL(x_);
  double *y = REAL(y_);
  
  int nverts = length(x_);
  int nsegs = length(v1_);
  
  SEXP fv1_   = PROTECT(allocVector(INTSXP, nsegs)); nprotect++;
  SEXP fv2_   = PROTECT(allocVector(INTSXP, nsegs)); nprotect++;
  SEXP fline_ = PROTECT(allocVector(INTSXP, nsegs)); nprotect++;
  
  // int *fv1   = INTEGER(fv1_);
  // int *fv2   = INTEGER(fv2_);
  // int *fline = INTEGER(fline_);
  
  // These are initially copies of the input vectors 
  memcpy(INTEGER(fv1_)  , INTEGER(v1_)  , nsegs * sizeof(int));
  memcpy(INTEGER(fv2_)  , INTEGER(v2_)  , nsegs * sizeof(int));
  memcpy(INTEGER(fline_), INTEGER(line_), nsegs * sizeof(int));
  
  // Convert them to C 0-indexing for processing
  convert_indexing_r_to_c(fv1_);
  convert_indexing_r_to_c(fv2_);
  convert_indexing_r_to_c(fline_);
  
  // for (int i = 0; i < length(fv2_); i++) {
  //   Rprintf("(%i) v2[%i] = %i    =   fv2[%i] = %i\n", nsegs, i, INTEGER(v2_)[i], i, INTEGER(fv2_)[i]);
  // }
  
  // This will contain the final number of vertices after merging  
  int nsegs_final = 0;
  
  merge_vertices_core_(asReal(tol_), 
                       nverts, x, y, 
                       nsegs, INTEGER(fline_), INTEGER(fv1_), INTEGER(fv2_), 
                       &nsegs_final, asInteger(verbosity_));
  
  trim_vec(fv1_  , nsegs_final, nsegs);
  trim_vec(fv2_  , nsegs_final, nsegs);
  trim_vec(fline_, nsegs_final, nsegs);
  
  convert_indexing_c_to_r(fv1_);
  convert_indexing_c_to_r(fv2_);
  convert_indexing_c_to_r(fline_);
  
  
  SEXP res_ = PROTECT(
    create_named_list(3, "line", fline_, "v1", fv1_, "v2", fv2_)
  ); nprotect++;
  set_df_attributes(res_, nsegs_final, nsegs_final);
  
  
  UNPROTECT(nprotect);
  return res_;
}



