


# extern SEXP bound_infinite_edges_(
# SEXP xmin_, SEXP ymin_, SEXP xmax_, SEXP ymax_,
# SEXP x_, SEXP y_,
# SEXP v1_, SEXP v2_);


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Find boundary intersectinons of infinite edges
#' 
#' @param x,y voronoi vertices
#' @param line,v1,v2 voronoi segment definition: line, vertex1, vertex2
#' @param xmin,ymin,xmax,ymax bounding box to trim to.  Must be larger than 
#'        bounding box of (x, y) points
#' @param a,b,c line defintion ax + by = c
#' @return list of new vertices and new edges
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bound_infinite_edges <- function(xmin, ymin, xmax, ymax,
                                 x, y, 
                                 a, b, c,
                                 line, v1, v2) {
  
  
  .Call(bound_infinite_edges_, 
        xmin, ymin, xmax, ymax,
        x, y,
        a, b, c,
        line, v1, v2)
}



if (FALSE) {
  
  set.seed(1)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  # x <- c(0, 0)
  # y <- c(0, 1)
  
  
  vor <- voronoi(x, y)
  
  dw <- 0.2
  lim <- c(0 - dw, 1 + dw)
  
  # Plotting all the finite segments
  buf <- 0.1
  plot(x, y, col = 'black', pch = '+',
       asp = 1, ann = FALSE, axes = FALSE, 
       xlim = c(vor$extents$xmin - buf, vor$extents$xmax + buf), 
       ylim = c(vor$extents$ymin - buf, vor$extents$ymax + buf)
  )
  
  # Finite segments
  fseg <- vor$segment
  fseg <- subset(fseg, v1 > 0 & v2 > 0)
  
  segments(
    vor$vertices$x[fseg$v1], vor$vertices$y[fseg$v1],
    vor$vertices$x[fseg$v2], vor$vertices$y[fseg$v2]
  )
  
  # Bounds
  with(vor$extents, rect(xmin, ymin, xmax, ymax, border = 'grey90'))
  
  with(vor$extents, rect(xmin - buf, ymin - buf, xmax + buf, ymax + buf, border = 'grey40'))
  
  
  # Add in the lines representing infinite segments
  # Extract segments which are unbounded. i.e. v1 or v2 is NA
  inf_seg <- subset(vor$segment, v1 < 0 | v2 < 0)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Find the intersections with the rectangular boundary
  #
  #  CORNER CASES
  #     * vertical line
  #     * horizontal line
  #     * unbounded edge in 2 directions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  seg_idx <- 1
  for (seg_idx in seq_len(nrow(inf_seg))) {
    
    seg <- inf_seg[seg_idx,]
    ll  <- vor$line[seg$line,]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Voronoi point
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    type <- 0
    x0   <- NA_real_
    y0   <- NA_real_
    if (seg$v1 > 0) {
      type <- 2
      x0   <- vor$vertices$x[seg$v1]
      y0   <- vor$vertices$y[seg$v1]
    } else if (seg$v2 > 0) {
      type <- 1
      x0   <- vor$vertices$x[seg$v2]
      y0   <- vor$vertices$y[seg$v2]
    }
    points(x0, y0, pch = 19, col = 'blue')
    
    
    # There can only be 2 intercept points which lie on the boundary
    intercepts <- list(NULL, NULL)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # intersection with left edge
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pleft <- list(
      x = vor$extents$xmin,
      y = with(ll, (c - a * vor$extents$xmin)/b)
    )
    
    # Is it within the uuppler/lower bounds
    if (pleft$y >= vor$extents$ymin && pleft$y <= vor$extents$ymax) {
      intercepts[[1]] <- pleft # This intercept is on the left
    }
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # intersection with right edge
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pright <- list(
      x = vor$extents$xmax,
      y = with(ll, (c - a * vor$extents$xmax)/b)
    )
    
    # Is it within the upper/lower bounds
    if (pright$y >= vor$extents$ymin && pright$y <= vor$extents$ymax) {
      intercepts[[2]] <- pright
    }
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # intersection with lower edge
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plower <- list(
      y = vor$extents$ymin,
      x = with(ll, (c - b * vor$extents$ymin)/a)
    )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # intersection with upper edge
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pupper <- list(
      y = vor$extents$ymax,
      x = with(ll, (c - b * vor$extents$ymax)/a)
    )
    
    if (is.null(intercepts[[1]])) {
      if (plower$x < pupper$x) 
        intercepts[[1]] <- plower
      else 
        intercepts[[1]] <- pupper
    }
    
    
    if (is.null(intercepts[[2]])) {
      if (pupper$x > plower$x) 
        intercepts[[2]] <- pupper
      else 
        intercepts[[2]] <- plower
    }
    
    
    
    
    if (type == 0) {
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Double ended unbounded edge
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      segments(intercepts[[1]]$x, intercepts[[1]]$y, intercepts[[2]]$x, intercepts[[2]]$y, col = 'orange')
    } else if (type == 1) {
      segments(x0, y0, intercepts[[1]]$x, intercepts[[1]]$y, col = 'red')
    } else {
      segments(x0, y0, intercepts[[2]]$x, intercepts[[2]]$y, col = 'red')
    }
    
  }
  
  
  
  res <- bound_infinite_edges(
    vor$extents$xmin, vor$extents$ymin, vor$extents$xmax, vor$extents$ymax,
    vor$vertices$x, vor$vertices$y,
    vor$line$a, vor$line$b, vor$line$c,
    vor$msegments$line, vor$msegments$v1, vor$msegments$v2
  )
  
  
  vor$vertices
  verts <- rbind(vor$vertices, res$vertex)  
  
  
  segments(verts$x[res$segment$v1], verts$y[res$segment$v1], verts$x[res$segment$v2], verts$y[res$segment$v2], col = 'blue', lwd = 3)
  
  
  
}




if (FALSE) {
  
  set.seed(2)
  N <- 100
  vor <- voronoi(runif(N), runif(N)) |> bench::mark()
  
  plot(vor$vertices, asp = 1, ann = F, axes = F, xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5))
  
  cols <- rainbow(length(vor$polygons))
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]], col = cols[i])  
  }
  
}









