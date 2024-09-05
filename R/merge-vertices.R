

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Merge vertices
#' 
#' @param x,y vertex coordinates
#' @param v1,v2 edges defined by two indices into the x, y coordinates
#' @param tol tolerance for merging close vertices. 
#' @param verbosity verbosity level. Default: 0
#' @return \code{list(v1, v2)}
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_vertices <- function(x, y, v1, v2, tol = 1e-5, verbosity = 0L) {
  .Call(merge_vertices_, x, y, v1, v2, tol, verbosity)
}





if (FALSE) {
  
  
  theta <- seq(0, 2*pi, length.out = 6)[-1]
  x <- c(0, cos(theta), 2 * cos(theta))
  y <- c(0, sin(theta), 2 * sin(theta))
  
  vor <- voronoi(x, y)
  
  vor$segment
  segs <- merge_vertices(vor$vertex$x, vor$vertex$y, vor$segment$v1, vor$segment$v2)

  polys <- extract_polygons(vor$vertex$x, vor$vertex$y, segs$v1, segs$v2)
  
  cols <- rainbow(length(polys))
  plot(vor$vertex, asp = 1, ann = F, axes = FALSE)
  for (i in seq_along(polys)) {
    polygon(polys[[i]], col = cols[i])
  }
  
}  

