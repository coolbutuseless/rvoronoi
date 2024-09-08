

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Merge vertices
#' 
#' @param x,y vertex coordinates
#' @param line edges defined by line two indices into the x, y coordinates
#' @param v1,v2 edges defined by line two indices into the x, y coordinates
#' @param tol tolerance for merging close vertices. 
#' @param verbosity verbosity level. Default: 0
#' @return \code{list(line, v1, v2)}
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_vertices <- function(x, y, line, v1, v2, tol = 1e-5, verbosity = 0L) {
  .Call(merge_vertices_, x, y, line, v1, v2, tol, verbosity)
}





if (FALSE) {
  
  
  theta <- seq(0, 2*pi, length.out = 6)[-1]
  x <- c(0, cos(theta), 2 * cos(theta))
  y <- c(0, sin(theta), 2 * sin(theta))
  
  vor <- voronoi(x, y)
  
  vor$segment
  segs <- merge_vertices(vor$vertex$x, vor$vertex$y, vor$segment$line, vor$segment$v1, vor$segment$v2)

  polys <- extract_polygons(vor$vertex$x, vor$vertex$y, segs$v1, segs$v2)
  
  cols <- rainbow(length(polys))
  plot(vor$vertex, asp = 1, ann = F, axes = FALSE)
  for (i in seq_along(polys)) {
    polygon(polys[[i]], col = cols[i])
  }
  
}  



if (FALSE) {
  
  set.seed(1)
  N <- 1000
  x <- runif(N)
  y <- runif(N)
  
  voronoi(x, y) |> bench::mark()
  
  vor <- voronoi(x, y)
  
  vor$segment
  segs <- merge_vertices(vor$vertex$x, vor$vertex$y, vor$segment$line, vor$segment$v1, vor$segment$v2)
  
  polys <- extract_polygons(vor$vertex$x, vor$vertex$y, segs$v1, segs$v2)
  extract_polygons(vor$vertex$x, vor$vertex$y, segs$v1, segs$v2) |> bench::mark()
  # Linear search: N = 3000,  24.6 itr/sec
  
  cols <- rainbow(length(polys))
  plot(vor$vertex, asp = 1, ann = F, axes = FALSE, xlim = c(-0.2, 1.2), ylim = c(-0.2, 1.2))
  for (i in seq_along(polys)) {
    polygon(polys[[i]], col = cols[i])
  }
  
}  

