

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Extract polygons
#' 
#' @inheritParams merge_vertices
#' @param xseed,yseed coorindates of the vornoi seed points.  If given then 
#'        the polygons will align with the seed points
#' @return \code{list(list(x = ..., y = ...), list(x = ..., y = ...), ...)}
#' @noRd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_polygons <- function(x, y, v1, v2, xseed = NULL, yseed = NULL, verbosity = 0L) {
  .Call(extract_polygons_, x, y, v1, v2, xseed, yseed, verbosity)
}




if (FALSE) {
  
  theta <- seq(0, 2*pi, length.out = 4)[-1]
  x <- c(0, cos(theta), 2 * cos(theta))
  y <- c(0, sin(theta), 2 * sin(theta))
  
  vor <- voronoi(x, y)
  
  vor$segment
  seg <- merge_vertices(vor$vertices$x, vor$vertices$y, vor$segment$line, vor$segment$v1, vor$segment$v2, verbosity = 0)
  seg
  
  polys <- extract_polygons(vor$vertices$x, vor$vertices$y, seg[[1]], seg[[2]])
  
  # extract_polygons  (vor$vertices$x, vor$vertices$y, seg[[1]], seg[[2]])
  # extract_polygons_r(vor$vertices$x, vor$vertices$y, seg[[1]], seg[[2]])
  
  
  cols <- rainbow(length(polys))
  plot(x, y, asp = 1, ann = FALSE, axes = FALSE)
  for (i in seq_along(polys)) {
    polygon(polys[[i]], col = cols[i])
  }
  
}


