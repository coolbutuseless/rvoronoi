

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot voronoi
#' 
#' @param vor voronoi object as returned by \code{voronoi()}
#' @param fill fill colours for polygons
#' @return None
#'
#' @import grid
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
voronoiGrob <- function(vor, fill = rainbow(length(vor$polygons))) {
  
  scale <- with(
    vor$extents, 
    max(xmax - xmin, ymax - ymin)
  )
  scale <- 1/scale
  
  
  if (length(fill) == 1) {
    fill <- rep(fill, length(vor$polygons))
  }
  if (length(fill) != length(vor$polygons)) {
    stop("Bad fill length: ", length(fill))
  }
  
  
  xfix <- function(x) { (x - vor$extents$xmin) * scale }
  yfix <- function(y) { (y - vor$extents$ymin) * scale }
  
  polys <- lapply(seq_along(vor$polygons), function(i) {
      grid::polygonGrob(
        x = xfix(vor$polygons[[i]]$x), 
        y = yfix(vor$polygons[[i]]$y), 
        gp = grid::gpar(fill = fill[i]),
        default.units = 'snpc'
      )
  })
  
  polys <- do.call(grid::grobTree, polys)
  
  
  # points <- pointsGrob(
  #   x = xfix(x), 
  #   y = yfix(y), 
  #   pch = '+', default.units = 'snpc'
  # )
  
  
  boundary <- with(
    vor$extents,
    grid::rectGrob(
      x      = xfix( (xmin + xmax) / 2 ),
      y      = yfix( (ymin + ymax) / 2 ),
      width  = (xmax - xmin) * scale,
      height = (ymax - ymin) * scale,
      gp = grid::gpar(fill = NA, col = 'black', lty = 2, lwd = 4),
      default.units = 'snpc'
    )
  )
  
  grobTree(polys, boundary)
}



if (FALSE) {
  
  set.seed(5)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  
  vor <- voronoi(x, y)
  
  grob <- voronoiGrob(vor)
  grid.newpage(); grid.draw(grob)
  
  grob$vp <- viewport(default.units = 'npc', width = 1, height = 0.25)
  grid.newpage(); grid.draw(grob)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Manual 'graphics' plottong
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plot(x, y, asp = 1, ann = FALSE, axes = FALSE, col = 'red', 
       xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1), mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
  
  cols <- rainbow(length(vor$polygons))
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]], col = cols[i])
  }
  
  points(x, y, pch = 19)
  
  
  
}