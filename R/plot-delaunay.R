

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Simple plot of Delaunay triangulation
#' 
#' @param x output from \code{delaunay()}
#' @param sites,tris,segments logical. draw geometric feature
#' @param site_pch,site_col graphics parameters for each input site
#' @param segment_col colour for segments (length = 1)
#' @param tri_col colour for polygons (length = 1 or N)
#' @param ... other arguments passed to \code{plot()}
#' 
#' @return None
#' @importFrom graphics segments points polygon
#' @importFrom grDevices rainbow
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.del <- function(x, 
                     sites = TRUE,
                     tris  = TRUE,
                     segments = FALSE,
                     site_pch = 19, site_col = 'black', 
                     segment_col = 'black', 
                     tri_col = rainbow(x$ntris), ...) {
  
  force(tri_col)
  del <- x
  x <- del$sites$x
  y <- del$sites$y
  
  v1 <- del$tris$v1
  v2 <- del$tris$v2
  v3 <- del$tris$v3
  
  # Plot the sites
  plot(
    x    = x, 
    y    = y, 
    asp  = 1, 
    ann  = F, 
    axes = F, 
    pch  = site_pch,
    col  = NA,
    ...
  )
  
  
  if (isTRUE(tris)) {
    if (length(tri_col) == 1) {
      tri_col <- rep(tri_col, del$ntris)
    }
    
    if (length(tri_col) != del$ntris) {
      stop("Bad length for tri_col");
    }
    
    tris <- split(del$polygons, del$polygons$idx)
    for (i in seq_along(tris)) {
      polygon(tris[[i]]$x, tris[[i]]$y, col = tri_col[i])
    }
  }

  if (isTRUE(segments)) {  
    # Plot all finite segments.  
    graphics::segments(x[v1], y[v1], x[v2], y[v2], col = segment_col)
    graphics::segments(x[v3], y[v3], x[v2], y[v2], col = segment_col)
    graphics::segments(x[v1], y[v1], x[v3], y[v3], col = segment_col)
  }
  
  if (isTRUE(sites)) {
    graphics::points(x, y, pch = site_pch, col = site_col)
  }
  
  
  invisible(del)
}


if (FALSE) {
  
  set.seed(1)
  N <- 99
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  plot(del, tris = FALSE, segments = TRUE)
  
  plot(x, y, asp = 1, ann = F, axes = F)
  with(del$segments, segments(x1, y1, x2, y2, col = grey(del$segments$dist), lwd = 2))
  
  
  library(grid)
  grid.polygon(x = del$polygons$x, y = del$polygons$y,
               id = del$polygons$idx, gp = gpar(fill = grey(del$areas$area ^ (0.5))))
}



if (FALSE) {
  set.seed(1)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  
}





