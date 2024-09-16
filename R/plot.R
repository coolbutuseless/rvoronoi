


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot a Voronoi tessellation
#' 
#' Barebones plotting of Voronoi
#' 
#' @param x object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param pch_site plot character for sites
#' @param pch_vertex plot character for Voronoi vertices
#' @param buffer buffer around extents. defualt: 0
#' @param col_polygons vector of colours for polygons
#' @return None.
#' 
#' @importFrom graphics par rect points segments
#' @importFrom grDevices rainbow
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.vor <- function(
    x, 
    pch_site   = '.', 
    pch_vertex = 19, 
    col_polygons = rainbow(length(vor$polygons)),
    buffer = 0,
    ...
) {
  
  vor <- x
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setup
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  xlim = c(vor$extents$xmin - buffer, vor$extents$xmax + buffer)
  ylim = c(vor$extents$ymin - buffer, vor$extents$ymax + buffer)
  
  oldpar <- graphics::par(mar = c(0, 0, 0, 0))
  on.exit(par(oldpar))
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialise the plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plot(
    x    = vor$vertices$x, 
    y    = vor$vertices$y, 
    asp  = 1, ann = FALSE, axes = FALSE, 
    xlim = xlim, 
    ylim = ylim, 
    pch  = pch_vertex,
    ...
  )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Finite segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    fseg <- subset(vor$segments, v1 > 0 & v2 > 0)
    
    graphics::segments(
      vor$vertices$x[fseg$v1],
      vor$vertices$y[fseg$v1],
      vor$vertices$x[fseg$v2],
      vor$vertices$y[fseg$v2]
    )
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Infinte segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    # Extract segments which are unbounded. i.e. v1 or v2 is NA
    segments <- vor$segments[vor$segments$v1 <=  0 & vor$segments$v2 <= 0, ]
    
    # Get the lines associated with these segments
    lines <- vor$lines[segments$line,]
    
    xlo <- vor$extents$xmin 
    xhi <- vor$extents$xmax
    ylo <- with(lines, (c - a * xlo)/b)
    yhi <- with(lines, (c - a * xhi)/b)
    
    for (i in seq_len(nrow(segments))) {
      graphics::segments(xlo, xlo[i], xhi, yhi[i])
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Semi-infinite segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    # Extract segments which are unbounded. i.e. v1 or v2 is NA
    segments <- vor$segments[xor(vor$segments$v1 <=  0, vor$segments$v2 <= 0), ]
    
    # Get the lines associated with these segments
    lines <- vor$lines[segments$line,]
    
    xlo <- vor$extents$xmin 
    xhi <- vor$extents$xmax
    ylo <- with(lines, (c - a * xlo)/b)
    yhi <- with(lines, (c - a * xhi)/b)
    
    for (i in seq_len(nrow(segments))) {
      seg <- segments[i, ]
      if (seg$v1 <= 0 && seg$v2 <= 0) {
        # warning("Double NA segment: ", i)
      } else if (seg$v1 > 0) {
        graphics::segments(vor$vertices$x[seg$v1], vor$vertices$y[seg$v1], xhi, yhi[i])
      } else if (seg$v2 > 0) {
        graphics::segments(xlo, ylo[i], vor$vertices$x[seg$v2], vor$vertices$y[seg$v2])
      } else {
        stop("Shouldn't get here")
      }
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Polygons
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    
    if (length(col_polygons) == 1) {
      col_polygons <- rep(col_polygons, length(vor$polygons))
    }
    
    if (length(col_polygons) != length(vor$polygons)) {
      stop("'col_polygons' must be length 1 or match length(vor$polygons)")
    }
    
    
    for (i in seq_along(vor$polygons)) {
      graphics::polygon(vor$polygons[[i]], col = col_polygons[i])
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot sites
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    graphics::points(
      x    = vor$vertices$x, 
      y    = vor$vertices$y, 
      pch  = pch_site
    )
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (TRUE) {
    graphics::points(
      x    = vor$sites$x, 
      y    = vor$sites$y, 
      pch  = pch_vertex
    )
  }

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Tight boundary  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  graphics::rect(
    vor$extents$xmin, vor$extents$ymin, 
    vor$extents$xmax, vor$extents$ymax, 
    lty = 3, border = 'grey80'
  )
  
  

  invisible(vor)
}



if (FALSE) {
  
  set.seed(1)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  vor <- voronoi(x, y)
  
  plot(vor) 
}








