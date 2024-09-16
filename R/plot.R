


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot a Voronoi tessellation
#' 
#' Barebones plotting of Voronoi
#' 
#' @param x object returned by \code{voronoi()}
#' @param sites,labels,verts,polys,fsegs,isegs,bounds logical values. Should this geometric feature 
#'        be plotted? Default: TRUE
#' @param ... other arguments passed to \code{plot()}
#' @param site_pch,site_cex,site_col parameters for sites
#' @param vert_pch,vert_cex,vert_col parameters for Voronoi vertices
#' @param fseg_col,iseg_col colours for finite and infinite segments
#' @param buffer buffer around extents. default: 0
#' @param poly_col vector of colours for polygons
#' @param label_cex size of text labels for sites
#' @return None.
#' 
#' @importFrom graphics par rect points segments
#' @importFrom grDevices rainbow
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.vor <- function(
    x, 
    sites     = TRUE, 
    labels    = !sites,
    verts     = TRUE,
    polys     = TRUE,
    fsegs     = !polys,
    isegs     = !polys,
    bounds    = TRUE,
    site_pch  = '.', site_cex = 1.0, site_col = 'black',
    vert_pch  =  19, vert_cex = 0.3, vert_col = 'black',
    label_cex = 0.5,
    fseg_col  = 'black',
    iseg_col  = 'red',
    poly_col  = rainbow(vor$polygons$npolygons),
    buffer    = 0,
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
    pch  = vert_pch,
    cex  = vert_cex,
    col  = vert_col,
    ...
  )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Finite segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(fsegs)) {
    fseg <- subset(vor$segments, v1 > 0 & v2 > 0)
    
    graphics::segments(
      vor$vertices$x[fseg$v1],
      vor$vertices$y[fseg$v1],
      vor$vertices$x[fseg$v2],
      vor$vertices$y[fseg$v2],
      col = fseg_col
    )
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Infinte segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(isegs)) {
    # Extract segments which are unbounded. i.e. v1 or v2 is NA
    segments <- vor$segments[vor$segments$v1 <=  0 & vor$segments$v2 <= 0, ]
    
    # Get the lines associated with these segments
    lines <- vor$lines[segments$line,]
    
    xlo <- vor$extents$xmin 
    xhi <- vor$extents$xmax
    ylo <- with(lines, (c - a * xlo)/b)
    yhi <- with(lines, (c - a * xhi)/b)
    
    for (i in seq_len(nrow(segments))) {
      graphics::segments(xlo, xlo[i], xhi, yhi[i], col = iseg_col)
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Semi-infinite segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(isegs)) {
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
        graphics::segments(vor$vertices$x[seg$v1], vor$vertices$y[seg$v1], xhi, yhi[i], 
                           col = iseg_col)
      } else if (seg$v2 > 0) {
        graphics::segments(xlo, ylo[i], vor$vertices$x[seg$v2], vor$vertices$y[seg$v2], 
                           col = iseg_col)
      } else {
        stop("Shouldn't get here")
      }
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Polygons
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(polys) && vor$polygons$npolygons > 0) {
    
    # Setup fill colours for polygon
    if (length(poly_col) == 1) {
      poly_col <- rep(poly_col, vor$polygons$npolygons)
    }
    
    if (length(poly_col) != vor$polygons$npolygons) {
      stop("'poly_col' must be length 1 or match vor$polygons$npolygons")
    }
    
    for (i in seq_len(vor$polygons$npolygons)) {
      graphics::polygon(vor$polygons$individual[[i]], col = poly_col[i])
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot sites
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(verts)) {
    graphics::points(
      x    = vor$vertices$x, 
      y    = vor$vertices$y, 
      pch  = vert_pch,
      cex  = vert_cex,
      col  = vert_col,
    )
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(sites)) {
    graphics::points(
      x    = vor$sites$x, 
      y    = vor$sites$y, 
      pch  = site_pch,
      cex  = site_cex,
      col  = site_col
    )
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Draw site labels
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(labels)) {
    graphics::text(
      x      = vor$sites$x, 
      y      = vor$sites$y,
      labels = seq_along(vor$sites$x),
      cex    = label_cex
    )
  }

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Tight boundary  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isTRUE(bounds)) {
    graphics::rect(
      vor$extents$xmin, vor$extents$ymin, 
      vor$extents$xmax, vor$extents$ymax, 
      lty = 3, border = 'grey80'
    )
  }
  
  

  invisible(vor)
}



if (FALSE) {
  
  set.seed(1)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  vor <- voronoi(x, y)
  
  plot(vor, polys = FALSE) 
  
  dim(vor$polygons$coords)
    
  
}








