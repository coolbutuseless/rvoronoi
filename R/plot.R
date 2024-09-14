# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi segments to an existing plot
#'
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#'
#' @return None.
#' @importFrom graphics segments
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_finite_segments <- function(vor, ...) {
  
  # finite segments
  fseg <- subset(vor$segments, v1 > 0 & v2 > 0)
  
  segments(
    vor$vertices$x[fseg$v1],
    vor$vertices$y[fseg$v1],
    vor$vertices$x[fseg$v2],
    vor$vertices$y[fseg$v2],
    ...
  )
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add semi-infinite voronoi segments to an existing plot
#'
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param col segment colour. Default: grey90
#' @return None.
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_semi_infinite_segments <- function(vor, col = 'grey90', ...) {
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
      segments(vor$vertices$x[seg$v1], vor$vertices$y[seg$v1], xhi, yhi[i], col = col, ...)
    } else if (seg$v2 > 0) {
      segments(xlo, ylo[i], vor$vertices$x[seg$v2], vor$vertices$y[seg$v2], col = col, ...)
    } else {
      stop("Shouldn't get here")
    }
  }
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi unbounded segments to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param col segment colour. Default: grey90
#' @return None.
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_infinite_segments <- function(vor, col = 'grey90', ...) {
  # Extract segments which are unbounded. i.e. v1 or v2 is NA
  segments <- vor$segments[vor$segments$v1 <=  0 & vor$segments$v2 <= 0, ]
  
  # Get the lines associated with these segments
  lines <- vor$lines[segments$line,]
  
  xlo <- vor$extents$xmin 
  xhi <- vor$extents$xmax
  ylo <- with(lines, (c - a * xlo)/b)
  yhi <- with(lines, (c - a * xhi)/b)
  
  for (i in seq_len(nrow(segments))) {
    segments(xlo, xlo[i], xhi, yhi[i], col = col, ...)
  }
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi bounded polygons to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param cols vector of colours
#' @param ... other arguments passed to \code{plot()}
#' @return None.
#' @importFrom graphics polygon
#' @importFrom grDevices heat.colors rainbow
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_polygons <- function(vor, cols = rainbow(length(vor$polygons)), ...) {
  
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]], col = cols[i], ...)
  }
  
  invisible(vor)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Draw vertex nodes
#' @param vor vor
#' @param pch plot character
#' @param nudge_x,nudge_y text adjustment
#' @param text show text
#' @param ... other arguments passed through
#' @importFrom graphics points
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_vertices <- function(vor, pch = 19, text = FALSE, nudge_x = 0.1, nudge_y = nudge_x, ...) {
  points(vor$vertices, pch = pch, ...)
  
  if (isTRUE(text)) {
    text(
      x = vor$vertices$x + nudge_x,
      y = vor$vertices$y + nudge_y, 
      labels = seq_along(vor$vertices$x)
    )
  }  
  invisible(vor)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Draw vertex nodes
#' @param vor vor
#' @param cex point size scale
#' @param pch plot character
#' @param text show text?
#' @param ... other arguments passed through
#' @importFrom graphics points
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_sites <- function(vor, pch = '+', cex = 0.5, text = FALSE, ...) {
  
  
  if (isTRUE(text)) {
    text(
      x = vor$sites$x,
      y = vor$sites$y, 
      labels = seq_along(vor$sites$x)
    )
    points(vor$sites, pch = '.', col = 'red', ...)
  } else {
    points(vor$sites, pch = pch, cex = cex, ...)
  }
   
  
  invisible(vor)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot the vertices of the voronoi tesselation
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param asp aspect ration. Default 1
#' @param ann annotation? Default: FALSE
#' @param pch plot character
#' @param axes draw axes? Default: FALSE
#' @param buffer buffer around extents. defualt: 0
#' @return None.
#' @importFrom graphics rect
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_vor <- function(vor, buffer = 0, asp = 1, ann = FALSE, axes = FALSE, pch = '.', ...) {
  
  xlim = c(vor$extents$xmin - buffer, vor$extents$xmax + buffer)
  ylim = c(vor$extents$ymin - buffer, vor$extents$ymax + buffer)
  
  par(mar = c(0, 0, 0, 0))
  plot(vor$vertices$x, vor$vertices$y, asp = asp, ann = ann, axes = axes, 
       xlim = xlim, ylim = ylim, ...)
  
  rect(vor$extents$xmin, vor$extents$ymin, vor$extents$xmax, vor$extents$ymax, 
       lty = 3, border = 'grey80')
  
  rect(vor$extents$xmin - buffer, vor$extents$ymin - buffer, 
       vor$extents$xmax + buffer, vor$extents$ymax + buffer, 
       lty = 3, border = 'grey80')
  
  
  invisible(vor)
}



if (FALSE) {
  
  set.seed(1)
  N <- 10
  x <- runif(N)
  y <- runif(N)
  vor <- voronoi(x, y)
  
  plot_vor(vor) |>
    draw_finite_segments() |>
    draw_semi_infinite_segments(col = 'hotpink') |>
    draw_infinite_segments(col = 'red') |>
    draw_polygons() |>
    draw_sites(text = TRUE) |>
    draw_vertices()
  
  
  
  
  
}








