
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi segments to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @return None.
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_segments <- function(vor, ...) {
  
  # finite segments
  fseg <- subset(vor$segment, v1 > 0 & v2 > 0)
  
  segments(
    vor$vertex$x[fseg$v1],
    vor$vertex$y[fseg$v1],
    vor$vertex$x[fseg$v2],
    vor$vertex$y[fseg$v2],
    ...
  )
  
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
draw_inf_segments <- function(vor, col = 'grey90', ...) {
  # Extract segments which are unbounded. i.e. v1 or v2 is NA
  inf_seg <- vor$segment[vor$segment$v1 <=  0 | vor$segment$v2 <= 0, ]
  
  # Get the lines associated with these segments
  inf_line <- vor$line[inf_seg$line,]
  
  xlo <- vor$extents$xmin 
  xhi <- vor$extents$xmax
  ylo <- with(inf_line, (c - a * xlo)/b)
  yhi <- with(inf_line, (c - a * xhi)/b)
  
  for (i in seq_len(nrow(inf_seg))) {
    seg <- inf_seg[i, ]
    if (seg$v1 <= 0 && seg$v2 <= 0) {
      # warning("Double NA segment: ", i)
    } else if (seg$v1 > 0) {
      segments(vor$vertex$x[seg$v1], vor$vertex$y[seg$v1], xhi, yhi[i], col = col, ...)
    } else if (seg$v2 > 0) {
      segments(xlo, ylo[i], vor$vertex$x[seg$v2], vor$vertex$y[seg$v2], col = col, ...)
    } else {
      stop("Shouldn't get here")
    }
  }
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi unbounded lines to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param lty linetype. default: 1
#' @param col colour. default: 'grey90'
#' @return None.
#' @importFrom graphics segments
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_inf_lines <- function(vor, lty = 1, col = 'grey90', ...) {
  # Extract segments which are unbounded. i.e. v1 or v2 is NA
  inf_seg <- vor$segment[vor$segment$v1 <= 0 | vor$segment$v2 <= 0, ]
  
  # Get the lines associated with these segments
  inf_line <- vor$line[inf_seg$line,]
  
  xlo <- vor$extents$xmin 
  xhi <- vor$extents$xmax
  ylo <- with(inf_line, (c - a * xlo)/b)
  yhi <- with(inf_line, (c - a * xhi)/b)
  
  segments(xlo, ylo, xhi, yhi, col = col, lty = lty, ...)
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi bounded polygons to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @return None.
#' @importFrom graphics polygon
#' @importFrom grDevices heat.colors rainbow
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_bounded_polygons <- function(vor, ...) {
  
  # segs <- merge_vertices(vor$vertex$x, vor$vertex$y, vor$segment$v1, vor$segment$v2)
  # ps <- extract_polygons(vor$vertex$x, vor$vertex$y, segs$v1, segs$v2)
  ps <- vor$polygons
  
  cols <- rainbow(length(ps))
  
  for (i in seq_along(ps)) {
    polygon(ps[[i]], col = cols[i], ...)
  }
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi unbounded polygons to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @return None.
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_unbounded_polygons <- function(vor, ...) {
  
  invisible(vor)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot the vertices of the voronoi tesselation
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @param asp aspect ration. Default 1
#' @param ann annotation? Default: FALSE
#' @param axes draw axes? Default: FALSE
#' @param buffer buffer around extents. defualt: 0
#' @return None.
#' @importFrom graphics rect
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_vor <- function(vor, buffer = 0, asp = 1, ann = FALSE, axes = FALSE, ...) {
  
  xlim = c(vor$extents$xmin - buffer, vor$extents$xmax + buffer)
  ylim = c(vor$extents$ymin - buffer, vor$extents$ymax + buffer)
  
  plot(vor$vertex$x, vor$vertex$y, asp = asp, ann = ann, axes = axes, 
       xlim = xlim, ylim = ylim, ...)
  
  rect(vor$extents$xmin, vor$extents$ymin, vor$extents$xmax, vor$extents$ymax, 
       lty = 3, border = 'grey80')
  
  rect(vor$extents$xmin - buffer, vor$extents$ymin - buffer, 
       vor$extents$xmax + buffer, vor$extents$ymax + buffer, 
       lty = 3, border = 'grey80')
  
  
  invisible(vor)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Draw vertex nodes
#' @param vor vor
#' @param nudge_x,nudge_y text adjustment
#' @param ... other arguments passed through
#' @importFrom graphics points
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_vertices <- function(vor, nudge_x = 0.1, nudge_y = nudge_x, ...) {
  points(vor$vertex, ...)
  # text(x = vor$vertex$x + nudge_x, 
  # y = vor$vertex$y + nudge_y, labels = seq_along(vor$vertex$x))
  
  
  invisible(vor)
}
