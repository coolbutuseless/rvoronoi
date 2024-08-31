
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Bounding box area of a polygon
#' 
#' Use this to find the overall large bounded polygon encompassing 
#' all the polygons - an artefact of how this algo works that is not needed
#' in the results
#' 
#' @param p list containing 'x' and 'y' members
#' @return area
#' @noRd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
polygon_bbox_area <- function(p) {
  diff(range(p$x)) * diff(range(p$y))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Extract all bounded polygons from a voronoi tessellation (R version)
#' 
#' @param vor result of \code{voronoi()}
#' @return list of polygons
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_polygons_r <- function(vor) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x <- vor$vertex$x
  y <- vor$vertex$y
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Edge list = pairs of vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  e <- vor$segment
  e$line <- NULL
  
  # Remove any unbounded vertices
  idx <- with(e, !is.na(v1) & !is.na(v2)) 
  e <- e[idx, ]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Phase 1: Find all the wedges
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # 1. Split undirected edge into two directed edges
  e <- rbind(e, data.frame(v1 = e$v2, v2 = e$v1))
  
  
  # 2 Add angle 
  theta <- atan2(y[e$v2] - y[e$v1], x[e$v2] - x[e$v1])
  theta <- ifelse(theta < 0, theta + 2 * pi, theta)
  
  e$theta <- theta
  e
  
  # 3. Sort by v1, theta
  e <- e[with(e, order(v1, theta)), ]
  
  # 4. Form wedges
  egroups <- split(e, e$v1)
  egroups <- lapply(egroups, function(df) {
    df$v0 <- with(df, c(v2[length(v2)], v2[-length(v2)])) # Roll
    df
  })
  wedge <- do.call(rbind, egroups)
  wedge <- wedge[, c('v0', 'v1', 'v2')]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Phase 2:  Group wedges into regions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # List of regions.
  # Each region is a list of wedge indices
  regions <- list()
  
  # 1. Sort the wedge list by v0, v1
  #    Sort needed as actual lookup to find matching wedge should be done 
  #    by binary search.
  wedge <- wedge[with(wedge, order(v0, v1)), ]
  
  # 2. Mark all wedges as unused
  wedge$used <- FALSE
  
  while(TRUE) {
    # 3. Find the next unused wedge - mark as 'used', add to region
    #    If no unused wedge, then algorithm complete
    orig <- which(!wedge$used)[1]
    if (length(orig) == 0 || is.na(orig)) {
      # message("All polygons extracted")
      break
    }
    region  <- c(orig); wedge$used[orig] <- TRUE
    current <- orig
    
    # 4. Search for matching continuing wedge.
    #    E.g. if initial wedge is c(1, 4, 7), look for wedge c(4, 7, x)
    #    Add to region
    while (TRUE) {
      match  <- with(wedge, which(!used & v0 == wedge$v1[current] & v1 == wedge$v2[current]))
      region <- c(region, match); wedge$used[match] <- TRUE
      current <- match
      
      # 5. If new wedge overlaps original wedge, 
      #       then: region is extracted. Go to 3
      #       else: go to 4
      if (wedge$v1[match] == wedge$v0[orig] & wedge$v2[match] == wedge$v1[orig]) {
        regions <- c(regions, list(region))
        # message("Finished region: ", deparse1(region))
        break
      }
    } # Return to step 4
  } # Return to step 3
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Convert the list of wedges for each region
  # 'vidxs' is a list of vectors
  #         each vector is a set of vertex indices into vor$vertex
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vidxs <- lapply(regions, function(region) {
    region_wedges <- wedge[region, ]
    region_verts  <- region_wedges$v1
    region_verts
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # List of polygons. Each polygon is a data.frame of x,y coordinates
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  polygons <- lapply(vidxs, function(idxs) {vor$vertex[idxs,]})
  class(polygons) <- union(class(polygons), "vorpolygons")
  
  areas <- vapply(polygons, polygon_bbox_area, double(1))
  biggest <- which.max(areas)  
  
  polygons[-biggest]
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Add voronoi segments to an existing plot
#' @param vor object returned by \code{voronoi()}
#' @param ... other arguments passed to \code{plot()}
#' @return None.
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_segments <- function(vor, ...) {
  segments(
    vor$vertex$x[vor$segment$v1],
    vor$vertex$y[vor$segment$v1],
    vor$vertex$x[vor$segment$v2],
    vor$vertex$y[vor$segment$v2],
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
  inf_seg <- vor$segment[is.na(vor$segment$v1) | is.na(vor$segment$v2),]
  
  # Get the lines associated with these segments
  inf_line <- vor$line[inf_seg$line,]
  
  xlo <- vor$extents$xmin 
  xhi <- vor$extents$xmax
  ylo <- with(inf_line, (c - a * xlo)/b)
  yhi <- with(inf_line, (c - a * xhi)/b)
  
  for (i in seq_len(nrow(inf_seg))) {
    seg <- inf_seg[i, ]
    if (is.na(seg$v1) && is.na(seg$v2)) {
      stop("Double NA segment: ", i)
    } else if (!is.na(seg$v1)) {
      segments(vor$vertex$x[seg$v1], vor$vertex$y[seg$v1], xhi, yhi[i], col = col, ...)
    } else if (!is.na(seg$v2)) {
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
  inf_seg <- vor$segment[is.na(vor$segment$v1) | is.na(vor$segment$v2),]
  
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
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_bounded_polygons <- function(vor, ...) {
  
  ps <- extract_polygons_r(vor)
  
  for (p in ps) {
    polygon(p$x, p$y, ...)
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




if (FALSE) {
  library(rvoronoi)
  
  
  set.seed(3)
  N <- 10
  x0 <- runif(N)
  y0 <- runif(N)
  
  vor <- voronoi(x0, y0)
  
  if (FALSE) {
    ps <- extract_polygons_r(vor)
    ps
  }
  
  plot_vor(vor, buffer = 0.3) |>
    draw_segments() |> 
    draw_inf_lines() |>
    draw_inf_segments(col = 'hotpink') |>
    draw_bounded_polygons(border = 'red')
  points(x0, y0, pch = 'o') 
  
}




if (FALSE) {
  library(rvoronoi)
  
  
  set.seed(3)
  N <- 10
  x0 <- runif(N)
  y0 <- runif(N)
  
  vor <- voronoi(x0, y0)
  
  if (FALSE) {
    ps <- extract_polygons_r(vor)
    ps
  }
  
  plot_vor(vor, buffer = 0.3) |>
    draw_segments() |> 
    draw_inf_lines() |>
    draw_inf_segments(col = 'hotpink') |>
    draw_bounded_polygons(border = 'red')
  points(x0, y0, pch = 'o') 
  
}



if (FALSE) {
  
  library(polyclip)
  library(rvoronoi)
  library(ggplot2)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # My original polygon
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  set.seed(3)
  N <- 20
  x0 <- runif(N)
  y0 <- runif(N)
  idx <- grDevices::chull(x0, y0)
  mypoly <- list(x = x0[idx], y = y0[idx])
  x0 <- x0[idx]
  y0 <- y0[idx]
  
  png("working/00-mypoly.png")
  plot(x0, y0, asp = 1, ann = FALSE, axes = FALSE, pch = '.')
  polygon(mypoly)
  invisible(dev.off())
  
  # bounding box of mypoly
  xr <- range(mypoly$x)
  yr <- range(mypoly$y)
  
  # Create many random seed points and find ones within mypoly
  Ns <- 5
  set.seed(1)
  xs <- runif(Ns * 2, xr[1], xr[2])
  ys <- runif(Ns * 2, yr[1], yr[2])
  
  idx <- polyclip::pointinpolygon(list(x = xs, y = ys), mypoly)
  idx <- which(idx == 1)
  idx <- sample(idx, Ns)
  
  xs <- xs[idx]
  ys <- ys[idx]
  
  # Expand our original polygon to use as seed points 
  expanded1 <- polyclip::polyoffset(mypoly, 0.25 * max(diff(xr), diff(yr)))[[1]]
  expanded2 <- polyclip::polyoffset(mypoly, 0.50 * max(diff(xr), diff(yr)))[[1]]
  
  
  
  xs <- c(xs, expanded1$x, expanded2$x)
  ys <- c(ys, expanded1$y, expanded2$y)
  
  
  plot(xs, ys, asp = 1, ann = FALSE, axes = FALSE)
  polygon(mypoly)
  
  
  
  vor <- voronoi(xs, ys)
  
  png("working/02-voronoi.png")
  plot_vor(vor, buffer = 0.3, pch = '.') |>
    # draw_segments() |>
    # draw_inf_lines() |>
    draw_inf_segments(col = 'hotpink') |>
    draw_bounded_polygons(border = 'red')
  points(xs, ys, pch = '+') 
  invisible(dev.off())
  
  
  
  polys <- extract_polygons_r(vor)
  
  ints <- lapply(polys, function(poly) {
    polyclip::polyclip(mypoly, poly) 
  })
  ints <- unlist(ints, recursive = FALSE)

  cols <- rainbow(length(ints))
  
  png("working/03-intersect.png")
  plot_vor(vor, buffer = 0.3, pch = '.') |>
    # draw_segments() |>
    # draw_inf_lines() |>
    draw_inf_segments(col = 'hotpink') |>
    draw_bounded_polygons(border = 'red')
  points(xs, ys, pch = '+') 
  
  for (i in seq_along(ints)) {
    polygon(ints[[i]], col = cols[i])
  }
  invisible(dev.off())
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Final segmentation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  png("working/04-final.png")
  plot(x0, y0, asp = 1, ann = FALSE, axes = FALSE, pch = '.')
  for (i in seq_along(ints)) {
    polygon(ints[[i]], col = cols[i])
  }
  invisible(dev.off())
  
  
}















