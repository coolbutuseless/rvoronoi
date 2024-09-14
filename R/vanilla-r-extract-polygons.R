
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
#' @param x,y vertex coordinates
#' @param v1,v2 edge definition
#' @param verbosity verbosity level. default: 0
#' @return list of polygons
#' @noRd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_polygons_r <- function(x, y, v1, v2, verbosity = 0) {
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Edge list = pairs of vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  e <- data.frame(v1 = v1, v2 = v2)
  
  # Remove any unbounded vertices
  idx <- with(e, v1 >= 0 & v2 >= 0) 
  e <- e[idx, ]
  
  if (nrow(e) == 0) {
    warning("No bounded edges at all!")
    return(NULL)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Phase 1: Find all the wedges
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # 1. Split undirected edge into two directed edges
  e <- rbind(e, data.frame(v1 = e$v2, v2 = e$v1))
  
  
  # 2 Add angle 
  theta <- atan2(y[e$v2] - y[e$v1], x[e$v2] - x[e$v1])
  theta <- ifelse(theta < 0, theta + 2 * pi, theta)
  
  e$theta <- theta
  
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
  
  if (verbosity > 0) {
    print(as.data.frame(wedge))
  }
  
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
  #         each vector is a set of vertex indices into vor$vertices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vidxs <- lapply(regions, function(region) {
    region_wedges <- wedge[region, ]
    region_verts  <- region_wedges$v1
    region_verts
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # List of polygons. Each polygon is a data.frame of x,y coordinates
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vert <- data.frame(x = x, y = y)
  polygons <- lapply(vidxs, function(idxs) {vert[idxs,]})
  
  class(polygons) <- union(class(polygons), "vorpolygons")
  
  areas <- vapply(polygons, polygon_bbox_area, double(1))
  biggest <- which.max(areas)  
  
  polygons[-biggest]
  # polygons
}
