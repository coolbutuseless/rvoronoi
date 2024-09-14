
globalVariables(c('v1', 'v2'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Merge close vertices and remove redundant segments
#' 
#' @param vor voronoi output
#' @return vor adjusted 'vor' object
#' @importFrom stats dist setNames
#' @noRd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_vertices_r <- function(vor) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate distance between all point pairs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vv <- vor$vertices
  d <- as.matrix(dist(vv))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Data frame of (v1, v2) index pairs which are close
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pairs <- which(d < 0.0001) |> 
    arrayInd(.dim = dim(d)) |>
    as.data.frame() |>
    setNames(c('v1', 'v2')) |> 
    subset(v1 < v2)
  
  print(pairs)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Renumber the indices in vor$segment.
  # Renumber backwards to ensure connection sanity
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vor$segment <- as.data.frame(vor$segment)
  ss <- vor$segment 
  for (i in rev(seq_len(nrow(pairs)))) {
    row <- pairs[i, ]
    ss$v1[ ss$v1 > 0 & ss$v1 == row$v2 ] <- row$v1
    ss$v2[ ss$v2 > 0 & ss$v2 == row$v2 ] <- row$v1
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Due to the way Fortune's works, when we've renumbered close vertices,
  # we end up with redundundant segments where v1 == v2
  # Remove these segments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ss <- ss[ss$v1 <= 0 | ss$v2 <= 0 | ss$v1 != ss$v2, ]
  vor$segment <- ss
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return the adjusted vor
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vor
}
