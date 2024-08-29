
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Voronoi tesselation on an infinite plane
#' @param x,y coordinates of seed points
#' @return named list of data.frames.
#' \describe{
#'   \item{vertex}{data.frame of vertex coordinates in the tesselation}
#'   \item{line}{data.frame of line equations in the tesselation of the form 'ax + by = c'}
#'   \item{segment}{data.frame of segments defined by 'line', 'v1' and 'v2'.  'line'
#'         is the row index into the 'line' data.frame.  'v1' and 'v2' are 
#'         the indices into 'vertex' which define the endpoints along the 
#'         specified 'line'.  If 'v1' or 'v2' are NA this indicates that the 
#'         segment continues to infinity}
#' }
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' voronoi(x, y)
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
voronoi <- function(x, y) {
  .Call(voronoi_, x, y)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate Voronoi
#' @inheritParams voronoi
#' @return data.frame. Each row specifies the indices of the three seed points
#'         which define a single delaunay triangle
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' delaunay(x, y)
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delaunay <- function(x, y) {
  .Call(delaunay_, x, y)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Normalise a delaunay triangulation so it can be conpared to other triangulations
#' 
#' @param del data.frame with three columns representing the output of a 
#'        delaunay triangulation.  Each column specifies the index of the
#'        original seed point which is the vertex of a triangle
#' @noRd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normalise_del <- function(del) {
  del <- del |> 
    apply(1, sort, simplify = F) |> 
    do.call(rbind, args = _) |>
    as.data.frame() 
  
  names(del) <- c('v1', 'v2', 'v3')
  
  del <- del[with(del, order(v1, v2, v3)),]
  rownames(del) <- NULL
  del
}


