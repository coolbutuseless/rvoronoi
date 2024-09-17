
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Voronoi Tessellation
#' 
#' @param x,y coordinates of seed sites. Duplicate points are not allowed.
#' @param calc_polygons Logical. Should voronoi polygons be calculated?
#'        Default: TRUE
#' @param match_sites Logical. Should the polygons be re-ordered to match 
#'        the seed points? Default: TRUE. This option only makes sense
#'        when \code{calc_polygons = TRUE}.
#' @param bound_segments logical. Default: TRUE.  If \code{calc_polygons = TRUE}
#'        then this option is always TRUE.
#' @return named list of data.frames.
#' \describe{
#'   \item{sites}{data.frame of original sites}
#'   \item{vertices}{data.frame of voronoi vertices. This is the raw output from Fortune's algorithm}
#'   \item{segments}{data.frame of segments defined by 'line', 'v1' and 'v2'.  'line'
#'         is the row index into the 'lines' data.frame.  'v1' and 'v2' are 
#'         the indices into 'vertex' which define the endpoints along the 
#'         specified 'line'.  If 'v1' or 'v2' are negative this indicates that the 
#'         segment continues to infinity}
#'   \item{polygons}{list of polygon information for each voronoi cell}
#'   \item{lines}{data.frame of line equations in the tessellation of the form 'ax + by = c'}
#'   \item{extents}{list of (xmin, ymin), (xmax, ymax) bounding box encompassing
#'   input sites and voronoi vertices}
#' }
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' voronoi(x, y)
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
voronoi <- function(x, y, calc_polygons = TRUE, match_sites = TRUE, bound_segments = TRUE) {
  .Call(voronoi_, x, y, calc_polygons, match_sites, bound_segments)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Delaunay Triangulation
#' 
#' @inheritParams voronoi
#' @param calc_polygons calculate polygon coordinates default: TRUE
#' @param calc_areas calculate the areas of each triangle. Default: FALSE
#' @param calc_segments calculate segments. Default: FALSE
#' @return named list of data
#' \describe{
#'  \item{segments}{data.frame. Each row specifies the indices of the three seed points
#'         which define a single delaunay triangle}
#'  \item{polygons}{data.frame of polygon coordinates}
#' }
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' delaunay(x, y)
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delaunay <- function(x, y, calc_polygons = TRUE, calc_areas = FALSE, 
                     calc_segments = FALSE) {
  .Call(delaunay_, x, y, calc_polygons, calc_areas, calc_segments)
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



