
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
#' @param merge_tolerance Limit of how close the ends of a segment must be
#'        before the segment is collapsed to a single vertex Default: 1e-10
#' @return An object of class "vor" which is a named list of data.frames.
#' \describe{
#'   \item{sites}{data.frame of original sites}
#'   \item{vertices}{data.frame of voronoi vertices. This is the raw output from Fortune's algorithm}
#'   \item{segments}{data.frame of segments defined by 'line', 'v1' and 'v2'.  'line'
#'         is the row index into the 'lines' data.frame.  'v1' and 'v2' are 
#'         the indices into 'vertices' which define the endpoints along the 
#'         specified 'line'.  If 'v1' or 'v2' are NA this indicates that the 
#'         segment continues to infinity}
#'   \item{polygons}{
#'      polygon information for each voronoi cell. Only calculated when
#'      \code{calc_polygons = TRUE}
#'      \describe{
#'        \item{npolygons}{Number of polygons}
#'        \item{coords}{data.frame of x,y coordiantes and polygon index for all polygons}
#'        \item{bbox}{data.frame of polygon bounding box information. Each row represents a polygon}
#'        \item{centrois}{data.frame of polygon centroids. Each row represents a polygon}
#'      }
#'   }
#'   \item{lines}{data.frame of line equations in the voronoi diagram of the form 'ax + by = c'}
#'   \item{extents}{list of (xmin, ymin), (xmax, ymax) bounding box which encompasses
#'         all input sites and voronoi vertices}
#'   \item{mvertices}{Merged and bounded vertices used to calculate polygons} 
#'   \item{msegments}{Merged and bounded segments used to calculate polygons} 
#' }
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' vor <- voronoi(x, y)
#' plot(vor)
#' vor
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
voronoi <- function(x, y, calc_polygons = TRUE, match_sites = TRUE, bound_segments = TRUE,
                    merge_tolerance = 1e-10) {
  
  if (anyNA(x) || anyNA(y)) {
    stop("voronoi(): 'x' and 'y' cannot contain NA values")
  }
  .Call(voronoi_, x, y, calc_polygons, match_sites, bound_segments, merge_tolerance)
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
#'  \item{sites}{coordinates of sites}
#'  \item{ntris}{Number of triangles in this delaunay triangulation} 
#'  \item{tris}{data.frame. Each row specifies the indices of the three sites
#'         which define a single delaunay triangle. If \code{calc_areas = TRUE} then
#'         this data.frame also includes an \code{area} column}
#'  \item{polygons}{data.frame of polygon coordinates with x,y and 
#'        polygon index. Only calculated when \code{calc_polygons = TRUE}}
#'  \item{segments}{The end-points and length of each unique undirected 
#'        edge in the triangulation. Only calculated when \code{calc_segments = TRUE}
#'     \describe{
#'        \item{v1,v2}{vertex indices referencing the original \code{sites}}
#'        \item{x1,y1,x2,y2}{coordinates of ends of segment}
#'        \item{len}{length of this segment}
#'     }      
#'  } 
#' }
#' @examples
#' set.seed(1)
#' x <- runif(10)
#' y <- runif(10)
#' del <- delaunay(x, y)
#' plot(del)
#' del
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delaunay <- function(x, y, calc_polygons = TRUE, calc_areas = FALSE, 
                     calc_segments = FALSE) {
  
  if (anyNA(x) || anyNA(y)) {
    stop("voronoi(): 'x' and 'y' cannot contain NA values")
  }
  
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

