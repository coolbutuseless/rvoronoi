

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Point in convex polygon?
#' 
#' @param x,y point
#' @param xp,yp polygon vertices
#' @return logical
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
point_in_polygon <- function(x, y, xp, yp) {
  .Call(point_in_polygon_, x, y, xp, yp)
}


if (FALSE) {
  
  n  <- 4
  px <- c(0, 1, 1, 0)
  py <- c(0, 0, 1, 1)
  
  plot(px, py, asp = 1, ann = F, axes = F)
  polygon(px, py)
  
  
  test <- function(x, y) {
    v <- numeric(4)
    for (i in seq(1, n - 1)) {
      v[i] <- (px[i + 1] - px[i]) * (y - py[i]) - (py[i + 1] - py[i]) * (x - px[i])
    }
    
    
    v[n] <- (px[n] - px[1]) * (y - py[n]) - (py[1] - py[n]) * (x - px[n])
    
    length(unique(sign(v))) == 1
  }
  
  N  <- 1000
  xs <- runif(N, -0.5, 1.5)
  ys <- runif(N, -0.5, 1.5)
  
  inside  <- mapply(test, x= xs, y = ys)
  inside2 <- mapply(point_in_polygon, x = xs, y = ys, MoreArgs = list(xp = px, yp = py)) 
  table(inside)
  table(inside2)
  
  points(xs, ys, col = c('red', 'blue')[inside2 + 1])
  
}