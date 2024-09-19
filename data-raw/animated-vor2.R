

library(nara)

set.seed(1)
N  <- 200
w <- h <- 400
cx <- runif(N, -0.25 * w, 1.25 * h)
cy <- runif(N, -0.25 * h, 1.25 * h)

r     <- runif(N, 1, w/4)
omega <- runif(N, 0.1, 2* pi/10)
# omega <- 0.2

t <- 1

nr <- nara::nr_new(w, h)

for (t in 1:30) {
  nr_fill(nr, 'black')
  x <- cx + r * cos(omega * t)
  y <- cy + r * sin(omega * t)

  del <- delaunay(x, y, calc_areas = TRUE)
  
  frac <- del$tris$area / max(del$tris$area)
  fills <- grey(frac)
  
  with(del$polygons, 
    nr_polygon(nr, x, y, id = idx, fill = fills)
  )
  # nr_point(nr, x, y, color = 'white')
  
  plot(nr)
  Sys.sleep(0.25)
}
