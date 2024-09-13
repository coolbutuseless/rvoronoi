

test_that("multiplication works", {
  
  set.seed(1)
  N <- 1
  x <- runif(N)
  y <- runif(N)
  
  theta <- seq(0, 2*pi, length.out = 6)[-1]
  x <- cos(theta)
  y <- sin(theta)
  
  vor <- voronoi(x, y) 
  
  # for (i in 1:1) {
  #   tryCatch({
  #     vor <- voronoi(x, y)
  #   }, error = function(e) {}
  #   )
  # }
  
  expect_true(TRUE)
  
  
})


if (FALSE) {
  
  seed <- 1024738599 # N = 1000
  seed <- 469561852
  seed <- as.integer(runif(1, 1, 2^30))
  print(seed)
  set.seed(seed)
  N <- 1000
  x <- runif(N)
  y <- runif(N)
  vor <- voronoi(x, y)
  
  cols <- rainbow(length(vor$polygons))
  plot(x, y, ann = F, asp = 1, axes = F, xlim = c(0, 1), ylim = c(0, 1), pch = '.')
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]]$x, vor$polygons[[i]]$y, col = cols[i])
  }
  points(x, y, pch = '.')
  
  polygon(vor$polygons[[5]], col = 'black')
  points(x[5], y[5], col = 'green', pch = 19)
  # points(x[6], y[6], col = 'black', pch = 19)
  
  
  
  xt <- x[2]
  yt <- y[2]
  poly <- vor$polygons[[5]]
  
  point_in_convex_polygon(xt, yt, poly$x, poly$y)  
  
  
  
  plot(poly, ann = FALSE, asp = 1, axes = F, pch = '.')
  polygon(poly)
  points(x[4], y[4], col = 'red')
  points(x[5], y[5], col = 'black')
  
  text(poly$x + 0.0, poly$y + 0.0, labels = seq_along(poly$x))
  
  
  points(x, y, pch = 19)
  # points(vor$vertex, col = 'hotpink', pch = 19)
  
}

