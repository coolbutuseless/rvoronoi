

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
  
  set.seed(1)
  
  
  # macOS
  # N = 5,  set.seed(1231197)
  set.seed(1231197)
  set.seed(7355344)
  set.seed(8077981)
  set.seed(8927592)
  
  # codespace
  N <- 6
  # Seed: 136721
  # Seed: 581684
  # Seed: 704818
  # Seed: 1735603
  
  start <- Sys.time()
  for (i in seq(1, 1000000)) {
    set.seed(i)
    N <- 6
    x <- runif(N) * 10
    y <- runif(N) * 10
    
    tryCatch(
      vor <- voronoi(x, y),
      warning = \(w) { message("Seed: ", i) }
    )
  }
  Sys.time() - start
  
  
  
  
  N <- 6
  set.seed(704818)
  x <- runif(N)
  y <- runif(N)
  
  min(dist(cbind(x, y))) ^ 2
  
  vor <- voronoi(x, y, match_polygons = TRUE)
  vor <- voronoi(x, y, match_polygons = FALSE)
  
  
  vor$segment
  vor$msegments
  vor$vertex
  
  
  cols <- rainbow(length(vor$polygons))
  plot(x/2, y/2, ann = F, asp = 1, axes = F, xlim = c(0, 1), ylim = c(0, 1))
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]]$x/2, vor$polygons[[i]]$y/2, col = cols[i])
  }
  points(x/2, y/2, pch = 19)
  
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

