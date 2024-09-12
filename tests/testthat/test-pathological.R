

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
  
  for (i in seq(100000)) {
    set.seed(i)
    N <- 4
    x <- runif(N)
    y <- runif(N)
  
    tryCatch(
      vor <- voronoi(x, y),
      warning = \(w) { message("Seed: ", i) }
    )
  }
  
  
  set.seed(327)
  N <- 400
  x <- runif(N)
  y <- runif(N)
  vor <- voronoi(x, y) 
  vor$segment
  vor$msegments
  vor$vertex
  
  
  cols <- rainbow(length(vor$polygons))
  plot(x, y, ann = F, asp = 1, axes = F, xlim = c(0, 1), ylim = c(0, 1))
  for (i in seq_along(vor$polygons)) {
    polygon(vor$polygons[[i]], col = cols[i])
  }
  points(x, y, pch = .)
  # points(vor$vertex, col = 'hotpink', pch = 19)
  
}


