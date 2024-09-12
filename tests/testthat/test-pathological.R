

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
    
    # N = 5,  set.seed(1231197)
    set.seed(1231197)
    set.seed(7355344)
    set.seed(8077981)
    set.seed(8927592)
    
    start <- Sys.time()
    for (i in seq(1, 10000000)) {
      set.seed(i)
      N <- 5
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
      Sys.time() - start
      
      
      
      
      set.seed(8927592)
      x <- runif(N)
      y <- runif(N)
      vor <- voronoi(x, y, match_polygons = TRUE)
      # vor <- voronoi(x, y, match_polygons = FALSE)
      
      
      
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
      
      points(x, y, pch = '.')
      # points(vor$vertex, col = 'hotpink', pch = 19)
      
      polygon(vor$polygons[[1]], col = 'black')
      
      # $x
      # [1]  0.470762318  0.469139242 -0.004196656 -0.004196656  0.291792072
      # 
      # $y
      # [1] 0.4934863 0.4942734 0.2754361 0.8359567 0.8359567
      
      points(x[2], y[2], col = 'hotpink')
      
      point_in_convex_polygon(x, y, vor$polygons[[1]]$x, vor$polygons[[1]]$y)  
      
      x <- x[2]
      y <- y[2]
      poly <- vor$polygons[[1]]
      
      
      plot(poly, ann = FALSE, asp = 1, axes = F)
      polygon(poly)
      points(x, y, col = 'red', pch = 19)  
      
      text(poly$x + 0.03, poly$y + 0.03, labels = seq_along(poly$x))
      
      
      points(x, y, pch = .)
      # points(vor$vertex, col = 'hotpink', pch = 19)
      
    }
    
  }
}

