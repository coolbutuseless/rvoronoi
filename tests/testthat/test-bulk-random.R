

# Just a bunch of random tests to stress test things a little

test_that("random point exploration", {
  
  set.seed(1)

  for (i in seq(100)) {
    
    N <- runif(1, min = 1, max = 1000)
    x <- runif(N)
    y <- runif(N)
    
    expect_no_error(
      delaunay(x, y, calc_polygons = TRUE, calc_areas = TRUE, calc_segments = TRUE)
    )
    
    expect_no_error(
      voronoi(x, y)
    )
    
  }  

})
