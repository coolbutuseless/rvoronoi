
test_that("Voronois with low point counts behave sanely", {

  set.seed(1)
  
  N <- 0
  x <- runif(N)
  y <- runif(N)
  expect_no_error({
    vor <- voronoi(x, y)
  })
  
  N <- 1
  x <- runif(N)
  y <- runif(N)
  expect_no_error({
    vor <- voronoi(x, y)
  })
  
  N <- 2
  x <- runif(N)
  y <- runif(N)
  expect_no_error({
    vor <- voronoi(x, y)
  })
  
  N <- 3
  x <- runif(N)
  y <- runif(N)
  expect_no_error({
    vor <- voronoi(x, y)
  })
  
  expect_true(TRUE)
})
