
test_that("Delaunays with low point counts behave sanely", {

  set.seed(1)
  
  N <- 0
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  expect_equal(del$ntris, 0)
  
  N <- 1
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  expect_equal(del$ntris, 0)
  
  N <- 2
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  expect_equal(del$ntris, 0)
  
  N <- 3
  x <- runif(N)
  y <- runif(N)
  del <- delaunay(x, y)
  expect_equal(del$ntris, 1)
})
