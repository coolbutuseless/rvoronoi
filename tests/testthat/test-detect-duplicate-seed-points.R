

test_that("duplicate seeds points are detected", {

  x <- c(0, 1, 2, 3, 3)
  y <- c(0, 1, 2, 3, 3)
  
  expect_error(
    voronoi(x, y),
    "duplicate"
  )
  
})
