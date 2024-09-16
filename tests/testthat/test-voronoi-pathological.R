

test_that("voronoi pathological1 works", {
  
  # 100 points in a cirlce
  
  theta <- seq(0, 2*pi, length.out = 101)[-1]
  x <- cos(theta)
  y <- sin(theta)
  
  expect_no_error({
    expect_no_warning({
      vor <- voronoi(x, y) 
    })
  })
  
})


test_that("voronoi pathological2 works", {
  
  # 100 points in circle. 
  # 1 point in centre
  
  theta <- seq(0, 2*pi, length.out = 101)[-1]
  x <- c(0, cos(theta))
  y <- c(0, sin(theta))
  
  expect_no_error({
    expect_no_warning({
      vor <- voronoi(x, y) 
    })
  })
  
})


test_that("voronoi pathological2a works", {
  
  # 2 concentric circles. 100 points each.
  # 1 point in centre
  
  theta <- seq(0, 2*pi, length.out = 101)[-1]
  x <- c(0, cos(theta), 2 * cos(theta))
  y <- c(0, sin(theta), 2 * sin(theta))
  
  vor <- voronoi(x, y)
  
  expect_no_error({
    expect_no_warning({
      vor <- voronoi(x, y) 
    })
  })
  
})


test_that("voronoi pathological3 works", {
  
  # 100 points in a line
  
  x <- seq(0, 2*pi, length.out = 100)
  y <- 0.5 * x  
  
  expect_no_error({
    expect_no_warning({
      vor <- voronoi(x, y) 
    })
  })
  
})


test_that("voronoi pathological4 works", {
  
  # * 2 concentric equalateral triangles 
  # * 1 point at the centre
  
  
  theta <- seq(0, 2*pi, length.out = 4)[-1]
  x <- c(0, cos(theta), 2 * cos(theta))
  y <- c(0, sin(theta), 2 * sin(theta))
  
  expect_no_error({
    expect_no_warning({
      vor <- voronoi(x, y) 
    })
  })
  
})


if (FALSE) {
  
  plot(vor)
}