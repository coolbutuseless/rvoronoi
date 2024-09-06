


test_that("Voronoi works", {
  vor <- voronoi(tpoints$x, tpoints$y)
  
  expect_identical(names(vor), c("vertex", "line", "segment", "extents", "polygons", "msegments"))
  
  expect_equal(nrow(vor$vertex ), 187)
  expect_equal(nrow(vor$line   ), 286)
  expect_equal(nrow(vor$segment), 286)
})

