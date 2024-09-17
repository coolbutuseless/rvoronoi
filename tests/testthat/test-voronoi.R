


test_that("Voronoi works", {
  vor <- voronoi(tpoints$x, tpoints$y)

  expect_identical(
    names(vor), 
    c("sites", "vertices", "segments", "polygons", "lines", "extents", "mvertices", "msegments")
  )

  expect_equal(nrow(vor$vertices), 187)
  expect_equal(nrow(vor$lines   ), 286)
  expect_equal(nrow(vor$segments), 286)
})

