

delaunay_ref <- structure(
  list(
    v1 = c(32L, 59L, 50L, 28L, 13L, 59L, 94L, 31L, 
           28L, 53L, 94L, 50L, 69L, 59L, 53L, 13L, 9L, 71L, 46L, 23L, 69L, 
           78L, 26L, 100L, 78L, 23L, 95L, 74L, 78L, 81L, 100L, 30L, 27L, 
           9L, 100L, 36L, 48L, 36L, 41L, 81L, 39L, 48L, 41L, 52L, 79L, 39L, 
           3L, 52L, 72L, 95L, 27L, 79L, 48L, 41L, 82L, 41L, 88L, 56L, 36L, 
           48L, 99L, 39L, 22L, 95L, 24L, 98L, 61L, 61L, 53L, 61L, 21L, 73L, 
           86L, 49L, 57L, 73L, 61L, 3L, 21L, 97L, 80L, 88L, 10L, 75L, 20L, 
           75L, 25L, 61L, 11L, 73L, 22L, 57L, 11L, 2L, 4L, 2L, 21L, 33L, 
           80L, 57L, 10L, 4L, 45L, 20L, 8L, 87L, 99L, 14L, 76L, 14L, 16L, 
           10L, 33L, 5L, 42L, 67L, 62L, 45L, 65L, 68L, 2L, 45L, 66L, 42L, 
           65L, 89L, 14L, 16L, 6L, 33L, 47L, 38L, 54L, 51L, 92L, 54L, 55L, 
           89L, 85L, 5L, 15L, 33L, 55L, 12L, 51L, 92L, 38L, 85L, 55L, 85L, 
           92L, 63L, 70L, 29L, 64L, 63L, 54L, 34L, 40L, 51L, 64L, 12L, 1L, 
           1L, 35L, 96L, 58L, 91L, 58L, 7L, 35L, 40L, 93L, 53L, 51L, 44L, 
           44L, 93L, 7L, 44L, 93L, 99L, 96L, 4L, 18L, 44L, 99L), 
    v2 = c(59L, 
           28L, 59L, 90L, 59L, 43L, 28L, 43L, 37L, 50L, 37L, 46L, 43L, 69L, 
           46L, 71L, 46L, 69L, 23L, 71L, 78L, 94L, 37L, 69L, 26L, 95L, 100L, 
           37L, 30L, 78L, 81L, 74L, 95L, 23L, 79L, 74L, 23L, 37L, 74L, 39L, 
           30L, 27L, 72L, 39L, 52L, 41L, 52L, 21L, 77L, 79L, 88L, 97L, 88L, 
           82L, 77L, 83L, 24L, 83L, 90L, 25L, 25L, 56L, 97L, 22L, 27L, 25L, 
           27L, 95L, 48L, 22L, 56L, 83L, 82L, 57L, 25L, 86L, 8L, 21L, 73L, 
           21L, 21L, 61L, 86L, 82L, 61L, 77L, 20L, 66L, 20L, 67L, 80L, 11L, 
           2L, 87L, 57L, 17L, 33L, 67L, 33L, 76L, 75L, 45L, 76L, 66L, 80L, 
           66L, 49L, 87L, 2L, 66L, 77L, 5L, 68L, 16L, 45L, 5L, 65L, 62L, 
           2L, 5L, 14L, 6L, 80L, 6L, 14L, 65L, 54L, 36L, 89L, 47L, 38L, 
           5L, 80L, 6L, 14L, 55L, 33L, 92L, 89L, 18L, 12L, 63L, 63L, 63L, 
           85L, 54L, 18L, 92L, 34L, 40L, 58L, 38L, 40L, 18L, 58L, 29L, 34L, 
           12L, 64L, 70L, 60L, 7L, 12L, 7L, 60L, 7L, 34L, 18L, 1L, 29L, 
           44L, 35L, 40L, 99L, 93L, 58L, 1L, 35L, 91L, 96L, 44L, 4L, 91L, 
           51L, 36L, 91L, 51L), 
    v3 = c(84L, 84L, 32L, 84L, 50L, 28L, 43L, 
           59L, 90L, 32L, 28L, 13L, 31L, 31L, 50L, 59L, 53L, 59L, 13L, 13L, 
           43L, 43L, 94L, 71L, 94L, 71L, 71L, 26L, 26L, 69L, 69L, 26L, 23L, 
           46L, 81L, 72L, 9L, 74L, 30L, 78L, 78L, 23L, 74L, 81L, 81L, 30L, 
           79L, 39L, 36L, 100L, 19L, 3L, 27L, 72L, 72L, 82L, 19L, 41L, 37L, 
           88L, 48L, 41L, 79L, 79L, 19L, 99L, 24L, 27L, 9L, 95L, 39L, 56L, 
           83L, 98L, 98L, 83L, 22L, 52L, 56L, 3L, 97L, 24L, 73L, 86L, 88L, 
           82L, 88L, 8L, 25L, 10L, 97L, 25L, 20L, 20L, 49L, 87L, 73L, 73L, 
           21L, 11L, 86L, 57L, 57L, 61L, 22L, 20L, 98L, 17L, 11L, 87L, 75L, 
           75L, 67L, 75L, 4L, 10L, 76L, 76L, 76L, 67L, 17L, 62L, 8L, 45L, 
           2L, 62L, 66L, 77L, 62L, 68L, 68L, 68L, 66L, 42L, 65L, 80L, 80L, 
           65L, 6L, 16L, 55L, 47L, 33L, 55L, 6L, 14L, 5L, 89L, 15L, 92L, 
           54L, 47L, 85L, 38L, 92L, 38L, 55L, 15L, 92L, 85L, 58L, 63L, 34L, 
           12L, 64L, 1L, 54L, 29L, 34L, 63L, 60L, 64L, 70L, 48L, 70L, 60L, 
           58L, 40L, 29L, 1L, 35L, 49L, 7L, 42L, 16L, 96L, 4L)), 
  row.names = c(NA, -187L), class = c("tbl_df", "tbl", "data.frame"))


test_that("Delaunay works", {
  
  del <- delaunay(tpoints$x, tpoints$y) 
  
  del <- rvoronoi:::normalise_del(del)
  ref <- rvoronoi:::normalise_del(delaunay_ref)
  
  expect_identical(del, ref)
})


test_that("Voronoi works", {
  vor <- voronoi(tpoints$x, tpoints$y)
  
  expect_identical(names(vor), c("vertex", "line", "segment", "extents"))
  
  expect_equal(nrow(vor$vertex ), 187)
  expect_equal(nrow(vor$line   ), 286)
  expect_equal(nrow(vor$segment), 286)
})

