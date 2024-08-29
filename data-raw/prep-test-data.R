
tpoints <- read.delim("data-raw/test-points.txt", sep = " ", header = FALSE) |>
  setNames(c('x', 'y'))


tdelaunay <- read.delim("data-raw/test-output-delaunay.txt", sep = " ", header = FALSE) |>
  setNames(c('v1', 'v2', 'v3'))

tdelaunay$v1 <- tdelaunay$v1 + 1L
tdelaunay$v2 <- tdelaunay$v2 + 1L
tdelaunay$v3 <- tdelaunay$v3 + 1L



tvoronoi <- readLines("data-raw/test-output-voronoi.txt")
# s a b
#   indicates that an input point at coordinates (a, b) was seen.
#
# l a b c
#   indicates a line with equation ax + by = c.
#
# v a b
#   indicates a vertex at (a b).
#
# e l v1 v2
#   indicates a Voronoi segment which is a subsegment of line number l
#   with endpoints numbered v1 and v2
#   If v1 or v2 is '-1' the line extends to infinity


parse_s <- function(txt) {
  res <- scan(textConnection(txt), what = list(character(), double(), double()), quiet = TRUE)
  data.frame(x = res[[2]], y = res[[3]])
}

parse_l <- function(txt) {
  res <- scan(textConnection(txt), what = list(character(), double(), double(), double()), quiet = TRUE)
  data.frame(a = res[[2]], b = res[[3]], c = res[[4]])
}

parse_v <- function(txt) {
  res <- scan(textConnection(txt), what = list(character(), double(), double()), quiet = TRUE)
  data.frame(x = res[[2]], y = res[[3]])
}

parse_e <- function(txt) {
  res <- scan(textConnection(txt), what = list(character(), integer(), integer(), integer()), quiet = TRUE)
  data.frame(
    line = res[[2]] + 1L, 
    v1   = ifelse(res[[3]] < 0, res[[3]], res[[3]] + 1L), 
    v2   = ifelse(res[[4]] < 0, res[[4]], res[[4]] + 1L)
  )
}


sites <- grep('^s', tvoronoi, value = TRUE) |>
  lapply(parse_s) |>
  do.call(rbind, args = _)

lines <- grep('^l', tvoronoi, value = TRUE) |>
  lapply(parse_l) |>
  do.call(rbind, args = _)

vertices <- grep('^v', tvoronoi, value = TRUE) |>
  lapply(parse_v) |>
  do.call(rbind, args = _)

segments <- grep('^e', tvoronoi, value = TRUE) |>
  lapply(parse_e) |>
  do.call(rbind, args = _)


tvoronoi <- list(
  site    = sites,
  line    = lines,
  vertex  = vertices,
  segment = segments
)






usethis::use_data(tpoints, tdelaunay, tvoronoi, overwrite = TRUE, internal = TRUE)
