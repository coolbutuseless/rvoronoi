
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rvoronoi

<!-- badges: start -->

![](https://img.shields.io/badge/cool-useless-green.svg)
![](https://img.shields.io/badge/status-dangerously_unstable-red.svg)
[![R-CMD-check](https://github.com/coolbutuseless/rvoronoi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/coolbutuseless/rvoronoi/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`rvoronoi` is a testing ground for some rendering ideas using fast
delaunay triangulation and voronoi tesselation.

For small sets of seed points (e.g. 20 points) this package can be 10x
faster than `RTriangle` package. For larger sets of seed points, the gap
closes up - e.g. for N = 1000, still 2x faster than `RTriangle` package.

The core of this package is Steven Fortune’s original C source code for
his sweep algorithm. This code has been updated and adapted to run
within R. [Original source code (packaged as a shell
archive)](https://netlib.sandia.gov/voronoi/sweep2)

## Installation

You can install from
[GitHub](https://github.com/coolbutuseless/rvoronoi) with:

``` r
# install.packages('remotes')
remotes::install_github('coolbutuseless/rvoronoi')
```

# In fair Voronoi, where we lay our scene

``` r
library(nara)
library(grid)

rj <- jpeg::readJPEG("man/figures/rj.jpg", native = TRUE)
dim(rj)
#> [1]  531 1278
# grid.raster(rj)


set.seed(1)
N <- 5000
x <- (runif(N, 1, ncol(rj))) |> sort() 
y <- (runif(N, 1, nrow(rj))) 

ind <- round(y) * ncol(rj) + round(x)

cols <- rj[ind]
cols <- nara::packed_cols_to_hex_cols(cols)

vor <- voronoi(as.double(x), as.double(y))
nr <- nara::nr_new(ncol(rj), nrow(rj))

for (i in seq_along(vor$polygons)) {
  nr_polygon(nr, vor$polygons[[i]]$x, vor$polygons[[i]]$y, fill = cols[i])
}
grid.raster(nr)
```

<img src="man/figures/README-romeo-1.png" width="100%" />

## Voronoi Tesselation

The following code calculates the voronoi tesselation on 20 random
points.

``` r
library(rvoronoi)

set.seed(2)
x <- runif(20)
y <- runif(20)

vor <- voronoi(x, y) 

# Plot the seed points
plot(x, y, asp = 1, ann = FALSE, axes = FALSE, col = 'red', 
     xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))

cols <- rainbow(length(vor$polygons))
for (i in seq_along(vor$polygons)) {
  polygon(vor$polygons[[i]], col = cols[i])
}

points(x, y, pch = 19)
```

<img src="man/figures/README-voronoi-1.png" width="100%" />

## Delaunay Triangulation

``` r
library(rvoronoi)

set.seed(2024)
x <- runif(10)
y <- runif(10)

del <- delaunay(x, y)$segment

# Plot the seed points
plot(x, y, asp = 1, col = 'red')

# Plot all finite segments.  
# This will not plot any of the segments which do not converge
segments(x[del$v1], y[del$v1], x[del$v2], y[del$v2])
segments(x[del$v3], y[del$v3], x[del$v2], y[del$v2])
segments(x[del$v1], y[del$v1], x[del$v3], y[del$v3])
```

<img src="man/figures/README-delaunay-1.png" width="100%" />

## Delaunay Benchmark

Compare Delaunay Triangulation using

- `{rvoronoi}`
- `{delone}`
  - `remotes::install_github("hypertidy/delone")`
- `{RTriangle}`
- `{deldir}`

``` r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delaunay
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rvoronoi)
# library(delone)
library(RTriangle)
library(deldir)
#> deldir 2.0-4      Nickname: "Idol Comparison"
#> 
#>      The syntax of deldir() has changed since version 
#>      0.0-10.  Read the help!!!.
#> deldir 2.0-4      Nickname: "Idol Comparison"
#> 
#>      The syntax of deldir() has changed since version 
#>      0.0-10.  Read the help!!!.

set.seed(1)
N <- 1000
x <- runif(N)
y <- runif(N)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# {rvoronoi}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
del_rvoronoi <- rvoronoi::delaunay(x, y)$segment

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# {delone} - normalise to 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# del_delone <- delone::xy_tri(x, y)
# del_delone <- matrix(del_delone, ncol = 3, byrow = TRUE) |> as.data.frame()
# 
# stopifnot(identical(
#   rvoronoi:::normalise_del(del_rvoronoi),
#   rvoronoi:::normalise_del(del_delone)
# ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RTriangle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
del_rtriangle <- RTriangle::triangulate(RTriangle::pslg(cbind(x, y)))$T
del_rtriangle <- as.data.frame(del_rtriangle)

stopifnot(identical(
  rvoronoi:::normalise_del(del_rvoronoi),
  rvoronoi:::normalise_del(del_rtriangle)
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# deldir
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
del_deldir <- deldir::deldir(x, y)
del_deldir <- deldir::triMat(del_deldir)
del_deldir <- as.data.frame(del_deldir)

stopifnot(identical(
  rvoronoi:::normalise_del(del_rvoronoi),
  rvoronoi:::normalise_del(del_deldir)
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Benchmark
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bench::mark(
  # delone    = xy_tri  (x, y),
  rvoronoi  = delaunay(x, y),
  rtriangle = triangulate(pslg(cbind(x, y))),
  deldir    = deldir(x, y),
  check = FALSE
)[,1:5]  |> knitr::kable()
```

| expression |     min |  median |    itr/sec | mem_alloc |
|:-----------|--------:|--------:|-----------:|----------:|
| rvoronoi   | 393.5µs | 419.9µs | 2371.91832 |  139.56KB |
| rtriangle  | 723.2µs | 787.2µs | 1242.62215 |  292.63KB |
| deldir     |  16.8ms |  16.9ms |   58.70806 |    5.67MB |

# Voronoi Tesselation Benchmark

Compare:

- `{rvoronoi}`
- `{deldir}`
- `{RTriangle}`

``` r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Voronoi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rvoronoi)
library(deldir)
library(RTriangle)

set.seed(1)
N <- 1000
x <- runif(N)
y <- runif(N)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rvoronoi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vor_rvoronoi <- rvoronoi::voronoi(x, y)

# plot(x, y, asp = 1, ann = F, axes = F, pch = '.')
# for (p in vor_rvoronoi$polygons) {
#   polygon(p)
# }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# deldir
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vor_deldir   <- deldir::cvt(deldir::deldir(x, y), stopcrit = 'maxit', maxit = 1)

# plot(x, y, asp = 1, ann = F, axes = F, pch = '.')
# for (p in vor_deldir$tiles) {
#   polygon(p)
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RTriangle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tri <- triangulate(pslg(P = cbind(x, y)))
ve <- tri$VE |> as.data.frame()
ve <- subset(ve, V1 > 0 & V2 > 0)

# plot(x, y, asp = 1, ann = F, axes = F, pch = '.')
# p1 <- tri$VP[ve$V1, ] |> as.data.frame()
# p2 <- tri$VP[ve$V2, ] |> as.data.frame()
# segments(p1$V1, p1$V2, p2$V1, p2$V2)


bench::mark(
  voronoi(x, y),
  cvt(deldir(x, y), stopcrit = 'maxit', maxit = 1),
  triangulate(pslg(P = cbind(x, y))),
  check = FALSE
)[,1:5] |> knitr::kable()
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
```

| expression | min | median | itr/sec | mem_alloc |
|:---|---:|---:|---:|---:|
| voronoi(x, y) | 3.38ms | 3.45ms | 287.813316 | 235KB |
| cvt(deldir(x, y), stopcrit = “maxit”, maxit = 1) | 157.6ms | 160.67ms | 6.113305 | 52.1MB |
| triangulate(pslg(P = cbind(x, y))) | 728.2µs | 758.5µs | 1229.584350 | 292.6KB |

# Pathological Test Cases

## Pathological 1

- 100 points in a circle

``` r
theta <- seq(0, 2*pi, length.out = 101)[-1]
x <- cos(theta)
y <- sin(theta)

vor <- voronoi(x, y)

plot_vor(vor) |>
  draw_segments() |>
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red')
points(x, y)
```

<img src="man/figures/README-path1-1.png" width="100%" />

## Pathological 2

- 100 points in a circle
- 1 point at the centre

``` r
theta <- seq(0, 2*pi, length.out = 101)[-1]
x <- c(0, cos(theta))
y <- c(0, sin(theta))

vor <- voronoi(x, y)

plot_vor(vor) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red')
points(x, y, pch = '+')
```

<img src="man/figures/README-path2-1.png" width="100%" />

## Pathological 2a

- 2 concentric circles (100 points each)
- 1 point at the centre

``` r
theta <- seq(0, 2*pi, length.out = 101)[-1]
x <- c(0, cos(theta), 2 * cos(theta))
y <- c(0, sin(theta), 2 * sin(theta))

vor <- voronoi(x, y)

plot_vor(vor) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red')
points(x, y, pch = '+')
```

<img src="man/figures/README-path2a-1.png" width="100%" />

## Pathological 3

- 100 points in a line

``` r
x <- seq(0, 2*pi, length.out = 100)
y <- 0.5 * x  

vor <- voronoi(x, y)

plot_vor(vor) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red')
points(x, y, pch = '+')
```

<img src="man/figures/README-path3-1.png" width="100%" />

## Pathological 4

- 2 concentric cirles (50 points each)
- 1 point at the centre

``` r
theta <- seq(0, 2*pi, length.out = 4)[-1]
x <- c(0, cos(theta), 2 * cos(theta))
y <- c(0, sin(theta), 2 * sin(theta))

vor <- voronoi(x, y)

plot_vor(vor) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red') |>
  draw_vertices()
points(x, y, pch = '+')
```

<img src="man/figures/README-path4-1.png" width="100%" />

``` r

if (FALSE) {
  plot(x, y, ann = F, axes = F, asp = 1)
  
  cols <- rainbow(length(vor$polygons))
  
  for (i in seq_along(vor$polygons)) {
  polygon(vor$polygons[[i]], col = cols[i])  
  text(x[i], y[i], labels = i)
  # points(x[i], y[i], col = 'black', pch = '+')
  }
  
}
```

# Algorithms

- Fortune’s Sweep Algorithm for Delaunay
  - Adapted original source code to call from R
- Polygon reconstruction
  - Use voronoi vertices and edge connectivity to reconstruct polygonal
    voronoi tiles
  - Using “An optimal algorithm for extracting the regions of a plane
    graph” by Jiang & Bunke, Pattern Recognition Letters 14 (1993)
    pp553-558.
  - Adapted algorithm from binary search to instead use pre-indexed
    search bounds.
- Point in polygon (to match seed points to voronoi polygons)
  - Search for bounding box collision to filter polygon candidates for
    exhaustive testing
  - Optimised point-in-polygon search with early termination possible
    (as voronoi polygons are **always** convex) Test edges using the
    “leftOf()” operator, and exit as soon as any vertices are right of a
    polygon edge.
  - No benefit in building an acceleration structure for
    point-in-bounding-box when there of the order of ~100-1000 bounding
    boxes. A linear search (with flagging of polygons already claimed)
    is fast enough. If use case if voronoi of \>\> 1000 points, then
    time taken to build acceleration structure (e.g. hierachical
    bounding boxes) would become worth it.

## Related

- [mdsumner’s helpful
  gist](https://gist.github.com/mdsumner/8db5ac01e47fa86f10e7ebc372e0ebda?permalink_comment_id=5171205)
- [`RTriangle::triangulate()`]()
- [`terra::voronoi()`](https://rspatial.github.io/terra/reference/voronoi.html)
- [`geos::geos_voronoi_polygons()`]()
