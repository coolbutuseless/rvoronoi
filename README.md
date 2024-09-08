
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

## Voronoi Tesselation

The following code calculates the voronoi tesselation on 20 random
points.

Only finite length segments are plotted. Segments which head to infinity
are not plotted here.

``` r
library(rvoronoi)

set.seed(2)
x <- runif(20)
y <- runif(20)

vor <- voronoi(x, y) 

# Plot the seed points
plot(x, y, asp = 1, col = 'red')

# Plot all finite segments.  
# This will not plot any of the segments which do not converge

fseg <- subset(vor$segment, v1 > 0 & v2 > 0)

segments(
  vor$vertex$x[ fseg$v1 ],
  vor$vertex$y[ fseg$v1 ],
  vor$vertex$x[ fseg$v2 ],
  vor$vertex$y[ fseg$v2 ],
)
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

Simple benchmark comparing this package with `{RTriangle}`

Benchmarking code for other packages is welcomed!

``` r
library(RTriangle)

set.seed(2024)
N <- 20
x <- runif(N)
y <- runif(N)

del_rtriangle <- RTriangle::triangulate(RTriangle::pslg(P = cbind(x, y)))$T
del_new       <- delaunay(x, y)$segment

identical(
  rvoronoi:::normalise_del(del_new),
  rvoronoi:::normalise_del(del_rtriangle)
)
#> [1] TRUE
```

| expression    |    min |  median |  itr/sec | mem_alloc |
|:--------------|-------:|--------:|---------:|----------:|
| del_rtriangle | 49.3µs | 56.64µs |  17456.8 |    5.76KB |
| del_new       |  4.8µs |  5.33µs | 170256.3 |    2.57KB |

## Debug plotting

- Used for debugging
- Not ready for general use

``` r
set.seed(3)
N <- 10
x0 <- runif(N)
y0 <- runif(N)

vor <- voronoi(x0, y0)

if (FALSE) {
  ps <- extract_polygons_r(vor)
  ps
}

plot_vor(vor) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red')
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Pathological 2a

- 2 concentric circles
- 1 point at the centre

``` r
theta <- seq(0, 2*pi, length.out = 100)[-1]
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

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

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

``` r
library(RTriangle)
theta <- seq(0, 2*pi, length.out = 6)[-1]
x <- c(0, cos(theta), 2 * cos(theta))
y <- c(0, sin(theta), 2 * sin(theta))

tri <- triangulate(pslg(P = cbind(x, y)))

v1 <- tri$VE[,1]
v2 <- tri$VE[,2]

seg <- merge_vertices(tri$VP[,1], tri$VP[,2], v1, v2)

polys <- extract_polygons(tri$VP[,1], tri$VP[,2], seg$v1, seg$v2)

cols <- rainbow(length(polys))
plot(tri$VP, asp = 1, ann = F, axes = FALSE)
for (i in seq_along(polys)) {
  polygon(polys[[i]], col = cols[i])
}
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

## Random

``` r
set.seed(1)
N <- 6
x <- runif(N)
y <- runif(N)

vor <- voronoi(x, y)

plot_vor(vor, buffer = 1) |>
  draw_segments() |> 
  draw_inf_lines() |>
  draw_inf_segments(col = 'hotpink') |>
  draw_bounded_polygons(border = 'red') |>
  draw_vertices()
points(x, y, pch = '+')
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

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
  - **TODO** Use STR to bulk load an R-tree for optimized
    point-in-bounding-box search
  - Search for bounding box collision to filter polygon candidates for
    exhaustive testing
  - Optimised point-in-polygon search with early termination as voronoi
    polygons are **always** convex. Test edges using the “leftOf()”
    operator, and exit as soon as any vertices are right of a polygon
    edge.

## Related

- [mdsumner’s helpful
  gist](https://gist.github.com/mdsumner/8db5ac01e47fa86f10e7ebc372e0ebda?permalink_comment_id=5171205)
- [`RTriangle::triangulate()`]()
- [`terra::voronoi()`](https://rspatial.github.io/terra/reference/voronoi.html)
- [`geos::geos_voronoi_polygons()`]()
