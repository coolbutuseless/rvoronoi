
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rvoronoi

<!-- badges: start -->

![](https://img.shields.io/badge/cool-useless-green.svg)
[![R-CMD-check](https://github.com/coolbutuseless/voronoi-dev/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/coolbutuseless/voronoi-dev/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of rvoronoi is to …

## Installation

You can install from
[GitHub](https://github.com/coolbutuseless/rvoronoi) with:

``` r
# install.packages('remotes')
remotes::install_github('coolbutuseless/rvoronoi')
```

## Original Source Code

- <https://www3.cs.stonybrook.edu/~algorith/implement/fortune/implement.shtml>
- <https://netlib.sandia.gov/voronoi/sweep2>
  - wget then ‘sh’ to unpack
- Another persons go at tidying up source:
  <https://github.com/stolk/forvor>
- tech info
  <https://www.cs.tufts.edu/comp/163/notes05/voronoi_handout.pdf>
  - max vertices = 2 \* n - 5
  - max edges = 3 \* n - 6
- Delaunay:
  - In the plane (d = 2), if there are b vertices on the convex hull,
    then any triangulation of the points has at most 2n – 2 – b triangle

## Voronoi Tesselation

``` r
library(rvoronoi)

set.seed(2024)
x <- runif(10)
y <- runif(10)

vor <- voronoi(x, y) 

# Plot the seed points
plot(x, y, asp = 1, col = 'red')

# Plot all finite segments.  
# This will not plot any of the segments which do not converge
segments(
  vor$vertex$x[vor$segment$v1],
  vor$vertex$y[vor$segment$v1],
  vor$vertex$x[vor$segment$v2],
  vor$vertex$y[vor$segment$v2],
)
```

<img src="man/figures/README-voronoi-1.png" width="100%" />

## Delaunay Triangulation

``` r
library(rvoronoi)

set.seed(2024)
x <- runif(10)
y <- runif(10)

del <- delaunay(x, y) 

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

``` r
library(RTriangle)

set.seed(2024)
N <- 100
x <- runif(N)
y <- runif(N)

del_rtriangle <- RTriangle::triangulate(RTriangle::pslg(P = cbind(x, y)))$T
del_new       <- delaunay(x, y)

identical(
  rvoronoi:::normalise_del(del_new),
  rvoronoi:::normalise_del(del_rtriangle)
)
#> [1] TRUE
```

| expression    |    min | median |   itr/sec | mem_alloc |
|:--------------|-------:|-------:|----------:|----------:|
| del_rtriangle | 91.2µs |  101µs |  9667.011 |   29.45KB |
| del_new       | 19.5µs | 20.1µs | 48257.578 |    2.48KB |

## Related Software
