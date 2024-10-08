% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-voronoi.R
\name{plot.vor}
\alias{plot.vor}
\title{Plot a Voronoi tessellation}
\usage{
\method{plot}{vor}(
  x,
  sites = TRUE,
  labels = !sites,
  verts = TRUE,
  polys = TRUE,
  fsegs = !polys,
  isegs = !polys,
  bounds = TRUE,
  site_pch = ".",
  site_cex = 1,
  site_col = "black",
  vert_pch = 19,
  vert_cex = 0.3,
  vert_col = "black",
  label_cex = 0.5,
  fseg_col = "black",
  iseg_col = "red",
  poly_col = rainbow(vor$polygons$npolygons),
  buffer = 0,
  ...
)
}
\arguments{
\item{x}{object returned by \code{voronoi()}}

\item{sites, labels, verts, polys, fsegs, isegs, bounds}{logical values. Should this geometric feature 
be plotted? Default: TRUE}

\item{site_pch, site_cex, site_col}{parameters for sites}

\item{vert_pch, vert_cex, vert_col}{parameters for Voronoi vertices}

\item{label_cex}{size of text labels for sites}

\item{fseg_col, iseg_col}{colours for finite and infinite segments}

\item{poly_col}{vector of colours for polygons}

\item{buffer}{buffer around extents. default: 0}

\item{...}{other arguments passed to \code{plot()}}
}
\value{
None.
}
\description{
Barebones plotting of Voronoi
}
