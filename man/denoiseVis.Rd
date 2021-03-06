% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{denoiseVis}
\alias{denoiseVis}
\title{Remove values that are too close to the average and are considered noise.}
\usage{
denoiseVis(data = NULL, noise_threshold = NULL, vis_bounds = "-1,1")
}
\arguments{
\item{noise_threshold}{the minimum difference a value can be from the
average reference in order for it not to be
removed as noise.}

\item{vis_bounds}{used as upper and lower bounds for values in the visualization.
Should be given in the form of '-1,1' (lower bound, upper bound).}

\item{smooth_matrix}{a matrix of values, smoothed, and with average
reference removed. Row = Genes, Col = Cells.}
}
\value{
Returns an 'infercnv' object including: denoised matrix.
}
\description{
Args:
}
