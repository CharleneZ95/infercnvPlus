% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infercnvPlus.R
\name{plotCNV}
\alias{plotCNV}
\title{Clustering and plotting cells based on cnv score matrix using ComplexHeatmap package.}
\usage{
plotCNV(
  data,
  assay,
  ref_lab,
  obs_lab,
  dist_method,
  clustering_method,
  cutree_k,
  colors,
  border,
  plot_dend,
  out_file,
  out_path
)
}
\arguments{
\item{data:}{an 'infercnv' object as produced by inferCNV.}

\item{assay:}{assay which assigned to use for plotting.}

\item{ref_lab:}{label for reference cells.}

\item{obs_lab:}{label for observations.}

\item{dist_method:}{distance measure used in clustering cells, possible
values are 'correlation' for Pearson correlation 
and all the distances supported by stats::dist.}

\item{clustering_method:}{clustering method used. Accepts the same values as stats::hclust.}

\item{cutree_k:}{an integer scalar or vector with the desired number of groups}

\item{colors:}{vector of colors used in heatmap.}

\item{border:}{whether draw border. The value can be logical or a string of color.}

\item{out_file:}{filename to save plot}

\item{out_path:}{output directory}
}
\value{
Returns an 'infercnv' object including:
    1. Clustering result.
    2. Dendlist of 'cutree'.
}
\description{
Args:
}
