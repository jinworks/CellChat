% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{runPCA}
\alias{runPCA}
\title{Dimension reduction using PCA}
\usage{
runPCA(data.use, do.fast = T, dimPC = 50, seed.use = 42, weight.by.var = T)
}
\arguments{
\item{data.use}{input data (samples in rows, features in columns)}

\item{do.fast}{whether do fast PCA}

\item{dimPC}{the number of components to keep}

\item{seed.use}{set a seed}

\item{weight.by.var}{whether use weighted pc.scores}
}
\description{
Dimension reduction using PCA
}
