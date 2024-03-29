% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nloglike.R
\name{internal_functions}
\alias{internal_functions}
\alias{nlogLike}
\alias{w2n}
\title{Internal function to compute the negative log likelihood}
\usage{
nlogLike(parvect, X.list, N, p, sld, tad, stationary)

w2n(parvect, N, p, sld, tad, stationary = FALSE)
}
\arguments{
\item{parvect}{Vector with parameters}

\item{X.list}{The data}

\item{N}{The number of sates}

\item{p}{The number of covariates}

\item{sld}{The step-length distribution}

\item{tad}{The turn-angle distribution}

\item{stationary}{Is Markov Chain assumed to be stationary?}
}
\description{
Internal function to compute the negative log likelihood

Transform parameters from the working scale to the natural scale
}
