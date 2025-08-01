% Generated by roxygen2: do not edit by hand
% Please edit documentation in _main.Rmd.
\name{rand}
\alias{rand}
\title{Rand index between responsibilities.}
\usage{
rand(resp_list1, resp_list2, times = NULL, smaller = TRUE, prop = 0.1)
}
\arguments{
\item{resp_list1}{One list.}

\item{resp_list2}{Another list.}

\item{times}{A subset of the time points (out of 1:length(resp_list1)) to
examine.}

\item{smaller}{If TRUE, use only a small (sampled) subset of the particles'
responsibilities for calculating RAND.}

\item{prop}{How much to downsample; defaults to 0.1.}
}
\description{
Rand index between responsibilities.
}