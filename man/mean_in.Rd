% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllHelperFunctions.R
\name{mean_in}
\alias{mean_in}
\title{Conditional Sampling from Fitted Poisson Gaussian Process Model}
\usage{
mean_in(SEED, coords_norm_sub, alpha.est, rho.est, mu.est, s.est)
}
\arguments{
\item{SEED}{Random seed for reproducibility}

\item{coords_norm_sub}{Normalized coordinates for spots within the domain}

\item{alpha.est}{Estimated alpha parameter}

\item{rho.est}{Estimated rho parameter}

\item{mu.est}{Estimated mean parameter}

\item{s.est}{Estimated s parameter}
}
\value{
A numeric vector representing the sampled expression levels within the domain
}
\description{
This function performs conditional sampling from a fitted Poisson Gaussian process model.
}
