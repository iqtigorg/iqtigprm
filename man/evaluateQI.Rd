% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_results.R
\name{evaluateQI}
\alias{evaluateQI}
\title{Computation of one QI with attached data and LE}
\usage{
evaluateQI(
  qi,
  meta_unit = NULL,
  conf_level = 0.95,
  a = 0.5,
  b = 0.5,
  nMC = 1e+05
)
}
\arguments{
\item{qi}{The QI to compute (object of class \code{qi}).}

\item{meta_unit}{The meta_unit (i.e. LE) to compute results for. Default (\code{NULL}) means that all data are used for computation of overall results.}

\item{conf_level}{Confidence level of the resulting uncertainty interval}

\item{a}{first prior parameter of Beta distribution for each Merkmal}

\item{b}{second prior parameter of Beta distribution for each Merkmal}

\item{nMC}{Number of Monte Carlo Samples (needed for QIs with >1 Merkmale)}
}
\value{
A list containing the point estimate and the uncertainty interval.
}
\description{
Apply Bayesian computation method to one QI with an attached data set for one LE.
Each individual merkmal of the QI is computed separately. Overall QI results are computed by MC-sampling
from the Merkmal-specific posteriors. The seed for MC-sampling is fixed such that repeated
evaluation of the QI yields the same results with respect to the credibility interval.
}
\section{Fields}{

\describe{
\item{\code{QI_hat}}{Posterior mean estimate for QI value.}

\item{\code{J}}{Number of patients in the data set for the LE.}

\item{\code{interval}}{A vector of length two containing the uncertainty interval.}
}}

