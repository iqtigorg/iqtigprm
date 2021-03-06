% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_results.R
\name{theta_bayes}
\alias{theta_bayes}
\title{Bayesian inference for theta for one Merkmal}
\usage{
theta_bayes(y, nClass, a, b, conf_level)
}
\arguments{
\item{y}{Observed input of all Items of the Merkmal with points}

\item{nClass}{Number of Item categories (supposed do be the same for each Item)}

\item{a}{Prior shape parameter of Beta distribution}

\item{b}{Prior scale parameter of Beta distribution}

\item{conf_level}{Confidence level of the resulting uncertainty interval.}
}
\value{
A \code{list} containing the prior and posterior parameters of the Beta distribution,
the corresponding posterior mean value for theta and a credibility interval.
}
\description{
This is a helper function for the Bayesian QI computation function.
The function computes the parameters of the Beta posterior for the underlying
merkmals parameter theta as well as a credibility interval for theta.
}
