% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/V_d.R
\name{V_d}
\alias{V_d}
\title{V_d}
\usage{
V_d(data, w = c(0, 1, 0), tau = 3, dec.fun, feat = NULL, SE = TRUE)
}
\arguments{
\item{data}{A data frame of the covariates, treatment, and observed outcome}

\item{w}{Weight Vector for patient preferences of state.}

\item{tau}{The total time of the trial}

\item{dec.fun}{Class of Decision Function ("linear", "rbf", or custom)}

\item{feat}{List of Column Names of Covariate}

\item{SE}{Standard Error calculation (bool)}
}
\description{
Estimator of the value function of an ITR
}
