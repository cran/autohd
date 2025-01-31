% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdcoxma.R
\name{hdcoxma}
\alias{hdcoxma}
\title{High dimensional multivariate cox proportional model with bayesian
       mediation analysis.}
\usage{
hdcoxma(m, n, survdur, event, ths, b, d, data)
}
\arguments{
\item{m}{Starting column number from where high dimensional variates to be selected.}

\item{n}{Ending column number till where high dimensional variates to be selected.}

\item{survdur}{"Column/Variable name" consisting duration of survival.}

\item{event}{"Column/Variable name" consisting survival event.}

\item{ths}{A numeric between 0 to 100.}

\item{b}{Number of MCMC iterations to burn.}

\item{d}{Number of draws for the iterations.}

\item{data}{High dimensional data containing survival observations and high dimensional covariates.}
}
\value{
Data frame containing the beta and alpha values of active variables among the significant variables.
}
\description{
Given the dimension of variables and survival information the function filters significant variables
by fitting multivariate Cox PH with 5 variables at a time. Further, it performs mediation analysis among the significant
variables and provides handful variables with their alpha.a values which are mediator model exposure coefficients
and beta.a coefficients.
}
\examples{
hdcoxma(m=8,n=106,survdur="os",event="death",ths=0.02,b=10,d=10,data=hnscc2)
}
