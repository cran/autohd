% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hidimcox.R
\name{hidimcox}
\alias{hidimcox}
\title{High dimensional univariate cox proportional hazard analysis.}
\usage{
hidimcox(m, n, survdur, event, sig, data)
}
\arguments{
\item{m}{Starting column number form where study variables of high dimensional data will get selected.}

\item{n}{Ending column number till where study variables of high dimensional data will get selected.}

\item{survdur}{Column name of survival duration event, a string value. i.e. "os"}

\item{event}{Column name of survival event, a string value. i.e "death"}

\item{sig}{Level of significance pre-determined by the user}

\item{data}{High dimensional data containing the survival, progression and genomic observations.}
}
\value{
Data set containing the list of selected variables with HR, LCL,UCL and Pvalues through survival analysis.
}
\description{
Given the dimension of variables and survival information the function performs univariate Cox PH.
}
\examples{
##
data(hnscc2)
hidimcox(m=8,n=50,survdur="os",event="death",sig=0.05,data=hnscc2)
##
}
