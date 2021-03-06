% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Usc.R
\name{Usc}
\alias{Usc}
\title{U-scores for (non-interval, sinle-column) Censored Data}
\usage{
Usc(y, ind, rnk = TRUE)
}
\arguments{
\item{y}{The column of data values plus detection limits}

\item{ind}{The column of indicators, where 1 (or \code{TRUE}) indicates a detection limit in the \code{y} column, and 0 (or \code{FALSE}) indicates a detected value in \code{y}.}

\item{rnk}{A \code{TRUE}/\code{FALSE} variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=\code{TRUE}, use the ranks. rnk = \code{FALSE} returns the uscores.}
}
\value{
Returns a single column of uscores or the ranks of uscores for a single pair of (concentration, indicator) censored data columns.
}
\description{
Computes the column of uscores from 2 columns of data in the indicator value format. Multiple detection limits allowed.  Called by the uscores function, Usc (this function) is not expected to be of much use to users on its own.
}
\examples{
data(Brumbaugh)
uscore(Brumbaugh$Hg,Brumbaugh$HgCen)
}
