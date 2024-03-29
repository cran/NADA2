% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cenQQ.R
\name{cenQQ}
\alias{cenQQ}
\title{Q-Q Plot censored data}
\usage{
cenQQ(x.var, cens.var, dist = "lnorm", Yname = yname)
}
\arguments{
\item{x.var}{The column of \code{x} (response variable) values plus detection limits}

\item{cens.var}{The column of indicators, where 1 (or \code{TRUE}) indicates a detection limit in the \code{y.var} column, and 0 (or \code{FALSE}) indicates a detected value in \code{y.var}.}

\item{dist}{One of three distributional shapes to fit to your data:  lognormal (\code{lnorm}), normal (\code{norm}) or gamma (\code{gamma}).}

\item{Yname}{Optional – input text in quotes to be used as the variable name on the Q-Q plot.  The default is the \code{Yname} name of the \code{y.var} input variable.}
}
\value{
A single Q-Q plot of data fitted by normal, lognormal or gamma distributions with Shapiro-Francia W value printed on plot.
}
\description{
Plots a quantile-quantile (Q-Q) plot of censored data versus a fitted data distribution
}
\examples{
\donttest{
data(Brumbaugh)
cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen)

# User defined distribution
cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen,dist="gamma")
}
}
\references{
Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.

Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
}
