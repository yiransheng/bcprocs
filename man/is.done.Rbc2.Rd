% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/2.2_Rbc2.R
\name{is.done.Rbc2}
\alias{is.done.Rbc2}
\title{Check if Rbc2 procedure has completed}
\usage{
\method{is.done}{Rbc2}(x)
}
\arguments{
\item{x}{procedure object of class Rbc2}
}
\value{
Boolean
}
\description{
Check if Rbc2 procedure has completed
}
\examples{
proc <- Rbc2(2, 56, 8)
is.done(proc)
# FALSE
}

