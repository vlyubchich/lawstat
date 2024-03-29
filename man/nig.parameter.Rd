% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nig.parameter.R
\name{nig.parameter}
\alias{nig.parameter}
\title{Generate Parameters for the Normal Inverse Gaussian (NIG) Distribution}
\usage{
nig.parameter(
  mean = mean,
  variance = variance,
  kurtosis = kurtosis,
  skewness = skewness
)
}
\arguments{
\item{mean}{mean of the NIG distribution.}

\item{variance}{variance of the NIG distribution.}

\item{kurtosis}{excess kurtosis of the NIG distribution.}

\item{skewness}{skewness of the NIG distribution.}
}
\value{
A list with the following numeric components:
\item{alpha}{tail-heavyness parameter of the NIG distribution.}
\item{beta}{asymmetry parameter of the NIG distribution.}
\item{delta}{scale parameter of the NIG distribution.}
\item{mu}{location parameter of the NIG distribution.}
}
\description{
Produce four parameters, alpha (tail heavyness), beta (asymmetry), 
delta (scale), and mu (location) from the four variables: mean, variance, 
kurtosis, and skewness.
}
\details{
The parameters are generated with three conditions: 
1) \eqn{3\times kurtosis > 5\times skewness^2}; 
2) \eqn{skewness > 0}, and
3) \eqn{variance > 0}.
See \insertCite{Atkinson_1982;textual}{lawstat}, 
\insertCite{BarndorffNielsen_Blaesild_1983;textual}{lawstat}, and 
\insertCite{Noguchi_Gel_2010;textual}{lawstat}.
}
\examples{
library(fBasics)
test <- nig.parameter(0, 2, 5, 1)
random <- rnig(1000000, alpha = test$alpha, beta = test$beta, 
               mu = test$mu, delta = test$delta)
mean(random)
var(random)
kurtosis(random)
skewness(random)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[fBasics]{rnig}}
}
\author{
Kimihiro Noguchi, Yulia R. Gel
}
\keyword{distribution}
