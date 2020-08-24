\name{msset}
\alias{msset}
\title{Testing and correcting for small study effects of multivariate meta-analysis}
\description{
Testing and correcting for small study effects of multivariate meta-analysis}
\usage{
msset(data, nm.y1, nm.s1, nm.y2, nm.s2, method, type, k)
}
\arguments{
  \item{data}{dataset}
  \item{nm.y1}{column name for outcome 1}
  \item{nm.s1}{column name for standard error of outcome 1}
  \item{nm.y2}{column name for outcome 2}
  \item{nm.s2}{column name for standard error of outcome 2}
  \item{method}{"nn.cl" indicating the score test for detecting small study effects of MMA}
  \item{type}{either "continuous" or "binary" indicating the type of outcomes}
  \item{k}{integer indicating the number of outcomes}
}
\details{
This function returns the test statistics for testing small study effects of multivariate
meta-analysis using regression method.} 
\section{A score test for detecting small study effects in multivariate meta-analysis}{
Small study effects occur when smaller studies show different, often larger, treatment effects than large ones, which may threaten the validity of systematic reviews and meta-analyses. Detecting small study effects in a multivariate meta-analysis setting remains an untouched research area. Hong et al. (2019) propose a pseudolikelihood-based score test for detecting small study effects in multivariate random-effects meta-analysis. 
This is the first test for detecting small study effects in multivariate meta-analysis setting.
}

\value{
\code{msset.TS} returns the test statistic and p value of the score test.
}
\references{Hong, C., Salanti, G., Morton, S., Riley, R., Chu, H., Kimmel, S.E. and Chen Y. (2019). Testing small study effects in multivariate meta-analysis (Biometrics).
}
\author{Chuan Hong}

\examples{
data(prostate)
fit.msset=msset(data=prostate, nm.y1="y1", nm.s1="s1", nm.y2="y2", nm.s2="s2", 
method = "nn.cl", type = "continuous", k=2)
summary(fit.msset)
}
\keyword{Multivariate meta-analysis}
\keyword{Small study effects}
\keyword{Score test}
