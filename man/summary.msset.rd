\name{summary.msset}
\alias{summary.msset}
\title{Summarize the objects \code{msset}}
\description{
    Summarize a model of class \code{msset} fitted by \code{msset}.
}
\usage{
      \method{summary}{msset}(object,...) 
}

\arguments{  
    \item{object}{an object inheriting from class \code{msset}.}
	 \item{...}{ additional arguments; currently none is used.}
}


\value{
  A list with the following components: test statistics (msset) and p-value.
}

\references{Hong, C., Salanti, G., Morton, S., Riley, R., Chu, H., Kimmel, S.E. and Chen Y. (2019). Testing small study effects in multivariate meta-analysis (Biometrics).
}



\seealso{\code{\link{msset}}}
\examples{
data(prostate)
fit.msset=msset(data=prostate, nm.y1="y1", nm.s1="s1", nm.y2="y2", nm.s2="s2", 
method = "nn.cl", type = "continuous", k=2)
summary(fit.msset)
}


\keyword{summary}
