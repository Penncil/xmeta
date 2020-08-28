\name{summary.mmeta}
\alias{summary.mmeta}
\title{Summarize the objects \code{mmeta}}
\description{
    Summarize a model of class \code{mmeta} fitted by \code{mmeta}.
}
\usage{
      \method{summary}{mmeta}(object,...) 
}

\arguments{  
    \item{object}{an object inheriting from class \code{mmeta}.}
	 \item{...}{ additional arguments; currently none is used.}
}


\value{
  A list with the following components: coefficients, covariance matrix.
}

\references{
Chen, Y., Hong, C. and Riley, R. D. (2015). An alternative pseudolikelihood method for multivariate random-effects meta-analysis. Statistics in medicine, 34(3), 361-380.

Chen, Y., Hong, C., Ning, Y. and Su, X. (2015). Meta-analysis of studies with bivariate binary outcomes: a marginal beta-binomial model approach, Statistics in Medicine (in press).

Hong, C., Riley, R. D. and Chen, Y. (2015). An improved method for multivariate random-effects meta-analysis (in preparation).

Chen, Y., Liu, Y., Ning, J., Nie, L., Zhu, H. and Chu, H. (2014). A composite likelihood method for bivariate meta-analysis in diagnostic systematic reviews. Statistical methods in medical research (in press).

Chen, Y., Cai, Y., Hong, C. and Jackson, D. (2015). Inference for correlated effect sizes using multiple univariate meta-analyses, Statistics in Medicine (provisional acceptance).

Chen, Y., Liu, Y., Ning, J., Cormier J. and Chu H. (2014). A hybrid model for combining case-control and cohort studies in systematic reviews of diagnostic tests, Journal
of the Royal Statistical Society: Series C (Applied Statistics) 64.3 (2015): 469-489.

Chen, Y., Liu, Y., Chu, H., Lee, M. and Schmid C. (2015). A simple and robust method for multivariate meta-analysis of diagnostic test accuracy, Statistics in Medicine (under revision). 
}



\seealso{\code{\link{mmeta}}}
\examples{
data(prostate)
fit.nn=mmeta(data=prostate, type="continuous", k=2, method="nn.cl") 
summary(fit.nn)
}


\keyword{summary}
