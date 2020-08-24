\name{ca125}
\alias{ca125}

\docType{data}

\title{Recurrent ovarian carcinoma study}

\description{A meta-analysis of 52 studies that were reported between January 1995 and November 2007.}

\format{
  The data frame contains the following columns:
  \describe{
    \item{n}{total number of subjects}
    \item{PiY}{disease prevalence}
    \item{SeY}{true positive}
    \item{n1}{subjects with disease}
    \item{SpY}{true negative}
    \item{n1}{health individuals}  
  }
}


\references{
Chen, Y., Liu, Y., Chu, H., Lee, M. and Schmid, C. (2017) A simple and robust method for multivariate meta-analysis of diagnostic test accuracy, Statistics in Medicine, 36, 105-121.

Gu P, Pan L, Wu S, Sun L, Huang G. CA 125, PET alone, PET-CT, CT and MRI in diagnosing recurrent ovarian carcinoma: a systematic review and meta-analysis. European journal of radiology 2009; 71(1):164-174.
}


\note{
The dataset \code{ca125} is used to conduct multivariate meta-analysis of diagnostic test accuracy.
}


\seealso{
  \code{\link{mmeta}},
  \code{\link{summary.mmeta}}
}

\examples{
data(ca125)
summary(ca125)
}

\keyword{datasets}
