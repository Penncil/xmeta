\name{sim_dat}
\alias{sim_dat}


\docType{data}

\title{Simulated data}

\description{A simulated dataset for galaxy function}

\format{
  The data frame contains the following columns:
  \describe{
    \item{study.id}{study id}
    \item{y1}{effect size for the first outcome}
    \item{s1}{within-study standard error for the first outcome}
    \item{y2}{effect size for the second outcome}
    \item{s2}{within-study stamdard error for the second outcome}
    \item{corr}{within-study correlation}
    \item{subgroup}{subgroup of the studies}
  }
}


\note{
The dataset \code{sim_dat} is used to illustrate the galaxy plot.
}


\seealso{
  \code{\link{galaxy}}
}

\examples{
data(sim_dat)
summary(sim_dat)
}

\keyword{datasets}
