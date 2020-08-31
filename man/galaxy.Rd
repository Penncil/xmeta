\name{galaxy}
\alias{galaxy}
\title{Galaxy Plot: A New Visualization Tool of Bivariate Meta-Analysis Studies}
\description{
A new visualization method that simultaneously presents the effect sizes of bivariate outcomes and their standard errors in a two-dimensional space.}
\usage{
galaxy(data, nm.y1, nm.s1, nm.y2, nm.s2, scale)
}
\arguments{
  \item{data}{dataset}
  \item{nm.y1}{column name for outcome 1}
  \item{nm.s1}{column name for standard error of outcome 1}
  \item{nm.y2}{column name for outcome 2}
  \item{nm.s2}{column name for standard error of outcome 2}
  \item{scale}{a pre-specified parameter to adjust for the length of the cross hair}
}
\details{
This function returns the galaxy plot to visualize bivariate meta-analysis data, which faithfully retains the information in two separate funnel plots, while providing useful insights into outcome correlations, between-study heterogeneity and joint asymmetry.}

\section{Galaxy plot: a new visualization tool of bivariate meta-analysis studies }{
Funnel plots have been widely used to detect small study effects in the results of univariate meta-analyses. However, there is no existing visualization tool that is the counterpart of the funnel plot in the multivariate setting. We propose a new visualization method, the galaxy plot, which can simultaneously present the effect sizes of bivariate outcomes and their standard errors in a two-dimensional space.
The galaxy plot is an intuitive visualization tool that can aid in interpretation of results of multivariate meta-analysis. It preserves all of the information presented by separate funnel plots for each outcome while elucidating more complex features that may only be revealed by examining the joint distribution of the bivariate outcomes.
}

\value{
\code{res} returns the galaxy plot.
}
\references{Chuan Hong, Rui Duan, Lingzhen Zeng, Rebecca A Hubbard, Thomas Lumley, Richard Riley, Haitao Chu, Stephen E Kimmel and Yong Chen. (2020). Galaxy Plot: A New Visualization Tool of Bivariate Meta-Analysis Studies, American Journal of Epidemiology, https://doi.org/10.1093/aje/kwz286.
}
\author{Chuan Hong, Yong Chen}

\examples{
data(sim_dat)
myplot=galaxy(data=sim_dat, nm.y1="y1", nm.s1="s1", nm.y2="y2", nm.s2="s2", scale=0.15)
}
\keyword{Multivariate meta-analysis}
\keyword{Small study effects}
\keyword{Visualization tool}
