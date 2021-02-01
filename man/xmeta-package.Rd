\name{xmeta-package}
\alias{xmeta-package}
\alias{xmeta}
\docType{package}
\title{
A Tool Box for Multivariate Meta-Analysis
}
\description{
The package \pkg{xmeta} consists of a collection of functions for making inference and detecting publication bias in multivariate meta-analysis (MMA).
}
\details{
\tabular{ll}{
Package: \tab xmeta\cr
Type: \tab Package\cr
Version: \tab 1.3-0\cr
Date: \tab 2021-01-31\cr
License: GPL>=2 \cr
}
}

\section{Inference}{The aim of the estimation methods is to estimate the coefficients 
\eqn{\boldsymbol{\beta}}{\beta} and the components of the between-study (co)variance matrix \eqn{\boldsymbol{\Psi}}{\Psi} 
for multivariate random-effects meta-analysis. 
One major challenge in MMA is the standard inference procedures, such as the maximum likelihood or maximum restricted likelihood inference,
require the within-study correlations, which are usually unavailable. 
Different estimators with and without the knowledge of within study correlation are implemented in the package \pkg{xmeta}. The estimation methods available in function \code{\link{mmeta}} are:
\itemize{
  \item \bold{Restricted maximum likelihood for MMA with continuous outcomes} 
  \item \bold{Composite likelihood method for MMA with continuous outcomes} 
  \item \bold{Method of moment for MMA with continuous outocmes} 
  \item \bold{Improved method for Riley model for MMA with continuous outcomes} 
  \item \bold{Marginal bivariate normal model for MMA with binary outcomes} 
  \item \bold{Marginal beta-binomial model for MMA with binary outcomes}
  \item \bold{Hybrid model for disease prevalence along with sensitivity and specificity for diagnostic test accuracy} 
  \item \bold{Trivariate model for multivariate meta-analysis of diagnostic test accuracy}
}
}

\section{Small study effects}{Detecting and accounting for small study effects are challenging in MMA setting. The multivariate nature is often not fully accounted for by the existing univariate methods. 
The score test for detecting small study effects in MMA when the within-study correlations are unknown is implemented in the function \code{\link{msset}}. 
}

\section{Galaxy Plot}{A New Visualization Tool of Bivariate Meta-Analysis Studies. This function \code{\link{galaxy}} returns the galaxy plot to visualize bivariate meta-analysis data, which faithfully retains the information in two separate funnel plots, while providing 
useful insights into outcome correlations, between-study heterogeneity and joint asymmetry. 
Galaxy plot is the counterpart of the funnel plot in the multivariate setting. 
The galaxy plot is an intuitive visualization tool that can aid in interpretation 
of results of multivariate meta-analysis. It preserves all of the information presented 
by separate funnel plots for each outcome while elucidating more complex features that 
may only be revealed by examining the joint distribution of the bivariate outcomes.
}

\section{Publication bias in bivariate meta-analysis}{The function \code{\link{galaxy.trimfill}} implements a bivariate T&F method accounting for publication bias in bivariate meta-analysis, based on symmetry of the galaxy plot. The bivariate T&F method assumes studies are suppressed based on a weighted sum of the two outcomes. 
We use a searching algorithm to find the optimal direction which gives the most trimmed studies. This is based on the observation that the closer a direction is 
to the truth, the more studies are expected to be trimmed along that direction.   
}

\author{
Author: Chuan Hong, Chongliang Luo, Jiayi Tong, Yong Chen
Maintainer: Jiayi Tong <Jiayi.Tong@pennmedicine.upenn.edu>
Contributor: Rui Duan, Haitao Chu, Yulun Liu 
}

\keyword{package}

