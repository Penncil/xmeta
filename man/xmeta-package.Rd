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
Version: \tab 1.1\cr
Date: \tab 2019-10-30\cr
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
\author{
Chuan Hong \cr
Chongliang Luo \cr
Rui Duan \cr
Jiayi Tong \cr
Yong Chen \cr
Maintainer: Jiayi Tong <Jiayi.Tong@pennmedicine.upenn.edu>
Contributor: Haitao Chu, Yulun Liu
}

\keyword{package}

