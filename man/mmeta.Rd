\name{mmeta}
\alias{mmeta}
\title{Methods for multiviarate random-effects meta-analysis}
\description{
Methods for multiviarate random-effects meta-analysis}
\usage{
mmeta(data, rhow, type, k, method)
}
\arguments{
  \item{data}{dataset}
  \item{rhow}{within-study correlation}
  \item{type}{either "continuous" or "binary", indicating the type of outcomes.}
  \item{k}{integer indicating the number of outcomes}
  \item{method}{either "nn.reml", "nn.cl", "nn.mom", "nn.rs", "bb.cl", "bn.cl", "tb.cl" or "tn.cl", indicating the estimation method. }
}
\details{
\bold{Inference on the multivariate random-effects meta-analysis for both continuous and binary outcomes}\cr
The function can be used in meta-analyses with continous outcomes and binary outcomes (e.g., mean differences, diagnostic test results in diagnostic accuracy studies, the exposure status of both cases and controls in case-control studies and so on). 
Different estimators with and without the knowledge of within-study correlations are implemented in this function. The estimation methods include \cr
\itemize{
  \item \bold{Restricted maximum likelihood for MMA with continuous outcomes}(nn.reml)
  \item \bold{Composite likelihood method for MMA with continuous outcomes} (nn.cl)
  \item \bold{Moment of method for MMA with continuous outocmes} (nn.mom)
  \item \bold{Improved method for Riley model for MMA with continuous outcomes} (nn.rs)
  \item \bold{Marginal bivariate normal model for MMA with binary outcomes} (bn.cl)
  \item \bold{Marginal beta-binomial model for MMA with binary outcomes}(bb.cl)
  \item \bold{Hybrid model for disease prevalence along with sensitivity and specificity for diagnostic test accuracy} (tb.cl)
  \item \bold{Trivariate model for multivariate meta-analysis of diagnostic test accuracy}(tn.cl)
}
}
\section{Multivariate random-effects meta analysis}{
We consider a meta-analysis with \eqn{m} studies where two outcomes in each study are of interest. 
For the \eqn{i}th study, denote \eqn{Y_{ij}} and \eqn{s_{ij}} the summary measure for the \eqn{j}th outcome of interest and associated standard error respectively, both assumed known, \eqn{i=1, \ldots, m}, and \eqn{j=1,2}. 
Each summary measure \eqn{Y_{ij}} is an estimate of the true effect size \eqn{\theta_{ij}}. 
To account for heterogeneity in effect size across studies, we assume \eqn{\theta_{ij}} to be independently drawn from a common distribution with overall effect size \eqn{\beta_j} and between study variance \eqn{\tau_j^2, j=1,2}. 
Under normal distribution assumption for \eqn{Y_{ij}} and \eqn{\theta_{ij}}, the general bivariate random-effects meta-analysis can be written as 

\deqn{
\left(\begin{array}{c}Y_{i1} \\ Y_{i2} \end{array} \right)\sim N\left( \left(\begin{array}{c}\theta_{i1} \\ \theta_{i2} \end{array} \right), {{\Delta_i}}\right),\quad
{{\Delta_i}}=\left(\begin{array}{cc} s_{i1}^2 & s_{i1}s_{i2}\rho_{\textrm{W}i}\\ s_{i1}s_{i2}\rho_{\textrm{W}i} & s_{i2}^2 \end{array} \right),
}
\deqn{
\left(\begin{array}{c}\theta_{i1} \\ \theta_{i2} \end{array} \right)\sim N\left(\left(\begin{array}{c}\beta_{1} \\ \beta_{2} \end{array} \right), {{\Omega}}\right),\quad
{{\Omega}}=\left(\begin{array}{cc} \tau_{1}^2 & \tau_{1}\tau_{2}\rho_{\textrm{B}}\\ \tau_{1}\tau_{2}\rho_{\textrm{B}} & \tau_{2}^2 \end{array} \right),
}
where \eqn{\Delta_i} and \eqn{\Omega} are the respective within-study and between-study covariance matrices, and \eqn{\rho_{\textrm{W}i}} and \eqn{\rho_{\textrm{B}}} are the respective within-study and between-study correlations.
}

\section{Restricted maximum likelihood for MMA}{
When the within-study correlations are known, inference on the overall effect sizes \eqn{\beta_1} and \eqn{\beta_2} or their comparative measures (e.g., \eqn{\beta_1-\beta_2}) 
can be based on the marginal distribution of \eqn{\left(Y_{i1}, Y_{i2}\right)}

\deqn{
\left(\begin{array}{c}Y_{i1} \\ Y_{i2} \end{array} \right)\sim N\left(\left(\begin{array}{c}\beta_{1} \\ \beta_{2} \end{array} \right), \bf{V_i}\right),
\bf{V_i}=\Delta_i+\Omega=\left(\begin{array}{cc} s_{i1}^2+\tau_{1}^2 & s_{i1}s_{i2}\rho_{wi}+\tau_{1}\tau_{2}\rho_{\textrm{B}}\\ s_{i1}s_{i2}\rho_{wi}+\tau_{1}\tau_{2}\rho_{\textrm{B}} & s_{i2}^2+\tau_{2}^2 \end{array} \right).
}

For simplicity of notation, denote \eqn{\bf{Y_i}=(Y_{i1}, Y_{i2})^{T}}, \eqn{{{\beta}}=(\beta_1, \beta_2)^T}, 
\eqn{\eta_1=(\beta_1,\tau_1^2)^{T}} and \eqn{\eta_2=(\beta_2,\tau_2^2)^{T}}. 
The restricted likelihood of \eqn{(\eta_1, \eta_2, \rho_{\textrm{B}})} can be written as
\deqn{
\log L({{\eta}}_1, {{\eta}}_2, \rho_{\textrm{B}})
=-{1\over 2} \left[\log \left( \Big{|}\sum_{i=1}^m \bf{V_i}^{-1}\Big{|}\right)+ \sum_{i=1}^m\left\{ \log |\bf{V_i}| + (\bf{Y_i}-{{\beta}})^T \bf{V_i}^{-1} (\bf{Y_i}-{{\beta}}) \right\}\right].
}
The parameters \eqn{(\eta_1, \eta_2, \rho_{\textrm{B}})} can be estimated by the restricted maximum likelihood (REML) approach as described in Van Houwelingen et al. (2002). The REML method for 
MMA is specified via \code{method} argument (\code{method="nn.reml"}). 

The standard inference procedures, such as the maximum likelihood or maximum restricted likelihood inference, require the within-study correlations, which are usually unavailable.
In case within-study correlations are unknown, then one can leave the \eqn{\rho_w} argument unspecified, and specify a method that does not require the within-study correlations via \code{method} argument. 
}
\section{Composite likelihood method for MMA with continuous outcomes}{
Chen et al. (2014) proposed a pseudolikelihood method for MMA with unknown within-study correlation. 
The pseudolikelihood method does not require within-study correlations, and is not prone to singular covariance matrix problem. 
In addition, it can properly estimate the covariance between pooled estimates for different outcomes, 
which enables valid inference on functions of pooled estimates, and can be applied to meta-analysis where some studies have outcomes MCAR. 
This composite likelihood method for MMA is specified via \code{method} argument (\code{method="nn.cl"}). 
}

\section{Moment of method for MMA with continuous outocmes}{
Chen et al. (2015) proposed a simple non-iterative method that can be used for the analysis of multivariate meta-analysis datasets 
that has no convergence problems and does not require the use of within-study correlations. 
The strategy is to use standard univariate methods for the marginal effects but also provides valid joint inference for multiple parameters. 
This method method can directly handle missing outcomes under missing completely at random assumption. 
This moment of method for MMA is specified via \code{method} argument (\code{method="nn.mom"})
}

\section{Improved method for Riley model for MMA with continuous outcomes}{
Riley et al.(2008) proposed a working model and an overall synthesis correlation parameter to account for the marginal correlation between outcomes, 
where the only data needed are those required for a separate univariate random-effects meta-analysis. As within-study correlations are not required, 
the Riley method is applicable to a wide variety of evidence synthesis situations. 
However, the standard variance estimator of the Riley method is not entirely correct under many important settings. 
As a consequence, the coverage of a function of pooled estimates may not reach the nominal level even when the number of studies in the multivariate meta-analysis is large. 
Hong et al. (2015) improved the Riley method by proposing a robust variance estimator, 
which is asymptotically correct even when the model is misspecified (i.e., when the likelihood function is incorrect).  
The improved method for Riley model MMA is specified via \code{method} argument (\code{method="nn.rs"})

}
\section{Marginal bivariate normal model for MMA with binary outcomes}{
Diagnostic systematic review is a vital step in the evaluation of diagnostic technologies. In many applications, it invovles pooling paris of sensitivity and specificity
of a dichotomized diagnostic test from multiple studies. 
Chen et al. (2014) proposed a composite likelihood method for bivariate meta-analysis in diagnostic systematic reviews.
The idea of marginal bivariate normal model for MMA with binary outcomes is to construct a composite likelihood (CL) funciton by using an independent working assumption between sensitivity and specificity.
There are three immediate advantages of using this CL method. First, the non-convergence or non positive definite covariance matrix problem is resolved since there is no correlation parameter involved in the CL. 
Secondly, because the two-dimensional integration involved in the standard likelihood is substituted by one-dimensional integrals, the approximation errors are substantially reduced. Thirdly, the inference based on the CL only relies on the marginal normality of logit sensitivity and specificity. 
Hence the proposed method can be more robust than the standard likelihood inference to mis-specifications of the joint distribution assumption. 
This method is specified via \code{method} argument (\code{method="bn.cl"})

}
\section{Marginal beta-binomial model for MMA with binary outcomes}{
When conducting a meta-analysis of studies with bivariate binary outcomes, challenges arise when the within-study correlation and between-study heterogeneity should be taken into account. 
Chen et al. (2015) proposed a marginal beta-binomial model for the meta-analysis of studies with binary outcomes. 
This model is based on the composite likelihood approach, and has several attractive features compared to the existing models such as bivariate generalized linear mixed model (Chu and Cole, 2006) and Sarmanov beta-binomial model (Chen et al., 2012). 
The advantages of the proposed marginal model include modeling the probabilities in the original scale, not requiring any transformation of probabilities or any link function, having closed-form expression of likelihood function, and no constraints on the correlation parameter. 
More importantly, since the marginal beta-binomial model is only based on the marginal distributions, 
it does not suffer from potential misspecification of the joint distribution of bivariate study-specific probabilities. 
Such misspecification is difficult to detect and can lead to biased inference using currents methods. 
This method is specified via \code{method} argument (\code{method="bb.cl"})
}

\section{Hybrid model for disease prevalence along with sensitivity and specificity for diagnostic test accuracy}{
Meta-analysis of diagnostic test accuracy often involves mixture of case-control and cohort studies. 
The existing bivariate random effects models, which jointly model bivariate accuracy indices (e.g., sensitivity and specificity), 
do not differentiate cohort studies from case-control studies, and thus do not utilize the prevalence information contained in the cohort studies. 
The trivariate generalized linear mixed models are only applicable to cohort studies, and more importantly, they assume the common correlation structure across studies, and the trivariate normality on disease prevalence, test sensitivity and specificity after transformation by some pre-specified link functions. 
In practice, very few studies provide justifications of these assumptions, and sometimes these assumptions are violated. 
Chen et al. (2015) evaluated the performance of the commonly used random effects model under violations of these assumptions and propose a simple and robust method to fully utilize 
the information contained in case-control and cohort studies. 
The proposed method avoids making the aforementioned assumptions and can provide valid joint inferences for any functions of overall summary measures of diagnostic accuracy. 
This method is specified via \code{method} argument (\code{method="tb.cl"})

}
\section{Trivariate model for multivariate meta-analysis of diagnostic test accuracy}{
The standard methods for evaluating diagnostic accuracy only focus on sensitivity and specificity and ignore the information on disease prevalence contained in cohort studies. 
Consequently, such methods cannot provide estimates of measures related to disease prevalence, such as population averaged or overall positive and negative predictive values, 
which reflect the clinical utility of a diagnostic test. 
Chen et al. (2014) proposed a hybrid approach that jointly models the disease prevalence along with the diagnostic test sensitivity and specificity in cohort studies, 
and the sensitivity and specificity in case-control studies. In order to overcome the potential computational difficulties in the standard full likelihood inference of the proposed hybrid model, 
an alternative inference procedure was proposed based on the composite likelihood. Such composite likelihood based inference does not suffer computational problems and maintains high relative efficiency. 
In addition, it is more robust to model mis-specifications compared to the standard full likelihood inference. 
This method is specified via \code{method} argument (\code{method="tn.cl"})
}

\value{
 An object of class \code{"mmeta"}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{beta.cov}{covariance matrix of the coefficients.}
}
\references{
Chen, Y., Hong, C. and Riley, R. D. (2015). An alternative pseudolikelihood method for multivariate random-effects meta-analysis. Statistics in medicine, 34(3), 361-380.

Chen, Y., Hong, C., Ning, Y. and Su, X. (2015). Meta-analysis of studies with bivariate binary outcomes: a marginal beta-binomial model approach, Statistics in Medicine (in press).

Hong, C., Riley, R. D. and Chen, Y. (2015). An improved method for multivariate random-effects meta-analysis (in preparation).

Chen, Y., Liu, Y., Ning, J., Nie, L., Zhu, H. and Chu, H. (2014). A composite likelihood method for bivariate meta-analysis in diagnostic systematic reviews. Statistical methods in medical research (in press).

Chen, Y., Cai, Y., Hong, C. and Jackson, D. (2015). Inference for correlated effect sizes using multiple univariate meta-analyses, Statistics in Medicine (provisional acceptance).

Chen, Y., Liu, Y., Ning, J., Cormier J. and Chu H. (2014). A hybrid model for combining case-control and cohort studies in systematic reviews of diagnostic tests, Journal
of the Royal Statistical Society: Series C (Applied Statistics) 64.3 (2015): 469-489.

Chen, Y., Liu, Y., Chu, H., Lee, M. and Schmid, C. (2017) A simple and robust method for multivariate meta-analysis of diagnostic test accuracy, Statistics in Medicine, 36, 105-121.

}
\author{Yong Chen, Yulun Liu}
	
\examples{
data(prostate)
fit.nn=mmeta(data=prostate, type="continuous", k=2, method="nn.cl") 
summary(fit.nn)

rhow=runif(dim(prostate)[1], -0.2, 0.8)
fit.reml=mmeta(data=prostate, rhow=rhow, type="continuous", k=2, method="nn.reml") 
print(fit.reml)

data(nat2)
fit.bb=mmeta(data=nat2, type="binary", k=2, method="bb.cl") 
summary(fit.bb)

data(ca125)
fit.tb=mmeta(data=ca125, type="binary", k=2, method="tb.cl") 
summary(fit.tb)
}

\keyword{random-effects}
\keyword{mutlviariate meta-analysis}
\keyword{meta-analysis of diagnostic accuracy study}

