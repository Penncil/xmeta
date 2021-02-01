
####' @useDynLib  
#' @import MASS
#' @import metafor
#' @import mvmeta 
#' @title Bivariate trim&fill method
#' @description  Bivariate T&F method accounting for small-study effects in bivariate meta-analysis, based on symmetry of the galaxy plot.
#' @author Chongliang Luo, Yong Chen
#' 
#' @param y1 vector of the effect size estimates of the first outcome
#' @param v1 estimated variance of \code{y1}
#' @param y2 vector of the effect size estimates of the second outcome
#' @param v2 estimated variance of \code{y2} 
#' @param n.grid number of grid (equally spaced) candidate directions that the optimal projection direction are searched among, see Details
#' @param angle angles of candidate projection directions not by grid, this will overwrite n.grid
#' @param estimator   estimator used for the number of trimmed studies in univariate T&F on the projected studies, one of c('R0', 'L0', 'Q0')
#' @param side    either "left" or "right", indicating on which side of the galaxy plot 
#'                   the missing studies should be imputed. If null determined by the univariate T&F 
#' @param rho    correlation between y1 and y2 when computing the variance of the projected studies. Default is the estimated cor(y1, y2)
#' @param method   method to estimate the center for the bivariate outcomes. Default is 'mm', i.e. random-effects model
#' @param method.uni   method to estimate the center for the univariate projected studies using a univariate T&F procedure. Default is 'DL', i.e. fixed-effects model
#' @param maxiter   max number of iterations used in the univariate T&F. Default is 20. 
#' @param var.names    names of the two outcomes used in the galaxy plot (if plotted). Default is c('y1', 'y2')
#' @param scale        constant scale for plotting the galaxy plot for the bivariate studies, Default is 0.02. 
#' @param verbose      plot the galaxy plot? Default is FALSE.

#' @return List with component:
#' @return res  a data.frame of 9 columns and n.grid rows. Each row is the result for projection along one candidate grid direction, and the columns are named: 
#'                'y1.c', 'y2.c' for projected bivariate center, 'y1.f', 'y2.f' for bivariate center using filled studies,
#'               'k0', 'se.k0' for estimated number of trimmed studies and its standard error, 'se.y1.f', 'se.y2.f' for standard errors of 'y1.f', 'y2.f', 
#'               'side.left' for the estimated side
#' @return ID.trim list of vectors of ids of studies been trimmed along each of the candidate direction. 
#' @details  The bivariate T&F method assumes studies are suppressed based on a weighted sum of the two outcomes, i.e. the studies with smallest values 
#'           of z_i = c_1 * y_1i + c_2 * y_2i, i=1,...,N  are suppressed. We use a searching algorithm to find the optimal ratio of c_1  and c_2 (i.e. a direction), 
#'           which gives the most trimmed studies. This is based on the observation that the closer a direction is to the truth, the more studies 
#'           are expected to be trimmed along that direction. We set a sequence of equally-spaced candidate directions with angle a_m = m*pi/M, 
#'           and (c_1, c_2) = (cos(a_m), sin(a_m)), m=1,...,M. 
#'          
#' @references Luo C, Marks-Anglin AK, Duan R, Lin L, Hong C, Chu H, Chen Y. Accounting for small-study effects 
#'                using a bivariate trim and fill meta-analysis procedure. medRxiv. 2020 Jan 1.
#' @examples  
#' require(MASS)
#' require(mvmeta)
#' require(metafor)
#' set.seed(123)
#' mydata <- dat.gen(m.o=50, m.m=20,       # # observed studies, # missing studies
#'                   s.m= c(0.5, 0.5),     #  c(mean(s1), mean(s2))
#'                   angle.LC = pi/4,      # suppress line direction 
#'                   mybeta=c(2,2),        # true effect size
#'                   tau.sq=c(0.1, 0.1),   # true between-study var
#'                   rho.w=0.5, rho.b=0.5, # true within-study and between-study corr
#'                   s.min = 0.1,          # s1i ~ Unif(s.min, 2*s.m[1]-s.min) 
#'                   verbose = TRUE)
#' 
#' y1 <- mydata$mydat.sps$y1
#' y2 <- mydata$mydat.sps$y2
#' v1 <- mydata$mydat.sps$s1^2
#' v2 <- mydata$mydat.sps$s2^2
#' 
#' ## unadjusted est 
#' mv_obs <- mvmeta(cbind(y1, y2), cbind(v1, v2), method='mm')
#' c(mv_obs$coef)
#' # 2.142687 2.237741
#' 
#' estimator <- 'R0' 
#' ## univariate T&F based on y1 or y2
#' y1.rma <- rma(y1, v1, method='FE')
#' y2.rma <- rma(y2, v2, method='FE')
#' y1.tf <- trimfill.rma(y1.rma, estimator = estimator, method.fill = 'DL') 
#' y2.tf <- trimfill.rma(y2.rma, estimator = estimator, method.fill = 'DL') 
#' c(y1.tf$beta, y2.tf$beta)
#' # 2.122231 2.181333
#' c(y1.tf$k0, y2.tf$k0)
#' # 2 8
#' 
#' ## bivariate T&F method (based on galaxy plot)
#' tf.grid <- galaxy.trimfill(y1, v1, y2, v2, n.grid = 12,  
#'                            estimator=estimator, side='left',
#'                            method.uni = 'FE',
#'                            method = 'mm', 
#'                            rho=0.5, maxiter=100, verbose=FALSE) 
#' tf.grid$res
#' tf.grid$res[which(tf.grid$res$k0==max(tf.grid$res$k0)),3:5] 
#' #     y1.f     y2.f k0
#' # 2.053306 2.162347 14  
#' 
#' ## less bias by the proposed bivariate T&F method
#' rbind(true = c(2,2),
#'       unadjusted=c(mv_obs$coef), 
#'       tf.uni = c(y1.tf$beta, y2.tf$beta),
#'       tf.biv = tf.grid$res[which(tf.grid$res$k0==max(tf.grid$res$k0)),3:4])
#' 
#' ## unlike the univariate T&Fs, biv T&F obtains one estimate of # missing studies
#' c(k0.true = 20,
#'   k0.tf.uni.y1 = y1.tf$k0, 
#'   k0.tf.uni.y2 = y2.tf$k0, 
#'   k0.tf.biv = tf.grid$res[which(tf.grid$res$k0==max(tf.grid$res$k0)),5])
#' # k0.true k0.tf.uni.y1 k0.tf.uni.y2    k0.tf.biv 
#' # 20            2            8           14
#' 
#' @export
galaxy.trimfill <- function(y1, v1, y2, v2, n.grid = 12, angle, estimator, side, rho=0, method='mm', method.uni='DL',
                            maxiter=20, var.names=c('y1', 'y2'), scale=0.02, verbose=FALSE){
  k <- length(y1)
  
  ## calculate initial mean est
  # y1.rma <- rma(y1, v1)$beta
  # y2.rma <- rma(y2, v2)$beta
  # y.center <- c(y1.rma$beta, y2.rma$beta)
  mv_obs = mvmeta(cbind(y1,y2), cbind(v1,v2), method=method)
  y.center <- mv_obs$coef
  
  ##
  if(verbose){
    ylim <- c(min(y2), (max(y2)-min(y2))*1.5+min(y2))
    plot(y2 ~ y1, pch='.', ylim=ylim, 
         xlab=var.names[1], ylab=var.names[2],
         main=paste('Trim-fill of galaxy plot, # study =', k))
    plotrix::draw.ellipse(x=y1, y=y2, 
                          a = sqrt(max(v1) / v1)*scale, 
                          b = sqrt(max(v2) / v2)*scale, angle = 0)
    text(y.center[1], y.center[2], '*', col=0+2, cex=3)
    text(min(y1)*0.3+max(y1)*0.7, ylim[1]+(ylim[2]-ylim[1])*(0.99-0*0.03), 
         paste(round(y.center[1],3), round(y.center[2],3)), col=0+2)  # 'iter =', 0, ', k0 =', 0, 
 
  }
  
  # the Linear Combination of (y1, y2) along the n.grid projection lines
  if(!missing(angle)){
    n.grid <- length(angle)
  } else{
    angle <- c(0:(n.grid-1)) * pi / n.grid
  }
  # cat(angle, '\n')
  res <- matrix(NA, n.grid, 2*4+1)  # add side, 1=left, 0=right
  ID.trim <- list()     # keep id trimmed for each direction
  for(ig in 1:n.grid){
    theta <- angle[ig]
    LC.grid <- y1 * cos(theta) + y2 * sin(theta)
    LC.v.grid <- v1*cos(theta)^2 + v2*sin(theta)^2 + rho*sin(theta*2)*sqrt(v1*v2)
    rma.grid <- rma.uni(LC.grid, LC.v.grid, method=method.uni)
    # cat(ig, rma.grid$beta, '\n')
    res[ig, 1:2] <- c(y.center) * c(cos(theta), sin(theta))  # c(rma.grid$beta)
    # LC.center <- sum( y.center * c(cos(theta), sin(theta)) )
    tf.grid <- trimfill.rma(rma.grid, side=side, estimator=estimator, maxiter=maxiter)
    res[ig, 5:6] <- c(tf.grid$k0, tf.grid$se.k0)
    # fill studies symmetric to the trimmed studies
    if(tf.grid$k0 > 0){
      id.trim <- order(LC.grid, decreasing = ifelse(tf.grid$side=='left',TRUE, FALSE))[1:tf.grid$k0]
      ID.trim[[ig]] <- id.trim
      y.c.grid <- c(mvmeta(cbind(y1,y2)[-id.trim,], cbind(v1,v2)[-id.trim,], method=method)$coef)
      y.fill <- rbind(cbind(y1,y2), t(2*y.c.grid - t(matrix(cbind(y1,y2)[id.trim,], ncol=2)) ))
      v.fill <- rbind(cbind(v1,v2), cbind(v1,v2)[id.trim,] ) 
      tmp <- mvmeta(y.fill, v.fill, method=method)
      res[ig, c(3:4,7:8)] <- c(tmp$coef, sqrt(diag(tmp$vcov)))
    } else {
      res[ig, c(3:4,7:8)] <- c(y.center, sqrt(diag(mv_obs$vcov)))
    }
    # add side, 1=left, 0=right
    res[ig, 9] <- ifelse(tf.grid$side=='left', 1, 0)
  } # end for grid 
  colnames(res) <- c('y1.c', 'y2.c', 'y1.f', 'y2.f', 'k0', 'se.k0', 'se.y1.f', 'se.y2.f', 'side.left')
  i.max <- which.max(res[,5])[1]
  angle.max <- angle[i.max]
  if(verbose){
    text(min(y1)*0.3+max(y1)*0.7, ylim[1]+(ylim[2]-ylim[1])*(0.95), paste(res[,5], collapse = ',') )
    text(min(y1)*0.3+max(y1)*0.7, ylim[1]+(ylim[2]-ylim[1])*(0.9), paste(estimator, ', max # trimmed= ', max(res[,5])) )
    abline(0, tan(angle.max+0.0001), lty='dashed' )
    text(res[i.max, 3], res[i.max, 4], '+', col=0+4, cex=3)
    # text(min(y1)*0.3+max(y1)*0.7, ylim[1]+(ylim[2]-ylim[1])*(0.8), paste(round(center.true[1],3), round(center.true[2],3)), col='green')
    # text(min(y1)*0.3+max(y1)*0.7, ylim[1]+(ylim[2]-ylim[1])*(0.75), paste(round(res[i.max,3],3), round(res[i.max,4],3)), col='blue')
    # if(is.null(center.true)){
      ll = c('* obs', '+ T&F')
      cc = c(2,4)
    # }else{
      # ll = c('* obs', '# true', '+ T&F')
      # cc = c(2:4)
    # }
    
    legend('topleft', legend=ll, col=cc)
  }
  return(list(res=data.frame(res), 
              ID.trim=ID.trim))
}


#' @import metafor
#' @title Trim&fill method for univariate meta analysis
#' 
#' @description  Modified metafor::trimfill.rma.uni to avoid the invalid sqrt in k0 calculation when estimator == "Q0"
#' @author Chongliang Luo, Yong Chen
#' 
#' @param x an object of class "rma.uni".  
#' @param side   the same as in metafor::trimfill
#' @param estimator the same as in metafor::trimfill
#' @param maxiter the same as in metafor::trimfill
#' @param method.trim the model used in rma.uni() for estimating the center when trimming studies, default is x$method
#' @param method.fill the model used in rma.uni() for estimating the center after filling studies, default is x$method
#' @param verbose  the same as in metafor::trimfill 
#' @param ilim limits for the imputed values as in metafor::trimfill. If unspecified, no limits are used. 
#' @return  the same as in metafor::trimfill
#' @details  It is recommend using fixed-effects for method.trim and random-effects for method.fill when heterogeneity exists.
#'          
#' @export
trimfill.rma <- function (x, side, estimator = "L0", maxiter = 100, method.trim=NULL, 
                          method.fill=NULL,verbose = FALSE, ilim) {
  if(is.null(method.trim)) method.trim <- x$method
  if(is.null(method.fill)) method.fill <- x$method
  Tn.coef=2
  
  if (missing(side)) 
    side <- NULL
  estimator <- match.arg(estimator, c("L0", "R0", "Q0"))
  if (x$k == 1) 
    stop("Stopped because k = 1.")
  
  yi <- x$yi
  vi <- x$vi
  wi <- x$weights
  ni <- x$ni
  res <- suppressWarnings(rma.uni(yi, vi, weights = wi, mods = sqrt(vi), method = 'DL', weighted = x$weighted))
  tau2 <- res$tau2
  if (is.null(side)) {
    if (res$beta[2] < 0) {
      side <- "right"
    }
    else {
      side <- "left"
    }
  }
  else {
    side <- match.arg(side, c("left", "right"))
  }
  if (side == "right") 
    yi <- -1 * yi
  ix <- sort(yi, index.return = TRUE)$ix
  yi <- yi[ix]
  vi <- vi[ix]
  wi <- wi[ix]
  ni <- ni[ix]
  k <- length(yi)
  k0.sav <- -1
  k0 <- 0
  iter <- 0
  if (verbose) 
    cat("\n")
  while (abs(k0 - k0.sav) > 0 & iter <= maxiter) {
    k0.sav <- k0
    iter <- iter + 1
    if (iter > maxiter) 
      warning("Trim and fill algorithm did not converge.") # output with warning
    
    yi.t <- yi[seq_len(k - k0)]
    vi.t <- vi[seq_len(k - k0)]
    wi.t <- wi[seq_len(k - k0)]
    res <- suppressWarnings(rma.uni(yi.t, vi.t, weights = wi.t, 
                                    method = method.trim, weighted = x$weighted))
    beta <- c(res$beta)
    yi.c <- yi - beta
 
    yi.cp <- yi.c
    
    yi.c.r <- rank(abs(yi.cp), ties.method = "first")
    yi.c.r.s <- sign(yi.cp) * yi.c.r
    
    if (estimator == "R0") {
      k0 <- (k - max(-1 * yi.c.r.s[yi.c.r.s < 0])) - 1
      se.k0 <- sqrt(2 * max(0, k0) + 2)
    }
    if (estimator == "L0") {
      Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
      k0 <- (4 * Sr - k * (k + 1))/(Tn.coef * k - 1)
      varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                         k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                         18 * k * k0 + 6 * k^2 * k0)
      se.k0 <- 4 * sqrt(varSr)/(Tn.coef * k - 1)
    }
    if (estimator == "Q0") {
      # Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
      Sr <- min(sum(yi.c.r.s[yi.c.r.s > 0]), k*(k-1)/2) # cap Sr to avoid sqrt in k0 invalid
      # cat(iter, Sr, '\n')
      k0 <- k - 1/2 - sqrt(2 * k^2 - 4 * Sr + 1/4)
      varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                         k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                         18 * k * k0 + 6 * k^2 * k0)
      se.k0 <- 2 * sqrt(varSr)/sqrt((k - 1/2)^2 - k0 * 
                                      (2 * k - k0 - 1))
    }
    k0 <- max(0, round(k0))
    se.k0 <- max(0, se.k0)
  }
  if (k0 > 0) {
    if (side == "right") {
      yi.c <- -1 * (yi.c - beta)
    }
    else {
      yi.c <- yi.c - beta
    }
    yi.fill <- c(x$yi.f, -1 * yi.c[(k - k0 + 1):k])
    if (!missing(ilim)) {
      ilim <- sort(ilim)
      yi.fill[yi.fill < ilim[1]] <- ilim[1]
      yi.fill[yi.fill > ilim[2]] <- ilim[2]
    }
    vi.fill <- c(x$vi.f, vi[(k - k0 + 1):k])
    wi.fill <- c(x$weights.f, wi[(k - k0 + 1):k])
    ni.fill <- c(x$ni.f, ni[(k - k0 + 1):k])
    attr(yi.fill, "measure") <- x$measure
    res <- suppressWarnings(rma.uni(yi.fill, vi.fill, weights = wi.fill, 
                                    ni = ni.fill, method = method.fill, weighted = x$weighted ))
    res$fill <- c(rep(FALSE, x$k.f), rep(TRUE, k0))
    res$ids <- c(x$ids, (max(x$ids) + 1):(max(x$ids) + k0))
    if (x$slab.null) {
      res$slab <- c(paste("Study", x$ids), paste("Filled", seq_len(k0)))
    }
    else {
      res$slab <- c(x$slab, paste("Filled", seq_len(k0)))
    }
    res$slab.null <- FALSE
  }
  else {
    res <- x
    res$fill <- rep(FALSE, k)
  }
  res$k0 <- k0
  res$se.k0 <- se.k0
  res$side <- side
  res$k0.est <- estimator
  res$k.all <- x$k.all + k0
  if (estimator == "R0") {
    m <- -1:(k0 - 1)
    res$p.k0 <- 1 - sum(choose(0 + m + 1, m + 1) * 0.5^(0 +  m + 2))
  }
  else {
    res$p.k0 <- NA
  }
  class(res) <- c("rma.uni.trimfill", class(res))
  return( res) 
}

 


#' @import MASS
#' @title Generate  bivariate meta analysis studies
#' @description  Generate  bivariate meta analysis studies based on random-effects model, 
#'               some studies with smallest weighted sum of the two outcomes are suppressed.
#' @author Chongliang Luo, Yong Chen
#' 
#' @param m.o number of observed studies 
#' @param m.m number of missing / suppressed studies 
#' @param s.m vector of the mean of the variances of the two outcomes
#' @param angle.LC  direction of suppressing line, default is pi/4, i.e. the studies on the left bottom corner are missing
#' @param mybeta  the true center of the effect sizes  
#' @param tau.sq  between-study variance, the larger it is the more heterogeneity.
#' @param rho.w  within-study correlation of the two outcomes
#' @param rho.b  between-study correlation of the two outcomes
#' @param s.min  minimum of the variances of the outcomes, default is 0.01
#' @param m.m.o  number of studies on one side of the suppressing line been observed, 
#'               i.e. non-deterministic suppressing, default is 0, i.e. deterministic suppressing
#' @param s2.dist options for generating the outcomes' variances. 1=runif, 2=runif^2, 3=runif^4, 4=rnorm 
#' @param verbose logical, galaxy plot the studies? Default FALSE 
#' @references Luo C, Marks-Anglin AK, Duan R, Lin L, Hong C, Chu H, Chen Y. Accounting for small-study effects 
#'                using a bivariate trim and fill meta-analysis procedure. medRxiv. 2020 Jan 1.
#' @export
dat.gen <- function(m.o, m.m, s.m, angle.LC = pi/4,
                    mybeta,             
                    tau.sq,             
                    rho.w,
                    rho.b,
                    s.min=0.01,
                    m.m.o=0,      
                    s2.dist = 2,
                    verbose=F ){  
  # total SS = observed + missing (suppressed)
  m <- m.o+m.m          
  n.m.o <- m.m * m.m.o
  # s1.sq <- runif(m, s1.m*0.1, s1.m*1.9); s2.sq <- runif(m, s2.m*0.1, s2.m*1.9);
  s1.m <- s.m[1]; s2.m <- s.m[2]; 
  
  tau1.sq <- tau.sq[1]; tau2.sq <- tau.sq[2]
  
  myCov.theta = matrix(c(tau1.sq, rho.b*sqrt(tau1.sq*tau2.sq), 
                         rho.b*sqrt(tau1.sq*tau2.sq), tau2.sq), nrow=2)
  myTheta = mvrnorm(n=m, mu=mybeta, Sigma=myCov.theta)
  
  # s1.sq <- rgamma(m, 3, 3/s1.m)^2 * sqrt(abs(myTheta[,1] - mybeta[1])); 
  # s2.sq <- rgamma(m, 3, 3/s2.m)^2 * sqrt(abs(myTheta[,2] - mybeta[2]));
  if(s2.dist==1){
    s1.sq <- runif(m, s.min, 2*s1.m-s.min)
    s2.sq <- runif(m, s.min, 2*s2.m-s.min)
  } else if(s2.dist==2){
    s1.sq <- runif(m, s.min, 2*s1.m-s.min)^2
    s2.sq <- runif(m, s.min, 2*s2.m-s.min)^2
  }else if(s2.dist==3){
    s1.sq <- runif(m, s.min, 2*s1.m-s.min)^4
    s2.sq <- runif(m, s.min, 2*s2.m-s.min)^4
  }else if(s2.dist==4){
    s1.sq <- rnorm(m, s1.m, s.min)^2
    s2.sq <- rnorm(m, s2.m, s.min)^2
  }
  
  y=array(NA, c(m,2))
  for (i in 1:m){
    mySigma = matrix(c(s1.sq[i], rho.w*sqrt(s1.sq[i]*s2.sq[i]), 
                       rho.w*sqrt(s1.sq[i]*s2.sq[i]), s2.sq[i]), nrow=2, byrow=TRUE)
    y[i,]=mvrnorm(n=1, myTheta[i,], Sigma=mySigma)
  }
  y1=y[,1]
  y2=y[,2]
  
  a.sps <- cos(angle.LC+0.0001)
  b.sps <- sin(angle.LC+0.0001)
  ##suppress by y1
  LC.sps=y1*a.sps + y2*b.sps 
  # mydat = data.frame(cbind(LC.sps, y1, y2, s1.sq, s2.sq, tau1.sq, tau2.sq, rho.w, rho.b))
  mydat = data.frame(LC.sps, y1, y2, s1=sqrt(s1.sq), s2=sqrt(s2.sq), tau1.sq, tau2.sq, rho.w, rho.b)
  mydat = mydat[order(mydat$LC.sps),]
  LC.cut <- (mydat$LC.sps[m.m] + mydat$LC.sps[m.m+1])/2
  
  id=c(sample(1:m.m, n.m.o), (m.m+1):nrow(mydat))
  mydat.sps=mydat[id,]
  mydat.sps=data.frame(cbind(id,na.omit(mydat.sps)))
  y.mvmeta <- list()
  y.sps.mvmeta <- list()
  
  if(verbose){
    # y.mvmeta <- mvmeta(cbind(y1,y2), cbind(s1.sq,s2.sq), data=mydat, method = "mm")
    # y.sps.mvmeta <- mvmeta(cbind(y1,y2), cbind(s1.sq,s2.sq), data=mydat.sps, method = "mm")
    # plot.galaxy(mydat$y1,mydat$y2,mydat$s1.sq,mydat$s2.sq, c(1:m)%in% mydat.sps$id, # mydat$LC.sps>LC.cut, 
    #             main=paste0('# studies = ', m.o, '+', m.m), xlab='y1', ylab='y2', scale=scale)
    # text(y.mvmeta$coef[1], y.mvmeta$coef[2], '#', col=1, cex=2)
    # text(y.sps.mvmeta$coef[1], y.sps.mvmeta$coef[2], '*', col=2, cex=2)
    # tmp_dat = data.frame(y1=y1, s1=sqrt(s1.sq), y2=y2, s2=sqrt(s2.sq))
    myplot = galaxy(data=mydat)
    abline(LC.cut / sin(angle.LC+0.0001), -1/tan(angle.LC+0.0001))
  }
  return(list(mydat=mydat, mydat.sps=mydat.sps, y.sps.mvmeta=y.sps.mvmeta, y.mvmeta=y.mvmeta, LC.cut=LC.cut))
}
