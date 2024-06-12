#' Methods for Composite Likelihood Network Meta-Analysis
#'
#' @import mvtnorm
#' @import MASS
#' @import Matrix
#' @title Methods for Composite Likelihood Network Meta-Analysis
#' @description  Methods for Composite Likelihood Network Meta-Analysis without knowledge of within-study variance and accounting for small sample effect sizes.
#' @author Yu-lun Liu, Bingyu Zhang, Haitao Chu, Yong Chen
#' @param dataout dataset
#' @param narm number of treatment arms in total
#' @param correct methods to correct small sample sandwich estimators
#' @return List with component:
#' @return mu estimated effect sizes
#' @return tau estimated between-study variance
#' @return diagV estimated variance of effect sizes
#' @return V estimated variance matrix of effect sizes
#' @examples
#' require(mvtnorm)
#' require(MASS)
#' require(Matrix)
#' 
#' data(sim_dat_clnma)
#' 3-arm example with REML
#' CLNMA.equal.tau(data, narm=3, correct='KC')
#' 3-arm example without REML
#' CLNMA.equal.tau(data, narm=3, correct='KC', reml_term=FALSE)


CLNMA.equal.tau = function(dataout, narm, correct='NONE', reml_term=TRUE){
  #get estimations
  tau=Iterate.REML.equal.tau(0.5,tol=10^(-5),maxiter=500,dataout,narm,reml_term)
  if (is.na(tau)) {
    return(list(mu=NA,tau=NA, diagV=NA,V=NA))
  }
  
  mu=get.mu.equal.tau(max(tau,0.01),dataout,narm)
  
  #getHessian
  H = Hessian.equal.tau(mu,max(tau,0.01),dataout,narm)
  H = as.matrix(bdiag(list(H)))
  
  #get SS
  SS = SS.equal.tau(mu,tau,dataout,narm,H,correct)
  
  #get variance
  V = solve(H)%*%SS%*%solve(H)
  
  return(list(mu=mu,tau=tau,diagV=diag(V),V=V))
}

######################### Iterative algorithm with equal variance ###########################
### Input: t0, tol, maxiter, dataout, narm
### Note: the reference drug has to be indicated by the largest number (narm)
### Output: t1

Iterate.REML.equal.tau = function (t0,tol=10^(-5),maxiter=50,dataout,narm,reml_term=TRUE){
  delta = 10
  m = 0
  while(delta>tol & m<maxiter){
    A = get.mu.equal.tau(t0,dataout,narm)
    re = function(t){RElik.equal.tau(t,a=A,dataout,narm, reml_term)}
    t1 = nlminb(t0, re, lower = 0, upper= Inf)$par
    delta = max(abs(t1-t0))
    t0 = t1
    m=m+1
  }
  {if(m>=maxiter){t1 =NA}}
  return(t1)
}


######################### Iterative algorithm with equal variance ###########################
### Input: t0, tol, maxiter, dataout, narm
### Note: the reference drug has to be indicated by the largest number (narm)
### Output: t1

Iterate.equal.tau = function (t0,tol=10^(-5),maxiter=50,dataout,narm){
  delta = 10
  m = 0
  while(delta>tol & m<maxiter){
    A = get.mu.equal.tau(t0,dataout,narm)
    re = function(t){lik.equal.tau(t,a=A,dataout,narm)}
    t1 = nlminb(t0, re, lower = 0, upper= Inf)$par
    delta = max(abs(t1-t0))
    t0 = t1
    m=m+1
  }
  {if(m>=maxiter){t1 =NA}}
  return(t1)
}

##################### Calculate mu with equal variance #########################
### Input: t, dataout, narm
### Output: mu

get.mu.equal.tau = function (t,dataout,narm){
  Y = W = rep(0,narm-1)
  for (i in 1:(narm-1)){
    datatemp = select(i,dataout)
    y = datatemp$outcome
    w = 1/(datatemp$sd^2+t)
    Y[i] = sum(y*w)
    W[i] = sum(w)
  }
  
  A = diag(W)
  
  for (j in 1:(narm-2)){
    for (k in (j+1):(narm-1)){
      #indext = i*(narm-1) - i*(i-1)/2+j-
      data = subset(dataout,dataout$t1==j&dataout$t2==k)
      A[j,k] = A[k,j] =-sum(1/(data$sd^2+t))
    }
  }
  mu = solve(-A, Y)
  return(mu)
}

########################### Data format function ##############################
### Input: i,dataout
### Output: data

select = function(i,dataout){
  #Function to select all compares to i
  data1 = subset(dataout,dataout$t2==i)
  data2 = subset(dataout,dataout$t1==i)
  data2$outcome = -data2$outcome
  data = rbind(data1,data2)
  return(data)
}

###################### REML function with equal variance #######################
### Input: t, a, dataout, narm
### Output: REML likelihood

RElik.equal.tau = function (t,a, dataout,narm, reml_term=TRUE){
  l=0
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      mu = a[i] - a[j]
      l.temp = -0.5*sum(log(t+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t+data$sd^2))
      l = l+l.temp
    }
  }
  for(k in 1:(narm-1)){
    data = subset(dataout,dataout$t1==i&dataout$t2==narm)
    mu =a[k]
    l.temp = -0.5*sum(log(t+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t+data$sd^2))
    l = l+l.temp
  }
  reml = l-0.5*log(sum(t+dataout$sd^2)^{-1})
  if (reml_term==TRUE){
    return(-reml)
  }else{
    return(-l)
  }
}


################ Calculate Hessian matrix with equal variance ##################
### Input: a, t, dataout, narm
### Output: H

Hessian.equal.tau = function (a,t,dataout,narm){
  H_mu = matrix(0,nrow = narm-1, ncol = narm-1)
  H_tm =matrix(0,nrow = narm-1, ncol = 1)
  T1 = T2 = temp_mu =temp_tau= temp_mutau= rep(0,narm*(narm-1)/2)
  for (i in 1:(narm-1)){
    data = subset(dataout,dataout$t1==i&dataout$t2==narm)
    y = data$outcome
    vv = data$sd^2
    mu = a[i]
    q1 = sum(1/(vv+t))
    q2 = sum(1/(vv+t)^2)
    q3 = sum((y-mu)/(vv+t)^2)
    q4 = sum((y-mu)^2/(vv+t)^3)
    temp_mu[i] = q1
    temp_mutau[i] = q3
    temp_tau[i] = q2/2-q4
    T1[i] = i
    T2[i] = narm
  }
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      indext = i*(narm-1) - i*(i-1)/2+j-i
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      y = data$outcome
      vv = data$sd^2
      mu = a[i]-a[j]
      q1 = sum(1/(vv+t))
      q2 = sum(1/(vv+t)^2)
      q3 = sum((y-mu)/(vv+t)^2)
      q4 = sum((y-mu)^2/(vv+t)^3)
      temp_mu[indext] = q1
      temp_mutau[indext] = q3
      temp_tau[indext] = q2/2-q4
      H_mu[i,j] = H_mu[j,i] = q1
      T1[indext] = i
      T2[indext] = j
    }}
  for (i in 1:(narm-1)){
    index1 = as.numeric(T1 ==i)
    index2 = as.numeric(T2 ==i)
    ind = index1+index2
    H_mu[i,i] = -sum(ind*temp_mu)
    H_tm[i,] = -sum(ind*temp_mutau)
  }
  
  return(H_mu)
}

################# Meat in sandwich estimator with equal variance ###############
### Input: a, t, dataout, narm, H, correct
### Output: SS

SS.equal.tau = function(a,t,dataout,narm, H, correct='NONE'){
  dd = dim(H)[1]
  SS = matrix(0,nrow = dd,ncol = dd)
  for (i in 1:dim(dataout)[1]){
    I = rep(0,dd)
    {if (dataout$t2[i] == narm){
      m = a[as.numeric(dataout$t1[i])]
      I[as.numeric(dataout$t1[i])] =1
    }
      else {
        m = a[as.numeric(dataout$t1[i])] - a[as.numeric(dataout$t2[i])]
        I[as.numeric(dataout$t1[i])] = 1
        I[as.numeric(dataout$t2[i])] = -1
      } }
    
    if (correct == 'NONE'){
      S = ((dataout$outcome[i]-m)/(t+dataout$sd[i]^2))*I
    }else if (correct == 'KC'){
      H.new = as.numeric(1- (t(I)%*%solve(-H)%*%(I)*((t + dataout$sd[i]^2)^(-1))))
      S = ((dataout$outcome[i]-m)/(t + dataout$sd[i]^2))*I*H.new^(-1/2)
    }else{
      H.new = as.numeric(1- (t(I)%*%solve(-H)%*%(I)*((t + dataout$sd[i]^2)^(-1))))
      S = ((dataout$outcome[i]-m)/(t + dataout$sd[i]^2))*I*H.new^(-1)
    }
    
    SS.temp = t(t(S))%*%S
    SS = SS.temp+SS
  }
  return(SS)}