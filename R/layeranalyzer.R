
##########################################################
#
# Prior representations
#
# layer.prior - Prior representation for layeranalyzer
#     Contains prior 95% credibility bands for core
#     parameters like mu, sigma and dt (characteristic
#     time). Prior 95% credibility bands for other
#     parameter types can also be specified.
#     There is also a method for loading,
#     layer.load.prior, and a couple of global
#     "standard priors", layer.standard.prior,
#     layer.standard.log.prior and layer.wide.prior
#     (even less informative prior).
#
#
# Trond Reitan, 06. June 2024
#
#####################################################


layer.prior=function(mu,dt,sigma,init=NULL,lin=NULL, beta=NULL, obs=NULL, islog=0)
{
  if(is.null(mu))
    stop("If prior is given, it must contain the 95% prior credibility for expected value, given as element 'mu'!")
  if(length(mu)!=2)
    stop("Prior for 'mu' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'mu'!")
  if(mu[1]>=mu[2])
    stop("First element of mu is the lower limit and must be smaller than the upper limit (second element)")
  mu=as.numeric(mu)
  
  if(is.null(dt))
    stop("If prior is given, it must contain the 95% prior credibility for the characteristic time, given as element 'dt'!")
  if(length(dt)!=2)
    stop("Prior for 'dt' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'dt'!")
  if(dt[1]>=dt[2])
    stop("First element of dt is the lower limit and must be smaller than the upper limit (second element)")
  if(dt[1]<=0)
    stop("Lower limit of dt must be greater than zero!")
  dt=as.numeric(dt)
  
  if(is.null(sigma))
    stop("If prior is given, it must contain the 95% prior credibility for the stochastic contributions, given as element 'sigma'!")
  if(length(sigma)!=2)
    stop("Prior for 'sigma' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'sigma'!")
  if(sigma[1]>=sigma[2])
    stop("First element of sigma is the lower limit and must be smaller than the upper limit (second element)")
  if(sigma[1]<=0)
      stop("Lower limit of sigma must be greater than zero!")
  sigma=as.numeric(sigma)
    
  ret=list(mu=mu,dt=dt,sigma=sigma)
  
  if(!is.null(init))
  {
    if(length(init)!=2)
     stop("Prior for 'init' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'init'!")
    if(init[1]>=init[2])
      stop("First element of init is the lower limit and must be smaller than the upper limit (second element)")
    init=as.numeric(init)
    ret$init=init
  }
  else
  {
    ret$init=NULL
  }
  
  if(!is.null(lin))
  {
    if(length(lin)!=2)
     stop("Prior for 'lin' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'lin'!")
  if(lin[1]>=lin[2])
    stop("First element of lin is the lower limit and must be smaller than the upper limit (second element)")
    lin=as.numeric(lin)
    ret$lin=lin
  }
  else
  {
   ret$lin=NULL
  }
  
  if(!is.null(beta))
  {
    if(length(beta)!=2)
     stop("Prior for 'beta' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'beta'!")
    if(beta[1]>=beta[2])
      stop("First element of beta is the lower limit and must be smaller than the upper limit (second element)")
    beta=as.numeric(beta)
    ret$beta=beta
  }
  else
  {
   ret$beta=NULL
  }

  if(!is.null(obs))
  {
    if(length(obs)!=2)
      stop("Prior for 'obs' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'obs'!")
    if(obs[1]>=obs[2])
      stop("First element of obs is the lower limit and must be smaller than the upper limit (second element)")
    if(obs[1]<=0)
      stop("Lower limit of obs must be greater than zero!")
    obs=as.numeric(obs)
    ret$obs=obs
  }
  else
  {
    ret$obs=NULL
  }
  
  if(is.null(islog))
    islog=0
  
  ret$islog=islog
  
  class(ret)="layer.prior"
  return(ret)
}


layer.load.prior=function(filename)
{
  tab=read.csv(filename, sep=";",header=TRUE)
  
  if(is.null(tab))
    stop("Unable to read file!")
  
  if(is.null(tab$is_log))
    stop("Element 'is_log' not found in prior file!")
  
  
  buffer=list(is_log=as.integer(tab$is_log))

  if(is.null(tab$mu1) | is.null(tab$mu2))
    stop("Element 'mu1' and/or 'mu2' not found in prior file!")
  buffer$mu=as.numeric(c(tab$mu1,tab$mu2))

  if(is.null(tab$dt1) | is.null(tab$dt2))
    stop("Element 'dt1' and/or 'dt2' not found in prior file!")
  buffer$dt=as.numeric(c(tab$dt1,tab$dt2))

  if(is.null(tab$s1) | is.null(tab$s2))
    stop("Element 's1' and/or 's2' not found in prior file!")
  buffer$s=as.numeric(c(tab$s1,tab$s2))

  if(!is.null(tab$lin1) & !is.null(tab$lin2))
    buffer$lin=as.numeric(c(tab$lin1,tab$lin2))

  if(!is.null(tab$beta1) & !is.null(tab$beta2))
    buffer$beta=as.numeric(c(tab$beta1,tab$beta2))

  if(!is.null(tab$init1) & !is.null(tab$init2))
    buffer$init=as.numeric(c(tab$init1,tab$init2))

  if(!is.null(tab$obs1) & !is.null(tab$obs2))
    buffer$obs=as.numeric(c(tab$obs1,tab$obs2))
  
  ret=layer.prior(mu=buffer$mu, dt=buffer$dt, sigma=buffer$s, lin=buffer$lin,
                  init=buffer$init, beta=buffer$beta, obs=buffer$obs, islog=buffer$is_log)

  return(ret)
}


layer.standard.prior=layer.prior(mu=c(-10,10),dt=c(0.001,1000),sigma=c(0.01,10),
	                         lin=c(-1,1),beta=c(-1,1),init=c(-100,100),
                                 obs=c(0.01,1))

layer.standard.log.prior=layer.prior(mu=c(-10,10),dt=c(0.001,1000),sigma=c(0.01,10),
                                     lin=c(-1,1),beta=c(-1,1),init=c(-100,100),
                                     obs=c(0.01,1),islog=as.integer(1))

layer.wide.prior=layer.prior(mu=c(-10000,10000),dt=c(0.000001,1000000),
                             sigma=c(0.001,1000),
                             lin=c(-1000,1000),beta=c(-1000,1000),
                             init=c(-100000,100000),
                             obs=c(0.00001,1000))



