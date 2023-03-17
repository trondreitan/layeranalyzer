
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


layer.series.total.number=as.integer(1)
layer.data.series=function(time.points, value.points, std.dev=NULL, 
                           num.meas.per.value=NULL, site=NULL, name=NULL)
{
  time=time.points
  value=value.points
  
  if(is.null(time))
    stop("No 'time' array in the incoming data list!")
    
  if(is.null(value))
    stop("No 'value' array in the incoming data list!")
    
  if(length(time)!=length(value))
    stop("The 'time.points' and 'value.points' part of the data list must correspond, so the lengths must be the same!")
    
  if(length(time)<2)
    stop("There must be at least two data points!")
     
  if(typeof(time)!="integer" & typeof(time)!="numeric" & 
     typeof(time)!="double" & typeof(time)!="POSIXct")
    stop("The 'time' array must be numeric or date-time (POSIXct)!")
  is.datetime=0
  if(typeof(time)=="POSIXct")
    is.datetime=1
  time=as.numeric(time)
  
  if(typeof(value)!="integer" & typeof(time)!="numeric" & typeof(time)!="double")
    stop("The 'value' array must be numeric!")
  if(sum(is.na(value))>0)
    value[is.na(value)]=-10000000
  value=as.numeric(value)
  
  if(!is.null(num.meas.per.value))
  {
    if(typeof(num.meas.per.value)!="integer" & 
      typeof(num.meas.per.value)!="numeric" & 
      typeof(num.meas.per.value)!="double")
      stop("The number of measurement array 'num.meas.per.value' array must be integer!")

    if(length(num.meas.per.value)!=length(value))
      stop("Number of measurements array 'num.meas.per.value' must correspond with the 'value' and 'time' arrays, so the lengths must be the same(1)!")

    if(sum(is.na(num.meas.per.value))>0)
      n[is.na(num.meas.per.value)]=-10000000

    num.meas.per.value=as.integer(num.meas.per.value)
  }
  
  if(!is.null(std.dev))
  {
    if(typeof(std.dev)!="integer" & typeof(std.dev)!="numeric" & 
      typeof(std.dev)!="double")
      stop("The standard deviation array 'std.dev' array must be numeric!")

    if(length(std.dev)!=length(value))
      stop("Standard deviation array 'std.dev' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")

    if(sum(is.na(std.dev))>0)
      std.dev[is.na(std.dev)]=-10000000

    std.dev=as.numeric(std.dev)
  }
  
  if(!is.null(site))
  {
    if(typeof(site)!="integer" & typeof(site)!="numeric" & 
      typeof(site)!="double")
      stop("The site array 'site' must be integer!")
      
    if(length(site)!=length(value))
      stop("Site indicator array 'site' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")
      
    if(sum(is.na(site))>0)
      stop("If given, sites cannot be missing!")
      
    if(sum(sort(unique(site))==(0:max(site)))!=length(unique(site)))
      stop("Sites must be numbered from 0 to #sites-1, with no sites in between missing!")

  site=as.integer(site)
  }
  
  if(!is.null(name) & name!="")
  {
    if(typeof(name)!="character")
      stop("The name must be string (called 'character' in R)!")
  }
  
  ret=list(time=time.points,value=value.points,
    name=sprintf("X%03d",layer.series.total.number))
  layer.series.total.number=layer.series.total.number+1
  if(!is.null(std.dev))
    ret$std.dev=std.dev
  if(!is.null(num.meas.per.value))
    ret$num.meas.per.value=num.meas.per.value
  if(!is.null(site))
    ret$site=site
  if(!is.null(name) & name!="")
    ret$name=name
  if(!is.null(is.datetime))
    ret$is.datetime=is.datetime
  class(ret)="layer.data.series"
  
  return(ret)
}  


read.layer.data.series=function(filename,column.type=c("time","value"), 
                name="",
                header=FALSE,dec=".",sep="", 
                quote="\"'", numerals = c("allow.loss", "warn.loss", "no.loss"),
                row.names=NULL, col.names=column.type, as.is = !stringsAsFactors,
                na.strings = "NA", colClasses = NA, nrows = -1,
                skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                strip.white = FALSE, blank.lines.skip = TRUE,
                comment.char = "#",
                allowEscapes = FALSE, flush = FALSE,
                stringsAsFactors = default.stringsAsFactors(),
                fileEncoding = "", encoding = "unknown", text, skipNul = FALSE)
{
  tab=read.table(filename,header=header,sep=sep,dec=dec,quote=quote,numerals=numerals, 
                 row.names=row.names, col.names=col.names, as.is=as.is,
                 na.strings=na.strings, colClasses=colClasses, nrows=nrows,
                 skip=skip, check.names=check.names, fill=fill,strip.white=strip.white,
                 blank.lines.skip=blank.lines.skip, comment.char=comment.char,
                 allowEscapes=allowEscapes, flush=flush, 
		 stringsAsFactors=stringsAsFactors, fileEncoding=fileEncoding,
                 encoding=encoding,text=text, skipNul=skipNul)
  
  if(is.null(tab))
    stop("File was not read!")

  if(dim(tab)[2]!=length(column.type))
    stop("File has different number of column.type as given in 'column.type' array!")

  for(i in 1:length(column.type))
  {
   if(column.type[i]!="time" & column.type[i]!="value" & 
      column.type[i]!="num.meas.per.value" &
      column.type[i]!="site" & column.type[i]!="std.dev")
     stop(sprintf("Uknown column type: \"%s\"! Available: \"time\",\"value\",\"std.dev\",\"num.meas.per.value\",\"site\"",column.type[i]))
   index=1:length(column.type)
   index=index[index!=i]
   if(sum(column.type[i]==column.type[index])>0)
     stop(sprintf("Column type \"%s\" is repeated!",column.type[i]))
  }

  if(sum(column.type=="time")==0)
   stop("One column type has to be \"time\"!") 
   
  if(sum(column.type=="value")==0)
   stop("One column type has to be \"value\"!") 

  std.dev=NULL
  if(sum(column.type=="std.dev")==1)
    std.dev=tab[,which(column.type=="std.dev")]

  num.meas.per.value=NULL
  if(sum(column.type=="num.meas.per.value")==1)
    num.meas.per.value=tab[,which(column.type=="num.meas.per.value")]

  site=NULL
  if(sum(column.type=="site")==1)
    site=tab[,which(column.type=="site")]
    
  ret=layer.data.series(time.points=tab[,which(column.type=="time")],
                        value.points=tab[,which(column.type=="value")],
                        std.dev=std.dev, num.meas.per.value=num.meas.per.value,
                        site=site,name=name)
}

# test:
# x1=read.layer.data.series("test_1layer.txt")
# x2=read.layer.data.series("brach_lspec.txt",column.type=c("time","value","std.dev"))

layer.series.structure=function(timeseries, numlayers=1, lin.time=FALSE,
                              time.integral=NULL, no.pull=FALSE, no.sigma=NULL,
                              regional.mu=FALSE, regional.lin.time=FALSE,
                              regional.pull=NULL, regional.sigma=NULL,
                              correlated.sigma=NULL, pairwise.correlated.sigma=NULL,
                              one.dim.sigma=NULL, grouping.sigma=NULL,
			      remove.sigma=NULL, differentiate.sigma=NULL,
			      differentiate.pull=NULL, differentiate.mu=FALSE,
                              differentiate.lin.time=FALSE,
                              init.0=FALSE, init.time=NULL, init.same.sites=FALSE,
			      init.same.layers=FALSE, init.specified=NULL,
			      allow.pos.pull=FALSE, period=NULL, 
                              prior=layer.standard.prior)
{
  if(is.null(timeseries$time))
    stop("No 'time' array in the incoming time series!")
    
  if(is.null(timeseries$value))
    stop("No 'value' array in the incoming time series!")
    
  if(length(timeseries$time)!=length(timeseries$value))
    stop("The 'time' and 'value' part of the time series must correspond, so the lengths must be the same!")
    
  if(length(timeseries$time)<2)
    stop("There must be at least two data points in time series!")
     
  if(typeof(timeseries$time)!="integer" & typeof(timeseries$time)!="numeric" & 
     typeof(timeseries$time)!="double" & typeof(timeseries$time)!="POSIXct")
    stop("The 'time' array in the time series must be numeric or date-time (POSIXct)!")
  timeseries$is.datetime=0
  if(typeof(timeseries$time)=="POSIXct")
    timeseries$is.datetime=1
  timeseries$time=as.numeric(timeseries$time)
  
  if(typeof(timeseries$value)!="integer" & typeof(timeseries$time)!="numeric" & typeof(timeseries$time)!="double")
    stop("The 'value' array in the time series must be numeric!")
  if(sum(is.na(timeseries$value))>0)
    timeseries$value[is.na(timeseries$value)]=-10000000
  timeseries$value=as.numeric(timeseries$value)
  
  if(!is.null(timeseries$num.meas.per.value))
  {
    if(typeof(timeseries$num.meas.per.value)!="integer" & 
      typeof(timeseries$num.meas.per.value)!="numeric" & 
      typeof(timeseries$num.meas.per.value)!="double")
      stop("The number of measurement array 'num.meas.per.value' array in the time series must be integer!")

    if(length(timeseries$num.meas.per.value)!=length(timeseries$value))
      stop("Number of measurements array 'num.meas.per.value' in the time series must correspond with the 'value' and 'time' arrays, so the lengths must be the same(2)!")

    if(sum(is.na(timeseries$num.meas.per.value))>0)
      timeseries$n[is.na(timeseries$num.meas.per.value)]=-10000000
    
    timeseries$num.meas.per.value=as.integer(timeseries$num.meas.per.value)
  }
  
  if(!is.null(timeseries$std.dev))
  {
    if(typeof(timeseries$std.dev)!="integer" & typeof(timeseries$std.dev)!="numeric" & 
      typeof(timeseries$std.dev)!="double")
      stop("The standard deviation array 'std.dev' array must be numeric!")

    if(length(timeseries$std.dev)!=length(timeseries$value))
      stop("Standard deviation array 'std.dev' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")

    if(sum(is.na(timeseries$std.dev))>0)
      timeseries$std.dev[is.na(timeseries$std.dev)]=-10000000

    timeseries$std.dev=as.numeric(timeseries$std.dev)
  }
  
  if(!is.null(timeseries$site))
  {
    if(typeof(timeseries$site)!="integer" & typeof(timeseries$site)!="numeric" & 
      typeof(timeseries$site)!="double")
      stop("The site array 'site' must be integer!")
      
    if(length(timeseries$site)!=length(timeseries$value))
      stop("Site indicator array 'site' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")
      
    if(sum(is.na(timeseries$site))>0)
      stop("If given, sites cannot be missing!")
      
    if(sum(sort(unique(timeseries$site))==(0:max(timeseries$site)))!=length(unique(timeseries$site)))
      stop("Sites must be numbered from 0 to #sites-1, with no sites in between missing!")
    
    timeseries$site=as.integer(timeseries$site)  
  }
  
  if(!is.null(timeseries$name))
  {
    if(typeof(timeseries$name)!="character")
      stop("The name must be string (called 'character' in R)!")
  }
  

  if(!is.null(numlayers))
  {
    if(typeof(numlayers)!="integer" & typeof(numlayers)!="double" & typeof(numlayers)!="numeric")
      stop("Number of layers, 'numlayers', must be an integer")
  }
  if(is.null(numlayers))
  {
    warning("Number of layers not given. Defaulting to 1 (Ornstein-Uhlenbeck)")
    numlayers=1
  }
  numlayers=as.integer(numlayers)
  if(numlayers<1)
  {
    stop("Number of layers must be 1 or more.")
  }
  if(numlayers>=100)
  {
    stop("Cannot handle 100 layers or more (at this time).")
  }
  
  
  ret=list(timeseries=timeseries, numlayers=numlayers)

  
  if(!is.null(lin.time))
  {
    if(typeof(lin.time)!="logical")
      stop("Linear time indicator 'lin.time' must be a logical")
    ret$lin.time=lin.time
  }
  
  if(!is.null(time.integral))
  {
    if(typeof(time.integral)!="integer" & typeof(time.integral)!="double")
      stop("Time integral layer specification, 'time.integral', must be an integer or a vector of such, specifying layers")
    time.integral=as.integer(time.integral)
    if(length(time.integral)>1 | time.integral[1]>0)
    {
      if(sum(time.integral>(numlayers-1) | time.integral<0)>0)
        stop("Time integral layer specification out of range!")
    }
    ret$time.integral=time.integral
  }
  
  if(!is.null(no.pull))
  {
    if(typeof(no.pull)!="logical")
      stop("No pull indicator, 'no.pull', must be a logical")
    ret$no.pull=no.pull
  }
  
  if(!is.null(no.sigma))
  {
    #show(typeof(no.sigma))
    if(typeof(no.sigma)!="integer" & typeof(no.sigma)!="double" & typeof(no.sigma)!="numeric")
      stop("No sigma indicator, 'no.sigma', must be an integer vector specifying layers")
    no.sigma=as.integer(no.sigma)
    if(length(no.sigma)>1 | no.sigma[1]>0)
    {
      if(sum(no.sigma>numlayers | no.sigma<0)>0)
        stop("No sigma specification out of range!")
    }
    ret$no.sigma=no.sigma
  }
  
  if(!is.null(regional.mu))
  {
    if(typeof(regional.mu)!="logical")
      stop("Regional expectation indicator, 'regional.mu', must be a logical")
    ret$regional.mu=regional.mu
  }

  if(!is.null(regional.lin.time))
  {
    if(typeof(regional.lin.time)!="logical")
      stop("Regional linear time indicator, 'regional.lin.time', must be a logical")
    ret$regional.lin.time=regional.lin.time
  }
  
  if(!is.null(regional.pull))
  {
    if(typeof(regional.pull)!="integer" & typeof(regional.pull)!="double" & typeof(regional.pull)!="numeric")
      stop("Regional pull specification, 'regional pull', must be a vector of integers specifying layers")
    regional.pull=as.integer(regional.pull)
    if(length(regional.pull)>1 | regional.pull[1]>0)
    {
      if(sum(regional.pull>numlayers | regional.pull<0)>0)
        stop("Regional pull layer specification out of range!")
    }
    ret$regional.pull=regional.pull
  }

  if(!is.null(regional.sigma))
  {
    if(typeof(regional.sigma)!="integer" & typeof(regional.sigma)!="double" & typeof(regional.sigma)!="numeric")
      stop("Regional sigma specification, 'regional sigma', must be a vector of integers specifying layers")
    regional.sigma=as.integer(regional.sigma)
    if(length(regional.sigma)>1 | regional.sigma[1]>0)
    {
      if(sum(regional.sigma>numlayers | regional.sigma<0)>0)
        stop("Rregional sigma layer specification out of range!")
    }
    ret$regional.sigma=regional.sigma
  }
  
  if(!is.null(correlated.sigma))
  {
    if(typeof(correlated.sigma)!="integer" & typeof(correlated.sigma)!="double" & typeof(correlated.sigma)!="numeric")
      stop("Correlated sigma indicator, 'correlated.sigma', must be an integer vector specifying layers")
    correlated.sigma=as.integer(correlated.sigma)
    if(length(correlated.sigma)>1 | correlated.sigma[1]>0)
    {
      if(sum(correlated.sigma>numlayers | correlated.sigma<0)>0)
        stop("Correlated sigma specification out of range!")
    }
    ret$correlated.sigma=correlated.sigma
  }

  if(!is.null(pairwise.correlated.sigma))
  {
    if(typeof(pairwise.correlated.sigma)!="integer" & typeof(pairwise.correlated.sigma)!="double" & typeof(pairwise.correlated.sigma)!="numeric")
      stop("Pairwise correlated sigma indicator, 'pairwise.correlated.sigma', must be an integer vector specifying layers")
    pairwise.correlated.sigma=as.integer(pairwise.correlated.sigma)
    if(length(pairwise.correlated.sigma)>1 | pairwise.correlated.sigma[1]>0)
    {
      if(sum(pairwise.correlated.sigma>numlayers | pairwise.correlated.sigma<0)>0)
        stop("Pairwise correlated sigma specification out of range!")
    }
    ret$pairwise.correlated.sigma=pairwise.correlated.sigma
  }
  
  if(!is.null(one.dim.sigma))
  {
    if(typeof(one.dim.sigma)!="integer" & typeof(one.dim.sigma)!="double" & typeof(one.dim.sigma)!="numeric")
      stop("1D sigma indicator, 'one.dim.sigma', must be an integer vector specifying layers")
    one.dim.sigma=as.integer(one.dim.sigma)
    if(length(one.dim.sigma)>1 | one.dim.sigma[1]>0)
    {
      if(sum(one.dim.sigma>numlayers | one.dim.sigma<0)>0)
        stop("1D sigma specification out of range!")
    }
    ret$one.dim.sigma=one.dim.sigma
  }

  if(!is.null(grouping.sigma))
  {
    if(typeof(grouping.sigma)!="integer" & typeof(grouping.sigma)!="double" & typeof(grouping.sigma)!="numeric")
      stop("Sigma grouping indicator, 'grouping.sigma', must be an integer vector specifying layers")
    grouping.sigma=as.integer(grouping.sigma)
    if(length(grouping.sigma)>1 | grouping.sigma[1]>0)
    {
      if(sum(grouping.sigma>numlayers | grouping.sigma<0)>0)
        stop("Sigma grouping specification out of range!")
    }
    ret$grouping.sigma=grouping.sigma
  }

  if(!is.null(remove.sigma))
  {
    if(typeof(remove.sigma)!="integer" & typeof(remove.sigma)!="double" & typeof(remove.sigma)!="numeric")
      stop("Sigma correlation removal indicator, 'remove.sigma', must be an integer vector specifying layers")
    remove.sigma=as.integer(remove.sigma)
    if(length(remove.sigma)>1 | remove.sigma[1]>0)
    {
      if(sum(remove.sigma>numlayers | remove.sigma<0)>0)
        stop("Sigma correlation removal specification out of range!")
    }
    ret$remove.sigma=remove.sigma
  }

  if(!is.null(differentiate.sigma))
  {
    if(typeof(differentiate.sigma)!="integer" & typeof(differentiate.sigma)!="double" & typeof(differentiate.sigma)!="numeric")
      stop("Differentiated sigma indicator, 'differentiate.sigma', must be an integer vector specifying layers")
    differentiate.sigma=as.integer(differentiate.sigma)
    if(length(differentiate.sigma)>1 | differentiate.sigma[1]>0)
    {
      if(sum(differentiate.sigma>numlayers | differentiate.sigma<0)>0)
        stop("Differentiated sigma specification out of range!")
    }
    ret$differentiate.sigma=differentiate.sigma
  }

  if(!is.null(differentiate.pull))
  {
    if(typeof(differentiate.pull)!="integer" & typeof(differentiate.pull)!="double" & typeof(differentiate.pull)!="numeric")
      stop("Differentiated pull indicator, 'differentiate.pull', must be an integer vector specifying layers")
    differentiate.pull=as.integer(differentiate.pull)
    if(length(differentiate.pull)>1 | differentiate.pull[1]>0)
    {
      if(sum(differentiate.pull>numlayers | differentiate.pull<0)>0)
        stop("Differentiated pull specification out of range!")
    }
    ret$differentiate.pull=differentiate.pull
  }

  if(!is.null(differentiate.mu))
  {
    if(typeof(differentiate.mu)!="logical")
      stop("Differentiated expectancy indicator, 'differentiate.mu', must be a logical")
    ret$differentiate.mu=differentiate.mu
  }
  
  
  if(!is.null(differentiate.lin.time))
  {
    if(typeof(differentiate.lin.time)!="logical")
      stop("Differentiated linear time indicator, 'differentiate.lin.time', must be a logical")
    ret$differentiate.lin.time=differentiate.lin.time
  }
  
  if(!is.null(init.0))
  {
    if(typeof(init.0)!="logical")
      stop("Initial value treatement indicator, 'init.0', must be a logical")
    ret$init.0=init.0
  }
  
  if(!is.null(init.time))
  {
    if(typeof(init.time)!="numeric" & typeof(init.time)!="double" & 
       typeof(init.time)!="POSIXct" )
      stop("Initial value at given time treatement indicator, 'init.time', must be a numeric or POSIXct")
    if(length(init.time)>1)
      stop("Initial value at given time treatement indicator, 'init.time', cannot be a vector")

    init.datetime=0
    if(typeof(init.time)=="POSIXct" )
      init.datetime=1

    init.time=as.numeric(init.time)
    init.datetime=as.integer(init.datetime)

    ret$init.time=init.time
    ret$init.datetime=init.datetime
  }
  
  if(!is.null(init.same.sites))
  {
    if(typeof(init.same.sites)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.sites', must be a logical")
    ret$init.same.sites=init.same.sites
  }
  
  if(!is.null(init.same.layers))
  {
    if(typeof(init.same.layers)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.layers', must be a logical")
    ret$init.same.layers=init.same.layers
  }
  
  if(!is.null(init.specified))
  {
    if(typeof(init.specified)!="numeric" & typeof(init.specified)!="double")
      stop("Specified initial value, 'init.specified', must be a numeric. If time has been specified with POSIXct, cast this as numeric here.")
    if(length(init.specified)!=2)
      stop("Specified initial value, 'init.specified', must be a vector with 2 values, 'time' and 'value'. If time has been specified with POSIXct, cast this as numeric")
    init.specified=as.numeric(init.specified)
    ret$init.specified=init.specified
  }
  
  if(!is.null(allow.pos.pull))
  {
    if(typeof(allow.pos.pull)!="logical")
      stop("Indicator for allowing positive (unstable) pulls, 'allow.pos.pull', must be a logical")
    ret$allow.pos.pull=allow.pos.pull
  }
  
  if(!is.null(period))
  {
    if(typeof(period)!="numeric" & typeof(period)!="double")
      stop("Periodial trigonometric regressor indicator, 'period', must be a numeric vector")
    period=as.numeric(period)
    ret$period=period
  }
  
  if(!is.null(prior))
  {
    if(is.null(prior$mu))
      stop("If prior is given, it must contain the 95% prior credibility for expected value, given as element 'mu'!")
    if(length(prior$mu)!=2)
      stop("Prior for 'mu' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'mu'!")

    if(is.null(prior$dt))
      stop("If prior is given, it must contain the 95% prior credibility for the characteristic time, given as element 'dt'!")
    if(length(prior$dt)!=2)
      stop("Prior for 'dt' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'dt'!")

    if(is.null(prior$s))
      stop("If prior is given, it must contain the 95% prior credibility for the stochastic contributions, given as element 's'!")
    if(length(prior$s)!=2)
      stop("Prior for 's' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 's'!")
    
    if(is.null(prior$init))
      prior$init=c(prior$mu[1]-3*(prior$mu[2]-prior$mu[1]),
                   prior$mu[2]+3*(prior$mu[2]-prior$mu[1]))
    if(length(prior$init)!=2)
      stop("Prior for 'init' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'init'!")
    
    if(is.null(prior$lin) & lin.time!=0)
      stop("If prior is given and linear time dependency is used, it must contain the 95% prior credibility for the linear time dependency, given as element 'lin'!")
    if(is.null(prior$lin))
      prior$lin=c(-1,1) # doesn't really matter since it's not going to be used
    if(length(prior$lin)!=2)
      stop("Prior for 'lin' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'lin'!")
    
    if(is.null(prior$beta))
      prior$beta=c(-1,1)
    if(length(prior$beta)!=2)
      stop("Prior for 'beta' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'beta'!")
    
    if(is.null(prior$obs) & is.null(timeseries$std.dev))
      stop("If prior is given and measurement-wise observational standard deviation is not given, the prior must contain the 95% prior credibility for the observational standard deviation, given as element 'obs'!")
    if(is.null(prior$obs))
      prior$obs=c(0.01,1) # doesn't really matter since it's not going to be used
    if(length(prior$obs)!=2)
      stop("Prior for 'obs' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'obs'!")
  }
  if(is.null(prior))
    prior=layer.standard.prior
  ret$prior=prior  
  
  class(ret)="layer.series.structure"
  return(ret)
}
		      


layer.analyzer=function(... ,
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.model.likelihood=TRUE,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE,causal=NULL,causal.symmetric=NULL,corr=NULL,
  smoothing.specs=
   list(do.smoothing=FALSE,smoothing.time.diff=0,smoothing.start=NULL,smoothing.end=NULL,
        num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE),
  realization.specs=list(do.realizations=FALSE,num.realizations=1000,strategy="N",
        realization.time.diff=0,realization.start=NULL,realization.end=NULL))
{
  data.structure=list(...)
  ret=layer.analyzer.timeseries.list(data.structure,
      num.MCMC=num.MCMC, spacing=spacing,burnin=burnin,num.temp=num.temp,
      do.model.likelihood=do.model.likelihood,
      do.maximum.likelihood=do.maximum.likelihood,
      maximum.likelihood.numstart=maximum.likelihood.numstart,
      silent.mode=silent.mode,talkative.burnin=talkative.burnin,
      talkative.likelihood=talkative.likelihood,
      id.strategy=id.strategy,use.stationary.stdev=use.stationary.stdev,
      T.ground=T.ground, # start.parameters=0,
      use.half.lives=use.half.lives, mcmc=mcmc,
      causal=causal,causal.symmetric=causal.symmetric,corr=corr,
      smoothing.specs=smoothing.specs, realization.specs=realization.specs)
  return(ret)
}


layer.analyzer.timeseries.list=function(data.structure ,
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.model.likelihood=TRUE,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE,causal=NULL,causal.symmetric=NULL,corr=NULL,
  smoothing.specs=
   list(do.smoothing=FALSE,smoothing.time.diff=0,smoothing.start=NULL,smoothing.end=NULL,
        num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE),
  realization.specs=list(do.realizations=FALSE,num.realizations=1000,strategy="N",
        realization.time.diff=0,realization.start=NULL,realization.end=NULL))
{
 n=length(data.structure)

 for(i in 1:n)
 {
  if(typeof(data.structure[[i]])!="list")
    stop("data.structure[[i]] must be a list, preferrably belonging to the 'layer.series.structure' class (as this class has the right elements).")
  
  if(is.null(data.structure[[i]]$timeseries))
    stop("data.structure[[i]] does not contain a time series (of type layer.data.series)!")
  
  if(is.null(data.structure[[i]]$timeseries$time))
    stop("No 'time' array in the incoming data list!")
    
  if(is.null(data.structure[[i]]$timeseries$value))
    stop("No 'value' array in the incoming data list!")
    
  if(length(data.structure[[i]]$timeseries$time)!=length(data.structure[[i]]$timeseries$value))
    stop("The 'time' and 'value' part of the data list must correspond, so the lengths must be the same!")
    
  if(length(data.structure[[i]]$timeseries$time)<2)
    stop("There must be at least two data points!")
     
  if(typeof(data.structure[[i]]$timeseries$time)!="integer" & 
     typeof(data.structure[[i]]$timeseries$time)!="numeric" & 
     typeof(data.structure[[i]]$timeseries$time)!="double" & 
     typeof(data.structure[[i]]$timeseries$time)!="POSIXct")
    stop("The 'time' array must be numeric or date-time (POSIXct)!")
  data.structure[[i]]$timeseries$is.datetime=0
  if(typeof(data.structure[[i]]$timeseries$time)=="POSIXct")
    data.structure[[i]]$timeseries$is.datetime=1
  data.structure[[i]]$timeseries$time=as.numeric(data.structure[[i]]$timeseries$time)
  
  if(typeof(data.structure[[i]]$timeseries$value)!="integer" & 
     typeof(data.structure[[i]]$timeseries$time)!="numeric" & 
     typeof(data.structure[[i]]$timeseries$time)!="double")
    stop("The 'value' array must be numeric!")
  if(sum(is.na(data.structure[[i]]$timeseries$value))>0)
    data.structure[[i]]$timeseries$value[is.na(data.structure[[i]]$timeseries$value)]=-10000000
  data.structure[[i]]$timeseries$value=as.numeric(data.structure[[i]]$timeseries$value)
  
  if(!is.null(data.structure[[i]]$timeseries$num.meas.per.value))
  {
    if(typeof(data.structure[[i]]$timeseries$num.meas.per.value)!="integer" & 
      typeof(data.structure[[i]]$timeseries$num.meas.per.value)!="numeric" & 
      typeof(data.structure[[i]]$timeseries$num.meas.per.value)!="double")
      stop("The number of measurement array 'num.meas.per.value' array must be integer!")
    
    if(length(data.structure[[i]]$timeseries$num.meas.per.value)!=
       length(data.structure[[i]]$timeseries$value))
      stop("Number of measurements array 'num.meas.per.value' must correspond with the 'value' and 'time' arrays, so the lengths must be the same(3)!")
    
    if(sum(is.na(data.structure[[i]]$timeseries$num.meas.per.value))>0)
      data.structure[[i]]$timeseries$num.meas.per.value[is.na(data.structure[[i]]$timeseries$num.meas.per.value)]=-10000000
  }
  if(is.null(data.structure[[i]]$timeseries$num.meas.per.value))
    data.structure[[i]]$timeseries$num.meas.per.value=as.integer(0)
  data.structure[[i]]$timeseries$num.meas.per.value=
    as.integer(data.structure[[i]]$timeseries$num.meas.per.value)
  
  if(!is.null(data.structure[[i]]$timeseries$std.dev))
  {
    if(typeof(data.structure[[i]]$timeseries$std.dev)!="integer" & 
      typeof(data.structure[[i]]$timeseries$std.dev)!="numeric" & 
      typeof(data.structure[[i]]$timeseries$std.dev)!="double")
      stop("The standard deviation array 'std.dev' array must be numeric!")

    if(length(data.structure[[i]]$timeseries$std.dev)!=length(data.structure[[i]]$timeseries$value))
      stop("Standard deviation array 'std.dev' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")

    if(sum(is.na(data.structure[[i]]$timeseries$std.dev))>0)
      data.structure[[i]]$timeseries$std.dev[is.na(data.structure[[i]]$timeseries$std.dev)]=-10000000
  }
  if(is.null(data.structure[[i]]$timeseries$std.dev))
    data.structure[[i]]$timeseries$std.dev=as.numeric(0)
  data.structure[[i]]$timeseries$std.dev=as.numeric(data.structure[[i]]$timeseries$std.dev)
  
  if(!is.null(data.structure[[i]]$timeseries$site))
  {
    if(typeof(data.structure[[i]]$timeseries$site)!="integer" & 
      typeof(data.structure[[i]]$timeseries$site)!="numeric" & 
      typeof(data.structure[[i]]$timeseries$site)!="double")
      stop("The site array 'site' must be integer!")
      
    if(length(data.structure[[i]]$timeseries$site)!=length(data.structure[[i]]$timeseries$value))
      stop("Site indicator array 'site' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")
      
    if(sum(is.na(data.structure[[i]]$timeseries$site))>0)
      stop("If given, sites cannot be missing!")
      
    if(sum(sort(unique(data.structure[[i]]$timeseries$site))==(0:max(data.structure[[i]]$timeseries$site)))!=length(unique(data.structure[[i]]$timeseries$site)))
      stop("Sites must be numbered from 0 to #sites-1, with no sites in between missing!")
  }
  if(is.null(data.structure[[i]]$timeseries$site))
    data.structure[[i]]$timeseries$site=as.integer(0)
  data.structure[[i]]$timeseries$site=as.integer(data.structure[[i]]$timeseries$site)
  
  if(!is.null(data.structure[[i]]$timeseries$name))
  {
    if(typeof(data.structure[[i]]$timeseries$name)!="character")
      stop("The name must be string (called 'character' in R)!")
  }
  if(is.null(data.structure[[i]]$timeseries$name))
  {
    data.structure[[i]]$timeseries$name="data.structure[[i]]"
  }
  
  if(!is.null(data.structure[[i]]$numlayers))
  {
    if(typeof(data.structure[[i]]$numlayers)!="integer" & typeof(data.structure[[i]]$numlayers)!="double" & typeof(data.structure[[i]]$numlayers)!="numeric")
      stop("Number of layers, 'numlayers', must be an integer")
  }
  if(is.null(data.structure[[i]]$numlayers))
  {
    warning("Number of layers not given. Defaulting to 1 (Ornstein-Uhlenbeck)")
    data.structure[[i]]$numlayers=1
  }
  data.structure[[i]]$numlayers=as.integer(data.structure[[i]]$numlayers)
  if(data.structure[[i]]$numlayers<1)
  {
    stop("Number of layers must be 1 or more.")
  }
  if(data.structure[[i]]$numlayers>=100)
  {
    stop("Cannot handle 100 layers or more (at this time).")
  }
  
  if(!is.null(data.structure[[i]]$lin.time))
  {
    if(typeof(data.structure[[i]]$lin.time)!="logical")
      stop("Linear time indicator 'lin.time' must be a logical")
  }
  if(is.null(data.structure[[i]]$lin.time))
  {
    data.structure[[i]]$lin.time=0
  }
  data.structure[[i]]$lin.time=as.integer(data.structure[[i]]$lin.time)
  
  if(!is.null(data.structure[[i]]$time.integral))
  {
    if(typeof(data.structure[[i]]$time.integral)!="integer" & typeof(data.structure[[i]]$time.integral)!="double")
      stop("Time integral layer specification, 'time.integral', must be an integer or a vector of such, specifying layers")
  }
  if(is.null(data.structure[[i]]$time.integral))
  {
    data.structure[[i]]$time.integral=0
  }
  data.structure[[i]]$time.integral=as.integer(data.structure[[i]]$time.integral)
  if(length(data.structure[[i]]$time.integral)>1 | data.structure[[i]]$time.integral[1]>0)
  {
    if(sum(data.structure[[i]]$time.integral>(data.structure[[i]]$numlayers-1) | data.structure[[i]]$time.integral<0)>0)
      stop("Time integral layer specification out of range!")
  }
  
  if(!is.null(data.structure[[i]]$no.pull))
  {
    if(typeof(data.structure[[i]]$no.pull)!="logical")
      stop("No pull indicator, 'no.pull', must be a logical")
  }
  if(is.null(data.structure[[i]]$no.pull))
  {
    data.structure[[i]]$no.pull=0
  }
  data.structure[[i]]$no.pull=as.integer(data.structure[[i]]$no.pull)
  
  if(!is.null(data.structure[[i]]$no.sigma))
  {
    #show(typeof(data.structure[[i]]$no.sigma))
    if(typeof(data.structure[[i]]$no.sigma)!="integer" & typeof(data.structure[[i]]$no.sigma)!="double" & typeof(data.structure[[i]]$no.sigma)!="numeric")
      stop("No sigma indicator, 'no.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$no.sigma))
  {
    data.structure[[i]]$no.sigma=0
  }
  data.structure[[i]]$no.sigma=as.integer(data.structure[[i]]$no.sigma)
  if(length(data.structure[[i]]$no.sigma)>1 | data.structure[[i]]$no.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$no.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$no.sigma<0)>0)
      stop("No sigma specification out of range!")
  }

  if(!is.null(data.structure[[i]]$regional.mu))
  {
    if(typeof(data.structure[[i]]$regional.mu)!="logical")
      stop("Regional expectation indicator, 'regional.mu', must be a logical")
  }
  if(is.null(data.structure[[i]]$regional.mu))
  {
    data.structure[[i]]$regional.mu=0
  }
  data.structure[[i]]$regional.mu=as.integer(data.structure[[i]]$regional.mu)

  if(!is.null(data.structure[[i]]$regional.lin.time))
  {
    if(typeof(data.structure[[i]]$regional.lin.time)!="logical")
      stop("Regional linear time indicator, 'regional.lin.time', must be a logical")
  }
  if(is.null(data.structure[[i]]$regional.lin.time))
  {
    data.structure[[i]]$regional.lin.time=0
  }
  data.structure[[i]]$regional.lin.time=as.integer(data.structure[[i]]$regional.lin.time)
  
  if(!is.null(data.structure[[i]]$regional.pull))
  {
    if(typeof(data.structure[[i]]$regional.pull)!="integer" & typeof(data.structure[[i]]$regional.pull)!="double" & typeof(data.structure[[i]]$regional.pull)!="numeric")
      stop("Regional pull specification, 'regional pull', must be a vector of integers specifying layers")
  }
  if(is.null(data.structure[[i]]$regional.pull))
  {
    data.structure[[i]]$regional.pull=0
  }
  data.structure[[i]]$regional.pull=as.integer(data.structure[[i]]$regional.pull)
  if(length(data.structure[[i]]$regional.pull)>1 | data.structure[[i]]$regional.pull[1]>0)
  {
    if(sum(data.structure[[i]]$regional.pull>data.structure[[i]]$numlayers | data.structure[[i]]$regional.pull<0)>0)
      stop("Regional pull layer specification out of range!")
  }

  if(!is.null(data.structure[[i]]$regional.sigma))
  {
    if(typeof(data.structure[[i]]$regional.sigma)!="integer" & typeof(data.structure[[i]]$regional.sigma)!="double" & typeof(data.structure[[i]]$regional.sigma)!="numeric")
      stop("Regional sigma specification, 'regional sigma', must be a vector of integers specifying layers")
  }
  if(is.null(data.structure[[i]]$regional.sigma))
  {
    data.structure[[i]]$regional.sigma=0
  }
  data.structure[[i]]$regional.sigma=as.integer(data.structure[[i]]$regional.sigma)
  if(length(data.structure[[i]]$regional.sigma)>1 | data.structure[[i]]$regional.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$regional.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$regional.sigma<0)>0)
      stop("Rregional sigma layer specification out of range!")
  }
  
  if(!is.null(data.structure[[i]]$correlated.sigma))
  {
    if(typeof(data.structure[[i]]$correlated.sigma)!="integer" & typeof(data.structure[[i]]$correlated.sigma)!="double" & typeof(data.structure[[i]]$correlated.sigma)!="numeric")
      stop("Correlated sigma indicator, 'correlated.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$correlated.sigma))
  {
    data.structure[[i]]$correlated.sigma=0
  }
  data.structure[[i]]$correlated.sigma=as.integer(data.structure[[i]]$correlated.sigma)
  if(length(data.structure[[i]]$correlated.sigma)>1 | data.structure[[i]]$correlated.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$correlated.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$correlated.sigma<0)>0)
      stop("Correlated sigma specification out of range!")
  }

  if(!is.null(data.structure[[i]]$pairwise.correlated.sigma))
  {
    if(typeof(data.structure[[i]]$pairwise.correlated.sigma)!="integer" & typeof(data.structure[[i]]$pairwise.correlated.sigma)!="double" & typeof(data.structure[[i]]$pairwise.correlated.sigma)!="numeric")
      stop("Pairwise correlated sigma indicator, 'pairwise.correlated.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$pairwise.correlated.sigma))
  {
    data.structure[[i]]$pairwise.correlated.sigma=0
  }
  data.structure[[i]]$pairwise.correlated.sigma=as.integer(data.structure[[i]]$pairwise.correlated.sigma)
  if(length(data.structure[[i]]$pairwise.correlated.sigma)>1 | data.structure[[i]]$pairwise.correlated.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$pairwise.correlated.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$pairwise.correlated.sigma<0)>0)
      stop("Pairwise correlated sigma specification out of range!")
  }
  
  if(!is.null(data.structure[[i]]$one.dim.sigma))
  {
    if(typeof(data.structure[[i]]$one.dim.sigma)!="integer" & typeof(data.structure[[i]]$one.dim.sigma)!="double" & typeof(data.structure[[i]]$one.dim.sigma)!="numeric")
      stop("1D sigma indicator, 'one.dim.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$one.dim.sigma))
  {
    data.structure[[i]]$one.dim.sigma=0
  }
  data.structure[[i]]$one.dim.sigma=as.integer(data.structure[[i]]$one.dim.sigma)
  if(length(data.structure[[i]]$one.dim.sigma)>1 | data.structure[[i]]$one.dim.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$one.dim.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$one.dim.sigma<0)>0)
      stop("1D sigma specification out of range!")
  }

  if(!is.null(data.structure[[i]]$grouping.sigma))
  {
    if(typeof(data.structure[[i]]$grouping.sigma)!="integer" & typeof(data.structure[[i]]$grouping.sigma)!="double" & typeof(data.structure[[i]]$grouping.sigma)!="numeric")
      stop("Sigma grouping indicator, 'grouping.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$grouping.sigma))
  {
    data.structure[[i]]$grouping.sigma=0
  }
  data.structure[[i]]$grouping.sigma=as.integer(data.structure[[i]]$grouping.sigma)
  if(length(data.structure[[i]]$grouping.sigma)>1 | data.structure[[i]]$grouping.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$grouping.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$grouping.sigma<0)>0)
      stop("Sigma grouping specification out of range!")
  }

  if(!is.null(data.structure[[i]]$remove.sigma))
  {
    if(typeof(data.structure[[i]]$remove.sigma)!="integer" & typeof(data.structure[[i]]$remove.sigma)!="double" & typeof(data.structure[[i]]$remove.sigma)!="numeric")
      stop("Sigma correlation removal indicator, 'remove.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$remove.sigma))
  {
    data.structure[[i]]$remove.sigma=0
  }
  data.structure[[i]]$remove.sigma=as.integer(data.structure[[i]]$remove.sigma)
  if(length(data.structure[[i]]$remove.sigma)>1 | data.structure[[i]]$remove.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$remove.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$remove.sigma<0)>0)
      stop("Sigma correlation removal specification out of range!")
  }

  if(!is.null(data.structure[[i]]$differentiate.sigma))
  {
    if(typeof(data.structure[[i]]$differentiate.sigma)!="integer" & typeof(data.structure[[i]]$differentiate.sigma)!="double" & typeof(data.structure[[i]]$differentiate.sigma)!="numeric")
      stop("Differentiated sigma indicator, 'differentiate.sigma', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$differentiate.sigma))
  {
    data.structure[[i]]$differentiate.sigma=0
  }
  data.structure[[i]]$differentiate.sigma=as.integer(data.structure[[i]]$differentiate.sigma)
  if(length(data.structure[[i]]$differentiate.sigma)>1 | data.structure[[i]]$differentiate.sigma[1]>0)
  {
    if(sum(data.structure[[i]]$differentiate.sigma>data.structure[[i]]$numlayers | data.structure[[i]]$differentiate.sigma<0)>0)
      stop("Differentiated sigma specification out of range!")
  }

  if(!is.null(data.structure[[i]]$differentiate.pull))
  {
    if(typeof(data.structure[[i]]$differentiate.pull)!="integer" & typeof(data.structure[[i]]$differentiate.pull)!="double" & typeof(data.structure[[i]]$differentiate.pull)!="numeric")
      stop("Differentiated pull indicator, 'differentiate.pull', must be an integer vector specifying layers")
  }
  if(is.null(data.structure[[i]]$differentiate.pull))
  {
    data.structure[[i]]$differentiate.pull=0
  }
  data.structure[[i]]$differentiate.pull=as.integer(data.structure[[i]]$differentiate.pull)
  if(length(data.structure[[i]]$differentiate.pull)>1 | data.structure[[i]]$differentiate.pull[1]>0)
  {
    if(sum(data.structure[[i]]$differentiate.pull>data.structure[[i]]$numlayers | data.structure[[i]]$differentiate.pull<0)>0)
      stop("Differentiated pull specification out of range!")
  }

  if(!is.null(data.structure[[i]]$differentiate.mu))
  {
    if(typeof(data.structure[[i]]$differentiate.mu)!="logical")
      stop("Differentiated expectancy indicator, 'differentiate.mu', must be a logical")
  }
  if(is.null(data.structure[[i]]$differentiate.mu))
  {
    data.structure[[i]]$differentiate.mu=0
  }
  data.structure[[i]]$differentiate.mu=as.integer(data.structure[[i]]$differentiate.mu)
  

  if(!is.null(data.structure[[i]]$differentiate.lin.time))
  {
    if(typeof(data.structure[[i]]$differentiate.lin.time)!="logical")
      stop("Differentiated linear time indicator, 'differentiate.lin.time', must be a logical")
  }
  if(is.null(data.structure[[i]]$differentiate.lin.time))
  {
    data.structure[[i]]$differentiate.lin.time=0
  }
  data.structure[[i]]$differentiate.lin.time=as.integer(data.structure[[i]]$differentiate.lin.time)
  
  if(!is.null(data.structure[[i]]$init.0))
  {
    if(typeof(data.structure[[i]]$init.0)!="logical")
      stop("Initial value treatement indicator, 'init.0', must be a logical")
  }
  if(is.null(data.structure[[i]]$init.0))
  {
    data.structure[[i]]$init.0=0
  }
  data.structure[[i]]$init.0=as.integer(data.structure[[i]]$init.0)

  if(!is.null(data.structure[[i]]$init.time))
  {
    if(typeof(data.structure[[i]]$init.time)!="numeric" & typeof(data.structure[[i]]$init.time)!="double" & 
       typeof(data.structure[[i]]$init.time)!="POSIXct" )
      stop("Initial value at given time treatement indicator, 'init.time', must be a numeric or POSIXct")
    if(length(data.structure[[i]]$init)>1)
      stop("Initial value at given time treatement indicator, 'init.time', cannot be a vector")

    data.structure[[i]]$init.datetime=0
    if(typeof(data.structure[[i]]$init.time)=="POSIXct" )
      data.structure[[i]]$init.datetime=1
  }
  if(is.null(data.structure[[i]]$init.time))
  {
    data.structure[[i]]$init.time=-10000000
    data.structure[[i]]$init.datetime=0
  }
  data.structure[[i]]$init.time=as.numeric(data.structure[[i]]$init.time)
  data.structure[[i]]$init.datetime=as.integer(data.structure[[i]]$init.datetime)
  
  if(!is.null(data.structure[[i]]$init.same.sites))
  {
    if(typeof(data.structure[[i]]$init.same.sites)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.sites', must be a logical")
  }
  if(is.null(data.structure[[i]]$init.same.sites))
  {
    data.structure[[i]]$init.same.sites=0
  }
  data.structure[[i]]$init.same.sites=as.integer(data.structure[[i]]$init.same.sites)
  
  if(!is.null(data.structure[[i]]$init.same.layers))
  {
    if(typeof(data.structure[[i]]$init.same.layers)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.layers', must be a logical")
  }
  if(is.null(data.structure[[i]]$init.same.layers))
  {
    data.structure[[i]]$init.same.layers=0
  }
  data.structure[[i]]$init.same.layers=as.integer(data.structure[[i]]$init.same.layers)
  
  if(!is.null(data.structure[[i]]$init.specified))
  {
    if(typeof(data.structure[[i]]$init.specified)!="numeric" & typeof(data.structure[[i]]$init.specified)!="double")
      stop("Specified initial value, 'init.specified', must be a numeric. If time has been specified with POSIXct, cast this as numeric here.")
    if(length(data.structure[[i]]$init.specified)!=2)
      stop("Specified initial value, 'init.specified', must be a vector with 2 values, 'time' and 'value'. If time has been specified with POSIXct, cast this as numeric")
  }
  if(is.null(data.structure[[i]]$init.specified))
  {
    data.structure[[i]]$init.specified=c(-10000000,-10000000)
  }
  data.structure[[i]]$init.specified=as.numeric(data.structure[[i]]$init.specified)
  
  if(!is.null(data.structure[[i]]$allow.pos.pull))
  {
    if(typeof(data.structure[[i]]$allow.pos.pull)!="logical")
      stop("Indicator for allowing positive (unstable) pulls, 'allow.pos.pull', must be a logical")
  }
  if(is.null(data.structure[[i]]$allow.pos.pull))
  {
    data.structure[[i]]$allow.pos.pull=0
  }
  data.structure[[i]]$allow.pos.pull=as.integer(data.structure[[i]]$allow.pos.pull)
  
  if(!is.null(data.structure[[i]]$period))
  {
    if(typeof(data.structure[[i]]$period)!="numeric" & typeof(data.structure[[i]]$period)!="double")
      stop("Periodial trigonometric regressor indicator, 'period', must be a numeric vector")
  }
  if(is.null(data.structure[[i]]$period))
  {
    data.structure[[i]]$period=0
  }
  data.structure[[i]]$period=as.numeric(data.structure[[i]]$period)
  
  if(!is.null(data.structure[[i]]$prior))
  {
    if(is.null(data.structure[[i]]$prior$mu))
      stop("If prior is given, it must contain the 95% prior credibility for expected value, given as element 'mu'!")
    if(length(data.structure[[i]]$prior$mu)!=2)
      stop("Prior for 'mu' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'mu'!")
    data.structure[[i]]$prior$mu=as.numeric(data.structure[[i]]$prior$mu)

    if(is.null(data.structure[[i]]$prior$dt))
      stop("If prior is given, it must contain the 95% prior credibility for the characteristic time, given as element 'dt'!")
    if(length(data.structure[[i]]$prior$dt)!=2)
      stop("Prior for 'dt' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'dt'!")
    data.structure[[i]]$prior$dt=as.numeric(data.structure[[i]]$prior$dt)
    if(use.half.lives)
      data.structure[[i]]$prior$dt=data.structure[[i]]$prior$dt/log(2.0)

    if(is.null(data.structure[[i]]$prior$s))
      stop("If prior is given, it must contain the 95% prior credibility for the stochastic contributions, given as element 's'!")
    if(length(data.structure[[i]]$prior$s)!=2)
      stop("Prior for 's' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 's'!")
    data.structure[[i]]$prior$s=as.numeric(data.structure[[i]]$prior$s)
    
    if(is.null(data.structure[[i]]$prior$init))
      data.structure[[i]]$prior$init=data.structure[[i]]$prior$mu
    if(length(data.structure[[i]]$prior$init)!=2)
      stop("Prior for 'init' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'init'!")
    data.structure[[i]]$prior$init=as.numeric(data.structure[[i]]$prior$init)
    
    if(is.null(data.structure[[i]]$prior$lin) & data.structure[[i]]$lin.time!=0)
      stop("If prior is given and linear time dependency is used, it must contain the 95% prior credibility for the linear time dependency, given as element 'lin'!")
    if(is.null(data.structure[[i]]$prior$lin))
      data.structure[[i]]$prior$lin=c(-1,1) # doesn't really matter since it's not going to be used
    if(length(data.structure[[i]]$prior$lin)!=2)
      stop("Prior for 'lin' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'lin'!")
    data.structure[[i]]$prior$lin=as.numeric(data.structure[[i]]$prior$lin)
    
    if(is.null(data.structure[[i]]$prior$beta))
      data.structure[[i]]$prior$beta=c(-1,1)
    if(length(data.structure[[i]]$prior$beta)!=2)
      stop("Prior for 'beta' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'beta'!")
    data.structure[[i]]$prior$beta=as.numeric(data.structure[[i]]$prior$beta)
    
    if(is.null(data.structure[[i]]$prior$obs) & is.null(data.structure[[i]]$std.dev))
      stop("If prior is given and measurement-wise observational standard deviation is not given, the prior must contain the 95% prior credibility for the observational standard deviation, given as element 'obs'!")
    if(is.null(data.structure[[i]]$prior$obs))
      data.structure[[i]]$prior$obs=c(0.01,1) # doesn't really matter since it's not going to be used
    if(length(data.structure[[i]]$prior$obs)!=2)
      stop("Prior for 'obs' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'obs'!")
    data.structure[[i]]$prior$obs=as.numeric(data.structure[[i]]$prior$obs)
    
    if(is.null(data.structure[[i]]$prior$islog))
      data.structure[[i]]$prior$islog=0
    data.structure[[i]]$prior$islog=as.integer(data.structure[[i]]$prior$islog)
  }
  if(is.null(data.structure[[i]]$prior))
    data.structure[[i]]$prior=layer.standard.prior
  
 }
 
 if(!is.null(causal))
 {
   if(class(causal)!="matrix")
     stop("The connection matrix 'causal' must be an integer matrix")
   if(dim(causal)[1]!=4)
     stop("The connection matrix 'causal' must have exactly four rows (representing cause and effect)!")
   if(typeof(causal)=="logical" & dim(causal)[2]>0)
     stop("The connection matrix 'causal' must have integer elements (represennting the number of each cause and effect series")
   if((typeof(causal)=="numeric" | typeof(causal)=="double") & 
      sum(round(causal)!=causal)>0)
     stop("The connection matrix 'causal' must have integer elements (represennting the number of each cause and effect series")
   if(typeof(causal)=="numeric" & typeof(causal)=="double" &  typeof(causal)=="integer")
     stop("The connection matrix 'causal' must have integer elements (representing the number of each cause and effect series")

   if(dim(causal)[2]>0)
    for(i in 1:dim(causal)[2])
     {
       if(causal[1,i]>n)
         stop(sprintf("Causal 'from' series %d > number of series, %d", 
                      causal[1,i],n))
       if(causal[1,i]<1)
         stop(sprintf("Causal 'from' series %d < 0!",causal[1,i]))
       if(causal[2,i]>data.structure[[causal[1,i]]]$numlayers)
         stop(sprintf("Causal 'from' series %d layer specification, %d > number of layers for that series, %d", 
              causal[1,i],causal[2,i],data.structure[[causal[1,i]]]$numlayers))
       if(causal[2,i]<0)
         stop(sprintf("Causal 'from' series %d layer specification, %d < 0!", 
              causal[1,i],causal[2,i]))
       if(causal[3,i]>n)
         stop(sprintf("Causal 'to' series %d > number of series, %d", 
                      causal[3,i],n))
       if(causal[3,i]<1)
         stop(sprintf("Causal 'to' series %d < 0!",causal[3,i]))
       if(causal[4,i]>data.structure[[causal[3,i]]]$numlayers)
         stop(sprintf("Causal 'to' series %d layer specification, %d > number of layers for that series, %d", 
              causal[3,i],causal[4,i],data.structure[[causal[3,i]]]$numlayers))
       if(causal[4,i]<0)
         stop(sprintf("Causal 'from' series %d layer specification, %d < 0!", 
              causal[3,i],causal[4,i]))
     }
 }
 if(is.null(causal))
   causal=as.matrix(array(NA,c(2,0)))
 
 if(!is.null(causal.symmetric))
 {
   if(class(causal.symmetric)!="matrix")
     stop("The connection matrix 'causal.symmetric' must be an integer matrix")
   if(dim(causal.symmetric)[1]!=4)
     stop("The connection matrix 'causal.symmetric' must have exactly four rows (representing cause and effect)!")
   if(typeof(causal.symmetric)=="logical" & dim(causal.symmetric)[2]>0)
     stop("The connection matrix 'causal.symmetric' must have integer elements (represennting the number of each cause and effect series")
   if((typeof(causal.symmetric)=="numeric" | typeof(causal.symmetric)=="double") & 
      sum(round(causal.symmetric)!=causal.symmetric)>0)
     stop("The connection matrix 'causal.symmetric' must have integer elements (represennting the number of each cause and effect series")
   if(typeof(causal.symmetric)=="numeric" & typeof(causal.symmetric)=="double" &  typeof(causal.symmetric)=="integer")
     stop("The connection matrix 'causal.symmetric' must have integer elements (representing the number of each cause and effect series")

   if(dim(causal.symmetric)[2]>0)
    for(i in 1:dim(causal.symmetric)[2])
     {
       if(causal.symmetric[1,i]>n)
         stop(sprintf("Causal.Symmetric 'from' series %d > number of series, %d", 
                      causal.symmetric[1,i],n))
       if(causal.symmetric[1,i]<1)
         stop(sprintf("Causal.Symmetric 'from' series %d < 0!",causal.symmetric[1,i]))
       if(causal.symmetric[2,i]>data.structure[[causal.symmetric[1,i]]]$numlayers)
         stop(sprintf("Causal.Symmetric 'from' series %d layer specification, %d > number of layers for that series, %d", 
              causal.symmetric[1,i],causal.symmetric[2,i],data.structure[[causal.symmetric[1,i]]]$numlayers))
       if(causal.symmetric[2,i]<0)
         stop(sprintf("Causal.Symmetric 'from' series %d layer specification, %d < 0!", 
              causal.symmetric[1,i],causal.symmetric[2,i]))
       if(causal.symmetric[3,i]>n)
         stop(sprintf("Causal.Symmetric 'to' series %d > number of series, %d", 
                      causal.symmetric[3,i],n))
       if(causal.symmetric[3,i]<1)
         stop(sprintf("Causal.Symmetric 'to' series %d < 0!",causal.symmetric[3,i]))
       if(causal.symmetric[4,i]>data.structure[[causal.symmetric[3,i]]]$numlayers)
         stop(sprintf("Causal.Symmetric 'to' series %d layer specification, %d > number of layers for that series, %d", 
              causal.symmetric[3,i],causal.symmetric[4,i],data.structure[[causal.symmetric[3,i]]]$numlayers))
       if(causal.symmetric[4,i]<0)
         stop(sprintf("Causal.Symmetric 'from' series %d layer specification, %d < 0!", 
              causal.symmetric[3,i],causal.symmetric[4,i]))
     }
 }
 if(is.null(causal.symmetric))
   causal.symmetric=as.matrix(array(NA,c(2,0)))
 
 
 if(!is.null(corr))
 {
   if(class(corr)!="matrix")
     stop("The connection matrix 'corr' must be an integer matrix")
   if(dim(corr)[1]!=4)
     stop("The connection matrix 'corr' must have exactly four rows (representing cause and effect)!")
   if(typeof(corr)=="logical" & dim(corr)[2]>0)
     stop("The connection matrix 'corr' must have integer elements (represennting the number of each cause and effect series")
   if((typeof(corr)=="numeric" | typeof(corr)=="double") & 
      sum(round(corr)!=corr)>0)
     stop("The connection matrix 'corr' must have integer elements (represennting the number of each cause and effect series")
   if(typeof(corr)=="numeric" & typeof(corr)=="double" &  typeof(corr)=="integer")
     stop("The connection matrix 'corr' must have integer elements (representing the number of each cause and effect series")

   if(dim(corr)[2]>0)
    for(i in 1:dim(corr)[2])
     {
       if(corr[1,i]>n)
         stop(sprintf("Corr 'from' series %d > number of series, %d", 
                      corr[1,i],n))
       if(corr[1,i]<1)
         stop(sprintf("Corr 'from' series %d < 0!",corr[1,i]))
       if(corr[2,i]>data.structure[[corr[1,i]]]$numlayers)
         stop(sprintf("Corr 'from' series %d layer specification, %d > number of layers for that series, %d", 
              corr[1,i],corr[2,i],data.structure[[corr[1,i]]]$numlayers))
       if(corr[2,i]<0)
         stop(sprintf("Corr 'from' series %d layer specification, %d < 0!", 
              corr[1,i],corr[2,i]))
       if(corr[3,i]>n)
         stop(sprintf("Corr 'to' series %d > number of series, %d", 
                      corr[3,i],n))
       if(corr[3,i]<1)
         stop(sprintf("Corr 'to' series %d < 0!",corr[3,i]))
       if(corr[4,i]>data.structure[[corr[3,i]]]$numlayers)
         stop(sprintf("Corr 'to' series %d layer specification, %d > number of layers for that series, %d", 
              corr[3,i],corr[4,i],data.structure[[corr[3,i]]]$numlayers))
       if(corr[4,i]<0)
         stop(sprintf("Corr 'from' series %d layer specification, %d < 0!", 
              corr[3,i],corr[4,i]))
     }
 }
 if(is.null(corr))
   corr=as.matrix(array(NA,c(2,0)))
 
 if(typeof(silent.mode)!="logical" & typeof(silent.mode)!="integer")
   stop("Option 'silent.mode' must be a logical!") 

 if(typeof(do.model.likelihood)!="logical" & typeof(do.model.likelihood)!="integer")
   stop("Option 'do.model.likelihood' must be a logical!") 

 if(typeof(do.maximum.likelihood)!="logical" & typeof(do.maximum.likelihood)!="integer")
   stop("Option 'do.maximum.likelihood' must be a logical!") 
   
 if(typeof(talkative.burnin)!="logical" & typeof(talkative.burnin)!="integer")
   stop("Option 'talkative.burnin' must be a logical!") 

 if(typeof(talkative.likelihood)!="logical" & typeof(talkative.likelihood)!="integer")
   stop("Option 'talkative.likelihood' must be a logical!") 

 if(typeof(use.stationary.stdev)!="logical" & typeof(use.stationary.stdev)!="integer")
   stop("Option 'use.stationary.stdev' must be a logical!") 
 
 if(use.stationary.stdev & do.maximum.likelihood)
   stop("Options 'use.stationary.stdev' and 'do.maximum.likelihood' in combination is not implemented, unfortunately!")
  
 if(typeof(use.half.lives)!="logical" & typeof(use.half.lives)!="integer")
   stop("Option 'use.half.lives' must be a logical!") 

 if(typeof(mcmc)!="logical" & typeof(mcmc)!="integer")
   stop("Option 'mcmc' must be a logical!") 

 if(typeof(num.MCMC)!="integer" & typeof(num.MCMC)!="numeric" & 
   typeof(num.MCMC)!="double" )
   stop("Option 'num.MCMC' must be an integer!") 
 num.MCMC=as.integer(num.MCMC)

 if(num.MCMC<1)
  stop("Option 'num.MCMC' must be 1 or larger  (preferrably much larger)!")

 if(typeof(spacing)!="integer" & typeof(spacing)!="numeric" & 
   typeof(spacing)!="double" )
   stop("Option 'spacing' must be an integer!") 
 spacing=as.integer(spacing)

 if(spacing<1)
  stop("Option 'spacing' must be 1 or larger!")

 if(typeof(burnin)!="integer" & typeof(burnin)!="numeric" & 
   typeof(burnin)!="double" )
   stop("Option 'burnin' must be an integer!")
 burnin=as.integer(burnin) 

 if(burnin<0)
  stop("Option 'num.temp' must be 0 or larger!")

 if(typeof(num.temp)!="integer" & typeof(num.temp)!="numeric" & 
   typeof(num.temp)!="double" )
   stop("Option 'num.temp' must be an integer!") 
 num.temp=as.integer(num.temp) 

 if(num.temp<1)
  stop("Option 'num.temp' must be 1 or larger!")

 if(typeof(id.strategy)!="integer" & typeof(id.strategy)!="numeric" & 
   typeof(id.strategy)!="double" )
   stop("Option 'id.strategy' must be an integer!") 
 id.strategy=as.integer(id.strategy) 
 
 if(id.strategy<0 | id.strategy>4)
  stop("Option 'id.strategy' must be an integer between 0 and 4!")
 
 if(typeof(T.ground)!="integer" & typeof(T.ground)!="numeric" & 
    typeof(T.ground)!="double")
  stop("Option 'T.ground' must be a numeric!")
 T.ground=as.numeric(T.ground)
 
 if(T.ground<=1.0)
  stop("Option 'T.ground' must be larger than 1.0!")
 
 if(typeof(maximum.likelihood.numstart)!="integer" & 
   typeof(maximum.likelihood.numstart)!="numeric" & 
   typeof(maximum.likelihood.numstart)!="double" )
   stop("Option 'maximum.likelihood.numstart' must be an integer!") 
 maximum.likelihood.numstart=as.integer(maximum.likelihood.numstart)
 
 if(maximum.likelihood.numstart<1)
  stop("Option 'maximum.likelihood.numstart' must be 1 or larger!")
 
 if(is.null(smoothing.specs))
 {
   smoothing.specs=
     list(do.smoothing=FALSE,smoothing.time.diff=0,smoothing.start=NULL,smoothing.end=NULL,
          num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE)
 }
 if(!is.null(smoothing.specs))
 {
   if(typeof(smoothing.specs)!="list")
   {	
     stop("\"smoothing.specs\" must be a list with minimum the content \"do.smoothing\" if \"do.smoothing\" is set to \"FALSE\"  and minimum \"do.smoothing\" and \"smoothing.time.diff\" if \"do.smoothing\" is set to \"TRUE\". \"smoothing.start\" and \"smoothing.end\" is optional.") 
   }
   
   if(is.null(smoothing.specs$do.smoothing))
   {	
     stop("\"smoothing.specs\" must be a list with minimum the content \"do.smoothing\" if \"do.smoothing\" is set to \"FALSE\" and minimum \"do.smoothing\", \"smoothing.time.diff\" and \"num.smooth.per.mcmc\" if \"do.smoothing\" is set to \"TRUE\". \"smoothing.start\" and \"smoothing.end\" is optional.") 
   }

   if(smoothing.specs$do.smoothing & (is.null(smoothing.specs$smoothing.time.diff) |
      is.null(smoothing.specs$num.smooth.per.mcmc)))
   {	
     stop("\"smoothing.specs\" must be a list with minimum the content \"do.smoothing\" if \"do.smoothing\" is set to \"FALSE\" and minimum \"do.smoothing\", \"smoothing.time.diff\" and \"num.smooth.per.mcmc\" if \"do.smoothing\" is set to \"TRUE\". \"smoothing.start\" and \"smoothing.end\" is optional.") 
   }
   
   if(smoothing.specs$do.smoothing & 
      (typeof(smoothing.specs$smoothing.time.diff)!="numeric" &
       typeof(smoothing.specs$smoothing.time.diff)!="integer" &
       typeof(smoothing.specs$smoothing.time.diff)!="double"))
   {	
     stop("\"smoothing.specs$smoothing.time.diff\" must a numeric!") 
   }   
   smoothing.specs$smoothing.time.diff=as.numeric(smoothing.specs$smoothing.time.diff)

   if(smoothing.specs$do.smoothing & 
      (typeof(smoothing.specs$num.smooth.per.mcmc)!="numeric" &
       typeof(smoothing.specs$num.smooth.per.mcmc)!="integer" &
       typeof(smoothing.specs$num.smooth.per.mcmc)!="double"))
   {	
     stop("\"smoothing.specs$num.smooth.per.mcmc\" must an integer!") 
   }
   smoothing.specs$num.smooth.per.mcmc=as.integer(smoothing.specs$num.smooth.per.mcmc) 
   if(smoothing.specs$num.smooth.per.mcmc<=0)
   {
     stop("\"smoothing.specs$num.smooth.per.mcmc\" must be a positive integer!")
   }
   
   if(is.null(smoothing.specs$do.return.smoothing.samples))
     smoothing.specs$do.return.smoothing.samples=smoothing.specs$do.smoothing
   
   if(smoothing.specs$do.smoothing &
     (typeof(smoothing.specs$do.return.smoothing.samples)!="logical"))
   {
     stop("\"smoothing.specs$do.return.smoothing.samples\" must be a logical!")
   }
   smoothing.specs$do.return.smoothing.samples=as.integer(smoothing.specs$do.return.smoothing.samples)
   
   if(smoothing.specs$do.smoothing & 
      is.null(smoothing.specs$smoothing.end)!=is.null(smoothing.specs$smoothing.start))
   {
     stop("If \"smoothing.specs$smoothing.start\" is given, then \"smoothing.specs$smoothing.end\" must also be given, and vice versa!")
   }

   smoothing.specs$start.end.given=0
   smoothing.specs$start.end.datetime=0
   if(smoothing.specs$do.smoothing & !is.null(smoothing.specs$smoothing.start) &
      !is.null(smoothing.specs$smoothing.end))
   {
     smoothing.specs$start.end.given=1
     if(typeof(smoothing.specs$smoothing.start)!=typeof(smoothing.specs$smoothing.end))
     {
       stop("Type mistmatch between \"smoothing.specs$smoothing.start\" and \"smoothing.specs$smoothing.end\"!")
     }
   }
   
   if(smoothing.specs$do.smoothing & !is.null(smoothing.specs$smoothing.start))
   {
     if(typeof(smoothing.specs$smoothing.start)!="integer" &
        typeof(smoothing.specs$smoothing.start)!="numeric" &
        typeof(smoothing.specs$smoothing.start)!="double" &
        typeof(smoothing.specs$smoothing.start)!="POSIXct")
	{
	  stop("\"smoothing.specs$smoothing.start\" must be a numeric or (if the time series contain date-times) a POSIXct")
	}

     if(typeof(smoothing.specs$smoothing.start)!="POSIXct")
     {
       smoothing.specs$smoothing.start=as.numeric(smoothing.specs$smoothing.start)
       smoothing.specs$smoothing.start.dt=NULL
     }
     if(typeof(smoothing.specs$smoothing.start)=="POSIXct")
     {
       smoothing.specs$start.end.datetime=1
       smoothing.specs$smoothing.start.dt=smoothing.specs$smoothing.start
       smoothing.specs$smoothing.start=NULL
     }
   }
      
   if(smoothing.specs$do.smoothing & !is.null(smoothing.specs$smoothing.end))
   {
     if(typeof(smoothing.specs$smoothing.end)!="integer" &
        typeof(smoothing.specs$smoothing.end)!="numeric" &
        typeof(smoothing.specs$smoothing.end)!="double" &
        typeof(smoothing.specs$smoothing.end)!="POSIXct")
	{
	  stop("\"smoothing.specs$smoothing.end\" must be a numeric or (if the time series contain date-times) a POSIXct")
	}

     if(typeof(smoothing.specs$smoothing.end)!="POSIXct")
     {
       smoothing.specs$smoothing.end=as.numeric(smoothing.specs$smoothing.end)
       smoothing.specs$smoothing.end.dt=NULL
     }
     if(typeof(smoothing.specs$smoothing.end)=="POSIXct")
     {
       smoothing.specs$start.end.datetime=1
       smoothing.specs$smoothing.end.dt=smoothing.specs$smoothing.end
       smoothing.specs$smoothing.end=NULL
     }
   }
 }


 if(is.null(realization.specs))
 {
   realization.specs=list(do.realizations=FALSE,num.realizations=100,strategy="N")
 }
 if(!is.null(realization.specs))
 {
   if(typeof(realization.specs)!="list")
   {	
     stop("\"realization.specs\" must be a list with minimum the content \"do.realizations\" if \"do.realizations\" is set to \"FALSE\" and minimum \"do.realizations\" and \"num.realizations\" if \"do.realizations\" is set to \"TRUE\". \"strategy\" is optional.") 
   }
   
   if(is.null(realization.specs$do.realizations) | 
     (typeof(realization.specs$do.realizations)!="logical"))
   {
     stop("\"realization.specs\" must be a list with minimum the content \"do.realizations\" if \"do.realizations\" is set to \"FALSE\" and minimum \"do.realizations\" and \"num.realizations\" if \"do.realizations\" is set to \"TRUE\". \"strategy\", \"realization.start\" and \"realization.end\" are optional.") 	
   }
   
   if(realization.specs$do.realizations & smoothing.specs$do.smoothing)
   {
     stop("Fetching smoothing results and realizations cannot be done on the same run!")
   }
   
   if(realization.specs$do.realizations & is.null(realization.specs$num.realizations))
   {	
     stop("\"realization.specs\" must be a list with minimum the content \"do.realizations\" if \"do.realizations\" is set to \"FALSE\" and minimum \"do.realizations\" and \"num.realizations\" if \"do.realizations\" is set to \"TRUE\". \"strategy\", \"realization.start\" and \"realization.end\" are optional.") 	
   }
   
   if(realization.specs$do.realizations & 
      (typeof(realization.specs$num.realizations)!="numeric" &
       typeof(realization.specs$num.realizations)!="integer" &
       typeof(realization.specs$num.realizations)!="double"))
   {	
     stop("\"realization.specs$num.realizations\" must an integer!") 
   }
   realization.specs$num.realizations=as.integer(realization.specs$num.realizations) 
   if(realization.specs$do.realizations & realization.specs$num.realizations<=0)
   {
     stop("\"realization.specs$num.realizations\" must be a positive integer!")
   }
   
   if(is.null(realization.specs$strategy))
     realization.specs$strategy="N"
   if(realization.specs$do.realizations & typeof(realization.specs$strategy)!="character")
   {
     stop("\"realization.specs$strategy\" must be a character array (string)")
   }
   if(realization.specs$do.realizations & 
      substr(realization.specs$strategy,1,1)!="N" &
      substr(realization.specs$strategy,1,1)!="A" &
      substr(realization.specs$strategy,1,1)!="D")
   {
     stop("Legal values for \"realization.specs$strategy\" are only \"N\" (None), \"A\" (Ascending) and \"D\" (Descending)")   
   }
   
   if(realization.specs$do.realizations & 
      (typeof(realization.specs$realization.time.diff)!="numeric" &
       typeof(realization.specs$realization.time.diff)!="integer" &
       typeof(realization.specs$realization.time.diff)!="double"))
   {	
     stop("\"realization.specs$realization.time.diff\" must a numeric!") 
   }   
   realization.specs$realization.time.diff=
         as.numeric(realization.specs$realization.time.diff)

   if(realization.specs$do.realizations & 
      is.null(realization.specs$realization.end)!=is.null(realization.specs$realization.start))
   {
     stop("If \"realization.specs$realization.start\" is given, then \"realization.specs$realization.end\" must also be given, and vice versa!")
   }
   
   realization.specs$start.end.given=0
   realization.specs$start.end.datetime=0
   if(realization.specs$do.realizations & !is.null(realization.specs$realization.start) &
      !is.null(realization.specs$realization.end))
   {
     realization.specs$start.end.given=1
     if(typeof(realization.specs$realization.start)!=typeof(realization.specs$realization.end))
     {
       stop("Type mistmatch between \"realization.specs$realization.start\" and \"realization.specs$realization.end\"!")
     }
   }
   
   if(realization.specs$do.realizations & !is.null(realization.specs$realization.start))
   {
     if(typeof(realization.specs$realization.start)!="integer" &
        typeof(realization.specs$realization.start)!="numeric" &
        typeof(realization.specs$realization.start)!="double" &
        typeof(realization.specs$realization.start)!="POSIXct")
	{
	  stop("\"realization.specs$realization.start\" must be a numeric or (if the time series contain date-times) a POSIXct")
	}
     
     if(typeof(realization.specs$realization.start)!="POSIXct")
     {
       realization.specs$realization.start=as.numeric(realization.specs$realization.start)
       realization.specs$realization.start.dt=NULL
     }
     if(typeof(realization.specs$realization.start)=="POSIXct")
     {
       realization.specs$start.end.datetime=1
       realization.specs$realization.start.dt=realization.specs$realization.start
       realization.specs$realization.start=NULL
     }
   }
      
   if(realization.specs$do.realizations & !is.null(realization.specs$realization.end))
   {
     if(typeof(realization.specs$realization.end)!="integer" &
        typeof(realization.specs$realization.end)!="numeric" &
        typeof(realization.specs$realization.end)!="double" &
        typeof(realization.specs$realization.end)!="POSIXct")
	{
	  stop("\"realization.specs$realization.end\" must be a numeric or (if the time series contain date-times) a POSIXct")
	}

     if(typeof(realization.specs$realization.end)!="POSIXct")
     {
       realization.specs$realization.end=as.numeric(realization.specs$realization.end)
       realization.specs$realization.end.dt=NULL
     }
     if(typeof(realization.specs$realization.end)=="POSIXct")
     {
       realization.specs$start.end.datetime=1
       realization.specs$realization.end.dt=realization.specs$realization.end
       realization.specs$realization.end=NULL
     }
   }

   realization.specs$do.realizations=as.integer(realization.specs$do.realizations)
 }
 
 out=.Call('layeranalyzer',data.structure,num.MCMC,burnin,spacing,num.temp,
      as.integer(do.model.likelihood),
      as.integer(do.maximum.likelihood),maximum.likelihood.numstart,
      as.integer(silent.mode),
      as.integer(talkative.burnin),as.integer(talkative.likelihood),
      id.strategy,as.integer(use.stationary.stdev),as.numeric(T.ground),
      #start.parameters=as.numeric(0),
      as.integer(use.half.lives),
      as.integer(mcmc), as.integer(causal), 
      as.integer(causal.symmetric),as.integer(corr), smoothing.specs, 
      realization.specs)
 
 out$data.structure=data.structure

 if(!is.null(out$mcmc))
 {
   out$mcmc=t(out$mcmc)
   colnames(out$mcmc)=out$parameter.names;
   out$mcmc=mcmc(out$mcmc, start=burnin, thin=1)
 }
 
 if(out$model.log.lik==-10000000)
  out$model.log.lik=NA
 	
 out$causal=causal
 out$causal.symmetric=causal.symmetric
 out$corr=corr

 class(out)="layered"
  
 return(out)
}
 
summary.layered=function(object, ...)
{
  input=list(...)
  
  if(is.null(object$ML.loglik))
  {
    mat=matrix(NA,nrow=length(object$parameter.names),ncol=4)
    for(i in 1:length(object$parameter.names))
     for(j in 1:4)
      mat[i,j]=as.numeric(object[[which(object$parameter.names[[i]]==names(object))]][j])
    dimnames(mat)=list(object$parameter.names, c("Mean","Median",
     "Lower 95%", 
     "Upper 95%"))
    
    ans=list(coefficients=round(mat,6), model.log.likelihood=object$model.log.lik)
    class(ans)="summary.layered"
    
    return(ans)
  }  
  if(!is.null(object$ML.loglik))
  {
    mat=matrix(NA,nrow=length(object$parameter.names),ncol=3)
    for(i in 1:length(object$parameter.names))
     for(j in 1:3)
      mat[i,j]=as.numeric(object[[i+10]][c(1,4,5)[j]])
    dimnames(mat)=list(object$parameter.names, c("ML estimate",
     "Bayesian Lower 95%", "Bayesian Upper 95%"))
    
    ans=list(coefficients=round(mat,6), ML.loglik=object$ML.loglik, 
      AIC=object$AIC,AICc=object$AICc,BIC=object$BIC)
    class(ans)="summary.layered"
    
    return(ans)
  }  
}

print.summary.layered=function(x, ...)
{
  input=list(...)
  
  if(is.null(x$ML.loglik))
  {
    print.srcref("Coefficients:")
    print(x$coefficients)
    if(!is.na(x$model.log.likelihood))
    {
      print.srcref("")
      print.srcref(sprintf("Model log-likelihood: %9.3f", x$model.log.likelihood))
    }
  }
  
  if(!is.null(x$ML.loglik))
  {
    print.srcref("Coefficients:")
    print(x$coefficients)
    
    print.srcref("")
    print.srcref(sprintf("Log-likelihood: %9.3f", x$ML.loglik))
    print.srcref(sprintf("AIC : %9.3f", x$AIC))
    print.srcref(sprintf("AICc: %9.3f", x$AICc))
    print.srcref(sprintf("BIC : %9.3f", x$BIC))
  }
}

compare.layered=function(...,p0=NULL,first.is.nullhypothesis=FALSE,
   ML.IC="AIC")
{
  input=list(...)

  if(class(input[[1]])!="layered")
    if(class(input[[1]][[1]])=="layered")
      input=input[[1]]
  
  n=length(input)
  
  all.ml=FALSE
  for(i in 1:n)
  {
   if(class(input[[i]])!="layered")
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
   
   if(i==1)
   {
     if(!is.null(input[[i]]$ML.loglik))
       all.ml=TRUE
   }
   
   if(i>1)
   {
     if((!all.ml && !is.null(input[[i]]$ML.loglik)) |
        (all.ml && is.null(input[[i]]$ML.loglik)))
         stop("Mix of Bayesian and ML objects!")
   } 
   
   if(!all.ml)
   {
     if(is.null(input[[i]]$model.log.lik) | is.na(input[[i]]$model.log.lik))
       stop(sprintf("Object number %d does not have a Bayesian model likelihood!",i))
   }
  }
  
  l=rep(0,n)
  for(i in 1:n)
  {
    if(!all.ml)
      l[i]=input[[i]]$model.log.lik
    if(all.ml)
    {
      if(ML.IC!="AIC" & ML.IC!="AICc" & ML.IC!="BIC")
       stop("ML information criterion weights, 'ML.IC' should be 'AIC', 'AICc' or 'BIC'")
      
      if(ML.IC=="AIC")
        l[i]=-0.5*input[[i]]$AIC
      if(ML.IC=="AICc")
        l[i]=-0.5*input[[i]]$AICc
      if(ML.IC=="BIC")
        l[i]=-0.5*input[[i]]$BIC
    }
  }
  
  if(!is.null(p0))
   if(length(p0)!=n)
     stop("p0 vector length does not match number of model inputs!")
  
  if(is.null(p0))
  {
   if(!first.is.nullhypothesis)
    p0=rep(1/n,n)
   if(first.is.nullhypothesis)
    p0=c(1/2,rep(1/2/(n-1),n-1))
  }
  
  maxl=max(l)
  p=p0*exp(l-maxl)/sum(p0*exp(l-maxl))*100
  
  outmatrix=cbind(l,p)
  if(is.null(names(input)))
    names(input)=sprintf("Model %3d",1:n)
  if(!all.ml)
    critname="log(lik)"
  if(all.ml)
  {
    if(ML.IC=="AIC")
      critname="weight=-0.5*AIC"
    if(ML.IC=="AICc")
      critname="weight=-0.5*AICc"
    if(ML.IC=="BIC")
      critname="weight=-0.5*BIC"
  }
  
  dimnames(outmatrix)=list(names(input),c(critname,"Post. Prob.(%)"))
  
  round(outmatrix,5)
}

nobs.layered=function(obj1,...)
{
  if(class(obj1)!="layered")
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  
 n=length(obj1$data.structure)
 len=0

 for(i in 1:n)
 {
  if(typeof(obj1$data.structure[[i]])!="list")
    stop("data.structure[[i]] must be a list, preferrably belonging to the 'layer.series.structure' class (as this class has the right elements).")
  
  if(is.null(obj1$data.structure[[i]]$timeseries))
    stop("data.structure[[i]] does not contain a time series (of type layer.data.series)!")
  
  if(is.null(obj1$data.structure[[i]]$timeseries$time))
    stop("No time points found in data material!")

  len=len+length(obj1$data.structure[[i]]$timeseries$time)
 }

 return(len)
}

logLik.layered=function(obj1,...)
{
  if(class(obj1)!="layered")
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  if(is.null(obj1$ML.loglik))
    stop("The 'layered' object need to be classically (ML) estimated!") 

  ll=obj1$ML.loglik

  df=length(obj1$parameter.names)
  attr(ll,"class")="logLik"
  attr(ll,"df")=df
  attr(ll,"nobs")=nobs.layered(obj1)

  return(ll)
}


anova.layered=function(obj1,...)
{
  if(class(obj1)!="layered")
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  
  input=list(...)
  n=length(input)
  
  if(is.null(obj1$ML.loglik))
    stop("All objects need to be classically (ML) estimated!") 
  for(i in 1:n)
  {
    if(class(input[[i]])!="layered")
     stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
    if(is.null(input[[i]]$ML.loglik))
      stop("All objects need to be classically (ML) estimated!") 
  }
  
  prev.df=length(obj1$parameter.names)
  prev.deviance=-2*obj1$ML.loglik
  
  res.df=rep(NA,n+1)
  res.dev=rep(NA,n+1)
  df=rep(NA,n+1)
  dev=rep(NA,n+1)
  pr=rep(NA,n+1)
  
  res.df[1]=nobs(obj1)-attr(logLik(obj1),"df")
  res.dev[1]=-2*as.numeric(logLik(obj1))
  
  if(n>0)
   for(i in 1:n)
   {
    if(class(input[[i]])!="layered")
     stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.") 
    
    res.df[1+i]=nobs(input[[i]])-attr(logLik(input[[i]]),"df")
    res.dev[1+i]=-2*as.numeric(logLik(input[[i]]))
    
    df[i+1]=res.df[i]-res.df[i+1]
    dev[i+1]=res.dev[i]-res.dev[i+1]
    if(df[i+1]>0)
      pr[i+1]=pchisq(dev[i+1],df[i+1],lower.tail=FALSE)
   }
  
  ret=list(res.df=res.df, res.dev=res.dev, df=df, dev=dev, pr=pr)
  names(ret)=c("Resid. Df","Resid. Dev", "Df", "Deviance", "Pr(>Chi)")
  class(ret)=c("anova","data.frame")
  attr(ret,"row.names")=1:(n+1)
  attr(ret,"heading")=c("Analysis of Deviance Table\n")
  
  return(ret)
}



n.combinations=function(length, base.number=2)
{
  mat=as.matrix(array(NA,c(base.number^length,length)))
  mat[1,]=rep(0,length)
  for(i in 2:(base.number^length))
  {
    mat[i,]=mat[i-1,]
    j=length
    mat[i,j]=(mat[i,j]+1)%%base.number
    while(mat[i,j]==0)
    {
      j=j-1
      mat[i,j]=(mat[i,j]+1)%%base.number
    }
  }
  return(mat)
}


traverse.standalone.layered=function(timeseries,
  max.layers=3, 
  talkative=FALSE, allow.one.feedback.loop=FALSE, 
  just.stationary=FALSE, no.rw=FALSE, regional.options=FALSE,    
  time.integrals.possible=FALSE, allow.deterministic.layers=TRUE,   
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, 
  use.half.lives=FALSE, mcmc=FALSE,
  prior=layer.standard.prior)
{
  if(is.null(timeseries$time))
    stop("No 'time' array in the incoming time series!")
    
  if(is.null(timeseries$value))
    stop("No 'value' array in the incoming time series!")
    
  if(length(timeseries$time)!=length(timeseries$value))
    stop("The 'time' and 'value' part of the time series must correspond, so the lengths must be the same!")
    
  if(length(timeseries$time)<2)
    stop("There must be at least two data points in time series!")
  
  if(typeof(timeseries$time)!="integer" & typeof(timeseries$time)!="numeric" & 
     typeof(timeseries$time)!="double" & typeof(timeseries$time)!="POSIXct")
    stop("The 'time' array in the time series must be numeric or date-time (POSIXct)!")
  timeseries$is.datetime=0
  if(typeof(timeseries$time)=="POSIXct")
    timeseries$is.datetime=1
  timeseries$time=as.numeric(timeseries$time)
  
  if(typeof(timeseries$value)!="integer" & typeof(timeseries$time)!="numeric" & typeof(timeseries$time)!="double")
    stop("The 'value' array in the time series must be numeric!")
  if(sum(is.na(timeseries$value))>0)
    timeseries$value[is.na(timeseries$value)]=-10000000
  timeseries$value=as.numeric(timeseries$value)
  
  if(!is.null(timeseries$num.meas.per.value))
  {
    if(typeof(timeseries$num.meas.per.value)!="integer" & 
      typeof(timeseries$num.meas.per.value)!="numeric" & 
      typeof(timeseries$num.meas.per.value)!="double")
      stop("The number of measurement array 'num.meas.per.value' array in the time series must be integer!")

    if(length(timeseries$num.meas.per.value)!=length(timeseries$value))
      stop("Number of measurements array 'num.meas.per.value' in the time series must correspond with the 'value' and 'time' arrays, so the lengths must be the same(4)!")

    if(sum(is.na(timeseries$num.meas.per.value))>0)
      timeseries$n[is.na(timeseries$num.meas.per.value)]=-10000000
    
    timeseries$num.meas.per.value=as.integer(timeseries$num.meas.per.value)
  }
  
  if(!is.null(timeseries$std.dev))
  {
    if(typeof(timeseries$std.dev)!="integer" & typeof(timeseries$std.dev)!="numeric" & 
      typeof(timeseries$std.dev)!="double")
      stop("The standard deviation array 'std.dev' array must be numeric!")

    if(length(timeseries$std.dev)!=length(timeseries$value))
      stop("Standard deviation array 'std.dev' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")

    if(sum(is.na(timeseries$std.dev))>0)
      timeseries$std.dev[is.na(timeseries$std.dev)]=-10000000

    timeseries$std.dev=as.numeric(timeseries$std.dev)
  }
  
  has.multiple.sites=F
  if(!is.null(timeseries$site))
  {
    if(typeof(timeseries$site)!="integer" & typeof(timeseries$site)!="numeric" & 
      typeof(timeseries$site)!="double")
      stop("The site array 'site' must be integer!")
      
    if(length(timeseries$site)!=length(timeseries$value))
      stop("Site indicator array 'site' must correspond with the 'value' and 'time' arrays, so the lengths must be the same!")
      
    if(sum(is.na(timeseries$site))>0)
      stop("If given, sites cannot be missing!")
      
    if(sum(sort(unique(timeseries$site))==(0:max(timeseries$site)))!=length(unique(timeseries$site)))
      stop("Sites must be numbered from 0 to #sites-1, with no sites in between missing!")
    
    timeseries$site=as.integer(timeseries$site) 
    
    if(max(timeseries$site)>0)
      has.multiple.sites=T
  }
  
  if(!is.null(timeseries$name))
  {
    if(typeof(timeseries$name)!="character")
      stop("The name must be string (called 'character' in R)!")
  }
  
  
  if(!is.null(max.layers))
  {
    if(typeof(max.layers)!="integer" & typeof(max.layers)!="double" & typeof(max.layers)!="numeric")
      stop("Number of layers, 'max.layers', must be an integer")
  }
  if(is.null(max.layers))
  {
    stop("Maximum number of layers not given!")
  }
  max.layers=as.integer(max.layers)
  if(max.layers<1)
  {
    stop("Maxmimum number of layers must be 1 or more.")
  }
  if(max.layers>=100)
  {
    stop("Cannot handle 100 layers or more (at this time).")
  }

  if(typeof(talkative)!="logical" & typeof(talkative)!="integer")
    stop("Option 'talkative' must be a logical!") 

  if(typeof(allow.one.feedback.loop)!="logical" & typeof(allow.one.feedback.loop)!="integer")
    stop("Option 'allow.one.feedback.loop' must be a logical!") 

  if(typeof(just.stationary)!="logical" & typeof(just.stationary)!="integer")
    stop("Option 'just.stationary' must be a logical!") 

  if(typeof(regional.options)!="logical" & typeof(regional.options)!="integer")
    stop("Option 'regional.options' must be a logical!") 
  
  if(typeof(time.integrals.possible)!="logical" & typeof(time.integrals.possible)!="integer")
    stop("Option 'time.integrals.possible' must be a logical!") 
  
  if(typeof(allow.deterministic.layers)!="logical" & typeof(allow.deterministic.layers)!="integer")
    stop("Option 'allow.deterministic.layers' must be a logical!") 
  
 if(typeof(use.stationary.stdev)!="logical" & typeof(use.stationary.stdev)!="integer")
   stop("Option 'use.stationary.stdev' must be a logical!") 
 
 if(typeof(do.maximum.likelihood)!="logical" & typeof(do.maximum.likelihood)!="integer")
   stop("Option 'do.maximum.likelihood' must be a logical!") 
   
 if(use.stationary.stdev & do.maximum.likelihood)
   stop("Options 'use.stationary.stdev' and 'do.maximum.likelihood' in combination is not implemented, unfortunately!")
  
 if(typeof(use.half.lives)!="logical" & typeof(use.half.lives)!="integer")
   stop("Option 'use.half.lives' must be a logical!") 

 if(typeof(mcmc)!="logical" & typeof(mcmc)!="integer")
   stop("Option 'mcmc' must be a logical!") 

 if(typeof(num.MCMC)!="integer" & typeof(num.MCMC)!="numeric" & 
   typeof(num.MCMC)!="double" )
   stop("Option 'num.MCMC' must be an integer!") 
 num.MCMC=as.integer(num.MCMC)
 
 if(num.MCMC<1)
   stop("Option 'num.MCMC' must be 1 or larger  (preferrably much larger)!")
 
 if(typeof(spacing)!="integer" & typeof(spacing)!="numeric" & 
   typeof(spacing)!="double" )
   stop("Option 'spacing' must be an integer!") 
 spacing=as.integer(spacing)
 
 if(spacing<1)
  stop("Option 'spacing' must be 1 or larger!")
 
 if(typeof(burnin)!="integer" & typeof(burnin)!="numeric" & 
   typeof(burnin)!="double" )
   stop("Option 'burnin' must be an integer!")
 burnin=as.integer(burnin) 
 
 if(burnin<0)
   stop("Option 'num.temp' must be 0 or larger!")
 
 if(typeof(num.temp)!="integer" & typeof(num.temp)!="numeric" & 
   typeof(num.temp)!="double" )
   stop("Option 'num.temp' must be an integer!") 
 num.temp=as.integer(num.temp) 

 if(num.temp<1)
   stop("Option 'num.temp' must be 1 or larger!")

 if(typeof(id.strategy)!="integer" & typeof(id.strategy)!="numeric" & 
   typeof(id.strategy)!="double" )
   stop("Option 'id.strategy' must be an integer!") 
 id.strategy=as.integer(id.strategy) 
 
 if(id.strategy<0 | id.strategy>4)
  stop("Option 'id.strategy' must be an integer between 0 and 4!")
 
 if(typeof(T.ground)!="integer" & typeof(T.ground)!="numeric" & 
    typeof(T.ground)!="double")
   stop("Option 'T.ground' must be a numeric!")
 T.ground=as.numeric(T.ground)
 
 if(T.ground<=1.0)
   stop("Option 'T.ground' must be larger than 1.0!")
 
 if(typeof(maximum.likelihood.numstart)!="integer" & 
   typeof(maximum.likelihood.numstart)!="numeric" & 
   typeof(maximum.likelihood.numstart)!="double" )
   stop("Option 'maximum.likelihood.numstart' must be an integer!") 
 maximum.likelihood.numstart=as.integer(maximum.likelihood.numstart)
 
 if(maximum.likelihood.numstart<1)
   stop("Option 'maximum.likelihood.numstart' must be 1 or larger!")
 
 if(!is.null(prior))
  {
    if(is.null(prior$mu))
      stop("If prior is given, it must contain the 95% prior credibility for expected value, given as element 'mu'!")
    if(length(prior$mu)!=2)
      stop("Prior for 'mu' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'mu'!")
    
    if(is.null(prior$dt))
      stop("If prior is given, it must contain the 95% prior credibility for the characteristic time, given as element 'dt'!")
    if(length(prior$dt)!=2)
      stop("Prior for 'dt' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'dt'!")
    
    if(is.null(prior$s))
      stop("If prior is given, it must contain the 95% prior credibility for the stochastic contributions, given as element 's'!")
    if(length(prior$s)!=2)
      stop("Prior for 's' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 's'!")
    
    if(is.null(prior$init))
      prior$init=c(prior$mu[1]-3*(prior$mu[2]-prior$mu[1]),
                   prior$mu[2]+3*(prior$mu[2]-prior$mu[1]))
    if(length(prior$init)!=2)
      stop("Prior for 'init' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'init'!")
    
    if(is.null(prior$lin))
      prior$lin=c(-1,1) # doesn't really matter since it's not going to be used
    if(length(prior$lin)!=2)
      stop("Prior for 'lin' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'lin'!")
    
    if(is.null(prior$beta))
      prior$beta=c(-1,1)
    if(length(prior$beta)!=2)
      stop("Prior for 'beta' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'beta'!")
    
    if(is.null(prior$obs) & is.null(timeseries$std.dev))
      stop("If prior is given and measurement-wise observational standard deviation is not given, the prior must contain the 95% prior credibility for the observational standard deviation, given as element 'obs'!")
    if(is.null(prior$obs))
      prior$obs=c(0.01,1) # doesn't really matter since it's not going to be used
    if(length(prior$obs)!=2)
      stop("Prior for 'obs' must contain exactly two values, namely lower and upper limit for a 95% prior credibility band for 'obs'!")
  }
 if(is.null(prior))
   prior=layer.standard.prior
 
 models=list()
 for(numlayers in 1:max.layers)
 {
   init.0=FALSE
   num.lowest=1
   lin.time=c(FALSE,FALSE,TRUE)
   no.pull=c(FALSE,TRUE,FALSE)
   
   if(!just.stationary)
   {
     init.0=TRUE
     num.lowest=3
   }
  
   num.no.sigma=1
   no.sigma=list(first=NULL)
   if(allow.deterministic.layers & numlayers>1)
   {
     #Find the ways to combine deterministic/stochastic layers:
     det.stoch=n.combinations(numlayers-1,base.number=2)
     det.stoch.list=list(first=det.stoch[1,])
     for(i in 2:(2^(numlayers-1)))
       det.stoch.list[[i]]=det.stoch[i,]
     whichone=function(x) which(x==1)
     no.sigma=lapply(det.stoch.list, whichone)
     num.no.sigma=length(no.sigma)
   }
   
   num.feedback=1
   causal=list(first=NULL)
   if(allow.one.feedback.loop & numlayers>1)
   {
     for(i in 1:(numlayers-1))
      for(j in (i+1):numlayers)
      {
        k=length(causal)
        causal[[k+1]]=c(1,i,1,j)
      }         
     num.feedback=length(causal)
   }
   
   num.time.int=1
   time.integral=list(first=NULL)
   if(time.integrals.possible & numlayers>1)
   {
     show(time.integral)
     #Find the ways to combine integral/no integral layers:
     timeint=bincombinations(numlayers-1)
     timeint.list=list(first=timeint[1,])
     for(i in 2:(2^(numlayers-1)))
       timeint.list[[i]]=timeint[i,]
     whichone=function(x) which(x==1)
     time.integral=lapply(timeint.list, whichone)
     show(time.integral)
     num.time.int=length(time.integral)
   }
   
   for(lowest in 1:num.lowest)
   if(lowest!=2 | !no.rw)
    for(sigma.strat in 1:num.no.sigma)
     for(feedback in 1:num.feedback)
      for(timeint in 1:num.time.int)
       if(length(intersect(no.sigma[[sigma.strat]],time.integral[[timeint]]))==0)
       {
         k=length(models)+1
         if(talkative)
         {
           print.srcref("")
           print.srcref(sprintf("Model %d:",k))
           print.srcref(sprintf("Number of layers: %d", numlayers))
           if(lowest==1)
             print.srcref("Lowest layer: OU")
           if(lowest==2)
             print.srcref("Lowest layer: RW")
           if(lowest==3)
             print.srcref("Lowest layer: Linear time trend OU")
           if(length(no.sigma[[sigma.strat]])>0)
             print.srcref(sprintf("No.sigma: %d", no.sigma[[sigma.strat]]))
           if(feedback!=1)
             print.srcref(sprintf("Feedback: layer %d -> layer %d", 
                          causal[[feedback]][2], causal[[feedback]][4]))
           if(timeint!=1)
             print.srcref(sprintf("Time.integral layer: %d", time.integral[[timeint]])) 
         }
         
         curr.causal=NULL
         if(feedback!=1)
	   curr.causal=as.matrix(causal[[feedback]],nrow=1)
         
	 curr.no.sigma=NULL
         if(sigma.strat>1)
           curr.no.sigma=no.sigma[[sigma.strat]]
         
	 curr.time.integral=NULL
         if(timeint>1)
           curr.time.integral=time.integral[[timeint]]
         
         struct=layer.series.structure(timeseries, numlayers=numlayers,
            lin.time=lin.time[lowest], time.integral=curr.time.integral,
            no.pull=no.pull[lowest], no.sigma=curr.no.sigma,
            init.0=init.0,prior=prior)
         
         res=layer.analyzer(struct, num.MCMC=num.MCMC, spacing=spacing, 
	     burnin=burnin, num.temp=num.temp,
             do.model.likelihood=TRUE, do.maximum.likelihood=do.maximum.likelihood, 
	     maximum.likelihood.numstart=maximum.likelihood.numstart,
             silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
             id.strategy=id.strategy,use.stationary.stdev=use.stationary.stdev,
             T.ground=T.ground, use.half.lives=use.half.lives, mcmc=mcmc,
             causal=curr.causal)
         
         models[[k]]=res
       }
 }

 return(models)
}


traverse.connections.layered=function(... ,  
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  talkative.traversal=TRUE, test.mode=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE, 
  allow.causal=TRUE, allow.correlation=TRUE, allow.direct.feedback=TRUE)
{
 data.structure=list(...)
 n=length(data.structure)
 
 if(n==0)
   stop("No time series structures given!")
 
 if(n==1)
   stop("Only one time series structure given! No connection model possible!")
 
 for(i in 1:n)
 {
  if(typeof(data.structure[[i]])!="list")
    stop("data.structure[[i]] must be a list, preferrably belonging to the 'layer.series.structure' class (as this class has the right elements).")
  
  if(is.null(data.structure[[i]]$timeseries))
    stop("data.structure[[i]] does not contain a time series (of type layer.data.series)!")
  
  if(is.null(data.structure[[i]]$timeseries$time))
    stop("No 'time' array in the incoming data list!")
    
  if(is.null(data.structure[[i]]$timeseries$value))
    stop("No 'value' array in the incoming data list!")
 }
 
  
 if(typeof(mcmc)!="logical" & typeof(mcmc)!="integer")
   stop("Option 'mcmc' must be a logical!") 

 if(typeof(silent.mode)!="logical" & typeof(silent.mode)!="integer")
   stop("Option 'silent.mode' must be a logical!") 

 if(typeof(do.maximum.likelihood)!="logical" & typeof(do.maximum.likelihood)!="integer")
   stop("Option 'do.maximum.likelihood' must be a logical!") 
   
 if(typeof(talkative.burnin)!="logical" & typeof(talkative.burnin)!="integer")
   stop("Option 'talkative.burnin' must be a logical!") 

 if(typeof(talkative.likelihood)!="logical" & typeof(talkative.likelihood)!="integer")
   stop("Option 'talkative.likelihood' must be a logical!") 
 
 if(typeof(talkative.traversal)!="logical" & typeof(talkative.traversal)!="integer")
   stop("Option 'talkative.traversal' must be a logical!") 
 
 if(typeof(test.mode)!="logical" & typeof(test.mode)!="integer")
   stop("Option 'test.mode' must be a logical!") 
 
 if(typeof(use.stationary.stdev)!="logical" & typeof(use.stationary.stdev)!="integer")
   stop("Option 'use.stationary.stdev' must be a logical!") 
 
 if(use.stationary.stdev & do.maximum.likelihood)
   stop("Options 'use.stationary.stdev' and 'do.maximum.likelihood' in combination is not implemented, unfortunately!")
  
 if(typeof(use.half.lives)!="logical" & typeof(use.half.lives)!="integer")
   stop("Option 'use.half.lives' must be a logical!") 

 if(typeof(num.MCMC)!="integer" & typeof(num.MCMC)!="numeric" & 
   typeof(num.MCMC)!="double" )
   stop("Option 'num.MCMC' must be an integer!") 
 num.MCMC=as.integer(num.MCMC)

 if(num.MCMC<1)
  stop("Option 'num.MCMC' must be 1 or larger  (preferrably much larger)!")

 if(typeof(spacing)!="integer" & typeof(spacing)!="numeric" & 
   typeof(spacing)!="double" )
   stop("Option 'spacing' must be an integer!") 
 spacing=as.integer(spacing)

 if(spacing<1)
  stop("Option 'spacing' must be 1 or larger!")

 if(typeof(burnin)!="integer" & typeof(burnin)!="numeric" & 
   typeof(burnin)!="double" )
   stop("Option 'burnin' must be an integer!")
 burnin=as.integer(burnin) 

 if(burnin<0)
  stop("Option 'num.temp' must be 0 or larger!")

 if(typeof(num.temp)!="integer" & typeof(num.temp)!="numeric" & 
   typeof(num.temp)!="double" )
   stop("Option 'num.temp' must be an integer!") 
 num.temp=as.integer(num.temp) 

 if(num.temp<1)
  stop("Option 'num.temp' must be 1 or larger!")

 if(typeof(id.strategy)!="integer" & typeof(id.strategy)!="numeric" & 
   typeof(id.strategy)!="double" )
   stop("Option 'id.strategy' must be an integer!") 
 id.strategy=as.integer(id.strategy) 
 
 if(id.strategy<0 | id.strategy>4)
  stop("Option 'id.strategy' must be an integer between 0 and 4!")
 
 if(typeof(T.ground)!="integer" & typeof(T.ground)!="numeric" & 
    typeof(T.ground)!="double")
  stop("Option 'T.ground' must be a numeric!")
 T.ground=as.numeric(T.ground)
 
 if(T.ground<=1.0)
  stop("Option 'T.ground' must be larger than 1.0!")
 
 if(typeof(maximum.likelihood.numstart)!="integer" & 
   typeof(maximum.likelihood.numstart)!="numeric" & 
   typeof(maximum.likelihood.numstart)!="double" )
   stop("Option 'maximum.likelihood.numstart' must be an integer!") 
 maximum.likelihood.numstart=as.integer(maximum.likelihood.numstart)
 
 if(maximum.likelihood.numstart<1)
  stop("Option 'maximum.likelihood.numstart' must be 1 or larger!")
 
 numlayers=rep(NA,n)
 for(i in 1:n)
   numlayers[i]=data.structure[[i]]$numlayers
 
 series1=NA
 layer1=NA
 series2=NA
 layer2=NA
 np=1 # number of pairings

 for(i in 1:(n-1))
  for(j in (i+1):n)
   {
     for(l in 1:numlayers[i])
      for(m in 1:numlayers[j])
       {   
         series1[np]=i
         layer1[np]=l
         series2[np]=j
         layer2[np]=m
         np=np+1
       }
   }
 np=length(series1)
 if(talkative.traversal)
 {
   print.srcref("Connection pairs:") 
   show(cbind(series1,layer1,series2,layer2))
   print.srcref("") 
  }

 types=c("none","F12","F21","F12,21","C12")
 if(!allow.causal)
 {
   if(!allow.correlation)
     stop("No causal nor correlative connections means no connections and thus no traversal to be done!")
   types=c("none","C12")
 }
 if(!allow.correlation)
 {
   if(allow.direct.feedback)
     types=c("none","F12","F21","F12,21")
   if(!allow.direct.feedback)
     types=c("none","F12","F21")
 }
 if(allow.causal & allow.correlation & !allow.direct.feedback)
   types=c("none","F12","F21","C12")
 numtypes=length(types) 

 combis=n.combinations(np, base.number=numtypes)
 models=list()
 
 nummodel=0
 for(i in 1:dim(combis)[1])
 {
   causal=matrix(nrow=4,ncol=0)
   corr=matrix(nrow=4,ncol=0)

   dont=F
   for(j in 1:dim(combis)[2])
   {
     #show("abc")
     #show(combis[i,j])
     #show(types[combis[i,j]+1])
     #show(layer2[j])
     #show(series2[j])
     #show(numlayers[series2[j]])
     #show(data.structure[[series2[j]]]$no.pull)
     
     if(types[combis[i,j]+1]=="F12")
     {
       if(layer2[j]==numlayers[series2[j]] &
          data.structure[[series2[j]]]$no.pull)
         dont=TRUE

       if(dim(causal)[2]>0)
         causal=cbind(causal,c(series1[j],layer1[j],series2[j],layer2[j]))
       if(dim(causal)[2]==0)
         causal=matrix(c(series1[j],layer1[j],series2[j],layer2[j]),ncol=1)
     }

     if(types[combis[i,j]+1]=="F21")
     {
       if(layer1[j]==numlayers[series1[j]] &
          data.structure[[series1[j]]]$no.pull)
         dont=TRUE

       if(dim(causal)[2]>0)
         causal=cbind(causal,c(series2[j],layer2[j],series1[j],layer1[j]))
       if(dim(causal)[2]==0)
         causal=matrix(c(series2[j],layer2[j],series1[j],layer1[j]),ncol=1)
     }

     if(types[combis[i,j]+1]=="F12,21")
     {
       if((layer1[j]==numlayers[series1[j]] &
           data.structure[[series1[j]]]$no.pull) |
          (layer2[j]==numlayers[series2[j]] &
           data.structure[[series2[j]]]$no.pull))
         dont=TRUE

       if(dim(causal)[2]>0)
         causal=cbind(causal,c(series1[j],layer1[j],series2[j],layer2[j]))
       if(dim(causal)[2]==0)
         causal=matrix(c(series1[j],layer1[j],series2[j],layer2[j]),ncol=1)
       causal=cbind(causal,c(series2[j],layer2[j],series1[j],layer1[j]))
     }
     
     if(types[combis[i,j]+1]=="C12")
     {
       if(dim(corr)[2]>0)
         corr=cbind(corr,c(series1[j],layer1[j],series2[j],layer2[j]))
       if(dim(corr)[2]==0)
         corr=matrix(c(series1[j],layer1[j],series2[j],layer2[j]),ncol=1)
     }
   }
   
   if(!dont)
   {
     nummodel=nummodel+1
     
     k=length(models)+1
     if(talkative.traversal | test.mode)
     {
       print.srcref("")
       print.srcref(sprintf("Model %d:",nummodel))
       print.srcref("Causal:")
       print.srcref(causal)
       print.srcref("Corr:")
       print.srcref(corr)
     }
     
     if(!test.mode)
     {
       res=layer.analyzer.timeseries.list(data.structure,
             num.MCMC=num.MCMC, spacing=spacing, 
	     burnin=burnin, num.temp=num.temp,
             do.model.likelihood=TRUE, do.maximum.likelihood=do.maximum.likelihood, 
	     maximum.likelihood.numstart=maximum.likelihood.numstart,
             silent.mode=silent.mode,talkative.burnin=talkative.burnin,
             talkative.likelihood=talkative.likelihood,
             id.strategy=id.strategy,use.stationary.stdev=use.stationary.stdev,
             T.ground=T.ground, use.half.lives=use.half.lives, mcmc=mcmc,
             causal=causal,corr=corr)
         
       models[[nummodel]]=res
     }
   }
 }
 
 return(models)
}
 




