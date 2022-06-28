
#####################################################
#
# Analysis methods.
#
# Requires that one or more 'layer.series.structure'
# objects have been made, representing time series
# and their internal model structure.
#
# Per default runs Bayesian MCMC inference, but can
# be convinced to run likelihood maximization
# afterwards.
#
# Trond Reitan, 17. Aug. 2018
#
#####################################################


#####################################################
# 'layer.analyzer' allows the user to input a set of
# time series+structures into the analysis.
# Options for how to run the analysis are available.
# This is however only a shell over the
# 'layer.analyzer.timeseries.list' function, which
# instead takes a list of time series+structures.
#####################################################

layer.analyzer=function(... ,
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.model.likelihood=TRUE,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE,causal=NULL,causal.symmetric=NULL,corr=NULL,
  smoothing.specs=
    list(do.smoothing=FALSE,smoothing.time.diff=0,
         smoothing.start=NULL,smoothing.end=NULL,
         num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE),
  realization.specs=
    list(do.realizations=FALSE,num.realizations=1000,strategy="N",
         realization.time.diff=0,realization.start=NULL,realization.end=NULL),
  return.residuals=FALSE)
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
      smoothing.specs=smoothing.specs, realization.specs=realization.specs,
      return.residuals=return.residuals, smooth.previous.run=FALSE,
      previous.run=NULL)
  return(ret)
}




###################################################################
# 'layer.analyzer.timeseries.list', is the real analysis function
# which takes a list of time series+structures.
###################################################################

layer.analyzer.timeseries.list=function(data.structure ,
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.model.likelihood=TRUE,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE,causal=NULL,causal.symmetric=NULL,corr=NULL,
  smoothing.specs=
    list(do.smoothing=FALSE,smoothing.time.diff=0,
         smoothing.start=NULL,smoothing.end=NULL,
         num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE),
  realization.specs=
    list(do.realizations=FALSE,num.realizations=1000,strategy="N",
         realization.time.diff=0,realization.start=NULL,realization.end=NULL),
  return.residuals=FALSE, smooth.previous.run=FALSE, previous.run=NULL
  )
{
 # Number of time series/structures:
 n=length(data.structure)

 ######################
 # User input checks:
 ######################

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
  if(data.structure[[i]]$numlayers<0)
  {
    stop("Number of layers must be 0 or more.")
  }
  if(data.structure[[i]]$numlayers>=100)
  {
    stop("Cannot handle 100 layers or more (at this time).")
  }
  
  if(!is.null(data.structure[[i]]$lin.time))
  {
    if(typeof(data.structure[[i]]$lin.time)!="logical" &
       (typeof(data.structure[[i]]$lin.time)!="integer" |
       (typeof(data.structure[[i]]$lin.time)=="integer" &
       (data.structure[[i]]$lin.time<0 | data.structure[[i]]$lin.time>1))))
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
    if(typeof(data.structure[[i]]$no.pull)!="logical" &
       (typeof(data.structure[[i]]$no.pull)!="integer" |
       (typeof(data.structure[[i]]$no.pull)=="integer" &
       (data.structure[[i]]$no.pull<0 | data.structure[[i]]$no.pull>1))))
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
    if(typeof(data.structure[[i]]$regional.mu)!="logical" &
      (typeof(data.structure[[i]]$regional.mu)!="integer" |
      (typeof(data.structure[[i]]$regional.mu)=="integer" & 
      (data.structure[[i]]$regional.mu<0 | data.structure[[i]]$regional.mu>1))))
      stop("Regional expectation indicator, 'regional.mu', must be a logical")
  }
  if(is.null(data.structure[[i]]$regional.mu))
  {
    data.structure[[i]]$regional.mu=0
  }
  data.structure[[i]]$regional.mu=as.integer(data.structure[[i]]$regional.mu)

  if(!is.null(data.structure[[i]]$regional.lin.time))
  {
    if(typeof(data.structure[[i]]$regional.lin.time)!="logical" &
       (typeof(data.structure[[i]]$regional.lin.time)!="integer" |
       (typeof(data.structure[[i]]$regional.lin.time)=="integer" &
       (data.structure[[i]]$regional.lin.time<0 |
        data.structure[[i]]$regional.lin.time>1))))
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
    if(typeof(data.structure[[i]]$differentiate.mu)!="logical" &
      (typeof(data.structure[[i]]$differentiate.mu)!="integer" |
      (typeof(data.structure[[i]]$differentiate.mu)=="integer" &
      (data.structure[[i]]$differentiate.mu<0 | 
       data.structure[[i]]$differentiate.mu>1))))
      stop("Differentiated expectancy indicator, 'differentiate.mu', must be a logical")
  }
  if(is.null(data.structure[[i]]$differentiate.mu))
  {
    data.structure[[i]]$differentiate.mu=0
  }
  data.structure[[i]]$differentiate.mu=as.integer(data.structure[[i]]$differentiate.mu)
  


  if(!is.null(data.structure[[i]]$differentiate.lin.time))
  {
    if(typeof(data.structure[[i]]$differentiate.lin.time)!="logical" &
      (typeof(data.structure[[i]]$differentiate.lin.time)!="integer" |
      (typeof(data.structure[[i]]$differentiate.lin.time)=="integer" &
      (data.structure[[i]]$differentiate.lin.time<0 |
       data.structure[[i]]$differentiate.lin.time>1))))
      stop("Differentiated linear time indicator, 'differentiate.lin.time', must be a logical")
  }
  if(is.null(data.structure[[i]]$differentiate.lin.time))
  {
    data.structure[[i]]$differentiate.lin.time=0
  }
  data.structure[[i]]$differentiate.lin.time=as.integer(data.structure[[i]]$differentiate.lin.time)
  
  if(!is.null(data.structure[[i]]$init.0))
  {
    if(typeof(data.structure[[i]]$init.0)!="logical" &
      (typeof(data.structure[[i]]$init.0)!="integer" |
      (typeof(data.structure[[i]]$init.0)=="integer" &
      (data.structure[[i]]$init.0<0 | data.structure[[i]]$init.0>1))))
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
    if(typeof(data.structure[[i]]$init.same.sites)!="logical" &
      (typeof(data.structure[[i]]$init.same.sites)!="integer" |
      (typeof(data.structure[[i]]$init.same.sites)=="integer" &
      (data.structure[[i]]$init.same.sites<0 | 
       data.structure[[i]]$init.same.sites>1))))
      stop("Initial value same for all sites indicator, 'init.same.sites', must be a logical")
  }
  if(is.null(data.structure[[i]]$init.same.sites))
  {
    data.structure[[i]]$init.same.sites=0
  }
  data.structure[[i]]$init.same.sites=as.integer(data.structure[[i]]$init.same.sites)
  
  if(!is.null(data.structure[[i]]$init.same.layers))
  {
    if(typeof(data.structure[[i]]$init.same.layers)!="logical" &
      (typeof(data.structure[[i]]$init.same.layers)!="integer" |
      (typeof(data.structure[[i]]$init.same.layers)=="integer" &
      (data.structure[[i]]$init.same.layers<0 | 
       data.structure[[i]]$init.same.layers>1))))
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
    if(typeof(data.structure[[i]]$allow.pos.pull)!="logical" &
      (typeof(data.structure[[i]]$allow.pos.pull)!="integer" |
      (typeof(data.structure[[i]]$allow.pos.pull)=="integer" &
      (data.structure[[i]]$allow.pos.pull<0 | 
       data.structure[[i]]$allow.pos.pull>1))))
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
   if(sum(class(causal)=="matrix")==0)
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
   causal=as.matrix(array(NA,c(4,0)))
 
 
 if(!is.null(causal.symmetric))
 {
   if(sum(class(causal.symmetric)=="matrix")==0)
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
   causal.symmetric=as.matrix(array(NA,c(4,0)))
 
 
 
 if(!is.null(corr))
 {
   if(sum(class(corr)=="matrix")==0)
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
   corr=as.matrix(array(NA,c(4,0)))
 
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

 if(typeof(return.residuals)!="logical" & typeof(return.residuals)!="integer")
   stop("Option 'return.residuals' must be a logical!") 
 ReturnResiduals=as.integer(return.residuals)

 if(typeof(smooth.previous.run)!="logical" & typeof(smooth.previous.run)!="integer")
   stop("Option 'smooth.previous.run' must be a logical!") 
 smooth.previous.run=as.logical(smooth.previous.run)


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

 if(smooth.previous.run==TRUE & smoothing.specs$do.smoothing==FALSE)
  {
    stop("smooth.previous.run==TRUE and smoothing.specs$do.smoothing==FALSE does not make any sense!")
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


  if(!is.null(previous.run))
  {
    if(sum(class(previous.run)=="layered")==0)
      stop("Previous run class must be a 'layered' object!")
  }

 #################################################
 # Call to C++ code for the actual analysis:
 #################################################
 
 input.mcmc=matrix(as.numeric(c(0)),nrow=1) # used for indicating 
      # that no input MCMC is used.
      # Input MCMC is only used for later smoothing.

 if(smooth.previous.run==TRUE)
 {
   if(is.null(previous.run))
     stop("previous.run not given!")
   
   if(is.null(previous.run$mcmc.origpar))
     stop("previous run has not MCMC samples (mcmc.origpar). Run again with the 'mcmc=TRUE' option!")

   if(sum(class(previous.run$mcmc.origpar)=="mcmc")==0)
   {
     stop("MCMC sample structure 'previous.run$mcmc.origpar' is not an MCMC object!")
   }

   input.mcmc=previous.run$mcmc.origpar

  #### todo: check if more needs to be done here
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
      realization.specs, ReturnResiduals,
      input.mcmc 
     )

 # Store the list of input time series/structures, in addition
 # to what the analysis returned:
 out$data.structure=data.structure

 # Store input options:
 out$input.options=
  list(num.MCMC=num.MCMC,burnin=burnin,
      spacing=spacing,num.temp=num.temp,
      do.model.likelihood=do.model.likelihood,
      do.maximum.likelihood=do.maximum.likelihood,
      maximum.likelihood.numstart=maximum.likelihood.numstart,
      silent.mode=silent.mode,
      talkative.burnin=talkative.burnin,
      talkative.likelihood=talkative.likelihood,
      id.strategy=id.strategy,
      use.stationary.stdev=use.stationary.stdev,
      T.ground=T.ground,
      use.half.lives=use.half.lives,
      mcmc=mcmc, causal=causal, 
      causal.symmetric=causal.symmetric,corr=corr, 
      smoothing.specs=smoothing.specs, 
      realization.specs=realization.specs, 
      ReturnResiduals=ReturnResiduals
      )

 # If MCMC samples were asked for, also return that:
 if(!is.null(out$mcmc))
 {
   out$mcmc=t(out$mcmc)
   colnames(out$mcmc)=out$parameter.names[1:dim(out$mcmc)[2]]
   out$mcmc=coda::mcmc(out$mcmc, start=burnin, thin=1)
 }
 if(!is.null(out$mcmc.origpar))
 {
   out$mcmc.origpar=t(out$mcmc.origpar)
   colnames(out$mcmc.origpar)=out$parameter.names[1:dim(out$mcmc.origpar)[2]]
   out$mcmc.origpar=coda::mcmc(out$mcmc.origpar, start=burnin, thin=1)
 }

 if(return.residuals & !is.null(out$standardized.residuals))
 {
   out$standardized.residuals[out$standardized.residuals==-10000000]=NA
 }
 if(return.residuals & is.null(out$standardized.residuals))
 {
   methods::show("Warning: Residuals called for but not returned!")
   methods::show(ReturnResiduals)
 }
 
 # If model log-likelihood was not (correctly) calculated, return
 # NA as the model log-likelihood.
 if(out$model.log.lik==-10000000)
  out$model.log.lik=NA


 # Also return the link specifications:
 out$causal=causal
 out$causal.symmetric=causal.symmetric
 out$corr=corr

 
 # Set class type of output:
 class(out)="layered"


 # Return output:
 return(out)
}
 





layer.predict.mcmc=function(... , analysis=NULL, 
   smoothing.time.diff=0,
   smoothing.start=NULL,smoothing.end=NULL,
   num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE)
{
  new.data.list=list(...)

  ret=layer.predict.mcmc.list(new.data.list, analysis, smoothing.time.diff,
     smoothing.start, smoothing.end, num.smooth.per.mcmc, 
     do.return.smoothing.samples)
  return(ret)
}


layer.predict.mcmc.list=function(new.data.list , analysis=NULL, 
   smoothing.time.diff=0,
   smoothing.start=NULL,smoothing.end=NULL,
   num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE,
  return.residuals=FALSE)
{
  if(is.null(analysis))
  {
    stop("Analysis must be given!")
  }

  if(sum(class(analysis)=="layered")==0)
  {
    stop("Analysis must be a result of running 'layer.analyzer'!")
  }
  
  n=length(new.data.list)
  for(i in 1:n)
  {
    if(sum(class(new.data.list[[i]])=="layer.data.series")==0)
    {
      stop("List of inputs not is not layer.data.series objects!")
    }
  }

  if(analysis$input.options$do.maximum.likelihood==TRUE & 
     is.null(analysis$mcmc.origpar))
  {
    stop("MCMC-based predictions are not available for maxmimum likelihood analyses!")
  }

  if(is.null(analysis$mcmc.origpar))
  {
    stop("MCMC-based predictions cannot be made unless the analysis returned MCMC samples! (Use the mcmc=TRUE option)")
  }
  
  if(sum(class(analysis$mcmc.origpar)=="mcmc")==0)
  {
    stop("MCMC sample structure 'mcmc.origpar' is not an MCMC object!")
  }

  ##### Todo: More checks may be needed!

  new.data.structure=analysis$data.structure
  n2=length(new.data.structure)
  if(n2!=n)
  {
    stop("Number of input series do not match those in the analysis!")
  }
 
  for(i in 1:n)
  {
     new.data.structure[[i]]$timeseries=new.data.list[[i]]
  }

  #####################################################
  #####################################################
  # Add input mcmc from previous analysis!
  #####################################################
  #####################################################

  ret=layer.analyzer.timeseries.list(new.data.structure,
      num.MCMC=analysis$input.options$num.MCMC, 
      spacing=analysis$input.options$spacing,
      burnin=analysis$input.options$burnin,
      num.temp=analysis$input.options$num.temp,
      do.model.likelihood=FALSE,
      do.maximum.likelihood=FALSE,
      maximum.likelihood.numstart=10,
      silent.mode=TRUE,talkative.burnin=FALSE,
      talkative.likelihood=FALSE,
      id.strategy=analysis$input.options$id.strategy,
      use.stationary.stdev=analysis$input.options$use.stationary.stdev,
      T.ground=analysis$input.options$T.ground, # start.parameters=0,
      use.half.lives=analysis$input.options$use.half.lives, mcmc=FALSE,
      causal=analysis$input.options$causal,
      causal.symmetric=analysis$input.options$causal.symmetric,
      corr=analysis$input.options$corr,
      smoothing.specs=list(do.smoothing=TRUE,
         smoothing.time.diff=smoothing.time.diff,
         smoothing.start=smoothing.start,smoothing.end=smoothing.end,
         num.smooth.per.mcmc=num.smooth.per.mcmc, 
         do.return.smoothing.samples=do.return.smoothing.samples),
      realization.specs=
         list(do.realizations=FALSE,num.realizations=1000,strategy="N",
         realization.time.diff=0,realization.start=NULL,realization.end=NULL),
  return.residuals=return.residuals,
  smooth.previous.run=TRUE, previous.run=analysis)
  
  return(ret)
}




layer.predict.estimate=function(... , analysis=NULL, 
   smoothing.time.diff=0,
   smoothing.start=NULL,smoothing.end=NULL,
  return.residuals=FALSE)
{
  new.data.list=list(...)
  n=length(new.data.list)
  for(i in 1:n)
  {
    if(sum(class(new.data.list[[i]])=="layer.data.series")==0)
    {
      stop("Estimate: List of inputs not is not layer.data.series objects!")
    }
  }

  analysis2=analysis
  analysis2$input.options$num.MCMC=1
  analysis2$mcmc=analysis2$mcmc.origpar=coda::mcmc(matrix(analysis$est.origpar,nrow=1))
  ret=layer.predict.mcmc.list(new.data.list, analysis2, smoothing.time.diff,
   smoothing.start,smoothing.end,num.smooth.per.mcmc=1,
   do.return.smoothing.samples=FALSE, return.residuals=return.residuals)

  return(ret)
}
