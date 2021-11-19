#####################################################
#
# Functions for traversing a set of models either
# for a checking multiple single time series model
# structures (number and nature of layers) or
# for checking for connections between multiple
# time series with pre-determined internal model
# structure. 
# 
# Trond Reitan, 17. Aug. 2018
#
#####################################################


# Gives a matrix of all combinations of 0's and 1'
# (or another base number), counting up from
# 0 to base number^length. Used by "traverse.standalone.layered"
# and "traverse.connections.layered".

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


# Traverse the set of model structures for a given time series,
# up to a given complexity (as set by the maximum number of layers).
# Traversal options can be given, that models with feedback-loops
# can be included or excluded, non-stationarity or deterministic layers
# can be allowed or disallowed etc. MCMC options can also be given.
# Returns the set of analyses examined, which can then directly be
# examined with 'compare.layered'.

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
 
 num.lowest=1
 init.0=FALSE
 if(!just.stationary)
   {
     init.0=TRUE
     num.lowest=3
   }

 models=list()
 for(numlayers in 1:max.layers)
 {
   lin.time=c(FALSE,FALSE,TRUE)
   no.pull=c(FALSE,TRUE,FALSE)
   
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
     methods::show(time.integral)
     #Find the ways to combine integral/no integral layers:
     timeint=n.combinations(numlayers-1,base.number=2)
     timeint.list=list(first=timeint[1,])
     for(i in 2:(2^(numlayers-1)))
       timeint.list[[i]]=timeint[i,]
     whichone=function(x) which(x==1)
     time.integral=lapply(timeint.list, whichone)
     methods::show(time.integral)
     num.time.int=length(time.integral)
   }
   
   for(lowest in 1:num.lowest)
   if(lowest!=2 | !no.rw)
    for(sigma.strat in 1:num.no.sigma)
     for(feedback in 1:num.feedback)
     {
      # Check if feedback is meaningful before proceeeding:
      dofeedback=FALSE
      if(feedback==1) # no feedback loops?
        dofeedback=TRUE # proceed
      if(feedback>1)  # feedback not to lowest layer or lowest layer isn't RW?
       if(causal[[feedback]][4]!=numlayers | lowest!=2)
        dofeedback=TRUE # proceed
 
      if(dofeedback) # proceed with feedback?
      {
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

layer.analyzer.timeseries.list(list(struct),num.MCMC=num.MCMC, spacing=spacing, 
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
     }
 }
 

 return(models)
}



# Traverse the set of connections between time series+model structure
# processes. The number of processes should be low, either just 3
# for single layered processes or 2 for multi-layered processes.
# Causal and correlative links and causal feedback loops
# can be allowed or disallowed separately. MCMC options can also be given.
# Returns the set of analyses examined, which can then directly be
# examined with 'compare.layered'.

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
   methods::show(cbind(series1,layer1,series2,layer2))
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
         
       ll=res$model.log.lik
       if(do.maximum.likelihood)
	 ll=res$input[[i]]$ML.loglik
       if(length(ll)==0)
         ll=NA
       iter=0
       while(iter<5 & (is.na(ll) | !(ll > -1e+199 & ll < 1e+199)))
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
         
            ll=res$model.log.lik
	    if(do.maximum.likelihood)
	      ll=res$input[[i]]$ML.loglik
            if(length(ll)==0)
              ll=NA
	    iter=iter+1
	 }

       models[[nummodel]]=res
     }
   }
 }
 
 return(models)
}
 





# Stepwise search through the set of connections between
# time series+model structure processes, when the number of
# such processes can be moderately large. Not well tested.
# use of parallel processing using the standalone program
# or multiple instances of R is any way recommended for speed.
# Causal and correlative links and causal feedback loops
# can be allowed or disallowed separately. MCMC options can also be given.
# Returns the set of analyses examined, which can then directly be
# examined with 'compare.layered', but note that not all connection
# models are examined in a step-wise search.

stepwise.connections.layered=function(... ,  
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  talkative.traversal=TRUE, test.mode=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE, 
  allow.causal=TRUE, allow.correlation=TRUE, allow.direct.feedback=TRUE,
  first.is.nullhypothesis=FALSE,ML.IC="AIC")
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
   methods::show(cbind(series1,layer1,series2,layer2))
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

 nullmodel.conn=rep("none",np) 
 nullmodel=layer.analyzer.timeseries.list(data.structure,
	        num.MCMC=num.MCMC, spacing=spacing, 
		burnin=burnin, num.temp=num.temp,
                do.model.likelihood=TRUE, 
                do.maximum.likelihood=do.maximum.likelihood,		
	        maximum.likelihood.numstart=maximum.likelihood.numstart,
                silent.mode=silent.mode,talkative.burnin=talkative.burnin,
                talkative.likelihood=talkative.likelihood,
                id.strategy=id.strategy,use.stationary.stdev=use.stationary.stdev,
                T.ground=T.ground, use.half.lives=use.half.lives, mcmc=mcmc)
 null.prob=0  
 best.prob=1 

 models=list()
 models[[1]]=nullmodel
 nummodel=1
 model.conn=list(nullmodelconn=nullmodel.conn)

 while(null.prob<best.prob)
 {
   # Try changing one connection at a time:
   for(i in 1:np)
   {
     for(j in 1:numtypes) # Traverse the possible connection types:
      if(types[j]!=nullmodel.conn[i])
       {
         conn=nullmodel.conn
         conn[i]=types[j]         

         if(talkative.traversal)
           print.srcref(sprintf("Try change connection %d.%d-%d.%d from %s to %s",
 	   	        series1[i],layer1[i],series2[i],layer2[i],
                        nullmodel.conn[i],types[j]))      
                    
         # Build connection descriptions:
         causal=matrix(nrow=4,ncol=0)
         corr=matrix(nrow=4,ncol=0)
         for(i2 in 1:np)
         {
            if(conn[i2]=="C12")
            {
              if(dim(corr)[2]>0)
                corr=cbind(corr,c(series1[i2],layer1[i2],series2[i2],layer2[i2]))
              if(dim(corr)[2]==0)
                corr=matrix(c(series1[i2],layer1[i2],series2[i2],layer2[i2]),ncol=1)
            }
            if(conn[i2]=="F12")
            {
              if(dim(causal)[2]>0)
                causal=cbind(causal,c(series1[i2],layer1[i2],series2[i2],layer2[i2]))
              if(dim(causal)[2]==0)
                causal=matrix(c(series1[i2],layer1[i2],series2[i2],layer2[i2]),ncol=1)
            }
            if(conn[i2]=="F21")
            {
              if(dim(causal)[2]>0)
                causal=cbind(causal,c(series2[i2],layer2[i2],series1[i2],layer1[i2]))
              if(dim(causal)[2]==0)
                causal=matrix(c(series2[i2],layer2[i2],series1[i2],layer1[i2]),ncol=1)
            }
            if(conn[i2]=="F12,F21")
            {
              if(dim(causal)[2]>0)
                causal=cbind(causal,c(series2[i2],layer2[i2],series1[i2],layer1[i2]))
              if(dim(causal)[2]==0)
                causal=matrix(c(series2[i2],layer2[i2],series1[i2],layer1[i2]),ncol=1)
              causal=matrix(c(series1[i2],layer1[i2],series2[i2],layer2[i2]),ncol=1)
            }
         }
         
         if(talkative.traversal | test.mode)     
         {
           print.srcref("")
           print.srcref(sprintf("Model %d:",nummodel))
           print.srcref("Causal:")
           print.srcref(causal)
           print.srcref("Corr:")
           print.srcref(corr)
         }

         nummodel=nummodel+1
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
	              
         ll=res$model.log.lik
	 if(do.maximum.likelihood)
	   ll=res$input[[i]]$ML.loglik
         if(length(ll)==0)
           ll=NA
	 iter=0
	 while(iter<5 & (is.na(ll) | !(ll > -1e+199 & ll < 1e+199)))
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
         
           ll=res$model.log.lik
	   if(do.maximum.likelihood)
	     ll=res$input[[i]]$ML.loglik
           if(length(ll)==0)
             ll=NA
	   iter=iter+1
	 }

         models[[nummodel]]=res
         model.conn[[nummodel]]=conn
       }        
   }   

   comp=compare.layered(models, ML.IC=ML.IC,
         first.is.nullhypothesis=first.is.nullhypothesis)
   if(talkative.traversal)
     methods::show(comp)
     
   best.prob=max(comp[,2])
   best.model.index=which(comp[,2]==max(comp[,2]))
   best.model=models[[best.model.index]]
   best.conn=model.conn[[best.model.index]]

   if(best.prob>null.prob)
   {
     nullmodel=best.model
     nullmodel.conn=best.conn
     models=list()
     models[[1]]=nullmodel
     nummodel=1
     model.conn=list(nullmodelconn=nullmodel.conn)
   }
 }

 return(nullmodel)
}
 




