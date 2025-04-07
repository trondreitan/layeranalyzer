###############################################################
#
# Time series and time series+model structure representation
#
# 'layer.data.series' represents a measured time series.
# Minimal specification consists of time points and
# corresponding measurements. Sites and standard deviations
# are optional (standard deviation will then be assumed
# constant and be a parameter to be estimated). Number
# of measurements in each time point can also be given, in
# which case the values will be interpreted as means.
# the mean standard deviation is then the individual
# measurement standard deviation divided by the square
# root of the number of measurements for each time point.
#
# 'layer.series.structure' represents a time series plus
# the internal model structure, i.e. number of layers (hidden
# layers plus measured layer) and the nature of each layer.
# The prior specification for the process can also be given
# and should be given if the user has better information
# concerning the process than that coded in the standard prior.
#
# Trond Reitan, 17. Aug. 2018
#
###############################################################


# 'layer.data.series' creates an object that contains
# time points and mean measurement values, as well as the possibility
# for sites, standard deviations and number of measurements per
# time point. A name for the time series must be given, in order
# to make sensible parameter names later.

layer.data.series=function(time.points, value.points, name, std.dev=NULL, 
                           num.meas.per.value=NULL, site=NULL)
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

  o=order(time)
  reorder=FALSE
  if(sum((1:length(time))!=o)>0)
  {
    reorder=TRUE
    print.srcref("Warning: time points are not ordered chronologiclaly (from smallest to highest value). Will reorder automatically.")
    print.srcref("(PS: This could be because you have fed the timeseries structure *age* instead of time. Set time=-age, if so, to avoid this warning.") 
  }	

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
      
    if(min(std.dev[!is.na(std.dev) & std.dev!=-10000000])<0)
      stop("Negative standard deviations found!")
      
    if(min(std.dev[!is.na(std.dev) & std.dev!=-10000000])==0)
      stop("Standard deviations=0 found!") 

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
  
  if(is.null(name))
    stop("Name must be given!")

  if(name=="")
    stop("Name must be given!")

  if(!reorder)
  {
    ret=list(time=time.points,value=value.points,name=name)
    if(!is.null(std.dev))
      ret$std.dev=std.dev
    if(!is.null(num.meas.per.value))
      ret$num.meas.per.value=num.meas.per.value
    if(!is.null(site))
      ret$site=site
    if(!is.null(is.datetime))
      ret$is.datetime=is.datetime
  }
  if(reorder)
  {
    ret=list(time=time.points[o],value=value.points[o],name=name)
    if(!is.null(std.dev))
      ret$std.dev=std.dev[o]
    if(!is.null(num.meas.per.value))
      ret$num.meas.per.value=num.meas.per.value[o]
    if(!is.null(site))
      ret$site=site[o]
    if(!is.null(is.datetime))
      ret$is.datetime=is.datetime[o]
  }
      
  class(ret)="layer.data.series"
  
  return(ret)
}  


# 'read.layer.data' reads a file having at least the columns 'time'
# and 'value' and tries to make a 'layer.data' object out of it.
# Input otherwise is the same as 'read.table', so that different
# file formats can be read.

read.layer.data.series=function(filename,name,column.type=c("time","value"), 
                header=FALSE,dec=".",sep="", 
                quote="\"'", numerals = c("allow.loss", "warn.loss", "no.loss"),
                row.names=NULL, col.names=column.type,
                na.strings = "NA", colClasses = NA, nrows = -1,
                skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                strip.white = FALSE, blank.lines.skip = TRUE,
                comment.char = "#",
                allowEscapes = FALSE, flush = FALSE,
		stringsAsFactors = FALSE, as.is = !stringsAsFactors,
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

  if(is.null(name))
    stop("Name must be given!")

  if(name=="")
    stop("Name must be given!")

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



# 'layer.series.structure' represents a time series plus internal
# model structure. Takes a 'layer.series' object as input, plus
# options for the model structure. Default is a one-layered
# Ornstein-Uhlenbeck process with standard prior specification.

layer.struct.name=function(struct)
{
 if(class(struct)!="layer.series.structure")
   stop("Input must be a 'layer.series.structure' object!")

 out=""
 if(struct$numlayers>1)
   out=sprintf("%d-layered:",struct$numlayers)

 if(struct$numlayers==0)
 {
   out=sprintf("%s%s",out,"No stochasticity")
   if(struct$lin.time)
     {
       out=paste(out," - lin.trend")
     }
   if(struct$lin.time)
     {
       out=paste(out," - constant")
     }
 }
 
 if(struct$numlayers>=1)
  for(layer in struct$numlayers:1)
 {
   if(struct$numlayers>1)
   {
     if(layer<struct$numlayers &
        (struct$numlayers>2 |
         (sum(names(struct$timeseries)=="site")>0 &
	  length(unique(struct$timeseries$site))>1)))
     {
       out=paste(out,sprintf("\n           Layer %d: ",layer))
     }
     if(layer<struct$numlayers &
        (struct$numlayers==2 &
         !(sum(names(struct$timeseries)=="site")>0 &
	   length(unique(struct$timeseries$site))>1)))
     {
       out=paste(out, sprintf(", Layer %d: ", layer))
     }
     if(layer==struct$numlayers)
     {
       out=paste(out, sprintf("Layer %d: ", layer))
     }
   }
   
   if(layer==struct$numlayers)
   {
     if(struct$no.pull)
     {
       out=sprintf("%s%s",out,"RW")
     }
     else
     {
       out=sprintf("%s%s",out,"OU")
     }
     if(struct$lin.time)
     {
       out=paste(out,", lin.trend")
     }
     if(struct$allow.pos.pull & !struct$no.pull)
     {
       out=paste(out,", positive pull allowed")
     }
   }
   
   if(layer<struct$numlayers)
   {
     ordinary=TRUE
     if(sum(struct$no.sigma==layer)>0)
     {
       out=sprintf("%s%s",out,"deterministic tracking")
       ordinary=FALSE
     }
     if(sum(struct$time.integral==layer)>0)
     {
       out=sprintf("%s%s",out,"integration layer")
       ordinary=FALSE
     }
     if(ordinary)
     {
       out=sprintf("%s%s",out,"OU-like tracking")
     }
   }

   # Site structure?
   if(sum(names(struct$timeseries)=="site")>0 &
      length(unique(struct$timeseries$site))>1)
   {
     if(layer==struct$numlayers)
     {
       if(struct$regional.mu)
       {
	 out=paste(out,", regional mu")
       }
       else
       {
         out=paste(out,", global mu")
       } 
    
       if(struct$lin.time)
       {
         if(struct$regional.lin.time)
	 {
             out=paste(out,", regional lin.time")
	   }
	 else
	 {
           out=paste(out,", global lin.time")
         }
       }
     }
     
     if(sum(struct$regional.pull==layer)==1)
     {
       out=paste(out,", regional pull")
     }
     else
     {
       out=paste(out,", global pull")
     }

     if(sum(struct$regional.sigma==layer)==1)
     {
       out=paste(out,", regional sigma")
     }
     else
     {
       out=paste(out,", global sigma")
     }

     has.corr=FALSE
     if(sum(struct$correlated.sigma==layer)==1)
     {
       out=paste(out,", correlation between sites")
       has.corr=TRUE
     }
     if(sum(struct$pairwise.correlated.sigma==layer)==1)
     {
       out=paste(out,", pairwise correlated")
       has.corr=TRUE
     }
     if(sum(struct$one.dim.sigma==layer)==1)
     {
       out=paste(out,", perfect correlation")
       has.corr=TRUE
     }
     if(sum(struct$grouping.sigma==layer)==1)
     {
       out=paste(out,", grouped sigma")
       has.corr=TRUE
     }
     if(sum(struct$remove.sigma)==1)
     {
       out=paste(out,", remove some corr.")
       has.corr=TRUE
     }
     if(sum(struct$differentiate.sigma)==1)
     {
       out=paste(out,", grouped corr. strength")
       has.corr=TRUE
     }
     if(!has.corr)
     {
       out=paste(out,", no site correlation")
     }

     if(sum(struct$differentiate.pull)==1)
     {
       out=paste(out,", grouped char.time")
     }
     if(sum(struct$differentiate.mu)==1)
     {
         out=paste(out,", grouped mu")
     }
     if(sum(struct$differentiate.lin.time)==1 &
            struct$lin.time)
     {
      out=paste(out,", grouped lin. time")
     }
   }
 }

 return(out)
}
   

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
########################
# Check user input:
########################

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
  if(numlayers<0)
  {
    stop("Number of layers must be 0 or more.")
  }
  if(numlayers>=100)
  {
    stop("Cannot handle 100 layers or more (at this time).")
  }

  ###################################
  # Make return object, containing
  # time series and number of layers:
  ###################################
  ret=list(timeseries=timeseries, numlayers=numlayers)


  # More user input checks:

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

  init.treatment=FALSE

  if(!is.null(init.0))
  {
    if(typeof(init.0)!="logical")
      stop("Initial value treatement indicator, 'init.0', must be a logical")
    ret$init.0=init.0
  init.treatment=TRUE
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
    init.treatment=TRUE
  }
  
  if(!is.null(init.same.sites))
  {
    if(typeof(init.same.sites)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.sites', must be a logical")
    ret$init.same.sites=init.same.sites

    if(!init.treatment && init.same.sites)
      print.srcref("Warning: Specifying init.same.sites=TRUE without any initial value treatment makes no sense!")
  }
  
  if(!is.null(init.same.layers))
  {
    if(typeof(init.same.layers)!="logical")
      stop("Initial value same for all sites indicator, 'init.same.layers', must be a logical")
    ret$init.same.layers=init.same.layers
    
    if(!init.treatment && init.same.layers)
      print.srcref("Warning: Specifying init.same.layers=TRUE without any initial value treatment makes no sense!")
  }
  
  if(!is.null(init.specified))
  {
    if(typeof(init.specified)!="numeric" & typeof(init.specified)!="double")
      stop("Specified initial value, 'init.specified', must be a numeric. If time has been specified with POSIXct, cast this as numeric here.")
    if(length(init.specified)!=2)
      stop("Specified initial value, 'init.specified', must be a vector with 2 values, 'time' and 'value'. If time has been specified with POSIXct, cast this as numeric")
    init.specified=as.numeric(init.specified)
    ret$init.specified=init.specified
    init.treatment=TRUE
  }

  if(no.pull && !init.treatment)
  {
    print.srcref("Warning: no.pull=TRUE (i.e. Brownian motion on the bottom layer) *should* have initial state treatment. If not, a very vague prior distribution is given to the initial state.")
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

  # Set output type and return:
  class(ret)="layer.series.structure"
  
  ret$description=layer.struct.name(ret)

  return(ret)
}

print.layer.series.structure=function(x)
  cat(sprintf("Timeseries: %s, structure:\n%s\n",x$timeseries$name, x$description))
  


