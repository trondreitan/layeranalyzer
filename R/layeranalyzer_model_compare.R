#####################################################
#
# Functions for comparing the output from multiple
# analysis, using Bayesian model testing or
# max-likelihood-based information criteria.
# 
# Trond Reitan, 17. Aug. 2018
#
#####################################################


# 'compare.layered' checks a set of specified analysis outputs
# for different models. Option 'p0' lets the user set all
# prior model proabilities. Without it, all models are assumed
# á prior equally probable. Option 'first.is.nullhypothesis'
# instead treats the first model as a null model and gives that
# á prior 50% model probability, while the rest of the models
# share equally the remaining 50% probability.
# Option 'ML.IC' allows the user to specify which information
# criterion is to be used, in case of classic analysis.
# Possibilities are "AIC", "BIC" or "AICc".
         
compare.layered=function(...,p0=NULL,first.is.nullhypothesis=FALSE,
   ML.IC="AIC")
{
  input=list(...)

  if(sum(class(input[[1]])=="layered")==0)
    if(sum(class(input[[1]][[1]])=="layered")>0)
      input=input[[1]]
  
  n=length(input)

  all.ml=FALSE
  for(i in 1:n)
  {
   if(sum(class(input[[i]])=="layered")==0)
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
  
  no.lineshift=TRUE
  for(i in 1:n)
    if(length(grep("\n",input[[i]]$description))>0)
      no.lineshift=FALSE
  maxdesc=0
  for(i in 1:n)
    maxdesc=max(maxdesc,nchar(input[[i]]$description))
  use.desc=no.lineshift & maxdesc<70
  
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

    if(is.na(l[i]) | !(l[i] > -1e+199 & l[i]< +1e+199))
      l[i]=-1e+200
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
  {
    names(input)=sprintf("Model %3d",1:n)
    if(use.desc)
    {
      for(i in 1:n)
        names(input)[i]=sprintf(sprintf("Model %%3d: %%%ds",maxdesc),
	                        i,input[[i]]$description)
    }
  }
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
  
  #  show("compare called")

  dimnames(outmatrix)=list(names(input),c(critname,"Post. Prob.(%)"))
  
  round(outmatrix,5)
}



# Returns number of observations for an output object from the analysis:

nobs.layered=function(object,...)
{
  if(sum(class(object)=="layered")==0)
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  
 n=length(object$data.structure)
 len=0

 for(i in 1:n)
 {
  if(typeof(object$data.structure[[i]])!="list")
    stop("data.structure[[i]] must be a list, preferrably belonging to the 'layer.series.structure' class (as this class has the right elements).")
  
  if(is.null(object$data.structure[[i]]$timeseries))
    stop("data.structure[[i]] does not contain a time series (of type layer.data.series)!")
  
  if(is.null(object$data.structure[[i]]$timeseries$time))
    stop("No time points found in data material!")

  len=len+length(object$data.structure[[i]]$timeseries$time)
 }

 return(len)
}

# Returns the ML log-likelihood for an object returned from an analysis:

logLik.layered=function(object,...)
{
  if(sum(class(object)=="layered")==0)
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  if(is.null(object$ML.loglik))
    stop("The 'layered' object need to be classically (ML) estimated!") 

  ll=object$ML.loglik

  df=length(object$parameter.names[substr(object$parameter.name,1,7)!="complex"])
  attr(ll,"class")="logLik"
  attr(ll,"df")=df
  attr(ll,"nobs")=nobs.layered(object)

  return(ll)
}


# Performs likelihood-ratio test (chi-squared) for a set of
# objects from analyses, representing different models.
# Returns p-value with more, as other 'anova' methods in R do.

anova.layered=function(object,...)
{
  if(sum(class(object)=="layered")==0)
    stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
  
  input=list(...)
  n=length(input)
  
  if(is.null(object$ML.loglik))
    stop("All objects need to be classically (ML) estimated!") 
  for(i in 1:n)
  {
    if(sum(class(input[[i]])=="layered")==0)
     stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.")
    if(is.null(input[[i]]$ML.loglik))
      stop("All objects need to be classically (ML) estimated!") 
  }
  
  prev.df=length(object$parameter.names)
  prev.deviance=-2*object$ML.loglik
  
  res.df=rep(NA,n+1)
  res.dev=rep(NA,n+1)
  df=rep(NA,n+1)
  dev=rep(NA,n+1)
  pr=rep(NA,n+1)
  
  res.df[1]=nobs.layered(object)-attr(logLik.layered(object),"df")
  res.dev[1]=-2*as.numeric(logLik.layered(object))
  
  if(n>0)
   for(i in 1:n)
   {
    if(sum(class(input[[i]])=="layered")==0)
     stop("Function can only be used on objects of type \"layered\",\n  i.e. objects returned from the 'layer.analysis' method.") 
    
    res.df[1+i]=nobs.layered(input[[i]])-attr(logLik.layered(input[[i]]),"df")
    res.dev[1+i]=-2*as.numeric(logLik.layered(input[[i]]))
    
    df[i+1]=res.df[i]-res.df[i+1]
    dev[i+1]=res.dev[i]-res.dev[i+1]
    if(df[i+1]>0)
      pr[i+1]=stats::pchisq(dev[i+1],df[i+1],lower.tail=FALSE)
   }
  
  ret=list(res.df=res.df, res.dev=res.dev, df=df, dev=dev, pr=pr)
  names(ret)=c("Resid. Df","Resid. Dev", "Df", "Deviance", "Pr(>Chi)")
  class(ret)=c("anova","data.frame")
  attr(ret,"row.names")=1:(n+1)
  attr(ret,"heading")=c("Analysis of Deviance Table\n")
  
  return(ret)
}



 




