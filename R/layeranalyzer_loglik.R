library(coda)


layer.param.loglik=function(analysis, new.param.values=NULL, silent.mode=TRUE,
		            num.optim=100, do.preanalysis.mcmc=FALSE,
			    num.MCMC=1000,spacing=10,
			    burnin=10000,num.temp=1)
{
 if(is.null(analysis))
 {
   stop("The 'analysis' object must exists (and be a 'layered' object)!")
 }
 if(sum(class(analysis)=="layered")==0)
 {
   stop("The 'analysis' object must be a 'layered' object!")
 }

 analysis2=analysis
 analysis2$input.options$num.MCMC=1
 if(is.null(new.param.values))
 {
  analysis2$mcmc=analysis2$mcmc.origpar=coda::mcmc(matrix(analysis2$est.par,nrow=1))
 }
 
 if(!is.null(new.param.values))
 {
  if(length(analysis2$est.origpar)!=length(new.param.values))
    stop("Number of new input parameter values has a length that does not match the models number of parameters (as specified in 'est.origpar')!")
  analysis2$mcmc=analysis2$mcmc.origpar=new.param.values
 }
 
 # need to transform back to characteristic times if half life
 numpar=length(analysis2$est.origpar)
 if(analysis2$input.options$use.half.lives) 
  for(j in 1:numpar)
  {
   if(substr(analysis2$parameter.names[j],1,3)=="dt_")
     analysis2$mcmc[j,]=analysis2$mcmc[j,]/log(2)
  }
 analysis2$mcmc[which(is.na(analysis2$mcmc), arr.ind=TRUE)] = -1e+7
 analysis2$mcmc.origpar[which(is.na(analysis2$mcmc.origpar), arr.ind=TRUE)] = -1e+7
 analysis2$mcmc=coda::as.mcmc(rbind(analysis2$mcmc))
 analysis2$mcmc.origpar=coda::as.mcmc(rbind(analysis2$mcmc.origpar))

 if(!do.preanalysis.mcmc)
   res=layer.analyzer.timeseries.list(analysis$data.structure,
     previous.run=analysis2,silent.mode=silent.mode,
      causal=analysis$causal,corr=analysis$corr,
      causal.symmetric=analysis$causal.symmetric,
     maximum.likelihood.numstart=num.optim,
     num.MCMC=0,spacing=0,burnin=0,num.temp=0,
     layer.analyzer.mode="Loglik-from-input")
 if(do.preanalysis.mcmc)
   res=layer.analyzer.timeseries.list(analysis$data.structure,
     previous.run=analysis2,silent.mode=silent.mode,
      causal=analysis$causal,corr=analysis$corr,
      causal.symmetric=analysis$causal.symmetric,
     maximum.likelihood.numstart=num.optim,
     num.MCMC=num.MCMC,spacing=spacing,burnin=burnin,num.temp=num.temp,
     layer.analyzer.mode="Loglik-from-input")
     
 return(res$loglik)
}


layer.param.logliks=function(analysis, new.param.value.sets, silent.mode=TRUE,
         num.optim=100, do.preanalysis.mcmc=FALSE,
	 num.MCMC=1000,spacing=10,burnin=10000,num.temp=1)
{
 if(is.null(analysis))
 {
   stop("The 'analysis' object must exists (and be a 'layered' object)!")
 }
 if(sum(class(analysis)=="layered")==0)
 {
   stop("The 'analysis' object must be a 'layered' object!")
 }
 
 analysis2=analysis
 analysis2$input.options$num.MCMC=1
 if(is.null(new.param.value.sets))
 {
   stop("Parameter sets must be given!")
 }
 
 if(length(analysis2$est.origpar)!=dim(new.param.value.sets)[2])
    stop("Number of new input parameter values has a length that does not match the models number of parameters (as specified in 'est.origpar')!")
  analysis2$mcmc=analysis2$mcmc.origpar=new.param.value.sets
 
 # need to transform back to characteristic times if half life
 numpar=length(analysis2$est.origpar)
 if(analysis2$input.options$use.half.lives) 
  for(j in 1:numpar)
  {
   if(substr(analysis2$parameter.names[j],1,3)=="dt_")
     analysis2$mcmc[j,]=analysis2$mcmc[j,]/log(2)
  }
 analysis2$mcmc[which(is.na(analysis2$mcmc), arr.ind=TRUE)] = -1e+7
 analysis2$mcmc.origpar[which(is.na(analysis2$mcmc.origpar), arr.ind=TRUE)] = -1e+7
 analysis2$mcmc=coda::as.mcmc(analysis2$mcmc)
 analysis2$mcmc.origpar=coda::as.mcmc(analysis2$mcmc.origpar)
     
 if(!do.preanalysis.mcmc)
   res=layer.analyzer.timeseries.list(analysis$data.structure,
      previous.run=analysis2,silent.mode=silent.mode,
      causal=analysis$causal,corr=analysis$corr,
      causal.symmetric=analysis$causal.symmetric,
      maximum.likelihood.numstart=num.optim,
      num.MCMC=0,spacing=0,burnin=0,num.temp=0,
      layer.analyzer.mode="Loglik-from-input")
 if(do.preanalysis.mcmc)
   res=layer.analyzer.timeseries.list(analysis$data.structure,
      previous.run=analysis2,silent.mode=silent.mode,
      causal=analysis$causal,corr=analysis$corr,
      causal.symmetric=analysis$causal.symmetric,
      maximum.likelihood.numstart=num.optim,
      num.MCMC=num.MCMC,spacing=spacing,burnin=burnin,num.temp=num.temp,
      layer.analyzer.mode="Loglik-from-input")
     
 return(res$loglik)
}


