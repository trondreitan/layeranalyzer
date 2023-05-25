install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/bleeding_edge/layeranalyzer_0.1.2.tar.gz",type="source",verbose=T)

# If this went well, load the package:
library(layeranalyzer)


# Define time series (time points, values and in this case also
# sample standard deviations and sample size, used for calculating
# standard errors for each measurement). 
X=layer.data.series(time.points=malta$Time.Year,
  value.points=malta$Mean..log.body.mass.,
  std.dev=sqrt(malta$Variance.calculated.from.the.data), 
  num.meas.per.value=malta$Sample.size,name="log.body.si1ze")

p=layer.prior(mu=c(log(1),log(1000)),init=c(log(0.1),log(10000)),dt=c(0.5,20),sigma=c(0.01,2),obs=c(0.01,1),lin=c(-0.1,0.1))

X.s=layer.series.structure(X,numlayers=2,prior=p)
r=layer.analyzer(X.s,num.MCMC=100,burnin=1000,mcmc=T, talkative.burnin=T,silent.mode=F)


summary.layered(r)
#$coefficients
#                           Mean   Median Lower 95% Upper 95%
#mu_log.body.si1ze      2.439053 2.419872  2.230325  2.828314
#dt_log.body.si1ze_1    4.539107 4.047118  0.100861 11.271348
#sigma_log.body.si1ze_1 0.032268 0.026536  0.006158  0.099078
#dt_log.body.si1ze_2    7.064895 6.581348  1.729561 16.736406
#sigma_log.body.si1ze_2 0.090366 0.079410  0.022892  0.215185
#
#$model.log.likelihood
#[1] 20.22657