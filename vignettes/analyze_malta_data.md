---
title: "analyze_malta_data"
output: rmarkdown::html_vignette
author: Trond Reitan
vignette: >
  %\VignetteIndexEntry{analyze_malta_data}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---
# Stuctural analyzis of Malta reed warbler data

First install the *layeranalyzer* package. If you have *devtools* installed,
this can be achived by the command:
> install_github(repo="trondreitan/layeranalyzer", build_vignettes=TRUE)

If devtools is not available, you need to reference the correct version:
> install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T, build_vignettes=TRUE)

The personal web page has newer "bleeding edge" versions of the package:
> install.packages("https://folk.uio.no/trondr/R/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T, build_vignettes=TRUE)

# If the package is installed, load it:
```r
library(layeranalyzer)
```

Usually, the next thing to happen would be for the data to be read:
> malta=read.table("malta.csv",sep=";",header=T)
However, this particular dataset is already in the package.

We will now defined the Bayesian prior for teh parameters. This
may seem a bit unnecessary, since the analysis is classic (ML-based),
but MCMC sampling is needed to provide the maximum likelihood
optimization with good starting points.

##Priors:
Initial values have a larger prior interval than the mean
in order to represent that the process might not start in equilibrium.
Wide priors for stochastic term and noise standard deviation.
characteristic time (dt) set to having a lower limit right below the
time resolution and upper limit comparable to the data length.
Priors are not so important when the analysis is classic, but
since the ML optimization starts from MCMC samples, a good prior
can give the optimization a good start.

```r
p=layer.prior(mu=c(log(1),log(1000)),init=c(log(0.1),log(10000)),dt=c(0.5,20),sigma=c(0.01,2),obs=c(0.01,1),lin=c(-0.1,0.1))
```

## Define the time series
Time points, values and in this case also sample standard deviations
and sample size, used for calculating standard errors for each measurement).

```r
X=layer.data.series(time.points=malta$Time.Year,
  value.points=malta$Mean..log.body.mass.,
  std.dev=sqrt(malta$Variance.calculated.from.the.data), 
  num.meas.per.value=malta$Sample.size,name="log.body.size")
```

We will now plot data, to see if it makes sense:

```r
plot(X$time, X$value,type="b",ylim=c(2.2,2.7))
for(i in 1:length(X$time))
  lines(c(X$time[i],X$time[i]),
        c(X$value[i]-1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i]),
          X$value[i]+1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i])))
```



## Call 'traverse.standalone.layered'
in order to automatically go through all process structures up to a
given complexity (max layers=2). do.maximum.likelihood = TRUE sets
ML-based estimation in stead of Bayesian inference, while
maximum.likelihood.numstart = 1000 sets the number of ML
optimizations performed (starting at different MCMC samples).

Note that initial value treatment is automatically used when
non-stationary models are among the possible candidates.
This is because the likelihood of the first measurement would
otherwise be determined by the stationary distribution of
the process, which is not available for non-stationary models.


PS: This will take a while:
```r
models.ml=traverse.standalone.layered(X, max.layers=2, 
  talkative=TRUE, allow.one.feedback.loop=TRUE, 
  just.stationary=FALSE, no.rw=FALSE,    
  time.integrals.possible=FALSE, 
  allow.deterministic.layers=TRUE,
  do.maximum.likelihood = TRUE, maximum.likelihood.numstart = 1000, 
  num.MCMC=1000,spacing=10,burnin=2000, num.temp = 4, prior=p)
```


## Compare using AICc:

```r
compare.layered(models.ml, ML.IC = "AICc")
```
This will output something like the following:  
```r
#weight=-0.5*AICc Post. Prob.(%)  
#Model   1         29.35569       16.57780  
#Model   2         23.78279        0.06299  
#Model   3         30.48248       51.15449  
#Model   4         26.14139        0.66615  
#Model   5         23.64615        0.05494  
#Model   6         29.46821       18.55215  
#Model   7         26.88072        1.39528  
#Model   8         25.92793        0.53811  
#Model   9         27.67578        3.08994  
#Model  10         25.07575        0.22949  
#Model  11         22.65920        0.02048  
#Model  12         28.51285        7.13650  
#Model  13         25.89695        0.52169  
```
Use "summary(models.ml[[1]])" (for instance) to look at the
parameter names and thus see what kind of model it is.

PS: There may be some variation in the results as the estimation
has some level of stochasticity.



## It turns out that the best model is linear trend OU.
Make a summary:

```r
summary(models.ml[[3]])
```
Output:
```r
#Coefficients:  
#                         ML estimate Bayesian Lower 95% Bayesian Upper 95%  
#mu_log.body.size            2.335158          -2.309826           2.483512  
#lin_t_log.body.size         0.007669          -0.044273           0.904957  
#dt_log.body.size_1          2.518253           2.211906         416.657806  
#sigma_log.body.size_1       0.000324           0.002351           0.052827  
#init_log.body.size_l1_s0    2.622032           2.490479           2.670063  
#  
#Log-likelihood:    36.937  
#AIC :   -63.874  
#AICc:   -61.147  
#BIC :   -59.708  
```



## Extra:
Take a further look at OU+linear trend:

Look at Bayesian analysis, in order to get process inference:
```r
X.mod3=layer.series.structure(X, numlayers=1,
		lin.time=T, prior=p, init.time=1996) 
mod3=layer.analyzer(X.mod3, num.MCMC=1000, burnin=10000,spacing=10,num.temp=6,
  do.model.likelihood = FALSE, 
  smoothing.specs=list(do.smoothing=TRUE, smoothing.time.diff=1, 
                       smoothing.start=1996, smoothing.end=2020.5,
		       num.smooth.per.mcmc=10))
```

Make plot of measurements, process inference estimate and
process inference uncertainty:

```r
png("malta.png",height=800,width=100)
par(cex=1.2)

## Plot data points:
plot(X$time, X$value,type="b",ylim=c(2.2,2.7),xlim=c(1996,2020),
    xlab="year",ylab="log(size)")

## Plot measurement uncertainties:
for(i in 1:length(X$time))
  lines(c(X$time[i],X$time[i]),
        c(X$value[i]-1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i]),
          X$value[i]+1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i])))
```

Plot process expectation value as a function of time:
PS: The 'mu' parameter is in this case the intercept of the
line at the mean measurement time.
PSS: The expected value will lag one characteristic time behind
the line that it is tracking, thus the subtraction of the characteristic 
time here.

```r
abline(c(models.ml[[3]]$mu_log.body.size$ML-models.ml[[3]]$lin_t_log.body.size$ML*((max(X$time)+min(X$time))/2 + models.ml[[3]]$dt_log.body.size_1$ML), models.ml[[3]]$lin_t_log.body.size$ML),lwd=3)

lines(mod3$process.time.points,mod3$process.mean, col="red",lwd=3)
lines(mod3$process.time.points,mod3$process.lower95, col="green",lwd=3)
lines(mod3$process.time.points,mod3$process.upper95, col="green",lwd=3)

dev.off()
```


## Look at residuals:

Define the linear trend OU structure:
```r
X.mod3=layer.series.structure(X, numlayers=1,
		lin.time=T, prior=p, init.time=1996)
```

Look at classic  analysis again:
```r
mod3=layer.analyzer(X.mod3, num.MCMC=1000, burnin=10000,spacing=10,num.temp=6,
  do.model.likelihood = TRUE,maximum.likelihood.numstart=100,
  return.residuals=TRUE)
```

Look for auto-correlation:
Raw  plot:
```r
plot(mod3$residuals.time, 
  mod3$standardized.residuals,type="b",xlab="Time", ylab="Residuals")
```

Partial auto-correlation plot:
```r
pacf(mod3$standardized.residuals)
```

Rest for autocorrelation as linear regression analysis:
```r
n=length(mod3$standardized.residuals)
summary(lm(mod3$standardized.residuals[2:n]~mod3$standardized.residuals[1:(n-1)]))
```
Output:  
```r
#Coefficients:  
#                                       Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                             0.10778    0.22092   0.488    0.633  
#mod3$standardized.residuals[1:(n - 1)] -0.07367    0.26731  -0.276    0.787  
```
  

Is there time trend for residuals?
```r
t=X$time-2000
summary(lm(mod3$standardized.residuals~t))
```
Output:  
```r
#t           -0.003547   0.035994  -0.099    0.923  
```
So, no linear trend detected.  
  
Look at non-linear trends using GAM:
```r
library(mgcv)
summary(gam(mod3$standardized.residuals~s(t)))
```
Output:  
```r
#Approximate significance of smooth terms:  
#       edf Ref.df     F p-value  
#s(t) 1.474  1.807 0.206   0.778  
```
  
No non-linear trend detected either


Look for whether the distribution is normal:
```r
qqnorm(mod3$standardized.residuals)
qqline(mod3$standardized.residuals)
shapiro.test(mod3$standardized.residuals)
```
Output: 
```r 
#  p-value = 0.2262  
```
Looks like normality is okay.  
  

Look on prior expected values vs residuals for each series:
```r
x1=mod3$prior.expected.values[,1]
plot(x1,mod3$standardized.residuals)
summary(lm(mod3$standardized.residuals~x1))
```
Output:  
```r
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    4.125      6.931   0.595    0.561  
# x1            -1.706      2.903  -0.588    0.565  
```
  
No discernable dependency between prior expected values and
standardized residuals.
