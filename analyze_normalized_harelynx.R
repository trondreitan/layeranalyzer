# Install the package. Check if this goes well before proceeding!
# Might require installation of Rcpp first.
#install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T)
#install.packages("https://folk.uio.no/trondr/R/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T)

# If using devtools:
# install_github(repo="trondreitan/layeranalyzer")


# Normalize original data (found in "hare.txt" and "lynx.txt")
# Send to files "hare_norm.txt" and "lynx_norm.txt":
source("normalize_harelynx.R")
# This also provides us with functions for transforming back:
# "invtrans.hare" and "invtrans.lynx"

# If this went well, load the package:
library(layeranalyzer)

##################################
# Read data and define prior:
##################################

hare.norm=read.layer.data.series("http://folk.uio.no/trondr/layered/hare_norm.txt",sep=" ",name="hare")
lynx.norm=read.layer.data.series("http://folk.uio.no/trondr/layered/lynx_norm.txt",sep=" ",name="lynx")

# Plot data:
#plot(hare.norm$time,hare.norm$value,type="b",ylim=c(-3,3))
#lines(lynx.norm$time,lynx.norm$value,type="b",col="red")

# Set priors:
p.hare=layer.prior(mu=c(-1,1), init=c(-5,5),
  dt=c(0.5,5),sigma=c(0.01,2),obs=c(0.01,1),lin=c(-0.1,0.1), beta=c(-2,2))
p.lynx=layer.prior(mu=c(-1,1), init=c(-5,5),
  dt=c(0.5,20),sigma=c(0.01,2),obs=c(0.01,1),lin=c(-0.1,0.1), beta=c(-2,2))

# Initial standalone structure analysis. Only mentioned in the main text.
models.hares=traverse.standalone.layered(hare.norm, max.layers=2, 
  talkative=TRUE, allow.one.feedback.loop=TRUE, 
  just.stationary=FALSE, no.rw=FALSE,    
  time.integrals.possible=FALSE, 
  allow.deterministic.layers=TRUE,
  use.stationary.stdev = FALSE, mcmc=TRUE,
  num.MCMC=1000,spacing=10,burnin=2000, num.temp = 4, prior=p.hare)
compare.layered(models.hares)

summary(models.hares[[5]])
#[1] "summary called"
#Coefficients:
#                        Mean    Median Lower 95% Upper 95%
#mu_hare             0.073010  0.083537 -0.378090  0.480828
#dt_hare_1           2.323119  2.151224  0.779086  4.912041
#sigma_hare_1        0.402682  0.408121  0.103110  0.696254
#dt_hare_2           4.086392  3.877787  1.981780  7.808524
#sigma_hare_2        0.981927  1.005536  0.180366  1.615923
#obs_sd_hare         0.211314  0.210086  0.071823  0.361114
#init_hare_l1_s0     1.737514  1.730336  1.308373  2.197310
#init_hare_l2_s0     1.341241  1.459415 -1.768511  3.844231
#beta_hare,1_hare,2 -1.419770 -1.391591 -3.027283  0.038752
#complex.eigen       0.934000  1.000000  0.000000  1.000000
#cycle01            18.857678 16.204408 10.588717 45.206529
#Model log-likelihood:   -71.079

# Initial standalone structure analysis. Only mentioned in the main text.
models.lynx=traverse.standalone.layered(lynx.norm, max.layers=2, 
  talkative=TRUE, allow.one.feedback.loop=TRUE, 
  just.stationary=FALSE, no.rw=FALSE,    
  time.integrals.possible=FALSE, 
  allow.deterministic.layers=TRUE,
  use.stationary.stdev = FALSE, mcmc=TRUE,
  num.MCMC=1000,spacing=10,burnin=2000, num.temp = 4, prior=p.lynx)
compare.layered(models.lynx)
summary(models.lynx[[13]])




# Set priors sufficiently for connection analysis:
p.hare=layer.prior(mu=c(-1,1), 
  dt=c(0.5,5),sigma=c(0.01,2),obs=c(0.01,1),beta=c(-2,2))
p.lynx=layer.prior(mu=c(-1,1), 
  dt=c(0.5,20),sigma=c(0.01,2),obs=c(0.01,1),beta=c(-2,2))


# Set structure to one-layered, as there is suspicition that the
# extra layer for hares are the lynx and the extra layer for lynx are
# the hares.
hare.struct=layer.series.structure(hare.norm, numlayers=1, prior=p.hare)
lynx.struct=layer.series.structure(lynx.norm, numlayers=1, prior=p.lynx)


# Traverse connection models between hares and lynx:
res=traverse.connections.layered(hare.struct,lynx.struct, 
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=10,T.ground=1.1,
 silent.mode=FALSE)

# Model comparison (using Bayesian model probabilities):
compare.layered(res)
#           log(lik) Post. Prob.(%)
#Model   1 -109.6261        0.02093
#Model   2 -103.8030        7.07560
#Model   3 -108.9128        0.04272
#Model   4 -101.2288       92.83724
#Model   5 -109.5099        0.02351

# Overwhelming support for model 4. Take a look at this:
summary(res[[4]])
#Coefficients:
#                        Mean    Median Lower 95% Upper 95%
#mu_hare            -0.026388 -0.026703 -0.426835  0.380152
#dt_hare_1           3.759316  3.581145  2.193782  6.320872
#sigma_hare_1        0.647815  0.644069  0.532804  0.800031
#obs_sd_hare         0.068970  0.050807  0.006233  0.224773
#mu_lynx            -0.021962 -0.019578 -0.386817  0.360592
#dt_lynx_1           5.138756  4.883584  2.888825  8.923214
#sigma_lynx_1        0.462834  0.456090  0.358721  0.592182
#obs_sd_lynx         0.049013  0.039922  0.006100  0.135660
#beta_hare,1_lynx,1  1.850127  1.770549  0.954505  3.046231
#beta_lynx,1_hare,1 -1.100573 -1.062133 -2.092366 -0.332004
#complex.eigen       0.998000  1.000000  1.000000  1.000000
#cycle01            21.576812 19.562330 13.378333 38.033484
#
#Model log-likelihood:  -101.229

# Two-way causal link, positive from hares to lynx,
# negative from lynx to hares. This seems reasonable!

# Run model with predictions forward in time:
res4=layer.analyzer(hare.struct,lynx.struct,
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=5,T.ground=1.1,
  causal=cbind(c(2,1,1,1),c(1,1,2,1)), 
  smoothing.specs=list(do.smoothing=TRUE, smoothing.time.diff = 0.1,
  smoothing.start = 1850, smoothing.end = 2018, 
  num.smooth.per.mcmc = 100, do.return.smoothing.samples = FALSE),
  return.residuals=TRUE)
 



# Make plot, using original scale:

# Plot hare data:

plot(hare.orig$year,hare.orig$val,xlim=c(1850,1960),
  log="y",ylim=c(500,200000),xlab="Year",ylab="#hares")
# Transform hare process inference back to original scale:
invtrans.hare.mean=res4$process.mean[1,]
for(i in 1:length(invtrans.hare.mean))
  invtrans.hare.mean[i]=exp(invtrans.hare(res4$process.mean[1,i]))
lines(res4$process.time.points, invtrans.hare.mean)
invtrans.hare.lower=res4$process.lower95[1,]
for(i in 1:length(invtrans.hare.mean))
  invtrans.hare.lower[i]=exp(invtrans.hare(res4$process.lower95[1,i]))
lines(res4$process.time.points, invtrans.hare.lower,col="red")
invtrans.hare.upper=res4$process.upper95[1,]
for(i in 1:length(invtrans.hare.mean))
  invtrans.hare.upper[i]=exp(invtrans.hare(res4$process.upper95[1,i]))
lines(res4$process.time.points, invtrans.hare.upper,col="red")

plot(hare.orig$year,hare.orig$val)
lines(res4$process.time.points, invtrans.hare.mean)
lines(res4$process.time.points, invtrans.hare.lower,col="red")
lines(res4$process.time.points, invtrans.hare.upper,col="red")

# Lynx:

plot(lynx.orig$year,lynx.orig$val,xlim=c(1850,1960),
  log="y",ylim=c(100,10000),xlab="Year",ylab="#lynx")
# Transform lynx process inference back to original scale:
invtrans.lynx.mean=res4$process.mean[2,]
for(i in 1:length(invtrans.lynx.mean))
  invtrans.lynx.mean[i]=exp(invtrans.lynx(res4$process.mean[2,i]))
lines(res4$process.time.points, invtrans.lynx.mean)
invtrans.lynx.lower=res4$process.lower95[2,]
for(i in 1:length(invtrans.lynx.mean))
  invtrans.lynx.lower[i]=exp(invtrans.lynx(res4$process.lower95[2,i]))
lines(res4$process.time.points, invtrans.lynx.lower,col="red")
invtrans.lynx.upper=res4$process.upper95[2,]
for(i in 1:length(invtrans.lynx.mean))
  invtrans.lynx.upper[i]=exp(invtrans.lynx(res4$process.upper95[2,i]))
lines(res4$process.time.points, invtrans.lynx.upper,col="red")

plot(lynx.orig$year,lynx.orig$val,log="y")
lines(res4$process.time.points, invtrans.lynx.mean)
lines(res4$process.time.points, invtrans.lynx.lower,col="red")
lines(res4$process.time.points, invtrans.lynx.upper,col="red")







##########################
# Look at residuals:
##########################


# Run model with predictions forward in time:
res4=layer.analyzer(hare.struct,lynx.struct,
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=1,T.ground=1.1,
  causal=cbind(c(2,1,1,1),c(1,1,2,1)),
  return.residuals=TRUE)
 
res.hare=res4$standardized.residuals[,1]
res.lynx=res4$standardized.residuals[,2]
t=res4$residuals.time

res.hare2=res.hare[!is.na(res.hare)]
t.hare=t[!is.na(res.hare)]-1900

res.lynx2=res.lynx[!is.na(res.lynx)]
t.lynx=t[!is.na(res.lynx)]-1900

# Look for auto-correlation:
# Raw plot:
plot(t, res.hare,type="b",xlab="Time", ylab="Hare residuals")
plot(t, res.lynx,type="b",xlab="Time", ylab="Lynx residuals")

# Partial auto-correlation plot:
pacf(res.hare2)
pacf(res.lynx2)
# Indication of auto-correlation in lynx residuals

# test for autocorrelation as linear regression analysis:
n.hare=length(t.hare)
summary(lm(res.hare2[2:n.hare]~res.hare2[1:(n.hare-1)]))
#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)
#(Intercept)               -0.01748    0.11350  -0.154    0.878
#res.hare2[1:(n.hare - 1)]  0.12511    0.11853   1.056    0.295

n.lynx=length(t.lynx)
summary(lm(res.lynx2[2:n.lynx]~res.lynx2[1:(n.lynx-1)]))
# res.lynx2[1:(n.lynx - 1)]  0.55608    0.13010   4.274  0.00012 ***

# time trend for residuals?
summary(lm(res.hare2~t.hare))
# t.hare      -0.005578   0.005552  -1.005    0.319

summary(lm(res.lynx2~t.lynx))
# t.lynx      -0.02514    0.01208  -2.082   0.0438 *


library(mgcv)
summary(gam(res.hare2~s(t.hare)))
#          edf Ref.df     F p-value
#s(t.hare)   1      1 1.009   0.319

summary(gam(res.lynx2~s(t.lynx)))
#          edf Ref.df     F p-value  
#s(t.lynx)   1      1 4.235   0.046 *

# Indication of autocorrelation or time trend in lynx residuals, but not
# hare residuals.

# look for whether the distribution is normal:
qqnorm(res.hare2)
qqline(res.hare2)
shapiro.test(res.hare2)
# p-value = 0.07761

qqnorm(res.lynx2)
qqline(res.lynx2)
shapiro.test(res.lynx2)
# p-value = 0.4673







# Try another model, two layers for lynx:
hare.struct2=layer.series.structure(hare.norm, numlayers=1, prior=p.hare)
lynx.struct2=layer.series.structure(lynx.norm, numlayers=2, prior=p.lynx)

# Run model with predictions forward in time:
res6=layer.analyzer(hare.struct2,lynx.struct2,
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=1,
  causal=cbind(c(2,1,1,1),c(1,1,2,1)), 
  return.residuals=TRUE)
  
compare.layered(res4,res6)
#[1] "compare called"
#            log(lik) Post. Prob.(%)
#Model   1 -101.24052        0.30896
#Model   2  -95.46389       99.69104
# 2-layered lynx process gives a better model for the same connections!

# Going through the same residual tests reveals no
# time trend/autocorrelation in the lynx residuals now.

res.hare=res6$standardized.residuals[,1]
res.lynx=res6$standardized.residuals[,2]
t=res6$residuals.time

res.hare2=res.hare[!is.na(res.hare)]
t.hare=t[!is.na(res.hare)]-1800

res.lynx2=res.lynx[!is.na(res.lynx)]
t.lynx=t[!is.na(res.lynx)]-1800

# Look for auto-correlation:
# Raw plot:
plot(t, res.hare,type="b",xlab="Time", ylab="Residuals")
plot(t, res.lynx,type="b",xlab="Time", ylab="Residuals")

# Partial auto-correlation plot:
pacf(res.hare2)
pacf(res.lynx2)
# Indication of auto-correlation in lynx residuals

# test for autocorrelation as linear regression analysis:
n.hare=length(t.hare)
summary(lm(res.hare2[2:n.hare]~res.hare2[1:(n.hare-1)]))
# res.hare2[1:(n.hare - 1)]  0.12812    0.11845   1.082    0.283

n.lynx=length(t.lynx)
summary(lm(res.lynx2[2:n.lynx]~res.lynx2[1:(n.lynx-1)]))
# res.lynx2[1:(n.lynx - 1)]  0.55725    0.13001   4.286 0.000115 ***

# time trend for residuals?
summary(lm(res.hare2~t.hare))
# t.hare      -0.005746   0.005568  -1.032    0.306

summary(lm(res.lynx2~t.lynx))
# t.lynx      -0.02533    0.01219  -2.079   0.0441 *


library(mgcv)
summary(gam(res.hare2~s(t.hare)))
#          edf Ref.df     F p-value
#s(t.hare)   1      1 1.065   0.306

summary(gam(res.lynx2~s(t.lynx)))
#          edf Ref.df     F p-value
#s(t.lynx)   1      1 4.322   0.044 *

# Indication of autocorrelation/time trend in lynx residuals, but not
# hare residuals.

# look for whether the distribution is normal:
qqnorm(res.hare2)
qqline(res.hare2)
shapiro.test(res.hare2)
# p-value = 0.07761

qqnorm(res.lynx2)
qqline(res.lynx2)
shapiro.test(res.lynx2)
# p-value = 0.4729






# Traverse connection models between hares and lynx:
res.new=traverse.connections.layered(hare.struct2,lynx.struct2, 
  num.MCMC=1000, burnin=10000, spacing=10,use.stationary.stdev=TRUE,
  silent.mode=FALSE)

# Model comparison (using Bayesian model probabilities):
compare.layered(res.new)
#            log(lik) Post. Prob.(%)
#Model   1 -100.76421        0.04633
#Model   2  -99.29643        0.20104
#Model   3 -101.30938        0.02686
#Model   4  -99.90567        0.10932
#Model   5  -98.38486        0.50023
#Model   6  -97.90837        0.80558
#Model   7  -98.01382        0.72496
#Model   8  -98.33383        0.52642
#Model   9  -98.21209        0.59457
#Model  10  -98.42398        0.48104
#Model  11 -100.12134        0.08811
#Model  12  -96.28704        4.07608
#Model  13  -96.99870        2.00066
#Model  14  -94.26391       30.82312
#Model  15  -98.01259        0.72585
#Model  16  -95.48237        9.11395
#Model  17  -94.98027       15.05799
#Model  18  -95.43251        9.57985
#Model  19  -94.67187       20.49762
#Model  20  -96.58559        3.02402
#Model  21 -100.88254        0.04116
#Model  22  -99.45364        0.17179
#Model  23 -100.43915        0.06412
#Model  24  -99.65768        0.14008
#Model  25  -98.23818        0.57926


# Best model, model 14 (followed by 17 and 19, then 18 and 16 (same as res6))
# Re-create it (in order that later runs do not need to re-do the
# connection traversal):
res14=layer.analyzer(hare.struct2,lynx.struct2, 
  num.MCMC=1000, burnin=10000, spacing=10,use.stationary.stdev=TRUE,
  causal=cbind(c(1,1,2,2),c(2,2,1,1),c(2,1,1,1)), return.residuals=TRUE)
  


summary(res14)
#Coefficients:
#                        Mean    Median Lower 95% Upper 95%
#mu_hare             0.067365  0.069041 -0.698753  0.856653
#dt_hare_1           2.373173  2.244645  1.425098  3.930496
#sd_hare_1           0.656256  0.646259  0.484738  0.892410
#obs_sd_hare         0.081879  0.064623  0.007445  0.249914
#mu_lynx            -0.013805 -0.015352 -0.397429  0.382947
#dt_lynx_1           1.448122  1.427379  0.812462  2.271165
#sd_lynx_1           0.073013  0.058676  0.005698  0.212808
#dt_lynx_2           2.728232  2.591797  1.721150  4.945554
#sd_lynx_2           1.181581  1.134244  0.765235  1.893839
#obs_sd_lynx         0.050966  0.042602  0.005949  0.135427
#beta_hare,1_lynx,2  1.097380  1.080033  0.289356  2.019172
#beta_lynx,2_hare,1  0.824386  0.816397  0.097219  1.516054
#beta_lynx,1_hare,1 -1.449238 -1.414978 -2.325471 -0.726898
#complex.eigen       0.997000  1.000000  1.000000  1.000000
#cycle01            19.026284 17.057398 12.609293 33.171597

# Interpretation. The second lynx layer is a mixture of everything
# that affects lynx reqruitment; environment including the number of hares
# available. The environment (excluding the hares) also affects the hare
# numbers in a similar way as it affects the lynx numbers.

# Make an even bigger model, with 3-layered lynx process,
# where we hope to separate environment
# (layer 3) and lynx requirement (layer 2):

hare.struct3=layer.series.structure(hare.norm, numlayers=1, prior=p.hare)
lynx.struct3=layer.series.structure(lynx.norm, numlayers=3, prior=p.lynx)

res26=layer.analyzer(hare.struct3,lynx.struct3,
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=1,
  causal=cbind(c(2,3,1,1),c(1,1,2,2),c(2,1,1,1)), 
  return.residuals=TRUE)

compare.layered(res.new[[14]], res26)
#           log(lik) Post. Prob.(%)
#Model   1 -94.28952         77.853
#Model   2 -95.54664         22.147
# Didn't work

# Also test 2 layers for hares:
hare.struct4=layer.series.structure(hare.norm, numlayers=2, prior=p.hare)
lynx.struct4=layer.series.structure(lynx.norm, numlayers=2, prior=p.lynx)

res27=layer.analyzer(hare.struct4,lynx.struct4,
  num.MCMC=1000, burnin=10000, spacing=10,num.temp=1,
  causal=cbind(c(2,2,1,1),c(1,1,2,2),c(2,1,1,1)), 
  return.residuals=TRUE)

compare.layered(res.new[[14]], res27)
#           log(lik) Post. Prob.(%)
#Model   1 -94.28952       80.54316
#Model   2 -95.71011       19.45684
# Didn't work


