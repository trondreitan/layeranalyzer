
# Install the package. Check if this goes well before proceeding!
# Might require installation of Rcpp first.
#install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T)
#install.packages("http://folk.uio.no/trondr/R/layeranalyzer_0.1.0.tar.gz", type="source",verbose=T)

# If this went well, load the package:
library(layeranalyzer)

##################################
# Read data:
##################################
bi.lext=read.layer.data.series("http://folk.uio.no/trondr/layered/bivalve_lext.txt",name="bivalve_lext",
	column.type=c("time","value","std.dev"))
br.lext=read.layer.data.series("http://folk.uio.no/trondr/layered/brach_lext.txt",name="brach_lext",
	column.type=c("time","value","std.dev"))
br.lorig=read.layer.data.series("http://folk.uio.no/trondr/layered/brach_lspec.txt",name="brach_lorig",
	column.type=c("time","value","std.dev"))

# Prior taken from Reitan&Liow 2017 work area:
lrate.pr=layer.load.prior("http://folk.uio.no/trondr/layered/lrate_prior.txt")


# 2-layers for extinction rates, 1 layer for brachipod origination rate:

bi.lext.struct2=layer.series.structure(bi.lext, numlayers=2,prior=lrate.pr)
br.lext.struct2=layer.series.structure(br.lext, numlayers=2,prior=lrate.pr)
br.lorig.struct1=layer.series.structure(br.lorig, numlayers=1,prior=lrate.pr)


# Best model accoring to Reitan&Liow 2017:
res1=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1)),
  use.half.lives = TRUE,use.stationary.stdev = TRUE,
  id.strategy=4, num.MCMC=1000, burnin=10000, spacing=20,num.temp=1)

# Show parameter estimates:
summary(res1)
#Coefficients:
#                                        Mean    Median Lower 95%   Upper 95%
#mu_bivalve_lext                    -3.931377 -3.925454 -4.251513   -3.618493
#dt_bivalve_lext_1                   0.945811  0.534617  0.000537    3.488160
#sd_bivalve_lext_1                   0.429385  0.306571  0.001653    1.066686
#dt_bivalve_lext_2                  41.081101  2.324459  0.803862  213.039004
#sd_bivalve_lext_2                   0.725230  0.867102  0.018118    1.272027
#mu_brach_lext                      -3.197016 -3.157338 -5.556519   -1.853327
#dt_brach_lext_1                     0.065006  0.028888  0.000229    0.284622
#sd_brach_lext_1                     0.082091  0.042955  0.000519    0.329365
#dt_brach_lext_2                   539.566741 47.900687  1.883028 3593.267783
#sd_brach_lext_2                     0.858529  0.534871  0.224426    3.580408
#mu_brach_lorig                     -3.161368 -3.167983 -3.472933   -2.830841
#dt_brach_lorig_1                    3.393715  3.179215  1.492998    6.559050
#sd_brach_lorig_1                    0.151878  0.119108  0.004710    0.443508
#beta_bivalve_lext,1_brach_lext,1    0.842496  0.844089  0.656697    1.010209
#beta_bivalve_lext,1_brach_lorig,1   1.197954  1.163697  0.716228    1.883574
#
#Model log-likelihood:  -394.375


# Traverse the 24 neighbouring models:
res2=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res3=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,2,2,1),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res4=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,2),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res5=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,2,2,2),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res6=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(2,1,1,1),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res7=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(2,1,1,2),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res8=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(2,2,1,1),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res9=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(2,2,1,2),c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res10=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,3,1)),corr=cbind(c(1,1,2,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res11=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,3,1)),corr=cbind(c(1,1,2,2)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res12=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,3,1)),corr=cbind(c(1,2,2,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res13=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,3,1)),corr=cbind(c(1,2,2,2)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res14=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res15=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,2,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res16=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(3,1,1,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res17=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(3,1,1,2)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res18=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1)),corr=cbind(c(1,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res19=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1)),corr=cbind(c(1,2,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res20=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1),c(2,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res21=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1),c(2,2,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res22=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1),c(3,1,2,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res23=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1),c(3,1,2,2)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res24=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1)),corr=cbind(c(2,1,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

res25=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
  causal=cbind(c(1,1,2,1),c(1,1,3,1)),corr=cbind(c(2,2,3,1)),
  id.strategy=4, num.MCMC=100, burnin=1000, spacing=10,num.temp=1)

date()
t2=Sys.time()
t2-t1


# Show comparison of models:
compare.layered(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,
  res12,res13,res14,res15,res16,res17,res18,res19,res20,res21,res22,res23,
  res24,res25)
#               log(lik) Post. Prob.(%)
#Model   1 -398.8481       24.41074
#Model   2 -420.7084        0.00000
#Model   3 -400.4024        5.15850
#Model   4 -421.9462        0.00000
#Model   5 -401.1084        2.54642
#Model   6 -400.2710        5.88341
#Model   7 -416.3637        0.00000
#Model   8 -400.9774        2.90287
#Model   9 -401.4153        1.87350
#Model  10 -399.3313       15.05546
#Model  11 -401.7427        1.35044
#Model  12 -400.6714        3.94179
#Model  13 -401.4561        1.79846
#Model  14 -421.8246        0.00000
#Model  15 -421.7367        0.00000
#Model  16 -411.6940        0.00006
#Model  17 -408.6291        0.00138
#Model  18 -405.1369        0.04533
#Model  19 -406.9415        0.00746
#Model  20 -401.6512        1.47971
#Model  21 -400.3943        5.20065
#Model  22 -402.3845        0.71075
#Model  23 -403.2051        0.31286
#Model  24 -399.0217       20.52077
#Model  25 -400.1263        6.79942

# Model 1 best among neighbouring models. Model 24 and model 10 second
# and third best respectively.


# Save:
save.image(file="2lay_local_id4.RData")


# Perform residual analysis
res1=layer.analyzer(bi.lext.struct2,br.lext.struct2,br.lorig.struct1,
   causal=cbind(c(1,1,2,1),c(1,1,3,1)),
   use.half.lives = TRUE,use.stationary.stdev = TRUE,
   id.strategy=4, num.MCMC=1000, burnin=10000, spacing=20,num.temp=1,
   return.residuals=TRUE)
   
resid1=res1$standardized.residuals[,1]
resid2=res1$standardized.residuals[,2]
resid3=res1$standardized.residuals[,3]
t=res1$residuals.time

resid1n=resid1[!is.na(resid1)]
t1=t[!is.na(resid1)]

resid2n=resid2[!is.na(resid2)]
t2=t[!is.na(resid2)]

resid3n=resid3[!is.na(resid3)]
t3=t[!is.na(resid3)]

# Look for auto-correlation:
# Raw plot:
plot(t, resid1,type="b",xlab="Time", ylab="Bivalve extinction resididuals")
plot(t, resid2,type="b",xlab="Time", ylab="Brachiopod extinction resididuals")
plot(t, resid3,type="b",xlab="Time", ylab="Brachiopod origination resididual")

# Partial auto-correlation plot:
pacf(resid1n)
pacf(resid2n)
pacf(resid3n)
# Indication of auto-correlation in lynx resididuals

# test for autocorrelation as linear regression analysis:
n1=length(t1)
summary(lm(resid1n[2:n1]~resid1n[1:(n1-1)]))
#resid1n[1:(n1 - 1)] -0.03567    0.12564  -0.284    0.777
n2=length(t2)
summary(lm(resid2n[2:n2]~resid2n[1:(n2-1)]))
#resid2n[1:(n2 - 1)] -0.03370    0.12987  -0.260    0.796
n3=length(t3)
summary(lm(resid3n[2:n3]~resid3n[1:(n3-1)]))
#resid3n[1:(n3 - 1)]  -0.1353     0.1131  -1.196  0.23527

summary(lm(resid1n~t1))
#t1          -0.0001138  0.0008549  -0.133    0.894
summary(lm(resid2n~t2))
#t2           0.0001776  0.0008439   0.210    0.834
summary(lm(resid3n~t3))
#t3          -0.0004969  0.0007264  -0.684   0.4960  

library(mgcv)
summary(gam(resid1n~s(t1)))
#        edf Ref.df    F p-value
# s(t1) 4.313  5.312 1.697   0.142
summary(gam(resid2n~s(t2)))
#        edf Ref.df    F p-value
#s(t2) 5.004  6.106 1.61   0.151
summary(gam(resid3n~s(t3)))
#        edf Ref.df    F p-value
#s(t3)   1      1 0.468   0.496


# No sign of auto-correlation or time trend

# look for whether the distribution is normal:
qqnorm(resid1n)
qqline(resid1n)
shapiro.test(resid1n)
# p-value = 8.49e-07
qqnorm(resid2n)
qqline(resid2n)
shapiro.test(resid2n)
#p-value = 1.203e-06
qqnorm(resid3n)
qqline(resid3n)
shapiro.test(resid3n)
#p-value = 0.0115

# extinction rates very clearly not normal,
# origination rate also not normal

# Non-normality of process or non-normality of measurement noise?

# Too few measurements 1.5 standard deviations below prediction,
# too many measurements 2 standard eviation or moe above prediction

# This could perhaps suggest sampling bias?


# Look on prior expected values vs residuals for each series:
x1=res1$prior.expected.values[,1]
x1=x1[!is.na(resid1)]

plot(x1, resid1n, xlab="Prior expected value, bivalve extinction",
  ylab="Residuals, bivalve extinction")
summary(lm(resid1n~x1))
#           Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.9596     1.3776  -0.697    0.488
#x1           -0.2037     0.3564  -0.572    0.569

summary(gam(resid1n~s(x1)))
#                edf Ref.df     F p-value  
#s(x1)         4.268  5.127 2.178  0.0648 .



x2=res1$prior.expected.values[,2]
x2=x2[!is.na(resid2)]

plot(x2, resid2n, xlab="Prior expected value, brachiopod extinction",
  ylab="Residuals, brachiopod extinction")
  
summary(lm(resid2n~x2))
#           Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.6857     0.7972  -0.860    0.393
#x2           -0.2024     0.2559  -0.791    0.431

summary(gam(resid2n~s(x2)))
#        edf Ref.df     F p-value
#s(x2) 1.001  1.002 0.625   0.432


x3=res1$prior.expected.values[,3]
x3=x3[!is.na(resid3)]

plot(x3, resid3n)
summary(lm(resid3n~x3))
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.8259     0.7783  -1.061    0.292
#x3           -0.1765     0.2533  -0.697    0.488

summary(gam(resid3n~s(x3)))
#        edf Ref.df     F p-value
#s(x3) 5.202  6.319 1.419    0.21

# No discernable dependency between prior expected values and residuals




# Time points for top extinction residual events
index1=order(resid1n,decreasing=TRUE)
t1[index1][1:5]
# -252.1700   -0.0117   -2.5800  -66.0000 -445.2000
resid1n[index1][1:5]
# 4.126867 3.737644 2.173538 2.032815 1.646106

index2=order(resid2n,decreasing=TRUE)
t2[index2][1:5]
# -252.1700   -0.0117 -168.3000 -201.3000 -443.4000
resid2n[index2][1:5]
# 4.268964 3.726002 1.685901 1.509380 1.472323

shapiro.test(resid1n[index1[6:length(resid1n)]])

shapiro.test(resid2n[index2[6:length(resid2n)]])


bi.lext2=bi.lext
index.use=rep(0,0)
t1.1=t1[index1]
for(j in 1:length(t1.1))
 {
   i=which(t1.1[j]==bi.lext$time)
   if(length(i)==1)
     index.use=c(index.use, i)
}
bi.lext2$time=bi.lext$time[index.use]
bi.lext2$value=bi.lext$value[index.use]
bi.lext2$std.dev=bi.lext$std.dev[index.use]

bi.lext2$std.dev




