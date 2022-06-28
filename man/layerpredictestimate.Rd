\name{layer.predict.estimate}
\alias{layer.predict.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Prediction based on previously performed analysis using single point estimate for the parameter set }
\description{
 When a linear SDE analysis has been performed, predictions for 
 the estimated parameter set can be made 
 afterwards, using 'layer.predict.estimate'. The data sent does not need to
 be the same as for what the parameter estimates in 'layer.analyzer' was
 trained on. It could be a new dataset or the old with some new data added
 or a few data points removed. The predictions are the mean and upper/lower 
 95% credibility interval for the process state using the Kalman smoother.
 The individual covariance matrices for each time point are also returned.
}
\usage{
layer.predict.estimate(... ,
   analysis=NULL, 
   smoothing.time.diff=0,
   smoothing.start=NULL,smoothing.end=NULL,
  return.residuals=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
 A set of 'layer.data.series' objects that represents the set of 
 time series that are to be predicted on. Note that these must match
 the ones used for running 'layer.analyzer' when it comes to number of
 time series and order. (See 'layer.data.series' for how they are defined). 
 }
 \item{analysis}{
 The 'analysis' variable must be a 'layer.analyzer' object. I.e. a 
 'layer.analyzer' analysis needs to be performed first. This variable
 is not an option. If it is NULL, the function will complain.
 }
 \item{smoothing.time.diff}{
  If 'smoothing.time.diff' is set to "0",
  only the states at the measurement times will be inferred. (PS:
  you can get the function to perform process state inference for 
  unmeasured time points by injecting them into the time series with'
  missing values). If 'smoothing.time.diff' is set to a positive real value,
  the state at the regularly spaced time points 
  (smoothing.start,
   smoothing.start+smoothing.time.diff, 
   smoothing.start+smoothing.time.diff, 
   ...
   smoothing.end) 
   will also be reported. Note that this means that 'smoothing.start'
   and 'smoothing.end' must also be given. 
  }
 \item{smoothing.start}{
  Start time for extra smoothing points.
 }
 \item{smoothing.end}{
  End time for extra smoothing points.
 }
 \item(return.residuals}{
  The difference between predicted and actual measurements are returned.
  If set to 'TRUE', the standardized residuals for all observed series
 as a matrix called 'standardized.residuals'
 (one column for each measured series) and the returns the time points for
 the residuals ('residuals.time'). Also returns the prior expected values from
 the Kalman filter ('prior.expected.values'), which should be
 independent of the residuals if the model assumptions are correct.
 Note that the residuals are calculated for the Bayesian median
 parameters if Bayesian analysis is used and maximum
 likelihood estimated parameters if classic (ML)
 analysis is used. Also note that the residuals for each observed
 series is the measurement minus the expected value divided by the
 uncertainty from the Kalman *filter* not the Kalman smoother (which
 is used for process inference), as it is the deviations in the
 Kalman filter that are assumed to be independent from measurement
 to measurement.
 }
}
\details{
 This function generates Kalman smoother-based process state samples
 for the estimated parameters in the previously performed analysis (returned
 in the 'est.origpar' variable. This is 
 performed in the same way as when the 'smoothing.specs' option
 is actived in 'layer.analyzer' and the same underlying C++ routines
 is utilized (though there is a fork in the code for when MCMC samples
 are predefined rather than in need of sampling). 

 PS: This uses the 'layer.predict.mcmc.list' function using
 'num.smooth.per.mcmc=1' and only one MCMC sample which is set ot the 
 'est.origpar. variable. This represents the parameter set in the 
 original parametrization, i.e. the parametrization used internally. 

}
\value{
Returns an object of type 'layered', just as 'layer.analyzer' does. 
However, note that no new parameter sampling or estimation is performed.
The variables of note are those that have to do with the process state
inference, so only those will be described here:

process.time.points, process.mean, process.lower95, process.upper95: Process time points, means, upper and lower credibility bands for each process. (PS: There can be multiple prcoesses per serie, since number of layers can be greater than one).  

P.k: A list of process state covariance matrices returned from the Kalman 
  smoother, one for each time point (see 'process.time.points').

standardized.residuals: Standardized residuals are measurements minus expected values from the Kalman filter, divided by the standard deviation also from the Kalman filter. Only returned if "return.residuals=TRUE" is set.

prior.expected.values: Expected values from the Kalman filter. Only returned if "return.residuals=TRUE" is set.

residuals.time: Time points for each residual. Only returned if "return.residuals=TRUE" is set.

}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.
}
\author{
Trond Reitan, trond.reitan@ibv.uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also layer.analyzer, layer.predict.mcmc.list and layer.data.series.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Define time series (time points, values and in this case also
# sample standard deviations and sample size, used for calculating
# standard errors for each measurement). 
X=layer.data.series(time.points=malta$Time.Year,
  value.points=malta$Mean..log.body.mass.,
  std.dev=sqrt(malta$Variance.calculated.from.the.data), 
  num.meas.per.value=malta$Sample.size,name="log.body.size")

# Define process structure for analysis purposes:
ser=layer.series.structure(X,numlayers=1)

# Perform Bayesian analysis:
res=layer.analyzer(ser,mcmc=TRUE)

# Now make up some new values much later:
t2=2100:2199
v2=rnorm(100,2+0.1*sin(2*pi*t2/50),0.02)
X2=layer.data.series(time.points=t2,
   value.points=v2, std.dev=rep(0.02,100), num.meas.per.value=rep(1,100),
   name="log.body.size")

# Perform predictions:
p3=layer.predict.estimate(X2, analysis=res,  smoothing.time.diff=0.25,
   smoothing.start=1500,smoothing.end=2800)

# Fetch the results
t3=p3$process.time.points
m3=p3$process.mean[1,]
l3=p3$process.lower95[1,]
u3=p3$process.upper95[1,]

# Plot the results
plot(t3,m3,type="l",ylim=c(min(l3),max(u3)))
points(t2,v2)
lines(t3,l3,col="red")
lines(t3,u3,col="red")


# Perform ML analysis:
# Do the same for ML estimate 
res.ml=layer.analyzer(ser,do.maximum.likelihood=TRUE,
  maximum.likelihood.numstart=100,do.model.likelihood=FALSE)

# Now make up some new values much later:
t2=2100:2199
v2=rnorm(100,2+0.1*sin(2*pi*t2/50),0.02)
X2=layer.data.series(time.points=t2,
   value.points=v2, std.dev=rep(0.02,100), num.meas.per.value=rep(1,100),
   name="log.body.size")

# Perform the prediction:
p4=layer.predict.estimate(X2, analysis=res.ml,  smoothing.time.diff=0.025,
   smoothing.start=1500,smoothing.end=2800)

# Fetch the results:
t4=p4$process.time.points
m4=p4$process.mean[1,]
l4=p4$process.lower95[1,]
u4=p4$process.upper95[1,]

# Plot the results:
plot(t4,m4,type="l",ylim=c(min(c(l4,m4)),max(c(m4,u4))))
points(t2,v2)
lines(t4,l4,col="red")
lines(t4,u4,col="red")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, prediction, predictions }


