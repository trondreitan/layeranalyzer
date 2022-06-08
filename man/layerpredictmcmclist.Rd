\name{layer.predict.mcmc.list}
\alias{layer.predict.mcmc.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Prediction based on previously performed MCMC analysis }
\description{
 When a linear SDE analysis has been performed using 'layer.analyzer' with
 the MCMC samples saved (the option 'mcmc=TRUE') or based on estimates, 
 predictions can be made 
 afterwards, using 'layer.predict.mcmc.list'. The data sent does not need to
 be the same as for what the parameter estimates in 'layer.analyzer' was
 trained on. It could be a new dataset or the old with some new data added
 or a few data points removed. The predictiong runs through the MCMC 
 samples or just a single estimate.
 This is a way for the functions 'layer.predict.mcmc' and 
 'layer.predict.estimate' to access the prediction module, and 
 is made for those rather than direct use.
}
\usage{
layer.predict.mcmc.list(new.data.list ,
   analysis=NULL, 
   smoothing.time.diff=0,
   smoothing.start=NULL,smoothing.end=NULL,
   num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new.data.list}{
 A list item of 'layer.data.series' objects that represents the set of 
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
 \item{num.smooth.per.mcmc}{
   This is the number of smoothing process state samples fetched 
   per MCMC sample.
 }
 \item{do.return.smoothing.samples}{
  'do.return.smoothing.sample' determines if the
  raw smoothing samples are to be returned or not, or only the 
  mean and lower/upper 95% credibility band of the processes. 
  These represents samples of the individual process states
  given the new data. Note that these are not process realizations,
  as the autocorrelation between one time point and the next is
  not taken into account when sampling. Process realizations
  may perhaps be added later.
 }
}
\details{

 This function generates Kalman smoother-based process state samples
 for each MCMC sample in the previously performed analysis. This is 
 performed in the same way as when the 'smoothing.specs' option
 is actived in 'layer.analyzer' and the same underlying c++ routines
 is utilized (though there is a fork in the code for when MCMC samples
 are predefined rather than in need of sampling). 

 PS: If 'num.smooth.per.mcmc=1' and there is only one MCMC sample, the
 routine will assume that a direct Kalman smoother result is instead required,
 based on the estimated parameter set.

}
\value{
Returns an object of type 'layered', just as 'layer.analyzer' does. 
However, note that no new parameter sampling or estimation is performed.
The variables of note are those that have to do with the process state
inference, so only those will be described here:

process.time.points, process.mean, process.lower95, process.upper95: Process time points, means, upper and lower credibility bands for each process. (PS: There can be multiple prcoesses per serie, since number of layers can be greater than one).  

smoothing.samples: Smoothing samples for each MCMC sample for each process. Only returned if option "do.return.smoothing.samples = TRUE)" is used.

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
See also layer.analyzer and layer.data.series.
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

# Perform analysis and make sure MCMC samples are returned:
res=layer.analyzer(ser,mcmc=TRUE)

# Now make up some new values much later:
t2=2100:2199
v2=rnorm(100,2+0.1*sin(2*pi*t2/50),0.02)

X2=layer.data.series(time.points=t2,
   value.points=v2, std.dev=rep(0.02,100), num.meas.per.value=rep(1,100),
   name="log.body.size")

# Perform predictions:
p2=layer.predict.mcmc.list(list(X2=X2), analysis=res,  
   smoothing.time.diff=0.25,
   smoothing.start=1500,smoothing.end=2800,
   num.smooth.per.mcmc=10, do.return.smoothing.samples=TRUE)

# Plot the results:
plot(p2$process.time.points,p2$process.mean[1,],type="l")
points(t2,v2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, prediction, predictions }


