\name{layer.analyzer}
\alias{layer.analyzer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Time series analysis tool using linear layered SDEs. }
\description{
 Linear SDE analysis allowing for hidden layers and multiple data series
 with correlative or causal links between them. Returns parameter estimates
 and Bayesian model likelihoods or information criteria, depending
 on whether the analysis is Bayesian (default) or maximum likelihood (ML)
 based. In case of causal feedback loops, the matrix operations (including 
 eigenvalue decompositions) allows for complex numbers. In this case,
 cyclic behavior can be expected.
}
\usage{
layer.analyzer(... ,
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.model.likelihood=TRUE,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  maximum.likelihood.strategy="MCMC-from-model",
  silent.mode=TRUE,talkative.burnin=FALSE,talkative.likelihood=FALSE,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, # start.parameters=0,
  use.half.lives=FALSE, mcmc=FALSE,causal=NULL,causal.symmetric=NULL,corr=NULL,
  smoothing.specs=
     list(do.smoothing=FALSE,smoothing.time.diff=0,
          smoothing.start=NULL,smoothing.end=NULL,
          num.smooth.per.mcmc=10, do.return.smoothing.samples=FALSE),
  realization.specs=
     list(do.realizations=FALSE,num.realizations=1000,strategy="N",
          realization.time.diff=0,realization.start=NULL,realization.end=NULL),
  return.residuals=FALSE,
  loglik.laxness="high")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
 A set of 'layer.series.structure' objects that represents the set of 
 time series that are to be analyzed and how the underlying processes of
 each such series is to be modelles. (See 'layer.series.structure'). At
 least one such object must be given. 
 }
 \item{num.MCMC}{
 num.MCMC is an integer specifying the number of Markov chain Monte Carlo 
 (MCMC) samples that the analysis will rest on. Default 1000. More MCMC 
 samples are run internally, depending on the spacing (default 10) and 
 burnin (default 1000).
 }
 \item{spacing}{
 Specifies the number of MCMC iteration between each sample that is used
 (default 10). A higher number means the samples will be less correlated
 but requires proptionally more computer time.
 }
 \item{burnin}{
 Markov chain Monte Carlo sampling takes some iterations to converge to 
 the posterior parameter distribution. The burn-in phase is the set of
 iterations performed where the samples are then discarded, so as to
 avoid including samples before convergence. Convergence can be check
 by for instance 'gelman.diag' in the 'coda' package, and the MCMC
 samples themselves can be gotten by using the option 'mcmc'. Still,
 stability of the results can be quite a good enough measure for most
 purposes. Default is 10000 iterations. PS: The burn-in algorithm here
 includes two steps where the random walk MCMC variance is adjusted
 in order to optimize the efficiency of the algorithm.
 }
 \item{num.temp}{
 The number of tempering chains (default 1, meaning it is just the MCMC
 chain itself that is beeing run). Parallel tempering (Geyer 1991) is a 
 method for overcoming local optima and can thus increase the stability of 
 the results. The way it is done is that in addition to the MCMC chain for 
 the posterior distribution, parallel chains sampling from the distribution 
 proportional to exp(-log(likelihood*prior)/T) where T is a "temperature" 
 larger than one (this is in practise a "smoothed" version of the posterior), 
 is also sampled from. Swaps between neighbouring chains are then allowed,
 making it possible for the MCMC chain for the posterior distribution to
 find a new optima.

 As the code is not parallel, the computer time will increase proportional 
 to the number of tempering chains. The set of "temperatures" is set to 
 (1, T.ground, T.ground^2, ..., T.ground^(numtemp-1)),
 where 'T.ground' is default set to 1.5. If no swapping is done (this can be seen
 if 'silent.mode' is set to 'FALSE'), this "ground temperature" must be
 lowered, in order to facilitate the swaps.
 }
 \item{do.model.likelihood}{
 Logical variable. Only makes sense for Bayesian analysis (which is default). 
 If set to 'TRUE' (default), the Bayesian model (marginal) likelihood is 
 estimated using an importance sampler. If set to 'FALSE', this will not happen, 
 which saves time but means that Bayesian model comparison will not be possible.
 }
 \item{do.maximum.likelihood}{
 Logical variable. If set to 'FALSE' (default), Bayesian analysis is performed.
 If set to 'TRUE', a classic maximum likelihood (ML) is performed, where the
 starting points of a set of optimizations is drawns from the MCMC samples.
 (This means that Bayesian MCMC sampling is performed anyway, but that the results
 from this is not focused on in the summary). The number of optimizations
 is determined by option 'maximum.likelihood.numstart', which is default set to 10. 
 Estimation of the Bayesian model (marginal) likelihood is switched off when
 this option is used.
 }
 \item{maximum.likelihood.numstart}{
 The number of ML optimizations performed, it applicable. (If 'do.maximum.likelihood'
 is set to 'TRUE'.) Default:10. 
 }
 \item{maximum.likelihood.strategy}{
 In layeranalyzer, ML optimizations are per default
 (maximum.likelihood.strategy="MCMC-from-model") started from
 MCMC samples,a s it was found that completely randomly drawn starting
 points worked poorly for large models. However, MCMC convergence
 is no guarantee either, especially for models with multiple datasets and
 connections between them. 'maximum.likelihood.strategy' allows
 for using alternative startegies for optimization in connection
 models by MCMC sampling from a model that is the same in structure
 but without connections and using these samples as starting points
 for the optimization instead, which is specified by setting
 maximum.likelihood.strategy="MCMC-from-nonconnected". The third
 possibility is to set maximum.likelihood.strategy="MCMC-from-hybrid",
 which alternates betwen using MCMC samples from the connection model
 and from MCMC samples for the non-connected model.
 } 
 \item{silent.mode}{
 If set to 'FALSE', the routine will show a lot debug information while running.
 This can be useful for tempering purposes, since the number of tempering
 swaps are shown in this debug information.
 }
 \item{talkative.burnin}{
 If set to 'TRUE', this triggers the printing of more specific debug information
 concerning the burn-in phase.
 }
 \item{talkative.likelihood}{
 If set to 'TRUE', this triggers the printing of more specific debug information
 concerning the likelihood calculation phase. (PS: This means a lot of printing.)
 }
 \item{id.strategy}{
 Determines how to handle the identifiability problem in Bayesian multi-layer 
 analysis. As mentioned in the supplementary of Reitan et al. (2012), one can
 switch the characteristic times of neighbouring layers and reorganize the 
 stochastic contributions in such a way that the same top layer process takes 
 place. This can be solved by requiring that the lower the layer, the larger
 the characteristic time. However, simply requiring this means the prior
 distribution of the characteristic times get shortened for multi-layered 
 models, as compared to one-layered model. As long as the data is within the
 range of these shortened prior distributions, this puts the one-layer model
 at a disadvantage when it comes to model comparison. There are various
 solutions to this problem, and which solution is used is determined by
 the variable 'id.strategy'.

 The variable 'id.strategy' should be an integer between 0 and 4, default 2. 
 0 - No identification treatment. (Default) PS: Can and
     even should yield multimodal characteristic times
 1 - Keep upper characteristic time prior. Add
     lognormally to beneath-lying characteristic times.
 2 - Keep lower characteristic time prior. Substract
     lognormally to above-lying characteristic times.
 3 - Keep lower characteristic time prior. Cut
     depending on that on the above-lying characteristic times. 
 4 - Keep upper characteristic time prior. Cut
     depending on that on the below-lying characteristic times.
 }
 \item{use.stationary.stdev}{
 If set to 'TRUE', instead of reporting the stochastic contribution size
 in the stochastic differential equations, the stationary standard deviation
 is reported. The stationary standard deviation is the standard deviation of
 the process state at any given time after convergence, unconditional on 
 previous states. In a multi-layered model, this is to be interpreted as the 
 stationary standard deviation if this layer when let alone (with no other 
 process affecting it). This is calculated as s*sqrt(characteristic time/2), 
 where s is the stochastic contribution size. PS: This option does not make 
 any sense if the bottom layer is a Wiener process, as that process is not 
 stationary.
 }
 \item{T.ground}{
 This determines how the parallel tempering chains are to be defined. Parallel 
 tempering (Geyer 1991) is a method for overcoming local optima and can thus 
 increase the stability of the results. The way it is done is that in addition 
 to the MCMC chain for the posterior distribution, parallel chains sampling 
 from the distribution proportional to exp(-log(likelihood*prior)/T) where T 
 is a "temperature" larger than one (this is in practise a "smoothed" version 
 of the posterior), is also sampled from. Swaps between neighbouring chains 
 are then allowed, making it possible for the MCMC chain for the posterior 
 distribution to find a new optima.

 As the code is not parallel, the computer time will increase proportional 
 to the number of tempering chains. The set of "temperatures" is set to 
 (1, T.ground, T.ground^2, ..., T.ground^(numtemp-1)),
 where 'T.ground' is default set to 1.5. If no swapping is done (this can be seen
 if 'silent.mode' is set to 'FALSE'), this "ground temperature" must be
 lowered, in order to facilitate the swaps.
 }
 \item{use.half.lives}{
 If set to 'TRUE', reports half-lives rather than characteristic times.
 Half-lives are the time the auto-correlation of an OU process drops to 1/2,
 or the time it takes for a perturbation from the expected value to drop to
 1/2 the original size in expected value. (As opposed to dropping to exp(-1)
 for characteristic times.) Half-life=log(2)*characteristic time. The
 auto-correlation for an OU process (and thus for the inner dynamics of any 
 linear layer) is exp(-diff.time/characteristic.time)=exp(-log(2)*diff.time/half.life)=
 (1/2)^(diff.time/half.life), where diff.time is the time difference between
 two process states. 
 }
 \item{mcmc}{
 If set to 'TRUE', returns the MCMC samples as an 'mcmc' object within the return
 object. The name of the list object will also be 'mcmc'. 
 }
 \item{causal}{
 If multiple time series are given, one can specify causal connections from
 underlying processes belonging to such time series to the underlying
 processes belonging to some other time series. This specification should come
 in the form of an 4 x n matrix, where n is the number of causal connections.
 For each column (i.e. for each connection), row 1 specifies the number of the
 cause time series (same as the ordering of the input 'layer.series.structure' input to
 the routine), row 2 specifies the layer of the cause process, row 3 specifies
 the number of the effect time series and row 4 specifies the layer of the
 effect process. So if one has two time series, X and Y, both being modelled as
 having 3 layers, one can specify a causal connection from layer 3 in X to
 layer 1 in Y by the matrix "causal<-matrix(c(1,3,2,1),ncol=1)". 

 Note that it is possible to specify a causal connection even with only one
 time series, if a multi-layered process is considered. So if X is two-layered,
 "causal<-matrix(c(1,1,1,2),ncol=1)" will work. As layers are defined so that there is
 always a causal connection from layer 'k' to layer 'k-1' (with strength=1 for
 identifiability purposes), the reverse specification 
 ("causal<-matrix(c(1,2,1,1),ncol=1)") will not make sense.
 }
 \item{causal.symmetric}{
 Specifies a set of causal connections going both ways with equal strength.
 For more on the usage, see option 'causal', but note that symmetric causal
 connection between one layer and the next does not make sense, as the
 connection strength from layer 'k' to layer 'k-1' is always set to "1" for
 identifiability purposes. Note that which series is first and which is second 
 does not matter, so the matrix "causal.symmetric<-matrix(c(2,1,1,3),ncol=1)" would 
 specify the same model as "causal.symmetric<-matrix(c(1,3,2,1),ncol=1)".
 }
 \item{corr}{
 A 4 x n matrix which specifies the correlative connections between the processes,
 where n is the number of such connections. For each column (i.e. for each 
 connection), row 1 specifies the number of the first time series (same as 
 the ordering of the input 'layer.series.structure' input to
 the routine), row 2 specifies the layer of the first process, row 3 specifies
 the number of the second time series and row 4 specifies the layer of the
 second process.  So if one has two time series, X and Y, both being modelled as
 having 3 layers, one can specify a correlative connection between layer 3 in X and
 layer 1 in Y by the matrix "corr<-matrix(c(1,3,2,1),ncol=1)". Note that 
 which series is first and which is second does not matter, so the
 matrix "corr<-matrix(c(2,1,1,3),ncol=1)" would specify the exact same model.
 }
 \item{smoothing.specs}{
 A list that specifies whether smoothing samples are to be returned and if so how.
 Smoothing samples are a posteriori samples of the process states at measurement
 times and possibly also other times. Smoothing samples can be
 turned on by setting smoothing.specs$do.smoothing=TRUE. If so, the 
 smoothing.specs$smoothing.time.diff and smoothing.specs$num.smooth.per.mcmc
 must also be given. If smoothing.specs$smoothing.time.diff is set to "0",
 only the states at the measurement times will be inferred. If 
 smoothing.specs$smoothing.time.diff is set to a positive real value,
 the state at the regularly spaced time points 
 (smoothing.specs$smoothing.start,
  smoothing.specs$smoothing.start+smoothing.specs$smoothing.time.diff, 
  smoothing.specs$smoothing.start+2*smoothing.specs$smoothing.time.diff, 
  ...
  smoothing.specs$smoothing.end) 
 will also be reported. Note that this means that smoothing.specs$smoothing.start
 and smoothing.specs$smoothing.end must also be given. Also, if smoothing
 samples are to be given, smoothing.specs$num.smooth.per.mcmc must be given.
 This is the number of smoothing samples fetched per MCMC sample. 
 Lastly, smoothing.specs$do.return.smoothing.sample determines if the
 smoothing samples are to be returned or not. If smoothing.specs$do.smoothing=TRUE,
 it rarely makes sense to set smoothing.specs$do.return.smoothing.sample<-FALSE, so
 if this is not given, the algorithm sets smoothing.specs$do.return.smoothing.sample<-
 smoothing.specs$do.smoothing.
 }
 \item{realization.specs}{
 A list that specifies whether realization samples are to be returned and if so how.
 Realization samples are a posteriori samples of the process itself (including
 the correlations between time points) at measurement
 times and possibly also other times. Realization samples can be
 turned on by setting realization.specs$do.realizations=TRUE. If so, the 
 realization.specs$realization.time.diff and realization.specs$num.realizations
 must also be given. If realization.specs$realization.time.diff is set to "0",
 only the states at the measurement times will be inferred. If 
 realization.specs$realization.time.diff is set to a positive real value,
 the state at the regularly spaced time points 
 (realization.specs$realization.start,
  realization.specs$realization.start+realization.specs$realization.time.diff, 
  realization.specs$realization.start+2*realization.specs$realization.time.diff, 
  ...
  realization.specs$realization.end) 
 will also be reported. Note that this means that realization.specs$realization.start
 and realization.specs$realization.end must also be given. Also, if realization
 samples are to be given, realization.specs$num.realizations must be given.
 This is the number of realization samples fetched in total. (This cannot be more
 than the number of MCMC samples. the reason it is not equal to that, is that
 when censoring is in use, not all MCMC samples will result in valid
 realizations.)
 
 Lastly, realization.specs$strategy determines whether censoring of
 realizations is to take place. Censoring has been used when the process in
 question represented sedimentation, in order to avoid realizations such
 that the process would have erased one or more of the measurements. 
 Thus if X is sedimentation depth and X_measurement < X_later_time_point
 for any combination of measurement or later time point, a censoring
 would take place. realization.specs$strategy="N" means no censoring.
 realization.specs$strategy="A" means ascending censoring (censors 
 realizations that anywhere ascends beyond the value of the previous 
 measurement point and realization.specs$strategy="D" means descending 
 censoring (censors realizations that anywhere descends beyond the
 value of the previous measurement point).
 }
 \item{return.residuals}{
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
 \item{loglik.laxness}{
 There are sanity checks for each measurement point handled during the
 log-likelihood calculation, in order to avoid errors in the calculation
 due to nuemrical errors. However, some nuemrical inaccuracies will
 always exist, so some laxness in these checks must be allowed.
 Too much laxness and erroneous log-likelihoods may be used.
 Too much rigidity, i.e. too little laxness and parameter values
 may needlessly be excluded and the optimal estiamtion is not
 achieved. Default laxness is "high", but options are also "low"
 (very rigid), "moderate" and "complete" (sanity checks are not performed).
 }
}
\details{

The function takes a set of time series along with their process
structure (represented by layer.series.structure objects),
together with specifications for causal and correlative links,
thus specifying data and the model.
Run conditions for the numerics can also be specified. Performs
parameter inference, either Bayesian using MCMC or classic (ML) using
optimization from MCMC samples. Bayesian model likelihoods are
calculated using an importance sampling method, while ML optimization gives
what is needed for information criterion based model comparison.
Process inference can also be performed, if specified. The numerical
work is performed by a compiled C++ routine, in order to minimize
computation time.

Each data set+process structure can also contain a prior hyperparameter
specification for the Bayesian analysis, see 'layer.prior'.

Analyses of different models can be compared using the 'compare.layered'
routine.

}
\value{
Returns an object of type 'layered', a list containing parameter estimates
and uncertainties, plus anything of extra output specified by the user.
Extra output includes log(BML) for bayesian analysis, AIC, AICc, BIC for classic
analysis, MCMC samples, process smoothing samples and process realization samples.

An object of the class 'layered' will typically have the following components:

mu_<series name>: A list containing "mean", "median", "lower02_5" (2.5%% quantile) and "upper97_5" (97.5%% quantile) for the expected value parameter of each series. (PS: Not given for Wiener processes or layered processes having the Wiener process as the bottom layer).

dt_<series_name>: A list containing "mean", "median", "lower02_5" (2.5%% quantile) and "upper97_5" (97.5%% quantile) for the expected characteristic of each series. (Alternatively, the half-life will be given if use.half.life=TRUE).

sigma_<series_name>: A list containing "mean", "median", "lower02_5" (2.5%% quantile) and "upper97_5" (97.5%% quantile) for the stochastic contribution parameter of each series. (Alternatively the stationary standard deviation will be given, if use.stationary.stdev=TRUE).

model.log.lik: Bayesian marginal model likelihood (unless do.model.likelihood=FALSE or ML analysis).

obs_sd_<series_name>: A list containing "mean", "median", "lower02_5" (2.5%% quantile) and "upper97_5" (97.5%% quantile) for the observational standard deviation (if not individual standard deviations are given in the data). 

computing.time: Computing time in seconds.

parameter.names: Array of parameter names.

data.structure: The data strucucture(s) sent to the analysis (including data, prior, number of layes etc.)

causal: Causal connections specified.

corr: Correlative connections specified.

causal.symmetric: Symmetric causal connections specified.

mcmc: Only returned if "mcmc=TRUE" in the input. Contains MCMC samples for each parameter.

mcmc.origpar: Only returned if "mcmc=TRUE" in the input. Contains MCMC samples for each parameter using the original parametrization, i.e. the parametrization used internally. (So log characteristic time instead of characteristic time, for instance.)
 
est.origpar: Estimates for original parametrization (see 'mcmc.origpar').

input.options: A full list of the options used in this analysis.

process.time.points, process.mean, process.lower95, process.upper95: Process time points, means, upper and lower credibility bands for each process. (PS: There can be multiple prcoesses per serie, since number of layers can be greater than one).  Only returned if option "smoothing.specs=list(do.smoothing=TRUE,...)" is used.

smoothing.samples: Smoothing samples for each MCMC sample for each process.  Only returned if option "smoothing.specs=list(do.smoothing=TRUE,...,do.return.smoothing.samples = TRUE)" is used.

standardized.residuals: Standardized residuals are measurements minus expected values from the Kalman filter, divided by the standard deviation also from the Kalman filter. Only returned if "return.residuals=TRUE" is set.

prior.expected.values: Expected values from the Kalman filter. Only returned if "return.residuals=TRUE" is set.

residuals.time: Time points for each residual. Only returned if "return.residuals=TRUE" is set.

realizations: Process realizations given measurements. Only returned if "realization.specs=list(do.realizations=TRUE,...)" is set.

realization.time.points: Time point for each realization. Only returned if "realization.specs=list(do.realizations=TRUE,...)" is set.

description: A description of the model can be found in the field 'description',containing structure descriptions and causal and correlative connections.

}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.

Geyer CJ (1991),
Comp Sci Stat Proc 23rd Symp Inter: 156.
Am Stat Ass, New York
}
\author{
Trond Reitan, trond.reitan@ibv.uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also layer.series.structure, layer.data.series, layer.prior, 
summary.layered and compare.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Create a sample from the Ornstein-Uhlenbeck (OU) process with mu=0, 
# characteristic time=50 and sigma=0.2:
n=1000
x=rep(rnorm(1),n)
for(i in 2:n)
 x[i]=exp(-1/50)*x[i-1]+0.2*rnorm(1)

# Create measurement set as a subsample of the process, plus 
# noise with standard deviation 0.1
t=sort(sample(1:n,300))
y=(x+rnorm(length(x),0,0.1))[t]

# Create the data series object from the time points and "measurements":
X=layer.data.series(time.points=t, value.points=y,name="X")

#res=layer.analyzer(layer.series.structure(X,numlayers=1))
#summary(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition }


