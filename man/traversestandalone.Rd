\name{traverse.standalone.layered}
\alias{traverse.standalone.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model traversal tool for standalone time series. }
\description{
Traverses all models for a single standalone time series, up to a 
given model complexity, as given by the maximum number of layers.
}
\usage{
traverse.standalone.layered(timeseries,
  max.layers=3, 
  talkative=FALSE, allow.one.feedback.loop=FALSE, 
  just.stationary=FALSE, no.rw=FALSE, regional.options=FALSE,    
  time.integrals.possible=FALSE, allow.deterministic.layers=TRUE,   
  num.MCMC=1000,spacing=10,burnin=10000,num.temp=1,
  do.maximum.likelihood=FALSE,maximum.likelihood.numstart=10,
  id.strategy=2,use.stationary.stdev=FALSE,T.ground=1.5, 
  use.half.lives=FALSE, mcmc=FALSE,
  prior=layer.standard.prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeseries}{
 A single time series in the form of a 'layer.data.series' object. This
 is the input time series for the analysis. 
 }
 \item{max.layers}{
 The maximal number of layers (hidden and observed) that is allowed in
 the model traversal. If any of the switches for how each layer is 
 handled (such as 'regional.options=TRUE', 'time.integrals.possible=TRUE',
 'allow.one.feedback.loop=TRUE' and 'allow.deterministic.layers=TRUE') are
 in use, the number of models to traverse will quickly increase with
 this maximal number of layers.
 }
 \item{talkative}{
 Gives a quick description of each model traversed and shows the progress.
 }
 \item{allow.one.feedback.loop}{
 If set to 'TRUE', allows for one feedback loop from a higher to a lower 
 layer (so 'against' the normal flow of causality in these models). This
 can provide an explanation for quasi-cyclical  behaviour.
 }
 \item{just.stationary}{
 If set to 'TRUE' allows for other models than OU (Ornstein-Uhlenbeck) at 
 the lowest layer, namely also the Wiener process (also called Random Walk) 
 and linear trend plus OU.
 }
 \item{no.rw}{
 Only in use when 'just.stationary=FALSE'. Switches off looking at the
 Wiener process (Random Walk) in the lowest layer.
 }
 \item{regional.options}{
 If set to 'TRUE' switches on traversal of the different regional options.
 This only makes sense if the number of sites is larger than one. Also,
 it can massively expland the number of models traversed. Not implemented yet!
 }
 \item{time.integrals.possible}{
 If set to 'TRUE' allows for time integral layers, i.e. the existence of 
 one or more layers which are the time integral of the layer below. This
 of course only makes sense if the number of layers is alrger than one.
 }
 \item{allow.deterministic.layers}{
 If set to 'TRUE', allows for traversing models where there is no
 stochasitcity (i.e. determinism) on some layers. This can easily be
 the case as long as the number of layers is larger than one. (The 
 lowest layer cannot be deterministic, as that just means a flat line.)
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
 \item{prior}{
 An object specifying the prior distribution, in the form of a 'layer.prior' object.
 }
}
\details{
 This routine traverses all models for a time series up to a certain model 
 complexity, as given by the maximal number of layers. All the traversed
 model are returned and this array of models can be sent directly to the
 'compare.layered' routine. Certain options can be explored, such as 
 stationarity, determinstic layers, time integral layers and the existance 
 of a feedback loop. Using these options means more models are traversed, though.
}
\value{
Returns a list object where each element will be of the type 'layered', 
representing a set of model analyses. If the turned object is called 'models', 
then each model analysis can be looked at through with 'models[[i]]' where
'i' is the model number. Using the 'talkative' option means a short description 
of each model will be shown during traversal. Or one can use 'summary(models[[i]])'
to see both the structure alongside the parameter estimates.
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
See also layer.data.series, layer.prior, 
summary.layered and compare.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

p=layer.prior(mu=c(-5,5),init=c(-100,100),dt=c(0.5,200),sigma=c(0.01,2),obs=c(0.01,1))
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_1layer.txt","X")

#models=traverse.standalone.layered(X, max.layers=2, 
#  talkative=TRUE, allow.one.feedback.loop=FALSE, 
#  just.stationary=FALSE, no.rw=FALSE,    
#  time.integrals.possible=FALSE, 
#  allow.deterministic.layers=TRUE,
#  num.MCMC=200,spacing=10,burnin=2000, prior=p)
#
#compare.layered(models)
#cat(models[[3]]$description)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition }


