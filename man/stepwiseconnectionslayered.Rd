\name{stepwise.connections.layered}
\alias{stepwise.connections.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Time series analysis tool using linear layered SDEs. }
\description{
Traverse locally around a given connection model in a step-wise fashion,
starting locally around the null model (no connections). Can perform
larger model selections, but may not be robust enough still. Parallel
runs using the standalone program 'layeranalyzer' should be preferred
for larger model selection tasks.
}
\usage{
stepwise.connections.layered(..., num.MCMC = 1000, spacing = 10, burnin = 10000,
                 num.temp = 1, do.maximum.likelihood = FALSE,
                 maximum.likelihood.numstart = 10, silent.mode = TRUE,
                 talkative.burnin = FALSE, talkative.likelihood =
                 FALSE, talkative.traversal = TRUE, test.mode = FALSE,
                 id.strategy = 2, use.stationary.stdev = FALSE,
                 T.ground = 1.5, use.half.lives = FALSE, mcmc = FALSE,
                 allow.causal = TRUE, allow.correlation = TRUE,
                 allow.direct.feedback = TRUE, first.is.nullhypothesis
                 = FALSE, ML.IC = "AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
 A set of 'layer.series.structure' objects that represents the set of 
 time series that are to be analyzed and how the underlying processes of
 each such series is to be modelles. (See 'layer.series.structure'). At
 least one such object must be given. 'traverse.standalone.layered' can be
 used for settling the structure of each time series. 
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
 \item{talkative.traversal}{
 If set to 'TRUE'(default), the connection pairs and the causal and correlative
 connections of each model.
 }
 \item{test.mode}{
 If set to 'TRUE', the analyses themselves are not performed. Only the 
 connection pairs and the causal and correlative connections of each model
 are shown. This can be useful if one wants to know how many models are
 are to be analyzed. 
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
 \item{allow.causal}{
 If set to 'FALSE', no causal connection models will be analyzed, only 
 correlative ones.}
 \item{allow.correlation}{
 If set to 'FALSE', no correlative conneciton models will be analyzed, only
 causal ones. (Incompatible with 'allow.causal=FALSE'). }
 \item{allow.direct.feedback}{
 If set to 'FALSE', no direct causal loops will be examined, i.e. causal
 links from series 1, layer 1 to series 2, layer 2 and another liunk back again.
 (Does not rule out more elaborate causal loops.)}
 \item{first.is.nullhypothesis}{
  If prior probabilities (item 'p0') is not given, the prior probabilities
  are default equal for all models. However, if first.is.nullhypothesis
  is set to 'TRUE', the first model is given 50\% of the prior probability, while
  the rest of the models shares the remaining 50\% equally.
  This variable is sent to 'compared.layered' for model comparison. 
 }
 \item{ML.IC}{
  If classic ML analysis is performed instead of Bayesian analysis, 
  "AIC","BIC" or "AICc" will be used for model comparison instead of 
  the Bayesian marginal likelihood. The item 'ML.IC' determines which of these 
  three model comparison criteria are to be used.
  This variable is sent to 'compared.layered' for model comparison. 
 }
}
\details{
Traverses all causal and correlative link combinations between a set of
time series. This can mean a lot of models. With even just 3 one layered
time series, that's 125 models, while 4 one layered time series means
15625 models! 
}
\value{
Returns an array of objects of type 'layered', each of which is a list 
containing parameter estimates and uncertainties, plus anything of extra 
output specified by the user.

This array can be inserted into 'compare.layered' or 'anova.layered' (if ML
analysis). 
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
See also layer.analyzer, layer.analyzer.timeseries.list, compare.layered
and anova.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Read data:
X1=read.layer.data.series("http://folk.uio.no/trondr/R/test_cause.txt",name="cause")
X2=read.layer.data.series("http://folk.uio.no/trondr/R/test_effect.txt",name="effect")

# Set structure:
X1struct=layer.series.structure(X1, numlayers=1)
X2struct=layer.series.structure(X2, numlayers=1)

# Perform traversal of analyses:
#res=stepwise.connections.layered(X1struct,X2struct, num.MCMC=100, burnin=1000)

# Compare analyses:
#compare.layered(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, model traversal, connections }
