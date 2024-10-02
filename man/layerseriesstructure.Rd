\name{layer.series.structure}
\alias{layer.series.structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Represents a measurement time series. }
\description{
This class function represents a measurement time series plus the way the underlying
process is modelled, not counting connections to other processes. Primarly,
that means the number of layers used for modelling that process and whether the
bottom layer is an Ornstein-Uhlenbeck (OU) process, a Wiener process or an 
OU with linear trend, but other options do also exists.
}
\usage{
layer.series.structure(timeseries, numlayers=1, lin.time=FALSE,
                       time.integral=NULL, no.pull=FALSE, no.sigma=NULL,
                       regional.mu=FALSE, regional.lin.time=FALSE,
                       regional.pull=NULL, regional.sigma=NULL,
                       correlated.sigma=NULL, pairwise.correlated.sigma=NULL,
                       one.dim.sigma=NULL, grouping.sigma=NULL,
		       remove.sigma=NULL, differentiate.sigma=NULL,
		       differentiate.pull=NULL, differentiate.mu=FALSE,
                       differentiate.lin.time=FALSE,
                       init.0=FALSE, init.time=NULL, init.same.sites=FALSE,
		       init.same.layers=FALSE, init.specified=NULL,
		       allow.pos.pull=FALSE, period=NULL, 
                       prior=layer.standard.prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeseries}{
  A 'layer.data.series' object representing a measurement time series.
 }
  \item{numlayers}{
 The number of layers in the model (hidden layers + measured layer). Setting
 numlayers=1 and nothing more means modellings the process behind the time series
 as an OU process. Note that increasing the number of layers can dramatically
 increase the computer resources needed (by as much as the cube of the number of
 layers), since matrix operations are needed in order to update knowledge of
 all layered processes.   
 }
  \item{lin.time}{
 If set to 'TRUE', this indicates that a linear time trend is to be included. 
 The bottom layer must be an OU process if so (i.e. 'no.pull' must be set to "False"). 
 }
  \item{time.integral}{
 An array of integers specifying which layers are the time integral of the layer below it.
 In differential equation terms, that means dX_n=X_(n+1) where 'n' is the layer number 
 in question. This option was introduced in order to model sedimention height by 
 letting it be the time integral of sedimentation rate and letting sedimentation rate
 be an OU process. Usually, only one time integral layer is needed, but the analysis
 allows for multiple such layers. Note that the bottom layer cannot be a time
 integral layer, as it is not tracking anything, so 'time.integral' must be less that
 'numlayers'. 
 }
  \item{no.pull}{
 If set to 'TRUE', it specifies that the bottom layers has no pull term, 
 i.e. that it is a Wiener process (random walk) rather than an OU process. 
 When such models are included, it pays to use the initial value options, see the
 'init' options in this routine. If not, a very wide generic prior for the initial 
 value is used, giving the Wiener process model an unfair disadvantage.
 }
  \item{no.sigma}{
 An array of integers indicating which layers are modelled to have no stochastic
 contributions. That means each layer specified will be a deterministic tracking
 of the layer below and will in practise be a smoothed version of the underlying 
 process. because of this, it does not make sense setting 'no.sigma' to
 'numlayers', as the bottom layers needs to have stoachasticity in order to 
 have any dynamics at all.
 }
  \item{prior}{
 Should be a 'layer.prior' object that specifies the prior distribution of the 
 parameters, see the 'layer.prior' documentation. If not given, a standard prior
 is used. (Warning: This will not be appropriate for all circumstances). 
 }
  \item{period}{
 If given, should specify an array of period lengths, in order to include in the 
 model periodic sine/cosine fluctuations (such as seasonality). Regression parameters
 for the sine/cosine terms are estimated.
 }
  \item{init.0}{
 If set to 'TRUE', the model includes a parameter for the initial values of
 each process (one for each combination of site and layer). The initial time
 is that of the first measurement.
 }
  \item{init.time}{
 If set, the model includes a parameter for the initial values of
 each process (one for each combination of site and layer). The initial time
 is that of 'init.time'.
 }
  \item{init.same.sites}{
 Only makes sense for time series with multiple sites. If set, an initial value 
 treatment will assume the same initial value for all sites (if applicable).
 }
  \item{init.same.layers}{
 If set, an initial value treatment will assume the same initial value for all 
 layers (if applicable, i.e. if 'numlayers'>1).
 }
  \item{init.specified}{
 If given, should be an array of two elements, initial time and initial value.
 }
  \item{regional.mu}{
 Only makes sense for time series with multiple sites. If set to 'TRUE', specifies 
 that the process of each site should have a separate expected value (mu). 
 }
  \item{regional.lin.time}{
 Only makes sense for time series with multiple sites. If set to 'TRUE', specifies 
 that the process of each site should have a separate linear time trend. 
 }
  \item{regional.pull}{
 Only makes sense for time series with multiple sites. Should be an integer array
 specifying which layers have separate pulls (and thus characteristic times/halflives
 and stationary standard deviations) for each site. 
 }
  \item{regional.sigma}{
 Only makes sense for time series with multiple sites. Should be an integer array
 specifying which layers have separate stochastic noise size (and thus stationary 
 standard deviations) for each site. 
 }
  \item{correlated.sigma}{
 Only makes sense for time series with multiple sites. Should be an array
 of integers specifying which layers have correlated stochastic contributions
 between the sites. The same correlation is assumed between any two sites.
 }
  \item{pairwise.correlated.sigma}{
 Only makes sense for time series with multiple sites. Should be an array
 of integers specifying which layers have correlated stochastic contributions
 between the sites. Different correlations are allowed between different pairs
 of sites.
 }
  \item{one.dim.sigma}{
 Only makes sense for time series with multiple sites. Should be an array
 of integers specifying which layers have perfect correlated stochastic contributions
 (meaning correlation=1, i.e. the same stochastic contribution for all sites) 
 between the sites. Can introduce numerical difficulties.
 }
  \item{grouping.sigma}{
 Only makes sense for time series with multiple sites. Should be an array
 of integers specifying which layers have correlated stochastic contributions
 between the sites. However, the sites will be put into two groups, with no
 correlation between sites beloning to one group and the other. The grouping 
 will be inferred. Not well tested.
 }
  \item{remove.sigma}{
  Only makes sense for time series with multiple sites. Should be an array
 of integers specifying which layers have a grouping between sites with and without
 correlation to the stochastic contributions of other sites. Not well tested.
 }
  \item{differentiate.sigma}{
 Only makes sense for time series with multiple sites. Attempts to find two
 groups of sites for the stochastic contribution sizes for the given set of layers,
 as specified by an array of integers.
 }
  \item{differentiate.pull}{
 Only makes sense for time series with multiple sites. Attempts to find two
 groups of sites for the pull (characteristic times/halftimes) for the given 
 set of layers, as specified by an array of integers.
 }
  \item{differentiate.mu}{
 Only makes sense for time series with multiple sites. If set to 'True',
 attempts to find two groups of sites for the expected value (mu) of the processes.
 }
  \item{differentiate.lin.time}{
 Only makes sense for time series with multiple sites. If set to 'True',
 attempts to find two groups of sites for the linear time trend.
 }
  \item{allow.pos.pull}{
 If set to 'TRUE', allows for positive pull parameters in the linear stochastic
 equations. Positive pulls means the process is prone to "explode", but for
 a short time (typically not much more than the characteristic time), such a process
 could none the less make sense.
 }
}
\details{
This class function is used as input to the analysis and represents how a 
single time series (which can include several sites) is to be modelled. This
includes the number of layers, the nature of the bottom layer and trends (linear
or period in time). Stochasticity can be switched off for given layers. if multiple
sites are included, there are options for controlling whether the parameters of 
the different sites are different or the same. Also, correlation between the processes
of different sites can be correlated with certain restrictions, or not correlated.
The objects of this class is used as input in the 'layer.analyzer' routine.
}
\value{
Returns the modelling specifications for a time series data set as a 
'layer.series.structure' object. Note the field 'description' which
contains a description of the object. Best outputted with
'cat(obj$description)
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
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","X")
# Calculate Bayesian model likelihood for an OU model:
#res1=layer.analyzer(layer.series.structure(X,numlayers=1,prior=pr))

# Calculate Bayesian model likelihood for a two-layered model:
#res2=layer.analyzer(layer.series.structure(X,numlayers=2,prior=pr))

# Calculate Bayesian model likelihood for a three-layered model:
#res3=layer.analyzer(layer.series.structure(X,numlayers=3,prior=pr))

# The models can then be compared, see 'compare.layered'.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, time series }
