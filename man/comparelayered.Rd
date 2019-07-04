\name{compare.layered}
\alias{compare.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform comparison of several models for the same dataset(s). }
\description{
The function compares a set of models for the same dataset(s). The comparison
is default done with Bayesian model probabilities. Prior model probababilities
can be given, but if not, each model is given the same prior probability. 
If the model analysis is classic, AIC is default used for model comparison, 
but BIC or AICc can also be specified.
}
\usage{
compare.layered(..., p0=NULL, first.is.nullhypothesis=FALSE, ML.IC="AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
 A set of 'layered' objects returned from the layer.analyzer function.
 This set represents the model analyses performed and should contain
 information used for model comparison (marginal likelihoods for Bayesian
 analyses, AIC/BIC/AICc for classic analyses).
 }
  \item{p0}{
 A numeric array describing the prior probabilities of each models
 sent to the function. The size of the array needs to match the number
 of models sent to the function. If prior probability array, 'p0', is not given
 each model is default given the same prior probabiliuty, but see item
 'first.is.nullhypothesis'.
 }
  \item{first.is.nullhypothesis}{
 If prior probabilities (item 'p0') is not given, the prior probabilities
 are default equal for all models. However, if first.is.nullhypothesis
 is set to 'TRUE', the first model is given 50\% of the prior probability, while
 the rest of the models shares the remaining 50\% equally.
 }
  \item{ML.IC}{
 If classic ML analysis is performed instead of Bayesian analysis, 
 "AIC","BIC" or "AICc" will be used for model comparison instead of 
 the Bayesian marginal likelihood. The item 'ML.IC' determines which of these 
 three model comparison criteria are to be used.
 }
}
\details{
The first arguments to this function should be a set of models represented
by 'layered' objects returned from the layer.analyzer function.

When Bayesian analyses have been performed (default), Bayesian
marginal likelihoods are used. (This means the option 'do.model.likelihood'
in the 'layer.analyzer' function must not be set to false in any of the analyses).
The model prior probabilities can be given, but if not, equal prior probability
will be assumed.

The item 'ML.IC' is used when classic analyses have been performed, in order
to determined which model comparison criteria will be used for that purpose.

The same dataset or set of datasets should be analysed for all models. 
This is something that can go wrong when comparing a hypothesis of no 
causal connection from process B to A versus a hypothesis of there being 
such a causal connection. Both time series A and B needs to be given, 
not only for the analysis of the connection hypothesis but also for the 
no connection hypothesis.

}
\value{
Returns a matrix containing the names of the model objects and the
posterior probabilties (or weights if classic analyses has been used).
}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.

Liow, L.H., Reitan, T., Harnik, P.G. (2015)
Ecological interactions on macroevolutionary time scales;
clams and brachiopods are more than ships that pass in the night,
Ecological Letters, Volume 18(10): 1030-1039. 
}
\author{
Trond Reitan, trond.reitan@ibv.uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also layer.analyzer.
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

#compare.layered(res1,res2,res3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, model comparison, AIC, BIC }
