\name{layer.prior}
\alias{layer.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specifying a prior distribution for the layer.analyzer. }
\description{
This class function represents a Bayesian prior distribution specified by a couple of hyper-parameters. 
}
\usage{
layer.prior(mu,dt,
  sigma, init,
  lin, beta,
  obs, islog=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
  An array with 2 elements that specifies the low and high end of
  the 95\% credibility band of mu, the expected value for the process.
 }
  \item{dt}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of dt, the characteristic time
 (or the half-life if that is specified).
 }
  \item{sigma}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of sigma, the stochastic contribution.
 }
  \item{init}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of init, the initical value of the process.
 }
  \item{lin}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of lin, the linear trend (if applicable).
 }
  \item{beta}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of beta, the strength of causal links
 (if applicable).
 }
  \item{obs}{
 An array with 2 elements that specifies the low and high end of
 the 95\% credibility band of obs, the standard deviation of the
 observational errors.
 }
  \item{islog}{
 If set to 1, the analysis will log-transform the data before analysis.
 If set to 2, the analysis will assume the data has already been
 log-transformed.
 }
}
\details{
Specifies the Bayesian prior distribution by specifying the boundries of
the 95\% credibiluity band of each class of parameters. Only 'mu', 'dt'
and 'sigma' are mandatory. The objects 'layer.standard.prior',
'layer.standard.log.prior', 'layer.wide.prior' are pre-specified
layer.prior objects.
}
\value{
Returns prior as a 'layer.prior' object.
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
See also layer.load.prior and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.prior(mu=c(-10,10),dt=c(0.001,1000),sigma=c(0.01,10),
	                         lin=c(-1,1),beta=c(-1,1),init=c(-10,10),
                                 obs=c(0.01,1),islog=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
