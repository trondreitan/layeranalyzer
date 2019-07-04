\name{print.summary.layered}
\alias{print.summary.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarizes the return object from 'layer.analyzer'. }
\description{
The function prints a summary, in the form of a 'summary.layered' object,
of a model analysis, as performed by 'layer.analyzer'. 
This comes in the form of a table with
parameter estimates and uncertainties, plus indicators of model fit.
}
\usage{
\method{print.summary}{layered}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 An 'summary.layered' object, as returned from the 'summary.layered' method. 
 }
  \item{\dots}{
 Not used. Only included for compatibility with other other 'print.summary' methods.
 }
}
\details{
Prints a 'summary.layered' object, which in turn performs a summary of a
'layered' object, as returned from 'layer.analyzer'. For 
each parameter, mean- and median-based estimates
and 95\% credibility intervals are given, when Bayesian analysis is performed.
In addition, the Bayesian model (marginal) likelihood is given. 
For classic maximum-likelihood (ML) analysis, the parameter estimates (ML) and various
model selection criteria are given. Bayesian credibility intervals are still given in
order to inform about parameter uncertainty.
}
\value{
No return value. The method only prints to screen. 
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
See also layer.analyzer and summary.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_1layer.txt","X")
# Calculate Bayesian model likelihood for an OU model:
#res1=layer.analyzer(layer.series.structure(X,numlayers=1,prior=pr))
#summary(res1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, AIC, BIC }
