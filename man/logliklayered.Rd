\name{logLik.layered}
\alias{logLik.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return maximum likelihood on the logarithmic scale for return object from 'layer.analyzer'.}
\description{
Retrieve maximum log-likelihood for one or more return object
from 'layer.analyzer'. (Requires classic analysis.)
}
\usage{
\method{logLik}{layered}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, \dots}{
 Objects of class 'layered', as returned from the 'layer.analyzer' function. 
 }
}
\details{
Used for extracting the maximum log-likelihood in order to do a
model comparison, using "anova.layered" (likelihood-ratio-testing).
}
\value{
An object of class "logLik".
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
See also layer.analyzer and anova.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","X")
# Calculate classic maximum likelihoods for an OU model, a 2 layered model
# and a 3 layered model. Look at log-likelihoods:
#res1=layer.analyzer(layer.series.structure(X,numlayers=1,prior=pr),do.maximum.likelihood=TRUE)
#res2=layer.analyzer(layer.series.structure(X,numlayers=2,prior=pr),do.maximum.likelihood=TRUE)
#res3=layer.analyzer(layer.series.structure(X,numlayers=2,prior=pr),do.maximum.likelihood=TRUE)
# c(logLik(res1), logLik(res2), logLik(res3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood, anova, deviance, linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, ML, classic, hypothesis testing }
