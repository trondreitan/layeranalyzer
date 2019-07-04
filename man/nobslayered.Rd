\name{nobs.layered}
\alias{nobs.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return number of observations for return object from 'layer.analyzer'.}
\description{
Retrieve number of observations for one or more return object
from 'layer.analyzer'. 
}
\usage{
\method{nobs}{layered}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, \dots}{
 Objects of class 'layered', as returned from the 'layer.analyzer' function. 
 }
}
\details{
Used for extracting the number of observations in order to do a
model comparison, using "anova.layered" (likelihood-ratio-testing).
}
\value{
An integer representing the number of observations.
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
See also layer.analyzer, logLik.layered and anova.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","X")
#res1=layer.analyzer(layer.series.structure(X,numlayers=1,prior=pr),do.maximum.likelihood=TRUE)
# nobs(res1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ number of observations, anova, deviance, linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, ML, classic, hypothesis testing }
