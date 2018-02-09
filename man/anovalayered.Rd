\name{anova.layered}
\alias{anova.layered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Analysis of Deviance for return object from 'layer.analyzer'. }
\description{
Compute an analysis of deviance table for one or more return object from 'layer.analyzer'.
}
\usage{
anova.layered(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,...}{
 Objects of class 'layered', as returned from the 'layer.analyzer' method. 
 }
}
\details{
Performs analysis of deviance for a set of return objects from 'layer.analyzer'.
These objects should be classically (ML) estimated, rather than Bayesian analysis
(option 'do.maximum.likelihood=TRUE' in 'layer.analyzer'). These objects should
represent nested models, with the simplest one first and the most complicated
one last. 

As long as the difference in degrees of freedom is larger than 0, a p-value
from the chi-squared-test, (D_0-D_1)~chisq(res.df.0-res.df.1), will also
be shown.
}
\value{
An object of class ‘"anova"’ inheriting from class ‘"data.frame"’.
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
See also layer.analyzer.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt")
# Calculate Bayesian model likelihood for an OU model, a 2 layered model
# and a 3 layered model. Test the number of layers. The dataset itself
# is simulated using 2 layers:
#res1=layer.analyzer(layer.series.structure(X,numlayers=1,prior=pr),do.maximum.likelihood=TRUE)
#res2=layer.analyzer(layer.series.structure(X,numlayers=2,prior=pr),do.maximum.likelihood=TRUE)
#res3=layer.analyzer(layer.series.structure(X,numlayers=2,prior=pr),do.maximum.likelihood=TRUE)
#anova.layered(res1,res2,res3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ anova, deviance, linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, ML, classic, hypothesis testing }
