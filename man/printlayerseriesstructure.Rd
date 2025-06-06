\name{print.layer.series.structure}
\alias{print.layer.series.structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Outputs the description of a 'layer'series'structure' object. }
\description{
This function (which should override 'print') outputs the description of the process structure of a series.
}
\usage{
\method{print}{layer.series.structure}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A 'layer.series.structure' object for which the description should be printed.
 }
}
\details{
Shows the assigned description of a 'layer.series.structure' object
(typically generated by 'layer.struct.name').
}
\value{
Prints a descriptive string.
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
See also layer.series.structure and layer.struct.name.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","X")
struct=layer.series.structure(X,numlayers=2)
print.layer.series.structure(struct)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
