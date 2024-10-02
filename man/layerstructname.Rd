\name{layer.struct.name}
\alias{layer.struct.name}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a description of the process structure of a time series. }
\description{
This function generates a string describing the structure specified for a time series.
}

\usage{
layer.struct.name(struct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{struct}{
  A 'layer.series.structure' object that should be described.
 }
}
\details{
Describes as a string the number of layers and the properties of
each layer, as specified in the input 'layer.series.structure' object.
For series with sites (regional data), the regional settings for
each layer is also described. Used when generating a 'layer.series.structure',
which stores the result in a field named 'description'. This is again
used when generating a description for an analysis in 'layer.analyzer'.
}
\value{
Returns a descriptive string.
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
See also layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","X")
struct=layer.series.structure(X,numlayers=2)
layer.struct.name(struct)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
