\name{layer.wide.prior}
\alias{layer.wide.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specifying a suggestion for a wide (little information) prior distribution. }
\description{
This 'layer.prior' object defines a Bayesian prior distribution suggested for use 
when the user wants a 'wide' prior, i.e. one specifying little prior knowledge.
}
\usage{
layer.wide.prior
}
%- maybe also 'usage' for other objects documented here.
\details{
Specifies a "standard" Bayesian prior distribution suitable for when very little
prior knowledge exists. Suitable when the user knows very little about the 
nature of the data and the underlying processes.
}
\value{
A 'layer.prior' object.
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
See also layer.prior, layer.standard.prior and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
layer.wide.prior
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
