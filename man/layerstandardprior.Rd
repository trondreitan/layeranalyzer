\name{layer.standard.prior}
\alias{layer.standard.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specifying a default prior distribution for the layer.analyzer. }
\description{
This 'layer.prior' object defines a Bayesian prior distribution used as the default
prior when the user does not give a prior specification. 
}
\usage{
layer.standard.prior
}
%- maybe also 'usage' for other objects documented here.
\details{
Specifies the standard Bayesian prior distribution used for default purposes.
Warning: Will not be suitable in every context.
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
See also layer.prior and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
layer.standard.prior
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
