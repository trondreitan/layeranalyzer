\name{layer.standard.log.prior}
\alias{layer.standard.log.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specifying a suggestion for prior distribution, 
for data that will be log-transformed. }
\description{
This 'layer.prior' object defines a Bayesian prior distribution suggested for use 
when the user wants readily defined prior for log-transformed data.
}
\usage{
layer.standard.log.prior
}
%- maybe also 'usage' for other objects documented here.
\details{
Specifies a "standard" Bayesian prior distribution suitable for log-transformed data.
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
See also layer.prior, layer.standard.prior and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
layer.standard.log.prior
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }
