\name{layer.load.prior}
\alias{layer.load.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reading prior specification for the layer.analyzer from a file. }
\description{
The function reads a prior specification from a csv file with semi-colon-separated instances. 
}
\usage{
layer.load.prior(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
 The name of the semi-colon-separated csv file that contains the
 prior specification. 
 }
 
}
\details{
This csv file should have a header row and one value row.
It should have the elements 
"is_log;mu1;mu2;dt1;dt2;s1;s2;lin1;lin2;beta1;beta2;init1;init2", where 
is_log is an indicator of whether the layer.analyzer should log-transform
the data or not (0=no, 1=yes, 2=already log-transformed). mu1;mu2 specifies
a 95% credibility interval for the expected value. Similarly, dt1;dt2, s1;s2,
lin1;lin2, beta1;beta2, init1;init2 specifies a 95% credbility interval
for the characteristic time, size of stochastic contributions, linear time trend,
causal connections and initial value, respectively. Only "is_log;mu1;mu2;dt1;dt2;s1;s2"
is mandatory for the simplest analysis. More options might make the rest necessary
also, see the documentation for 'layer.analyzer'.

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
See also layer.prior.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pr=layer.load.prior("http://folk.uio.no/trondr/R/test_prior.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, Bayesian, prior }


