\name{bi.lext}
\alias{bi.lext}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A dataset of reed warbler body mass over several years. }
\description{
This example dataset was collected on bivalve/brachipod diversification
rates in Liow et al. (2015) and re-analyzed by Reitan & Liow (2017). The
specific dataset contains log-transformed bivalve extinction rates.
}
\usage{
}
%- maybe also 'usage' for other objects documented here.
\details{
The data is a list of class layer.data.series having five elements:
time (millions of years (MY) since fozzilisation, negative value so that
events are ordered chronologically).
value (log-transformed extinction rates)
name ("bivalve_lext")
std.dev (standard error estimates of log-rates).
is.datetime (0=false)
}
\value{
}
\references{
Liow, L. H., Reitan, T. and P. G. Harnik (2015). Ecological interactions on
macroevolutionary time scales: clams and brachiopods are more than ships
passing in the night. Ecological Letters 18(10): 1030-1039

Reitan, T. and L. H. Liow (2017). An unknown Phanerozoic driver of
brachiopod extinction rates unveiled by multivariate linear stochastic
differential equations. Paleobiology 43(4): 1-13.
}
\source{https://datadryad.org/resource/doi:10.5061/dryad.hj30r}
\author{
Trond Reitan, trond.reitan@ibv.uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also other datasets: hare, lynx, br.lorig, br.lext.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ macroevolution, bivalve, brachiopod, diversification, extinction rate, linear stochastic differential equations }
