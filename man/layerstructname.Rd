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
library(layeranalyzer)
X<-read.layer.data.series(
 "https://github.com/trondreitan/layeranalyzer/raw/master/doc/test_2layer.txt",
 name="X")

# Insert the data intro a structure specificatiion
# with the default structure (OU):
struct1=layer.series.structure(X) 
layer.struct.name(struct1)
# Output should be "OU"

# Insert the data into a structure with 2 layers, where the
# lowest layer has no pull (i.e. Brownian motion, here also called
# "Random Walk" which is abbreviated to "RW"). Note that when the
# lowest layer is Brownian, the processes are nonstationary and
# need initial value treatment (init.0=TRUE).
# The top layer is specified to have no stochasticity
# (i.e. deterministic given the layer below it). 
struct2=layer.series.structure(X,numlayers=2, no.pull=TRUE, no.sigma=1,
   init.0=TRUE) 
layer.struct.name(struct2)
# Output should be "2-layered: Layer 2: RW , Layer 1: deterministic tracking"

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{linear stochastic differential equations}
\concept{layered}
\concept{hidden layers}
\concept{causal}
\concept{correlative}
\concept{matrix inverse}
\concept{Lapack}
\concept{eigenvalue decomposition}
\concept{Bayesian}
\concept{prior}
