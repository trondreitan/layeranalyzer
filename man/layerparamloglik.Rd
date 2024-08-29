\name{layer.param.loglik}
\alias{layer.param.loglik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returning log-likelihood for a set of parameter values. }
\description{
  Calculates the log-likelihood for a set of parameter values (or
  estimates if not given) for a given model specified by a previous
  analysis result. If silent.mode is set off, extra debug info is given.
  Meant primarily for debug and for exploration of the likelihood
  surface, post-analysis.
}
\usage{
layer.param.loglik(analysis, new.param.values=NULL, silent.mode=TRUE, num.optim=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{analysis}{
 The 'analysis' variable must be a 'layered' object. I.e. a 
 'layer.analyzer' analysis needs to be performed first. This variable
 is not an option. If it is NULL, the function will complain. it is
 from here the function gets the model specification (including
 conneciton elements in the case of multiple time series). Thus
 a small analysis should have been performed alread. Alternatively,
 one can try to fill in a 'layered' object by hand, but that takes
 effort.
 }
 \item{new.param.values}{
  Specifies the parameter values for the log-likelihood calculation.
  If this is not specified (NULL), then the parameter estimates from
  the analysis is used instead. For the order of the parameters,
  see for instance analysis$parameter.names. Note that the scale
  is the original scale, not the reparametrized scale used internally
  (where for instance stochastic contributions, characteristic times and
  observational noise level is log-trnasformed, and correlation is
  logit-transformed.) It is possible to set some of the parameter
  values to NA, in which case these will be optimized over.
 }
 \item{silent.mode}{
  If set to FALSE, shows debug info for each iteration in the Kalman filter.
 }
 \item{num.optim}{
  If some parameter values are set to NA, optimization is performed
  on those parameters. 'num.optim' specifies the number of optimizations
  performed on each such parameter (since optimization is performed
  by hill-climbing and can thus converge to different optima from
  different starting positions).
 }
}
\details{
 This function calculated log-likelihood for a specified (or under-specified)
 parameter set. If under-specified (having NA values) these are optimized
 over. Meant mostly as a debug tool (if silent.mode=FALSE) but can also be
 used for exploring the likelihood surface. Note that the specified
 parameter set is on the original scale, not the reparametrized version
 used internally.
}
\value{
Returns an numeric value for the log-likelihood.
}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.
}
\author{
Trond Reitan, trond.reitan@geo.uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also layer.analyzer.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Define time series (time points, values and in this case also
# sample standard deviations and sample size, used for calculating
# standard errors for each measurement).
# install.packages("../www_docs/R/layeranalyzer_0.2.1.tar.gz",type="source")
library(layeranalyzer)
X=layer.data.series(time.points=malta$Time.Year,
  value.points=malta$Mean..log.body.mass.,
  std.dev=sqrt(malta$Variance.calculated.from.the.data), 
  num.meas.per.value=malta$Sample.size,name="log.body.size")

# Define process structure for analysis purposes:
ser=layer.series.structure(X,numlayers=1)

# Perform Bayesian analysis:
res=layer.analyzer(ser)

layer.param.loglik(res, silent.mode=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, likelihood, debug }


