\name{layer.data.series}
\alias{layer.data.series}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Represents a measurement time series. }
\description{
This class function represents a measurement time series. Can include multiple sites. Individual measurement error can be added. 
}
\usage{
layer.data.series(time.points, value.points, name, std.dev=NULL, 
                  num.meas.per.value=NULL, site=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time.points}{
  The time points of the time series in chronological (ascending) order. Does not need 
 }
  \item{value.points}{
 The corresponding measurements of the time series (or the average measurement if more than one goes into the value, if 'num.meas.per.value' is given and larger than 1). Must necessarily be the same size as 'time.points'.
 }
  \item{name}{
 The name of the time series. Will be used in assigning parameter names for parameters associated with this series.
 }
  \item{std.dev}{
 The individual measurement standard deviations of each measurement (if applicable). If this does not exist, the model will estimate an overall measurement standard deviation belonging to all measurements. If 'num.meas.per.value', these are the standard deviations of of the individual measurements going into each value. Must necessarily be the same size as 'time.points', if given.
 }
  \item{num.meas.per.value}{
 The number of each individual measurement that goes into each value in the time series. Must necessarily be the same size as 'time.points', if given.
 }
  \item{site}{
 The site number of each measurement, if applicable. Should be labelled from 0 to <number of sites>-1. Must necessarily be the same size as 'time.points', if given.
 }
}
\details{
This class function is used as input to the analysis and represents a single time series. It should contain time points and values. The values are the measurements or measurement averages, if many individual at a given time point goes into each value. Standard deviation of the individual values can also be given and the number of individual measurements per value, if applicable. Multiple sites are also possible. PS: There are extra options for the process structure if sites are involved.
}
\value{
Returns a time series data set as a 'layer.data.series' object.
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
See also read.layer.data.series and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Create a sample from the Ornstein-Uhlenbeck (OU) process with mu=0, 
# characteristic time=50 and sigma=0.2:
n=1000
x=rep(rnorm(1),n)
for(i in 2:n)
 x[i]=exp(-1/50)*x[i-1]+0.2*rnorm(1)

# Create measurement set as a subsample of the process, plus 
# noise with standard deviation 0.1
t=sort(sample(1:n,300))
y=(x+rnorm(length(x),0,0.1))[t]

# Create the data series object from the time points and "measurements":
X=layer.data.series(time.points=t, value.points=y, name="X")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, time series }
