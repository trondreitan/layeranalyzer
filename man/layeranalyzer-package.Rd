\name{layeranalyzer}
\alias{layeranalyzer}
\docType{package}
\title{
\packageTitle{layeranalyzer}
}
\description{
\packageDescription{layeranalyzer}
}
\details{

The DESCRIPTION file:
\packageIndices{layeranalyzer}
Time series analysis tool using linear layered SDEs. 
 Allows for multiple time series with correlative or causal links between them.
 In case of causal feedback loops, the matrix operations (including 
 eigenvalue decompositions) allows for complex numbers. In this case,
 cyclic behavior can be expected.
}
\author{
\packageAuthor{layeranalyzer}

Maintainer: \packageMaintainer{layeranalyzer}
}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.
}
\keyword{ package }
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
X=layer.data.series(time.points=t, value.points=y,name="X")

# Specify layered structure. (In this case default structure, one layered=OU):
X.struct=layer.series.structure(X)

# Perform analysis:
#res=layer.analyzer(X.struct)

# Show the results:
#summary(res)

}
