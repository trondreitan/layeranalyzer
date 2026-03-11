\name{layer.connections.list}
\alias{layer.connections.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Alternative connection specification, taking a structure list. }
\description{
 This is a version of "layer.connections" that takes a list
 of "layer.series.structure" objects as input together with
 an alternative connecton specification and translates
 that into the default connection specification.
}
\usage{
layer.connections.list(data.structure , conn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.structure}{
 A list containing 'layer.series.structure' objects that represents the set of 
 time series that are to be analyzed and how the underlying processes of
 each such series is to be modelles. (See 'layer.series.structure'). At
 least two such objects must be given. 
 }
 \item{conn}{
 An array of integers between 0 and 4, one for each possible pairing
 of two processes (series+layer) belonging to different series
 (see 'layer.pairs').
 For each pair, 0 means no connection between the two,
 1 means correlative connection between the two,
 2 means causal connection from the first process in the pair
 to the second, 3 means causal connection from the second
 process in the pair to the first and 4 means causal connection
 going both ways between the pair of processes. 
 }
}
\details{

The function takes a set of time series along with their process
structure (represented by layer.series.structure objects),
and an array of integers specifying the nature of the connection
between each possible pair of processes. It then outputs
a specification of what the corresponding specification for
the 'causal' and 'corr' input for layer.analyzer should be.

This function is used in 'layer.analzyer' but is also provided
separately for a user that wants to see how this alternative
connection specification works.

The alternative way of specifying a connection was introduced in order
to make it easy to produce different connection models by
simply traversing the possible connection ebtween each pair of
processes. If one use a script to traverse all possible
connections models, reads a script and corrects the connection
specification, one has a tool for outputing one script for each
connection models. These scripts can then possibly be run in parallel
on a computer with many cores or a computing cluster. 

}
\value{
Returns an object of type 'layer.connections', a list
two specifications called 'causal' and 'corr'. In 'layer.analyzer',
these specifications are then sent further into the analysis system.
But the output can also be used for checking that the user understands
the alternative connection specification. 
}
\references{
Reitan, T., Schweder, T., Henderiks, J. (2012),
Phenotypic Evolution studied by Layered Stochastic Differential Equations,
Annals of Applied Statistics, Volume 6 (4): 1531-1551.
}
\author{
Trond Reitan, trondr@uio.no
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also layer.connections, layer.pairs, layer.pairs.list and
layer.analyzer.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Creates two 2-layered structures based on the two hare/lynx datasets.
# outputs the connection structure specified by the
# alternative conenction specification array c(0,1,2,4)

hare.st=layer.series.structure(hare.norm,numlayers=2)
lynx.st=layer.series.structure(lynx.norm,numlayers=2)

layer.connections.list(list(hare.st, lynx.st), conn=c(0,1,2,4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ process pairs, alternative connection specification, linear stochastic differential equations, layered, hidden layers, causal, correlative}


