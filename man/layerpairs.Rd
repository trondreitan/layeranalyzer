\name{layer.pairs}
\alias{layer.pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Alternative connection specification, taking a structure list. }
\description{
 This function takes a set of
 of "layer.series.structure" objects as input and
 then shows the possible pairings in the order that an
 alternative connection specification (using 'conn' in 'layer.analzyer')
 will refer to them.
}
\usage{
layer.pairs(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{...}{
 A set of 'layer.series.structure' objects that represents the set of 
 time series that are to be analyzed and how the underlying processes of
 each such series is to be modelles. (See 'layer.series.structure'). At
 least two such objects must be given. 
 }
}
\details{

The function takes a set of time series along with their process
structure (represented by layer.series.structure objects),
and outputs an array describing the pairings and also a
set of strings containing a verbal description.

Used by
'layer.analyzer' via the 'layer.connections' function, for assigning
and referring to pairs of processes when using alternative
connection specification (the 'conn' option in 'layer.analyzer').

This function  is also provided
separately for a user that wants to see what the different paris are
and how they are ordered.

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
Returns an object of type 'layer.pairs', a list with
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
See also layer.connections, layer.connections.list, layer.pairs.list and
layer.analyzer.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Creates two 2-layered structures based on the two hare/lynx datasets.
# outputs the pairs and a description

hare.st=layer.series.structure(hare.norm,numlayers=2)
lynx.st=layer.series.structure(lynx.norm,numlayers=2)

layer.pairs(hare.st, lynx.st)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ process pairs, alternative connection specification, linear stochastic differential equations, layered, hidden layers, causal, correlative}


