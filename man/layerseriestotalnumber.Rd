\name{layer.series.total.number}
\alias{layer.series.total.number}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reading prior specification for the layer.analyzer from a file. }
\description{
A variable that counts the total number of series (of class 'layer.data.series')
that has been read in this session.
}
\usage{
layer.series.total.number
}
%- maybe also 'usage' for other objects documented here.
\details{
This variable keeps track of the number of 'layer.data.series' objects made.
Used for making a unique series name when no series name has been given  by the user.
}
\value{
Integer.
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
layer.series.total.number
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers }


