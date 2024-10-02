\name{n.combinations}
\alias{n.combinations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns all ways of combining 'n' objects going form 0 to <base>. }
\description{
Gives all ways of combining 'n' objects going form 0 to <base>. Used for
traversing all model combinations in "traverse.standalone.layered"
and "traverse.connections.layered".

}
\usage{
n.combinations(length, base.number=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{length}{
  The number of objects to be combined.
  }
  \item{base.number}{
  The span of the possible object values, going from 0,1,...,base.number.
  Thus, base.number=2 (default) only allows "0" and "1".
 }
}
\details{
Used for going through all combination of options in
"traverse.standalone.layered"
and "traverse.connections.layered".
}
\value{
A matrix where each row represents a unique sequence of 
\{0,1,...,base.number\}.
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
See also traverse.standalone.layered
and traverse.connections.layered.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n.combinations(3,2)
#    [,1] [,2] [,3]
#[1,]    0    0    0
#[2,]    0    0    1
#[3,]    0    1    0
#[4,]    0    1    1
#[5,]    1    0    0
#[6,]    1    0    1
#[7,]    1    1    0
#[8,]    1    1    1

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ number of combinations, model selection }
