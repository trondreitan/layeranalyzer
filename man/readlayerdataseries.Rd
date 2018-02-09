\name{read.layer.data.series}
\alias{read.layer.data.series}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Represents a measurement time series. }
\description{
This class function reads a file representing a measurement time series. Can include multiple sites. Individual measurement error can be added. 
}
\usage{
read.layer.data.series(filename,column.type=c("time","value"), 
                name="",
                header=FALSE,dec=".",sep="", 
                quote="\"'", numerals = c("allow.loss", "warn.loss", "no.loss"),
                row.names=NULL, col.names=column.type, as.is = !stringsAsFactors,
                na.strings = "NA", colClasses = NA, nrows = -1,
                skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                strip.white = FALSE, blank.lines.skip = TRUE,
                comment.char = "#",
                allowEscapes = FALSE, flush = FALSE,
                stringsAsFactors = default.stringsAsFactors(),
                fileEncoding = "", encoding = "unknown", text, skipNul = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  Name of the input file.
 }
  \item{column.type}{
 This should specify the content of each column in the file. The file must contain a "time" and a "value" column, as that is the bare minimum in order to specify a time series. In addition, the a measurement standard deviation column can be specified by "std.dev". If in addition a "num.meas.per.value" column exists, the standard error of each mueasurement will be taken to be the content of the "std.dev" column divied by the square root of the content of the "num.meas.per.value" column. A "site" column can also be added, if the time series contains the same typ of measurements from different sites (note that regional options then becomes a possibility). The site index must run from 0 to the number of sites minus one.
 }
  \item{name}{
 The name given to the series. Will be used in assigning parameter names for parameters associated with this series.
 }
  \item{The rest of the input}{
 The rest of the input is sent directly to the 'read.table' routine.
 }
}
\details{
This routine reads from a file representing a time series and creates a 'layer.data.series' object which can be used in the analysis. The file should at minimum contain two columns containing the time points and the corresponding measurements. Which column is value and which is time should be specified in the input. (If not, it is assumed that time points comes first, then values). The file can also have additional columns descibing the standard deviation, the number of individual measurements per value and site indexes. 
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
See also layer.data.series and layer.series.structure.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Read a simulated 2 layer data set from an external web site:
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, time series }
