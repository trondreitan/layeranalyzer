\name{read.layer.data.series}
\alias{read.layer.data.series}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read a series representing a measurement time series. }
\description{
This class function reads a file representing a measurement time series. Can include multiple sites. Individual measurement error can be added. 
}
\usage{
read.layer.data.series(filename,name,column.type=c("time","value"), 
                header=FALSE,dec=".",sep="", 
                quote="\"'", numerals = c("allow.loss", "warn.loss", "no.loss"),
                row.names=NULL, col.names=column.type, 
                na.strings = "NA", colClasses = NA, nrows = -1,
                skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                strip.white = FALSE, blank.lines.skip = TRUE,
                comment.char = "#",
                allowEscapes = FALSE, flush = FALSE,
                stringsAsFactors = FALSE, as.is = !stringsAsFactors,
                fileEncoding = "", encoding = "unknown", text, skipNul = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  Name of the input file.
 }
  \item{name}{
 The name given to the series. Will be used in assigning parameter names for parameters associated with this series.
 }
  \item{column.type}{
 This should specify the content of each column in the file. The file must contain a "time" and a "value" column, as that is the bare minimum in order to specify a time series. In addition, the a measurement standard deviation column can be specified by "std.dev". If in addition a "num.meas.per.value" column exists, the standard error of each mueasurement will be taken to be the content of the "std.dev" column divied by the square root of the content of the "num.meas.per.value" column. A "site" column can also be added, if the time series contains the same typ of measurements from different sites (note that regional options then becomes a possibility). The site index must run from 0 to the number of sites minus one.
 }

 \item{header}{ This and the rest of the options are just sent to
          "read.table", so the help text has just been lifted from the
	  manual of that function. You should check that manual for
	  an updated version of this help text.
	  
          A logical value indicating whether the file contains the
          names of the variables as its first line.  If missing, the
          value is determined from the file format: ‘header’ is set to
          ‘TRUE’ if and only if the first row contains one fewer field
          than the number of columns.
 }
 \item{dec}{ the character used in the file for decimal points.
  }
 \item{sep}{ the field separator character.  Values on each line of the
          file are separated by this character.  If ‘sep = ""’ (the
          default for ‘read.table’) the separator is ‘white space’,
          that is one or more spaces, tabs, newlines or carriage
          returns.
  }
  \item{quote}{ the set of quoting characters. To disable quoting altogether,
          use ‘quote = ""’.  See ‘scan’ for the behaviour on quotes
          embedded in quotes.  Quoting is only considered for columns
          read as character, which is all of them unless ‘colClasses’
          is specified.
  }
  \item{numerals}{ string indicating how to convert numbers whose conversion to
          double precision would lose accuracy, see ‘type.convert’.
          Can be abbreviated.  (Applies also to complex-number inputs.)
  }
  \item{row.names}{ a vector of row names.  This can be a vector giving the
          actual row names, or a single number giving the column of the
          table which contains the row names, or character string
          giving the name of the table column containing the row names.
 
          If there is a header and the first row contains one fewer
          field than the number of columns, the first column in the
          input is used for the row names.  Otherwise if ‘row.names’ is
          missing, the rows are numbered.

          Using ‘row.names = NULL’ forces row numbering. Missing or
          ‘NULL’ ‘row.names’ generate row names that are considered to
          be ‘automatic’ (and not preserved by ‘as.matrix’).
  }
  \item{col.names}{ a vector of optional names for the variables.  The default
          is to use ‘"V"’ followed by the column number.
  }
  \item{na.strings}{ a character vector of strings which are to be interpreted
          as ‘NA’ values.  Blank fields are also considered to be
          missing values in logical, integer, numeric and complex
          fields.  Note that the test happens _after_ white space is
          stripped from the input, so ‘na.strings’ values may need
          their own white space stripped in advance.
  }
  \item{colClasses}{ character.  A vector of classes to be assumed for the
          columns.  If unnamed, recycled as necessary.  If named, names
          are matched with unspecified values being taken to be ‘NA’.

          Possible values are ‘NA’ (the default, when ‘type.convert’ is
         used), ‘"NULL"’ (when the column is skipped), one of the
          atomic vector classes (logical, integer, numeric, complex,
          character, raw), or ‘"factor"’, ‘"Date"’ or ‘"POSIXct"’.
          Otherwise there needs to be an ‘as’ method (from package
          ‘methods’) for conversion from ‘"character"’ to the specified
          formal class.

          Note that ‘colClasses’ is specified per column (not per
          variable) and so includes the column of row names (if any).
  }
  \item{nrows}{ integer: the maximum number of rows to read in.  Negative and
          other invalid values are ignored.
  }
  \item{skip}{ integer: the number of lines of the data file to skip before
          beginning to read data.
  }
  \item{check.names}{ logical.  If ‘TRUE’ then the names of the variables in the
          data frame are checked to ensure that they are syntactically
          valid variable names.  If necessary they are adjusted (by
          ‘make.names’) so that they are, and also to ensure that there
          are no duplicates.
  }
  \item{fill}{ logical. If ‘TRUE’ then in case the rows have unequal length,
         blank fields are implicitly added.  See ‘Details’.
  }
  \item{strip.white}{ logical. Used only when ‘sep’ has been specified, and
          allows the stripping of leading and trailing white space from
          unquoted ‘character’ fields (‘numeric’ fields are always
          stripped).  See ‘scan’ for further details (including the
          exact meaning of ‘white space’), remembering that the columns
          may include the row names.
  }
  \item{blank.lines.skip}{ logical: if ‘TRUE’ blank lines in the input are
          ignored.
  }
  \item{comment.char}{ character: a character vector of length one containing a
          single character or an empty string.  Use "" to turn off
          the interpretation of comments altogether.
  }
  \item{allowEscapes}{ logical.  Should C-style escapes such as "\\n" be
          processed or read verbatim (the default)?  Note that if not
          within quotes these could be interpreted as a delimiter (but
          not as a comment character).  For more details see ‘scan’.
  }
  \item{flush}{ logical: if "TRUE", "scan" will flush to the end of the line
          after reading the last of the fields requested.  This allows
         putting comments after the last field.
  }
  \item{stringsAsFactors}{ logical: should character vectors be converted to
          factors?  Note that this is overridden by ‘as.is’ and
          ‘colClasses’, both of which allow finer control.
  }
  \item{as.is}{ the default behavior of ‘read.table’ is to convert character
          variables (which are not converted to logical, numeric or
          complex) to factors.  The variable ‘as.is’ controls the
          conversion of columns not otherwise specified by
          ‘colClasses’.  Its value is either a vector of logicals
         (values are recycled if necessary), or a vector of numeric or
          character indices which specify which columns should not be
          converted to factors.

          Note: to suppress all conversions including those of numeric
          columns, set ‘colClasses = "character"’.

          Note that ‘as.is’ is specified per column (not per variable)
          and so includes the column of row names (if any) and any
          columns to be skipped.
  }
  \item{fileEncoding}{ character string: if non-empty declares the encoding used
          on a file (not a connection) so the character data can be
          re-encoded.  See the ‘Encoding’ section of the help for
          ‘file’, the ‘R Data Import/Export Manual’ and ‘Note’.
  }
  \item{encoding}{ encoding to be assumed for input strings.  It is used to mark
          character strings as known to be in Latin-1 or UTF-8 (see
          ‘Encoding’): it is not used to re-encode the input, but
          allows R to handle encoded strings in their native encoding
          (if one of those two).  See ‘Value’ and ‘Note’.
  }
  \item{text}{ character string: if ‘file’ is not supplied and this is, then
          data are read from the value of ‘text’ via a text connection.
          Notice that a literal string can be used to include (small)
          data sets within R code.
  }
  \item{skipNul}{logical: should nuls be skipped?
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

# Read a simulated 2 layer data set from an external web site, call it "Xseries":
X=read.layer.data.series("http://folk.uio.no/trondr/R/test_2layer.txt","Xseries")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear stochastic differential equations, layered, hidden layers, causal, correlative, matrix inverse, Lapack, eigenvalue decomposition, time series }
