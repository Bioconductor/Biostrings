\name{readFASTA}
\alias{readFASTA}
\title{ A function to read FASTA formatted files}
\description{
 FASTA is a simple file format for biological sequence data. A file may
contain one or more sequences, for each sequence there is a description 
line which begins with a \code{>}.
}
\usage{
readFASTA(file, checkComments = TRUE)
}
\arguments{
  \item{file}{ The name of the file to be read from.}
  \item{checkComments}{ Whether or not comments, lines beginning with
a semi-colon should be found and removed. }
}
\details{
 FASTA is a widely used format in biology. It is a relatively simple markup.
 I am not aware of a standard. It might be nice to check to see if the 
 data that were parsed are sequences of some appropriate type, but without
 a standard that does not seem possible.

 There are many other packages that provide similar, but different 
 capabilities.  The one in the package seqinr seems most similar but they
 separate the biological sequence into single character strings, which
 is too inefficient for large problems.
}
\value{
 A list with one element for each sequence in the file. The elements 
are in two parts, one the description and the second a character string
of the biological sequence.
}
\author{ R. Gentleman }

\seealso{\code{\link{read.table}}, \code{\link{scan}}}
\examples{
f1 = system.file("Exfiles/fastaEx", package="Biostrings")
ff = readFASTA(f1)
}
\keyword{ manip}
