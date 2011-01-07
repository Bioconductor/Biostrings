\name{matchprobes}
\alias{matchprobes}
\title{A function to match a query sequence to the sequences of a set of
  probes.}
\description{
  The \code{query} sequence, a character string (probably representing
  a transcript of interest), is scanned for the presence of exact
  matches to the sequences in the character vector \code{records}.
  The indices of the set of matches are returned.

  The function is inefficient: it works on R's character vectors, and
  the actual matching algorithm is of time complexity \code{length(query)}
  times \code{length(records)}!

  See \code{\link{matchPattern}}, \code{\link{vmatchPattern}} and
  \code{\link{matchPDict}} for more efficient sequence matching functions.
}
\usage{
matchprobes(query, records, probepos=FALSE)
}
\arguments{
  \item{query}{A character vector. For example, each element may represent
    a gene (transcript) of interest. See Details.}
  \item{records}{A character vector. For example, each element may represent
    the probes on a DNA array.}
  \item{probepos}{A logical value. If TRUE, return also the start
    positions of the matches in the query sequence.}
}
\details{
  \code{\link{toupper}} is applied to the arguments \code{query} and
  \code{records} before matching. The intention of this is to make
    the matching case-insensitive.
    The function is embarrassingly naive.
    The matching is done using the C library function \code{strstr}. 
}
\value{
  A list.
  Its first element is a list of the same length as the input vector.
  Each element of the list is a numeric vector containing the indices of
  the probes that have a perfect match in the query
  sequence.

  If \code{probepos} is TRUE,
  the returned list has a second element: it is of the same shape
  as described above, and gives the respective positions of the
  matches.
}
\author{R. Gentleman, Laurent Gautier, Wolfgang Huber}

\seealso{\code{\link{matchPattern}}, \code{\link{vmatchPattern}}, \code{\link{matchPDict}}}

\examples{
  if(require("hgu95av2probe")){
    data("hgu95av2probe")
    seq <- hgu95av2probe$sequence[1:20]
    target <- paste(seq, collapse="")
    matchprobes(target, seq, probepos=TRUE)
  }
}
\keyword{manip}
