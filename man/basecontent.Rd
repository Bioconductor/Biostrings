\name{basecontent}
\alias{basecontent}
\alias{countbases}

\title{Obtain the ATCG content of a gene}
\description{
  WARNING: Both \code{basecontent} and \code{countbases} are now defunct
  and will be removed soon together with this man page.
  See the examples at the bottom of the man page for how to use
  \code{\link{alphabetFrequency}} instead.

  These functions accept a character vector representing the nucleotide
  sequences and compute the frequencies of each base (A, C, G, T).
}
\usage{
basecontent(seq)
countbases(seq, dna = TRUE)
}
\arguments{
  \item{seq}{Character vector.}
  \item{dna}{Logical value indicating whether the sequence is DNA
    (\code{TRUE}) or RNA (\code{FALSE})}
}
\details{
  The base frequencies are calculated separately for
  each element of \code{x}.
  The elements of \code{x} can be in upper case, lower case
  or mixed. 
}

\value{
  A matrix with 4 columns and \code{length(x)} rows.  The columns are
  named \code{A}, \code{C}, \code{T}, \code{G}, and the values in each
  column are the counts of the corresponding bases in the elements of
  \code{x}.  When \code{dna=FALSE}, the \code{T} column is replaced with
  a \code{U} column.
}

\author{R. Gentleman, W. Huber, S. Falcon}

\seealso{\code{\link{alphabetFrequency}}, \code{\link{reverseComplement}}}

\examples{
 v<-c("AAACT", "GGGTT", "ggAtT")

 ## You can't use these functions anymore (defunct):
 if (interactive()) {
   basecontent(v)
   countbases(v)
 }

 ## But you can do this instead:
 v <- DNAStringSet(v)
 alphabetFrequency(v, baseOnly=TRUE)
}
\keyword{manip}
