\name{reverseSeq}
\alias{reverseSeq}
\alias{revcompDNA}
\alias{revcompRNA}

\title{Reverse Sequence}
\description{
  WARNING: The functions described in this man page are now defunct
  and will be removed soon together with this man page.
  See the examples at the bottom of the man page for how to use
  \code{\link{reverse,XString-method}} and \code{\link{reverseComplement}}
  instead.

  Functions to obtain the reverse and reverse complement of a sequence}
\usage{
reverseSeq(seq)
revcompDNA(seq)
revcompRNA(seq)
}
\arguments{
  \item{seq}{Character vector.  For \code{revcompRNA} and
  \code{revcompDNA} the sequence should consist of appropriate letter
  codes: \code{[ACGUN]} and \code{ACGTN}, respectively.}
}
\details{
  The function reverses the order of the constituent
  character strings of its argument.
}
\value{
  A character vector of the same length as \code{seq}.
}

\author{R. Gentleman, W. Huber, S. Falcon}

\seealso{\code{\link{alphabetFrequency}}, \code{\link{reverseComplement}}}

\examples{
 w <-  c("hey there", "you silly fool")
 ## You can't do this anymore (defunct):
 if (interactive()) {
   reverseSeq(w)  # was inefficient on large vectors
 }
 ## But you can do this instead:
 reverse(BStringSet(w))

 w <- "able was I ere I saw Elba"
 ## You can't do this anymore (defunct):
 if (interactive()) {
   reverseSeq(w)  # was inefficient on large vectors
 }
 ## But you can do this instead:
 reverse(BStringSet(w))

 rna1 <- "UGCA"
 ## You can't do this anymore (defunct):
 if (interactive()) {
   revcompRNA(rna1)  # was inefficient on large vectors
 }
 ## But you can do this instead:
 reverseComplement(RNAString(rna1))

 dna1 <- "TGCA"
 ## You can't do this anymore (defunct):
 if (interactive()) {
   revcompDNA(dna1)  # was inefficient on large vectors
 }
 ## But you can do this instead:
 reverseComplement(DNAString(dna1))
}
\keyword{manip}
