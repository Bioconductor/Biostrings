\name{reverseSeq}
\alias{reverseSeq}
\alias{revcompDNA}
\alias{revcompRNA}

\title{Reverse Sequence}
\description{
  WARNING: The functions described in this man page have been deprecated
  in favor of \code{\link{reverse,XString-method}} and \code{\link{reverseComplement}}.

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
 if (interactive()) {
   reverseSeq(w)  # deprecated (inefficient on large vectors)
 }
 reverse(BStringSet(w))  # more efficient

 w <- "able was I ere I saw Elba"
 if (interactive()) {
   reverseSeq(w)  # deprecated (inefficient on large vectors)
 }
 reverse(BStringSet(w))  # more efficient

 rna1 <- "UGCA"
 if (interactive()) {
   revcompRNA(rna1)  # deprecated (inefficient on large vectors)
 }
 reverseComplement(RNAString(rna1))  # more efficient

 dna1 <- "TGCA"
 if (interactive()) {
   revcompDNA(dna1)  # deprecated (inefficient on large vectors)
 }
 reverseComplement(DNAString(dna1))  # more efficient

 ## Comparing efficiencies:
 if (interactive()) {
   library(hgu95av2probe)
   system.time(y1 <- reverseSeq(hgu95av2probe$sequence))
   x <- DNAStringSet(hgu95av2probe)
   system.time(y2 <- reverse(x))

   system.time(y3 <- revcompDNA(hgu95av2probe$sequence))
   system.time(y4 <- reverseComplement(x))
 }
}
\keyword{manip}
