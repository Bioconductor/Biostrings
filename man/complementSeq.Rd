\name{complementSeq}
\alias{complementSeq}

\title{Complementary sequence.}
\description{
  WARNING: \code{complementSeq} has been
  deprecated in favor of \code{\link{complement}}.

  Function to obtain the complementary sequence.
}
\usage{
complementSeq(seq, start=1, stop=0)
}
\arguments{
  \item{seq}{Character vector consisting of the letters A, C, G and T.}
  \item{start}{Numeric scalar: the sequence position at which to start
    complementing. If 1, start from the beginning.}
  \item{stop}{Numeric scalar: the sequence position at which to stop
    complementing. If 0, go until the end.}
}
\details{
  The complemented sequence for each element of the input is computed and
  returned.  The complement is given by the mapping:
  A -> T, C -> G, G -> C, T -> A.

  An important special case is \code{start=13}, \code{stop=13}:
  If \code{seq} is a vector of 25mer sequences on an Affymetrix
  GeneChip, \code{complementSeq(seq, start=13, stop=13)}
  calculates the so-called \emph{mismatch} sequences.

  The function deals only with sequences that represent DNA.
  These can consist only of the letters  \code{A}, \code{C}, \code{T}
  or \code{G}. Upper, lower or mixed case is allowed and honored.
  
}
\value{
  A character vector of the same length as \code{seq} is
  returned. Each component represents the transformed sequence for the
  input value.
}

\author{R. Gentleman, W. Huber}

\seealso{\code{\link{alphabetFrequency}}, \code{\link{reverseComplement}}}

\examples{
  ## ---------------------------------------------------------------------
  ## EXAMPLE 1
  ## ---------------------------------------------------------------------
  seq <- c("AAACT", "GGGTT")

  ## Don't do this anymore (deprecated):
  if (interactive()) {
    complementSeq(seq)  # inefficient on large vectors
  }
  ## But do this instead:
  complement(DNAStringSet(seq))  # more efficient

  ## ---------------------------------------------------------------------
  ## EXAMPLE 2
  ## ---------------------------------------------------------------------
  seq <- c("CGACTGAGACCAAGACCTACAACAG", "CCCGCATCATCTTTCCTGTGCTCTT")

  ## Don't do this anymore (deprecated):
  if (interactive()) {
    complementSeq(seq, start=13, stop=13)
  }
  ## But do this instead:
  pm2mm <- function(probes)
  {
      parts <- threebands(DNAStringSet(probes), start=13, end=13)
      xscat(parts$left, complement(parts$middle), parts$right)
  }
  pm2mm(seq) 

  ## ---------------------------------------------------------------------
  ## SPEED OF complementSeq() VS complement()
  ## ---------------------------------------------------------------------
  if (interactive()) {
    library(hgu95av2probe)
    system.time(y1 <- complementSeq(hgu95av2probe$sequence))
    probes <- DNAStringSet(hgu95av2probe$sequence)
    system.time(y2 <- complement(probes))
  }
}
\keyword{manip}
