\name{matchProbePair}

\alias{matchProbePair}
\alias{matchProbePair,DNAString-method}
\alias{matchProbePair,XStringViews-method}
\alias{matchProbePair,MaskedDNAString-method}


\title{Find "theoretical amplicons" mapped to a probe pair}

\description{
  In the context of a computer-simulated PCR experiment, one wants to find
  the amplicons mapped to a given primer pair.
  The \code{matchProbePair} function can be used for this: given a forward and a
  reverse probe (i.e. the chromosome-specific sequences of the forward and
  reverse primers used for the experiment) and a target sequence (generally a
  chromosome sequence), the \code{matchProbePair} function will return all
  the "theoretical amplicons" mapped to this probe pair.
}

\usage{
matchProbePair(Fprobe, Rprobe, subject, algorithm="auto", logfile=NULL, verbose=FALSE)
}

\arguments{
  \item{Fprobe}{
    The forward probe.
  }
  \item{Rprobe}{
    The reverse probe.
  }
  \item{subject}{
    A \link{DNAString} object (or an \link{XStringViews} object
    with a \link{DNAString} subject) containing the target sequence.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"} or \code{"shift-or"}.
    See \code{\link{matchPattern}} for more information.
  }
  \item{logfile}{
    A file used for logging.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
}

\details{
  The \code{matchProbePair} function does the following: (1) find all
  the "plus hits" i.e. the Fprobe and Rprobe matches on the "plus" strand,
  (2) find all the "minus hits" i.e. the Fprobe and Rprobe matches on the
  "minus" strand and (3) from the set of all (plus_hit, minus_hit) pairs,
  extract and return the subset of "reduced matches" i.e. the (plus_hit,
  minus_hit) pairs such that (a) plus_hit <= minus_hit and (b) there are
  no hits (plus or minus) between plus_hit and minus_hit.
  This set of "reduced matches" is the set of "theoretical amplicons".
}

\value{
  An \link{XStringViews} object containing the set of "theoretical amplicons".
}

\author{H. Pagès}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \code{\link{findPalindromes}},
  \code{\link{reverseComplement}},
  \link{XStringViews-class}
}

\examples{
library(BSgenome.Dmelanogaster.UCSC.dm3)
subject <- Dmelanogaster$chr3R

## With 20-nucleotide forward and reverse probes:
Fprobe <- "AGCTCCGAGTTCCTGCAATA"
Rprobe <- "CGTTGTTCACAAATATGCGG"
matchProbePair(Fprobe, Rprobe, subject) # 1 "theoretical amplicon"

## With shorter forward and reverse probes, the risk of having multiple
## "theoretical amplicons" increases:
Fprobe <- "AGCTCCGAGTTCC"
Rprobe <- "CGTTGTTCACAA"
matchProbePair(Fprobe, Rprobe, subject) # 2 "theoretical amplicons"
Fprobe <- "AGCTCCGAGTT"
Rprobe <- "CGTTGTTCACA"
matchProbePair(Fprobe, Rprobe, subject) # 9 "theoretical amplicons"
}

\keyword{methods}
