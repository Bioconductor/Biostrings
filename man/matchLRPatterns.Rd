\name{matchLRPatterns}

\alias{matchLRPatterns}
\alias{matchLRPatterns,BString-method}


\title{Find patterns defined by a left and a right part in a sequence}

\description{
}

\usage{
  matchLRPatterns(Lpattern, Rpattern, max.ngaps, subject, Lmismatch=0, Rmismatch=0,
                  Lfixed=TRUE, Rfixed=TRUE)
}

\arguments{
  \item{Lpattern}{
    The left part of the pattern.
  }
  \item{Rpattern}{
    The right part of the pattern.
  }
  \item{max.ngaps}{
    The max number of gaps in the middle i.e the max distance between the left
    and right parts of the pattern.
  }
  \item{subject}{
    A \link{BString} (or derived) object containing the target sequence.
  }
  \item{Lmismatch}{
    The number of mismatches allowed for the left part of the pattern.
    If non-zero, an inexact matching algorithm is used (see the
    \code{\link{matchPattern}} function for more information).
  }
  \item{Rmismatch}{
    Same as \code{Lmismatch} but for the right part of the pattern.
  }
  \item{Lfixed}{
    Only with a \link{DNAString} or \link{RNAString} subject can a \code{Lfixed}
    value other than the default (\code{TRUE}) be used.

    With \code{Lfixed=FALSE}, ambiguities (i.e. letters from the IUPAC Extended
    Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}}) that are not from the
    base alphabet) in the left pattern \_and\_ in the subject are interpreted as
    wildcards i.e. they match any letter that they stand for.

    See the \code{fixed} argument of the \code{\link{matchPattern}} function
    for more information.
  }
  \item{Rfixed}{
    Same as \code{Lfixed} but for the right part of the pattern.
  }
}

\details{
}

\value{
  A \link{BStringViews} object containing all the matches, even when they are
  overlapping (see the examples below), and where the matches are ordered
  from left to right (i.e. by ascending starting position).
}

\author{H. Pages}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchProbePair}},
  \code{\link{findPalindromes}},
  \code{\link{reverseComplement}},
  \code{\link{BStringViews}}
}

\examples{
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  subject <- Dmelanogaster[["3R"]]
  Lpattern <- "AGCTCCGAG"
  Rpattern <- "TTGTTCACA"
  matchLRPatterns(Lpattern, Rpattern, 500, subject) # 1 match

  ## Note that matchLRPatterns() will return all matches, even when they are
  ## overlapping:
  subject <- DNAString("AAATTAACCCTT")
  matchLRPatterns("AA", "TT", 0, subject) # 1 match
  matchLRPatterns("AA", "TT", 1, subject) # 2 matches
  matchLRPatterns("AA", "TT", 3, subject) # 3 matches
  matchLRPatterns("AA", "TT", 7, subject) # 4 matches
}

\keyword{methods}

