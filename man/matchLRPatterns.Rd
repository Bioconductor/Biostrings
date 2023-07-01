\name{matchLRPatterns}

\alias{matchLRPatterns}
\alias{matchLRPatterns,XString-method}
\alias{matchLRPatterns,XStringViews-method}
\alias{matchLRPatterns,MaskedXString-method}


\title{Find paired matches in a sequence}

\description{
  The \code{matchLRPatterns} function finds paired matches in a sequence
  i.e. matches specified by a left pattern, a right pattern and a maximum
  distance between the left pattern and the right pattern.
}

\usage{
matchLRPatterns(Lpattern, Rpattern, max.gaplength, subject,
                max.Lmismatch=0, max.Rmismatch=0,
                with.Lindels=FALSE, with.Rindels=FALSE,
                Lfixed=TRUE, Rfixed=TRUE)
}

\arguments{
  \item{Lpattern}{
    The left part of the pattern.
  }
  \item{Rpattern}{
    The right part of the pattern.
  }
  \item{max.gaplength}{
    The max length of the gap in the middle i.e the max distance between
    the left and right parts of the pattern.
  }
  \item{subject}{
    An \link{XString}, \link{XStringViews} or \link{MaskedXString} object
    containing the target sequence.
  }
  \item{max.Lmismatch}{
    The maximum number of mismatching letters allowed in the left part of the
    pattern.
    If non-zero, an inexact matching algorithm is used (see the
    \code{\link{matchPattern}} function for more information).
  }
  \item{max.Rmismatch}{
    Same as \code{max.Lmismatch} but for the right part of the pattern.
  }
  \item{with.Lindels}{
    If \code{TRUE} then indels are allowed in the left part of the pattern.
    In that case \code{max.Lmismatch} is interpreted as the maximum "edit
    distance" allowed in the left part of the pattern.

    See the \code{with.indels} argument of the \code{\link{matchPattern}}
    function for more information.
  }
  \item{with.Rindels}{
    Same as \code{with.Lindels} but for the right part of the pattern.
  }
  \item{Lfixed}{
    Only with a \link{DNAString} or \link{RNAString} subject can a
    \code{Lfixed} value other than the default (\code{TRUE}) be used.

    With \code{Lfixed=FALSE}, ambiguities (i.e. letters from the IUPAC
    Extended Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}}) that
    are not from the base alphabet) in the left pattern \emph{and} in
    the subject are interpreted as wildcards i.e. they match any letter
    that they stand for.

    \code{Lfixed} can also be a character vector, a subset of
    \code{c("pattern", "subject")}.
    \code{Lfixed=c("pattern", "subject")} is equivalent to
    \code{Lfixed=TRUE} (the default).
    An empty vector is equivalent to \code{Lfixed=FALSE}.
    With \code{Lfixed="subject"}, ambiguities in the pattern only are
    interpreted as wildcards.
    With \code{Lfixed="pattern"}, ambiguities in the subject only are
    interpreted as wildcards.
  }
  \item{Rfixed}{
    Same as \code{Lfixed} but for the right part of the pattern.
  }
}

\value{
  An \link{XStringViews} object containing all the matches, even when they are
  overlapping (see the examples below), and where the matches are ordered
  from left to right (i.e. by ascending starting position).
}

\author{H. Pagès}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchProbePair}},
  \code{\link{trimLRPatterns}},
  \code{\link{findPalindromes}},
  \code{\link{reverseComplement}},
  \link{XString-class},
  \link{XStringViews-class},
  \link{MaskedXString-class}
}

\examples{
library(BSgenome.Dmelanogaster.UCSC.dm3)
subject <- Dmelanogaster$chr3R
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
