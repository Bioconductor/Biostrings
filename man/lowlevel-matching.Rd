\name{lowlevel-matching}

\alias{lowlevel-matching}

\alias{neditStartingAt}
\alias{neditEndingAt}
\alias{neditAt}
\alias{neditStartingAt,character-method}
\alias{neditStartingAt,XString-method}
\alias{neditStartingAt,XStringSet-method}
\alias{neditEndingAt,character-method}
\alias{neditEndingAt,XString-method}
\alias{neditEndingAt,XStringSet-method}

\alias{isMatchingStartingAt}
\alias{isMatchingEndingAt}
\alias{isMatchingAt}
\alias{isMatchingStartingAt,character-method}
\alias{isMatchingStartingAt,XString-method}
\alias{isMatchingStartingAt,XStringSet-method}
\alias{isMatchingEndingAt,character-method}
\alias{isMatchingEndingAt,XString-method}
\alias{isMatchingEndingAt,XStringSet-method}

\alias{which.isMatchingStartingAt}
\alias{which.isMatchingEndingAt}
\alias{which.isMatchingAt}
\alias{which.isMatchingStartingAt,character-method}
\alias{which.isMatchingStartingAt,XString-method}
\alias{which.isMatchingStartingAt,XStringSet-method}
\alias{which.isMatchingEndingAt,character-method}
\alias{which.isMatchingEndingAt,XString-method}
\alias{which.isMatchingEndingAt,XStringSet-method}

\alias{hasLetterAt}

% Old stuff:
\alias{nmismatchStartingAt}
\alias{nmismatchEndingAt}
\alias{isMatching}


\title{Low-level matching functions}

\description{
  In this man page we define precisely and illustrate what a "match" of a
  pattern P in a subject S is in the context of the Biostrings package.
  This definition of a "match" is central to most pattern matching functions
  available in this package: unless specified otherwise, most of them will
  adhere to the definition provided here.

  \code{hasLetterAt} checks whether a sequence or set of sequences has the
  specified letters at the specified positions.
  
  \code{neditAt}, \code{isMatchingAt} and \code{which.isMatchingAt} are
  low-level matching functions that only look for matches at the specified
  positions in the subject.
}

\usage{
hasLetterAt(x, letter, at, fixed=TRUE)

## neditAt() and related utils:
neditAt(pattern, subject, at=1,
        with.indels=FALSE, fixed=TRUE)
neditStartingAt(pattern, subject, starting.at=1,
        with.indels=FALSE, fixed=TRUE)
neditEndingAt(pattern, subject, ending.at=1,
        with.indels=FALSE, fixed=TRUE)

## isMatchingAt() and related utils:
isMatchingAt(pattern, subject, at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
isMatchingStartingAt(pattern, subject, starting.at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
isMatchingEndingAt(pattern, subject, ending.at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)

## which.isMatchingAt() and related utils:
which.isMatchingAt(pattern, subject, at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
        follow.index=FALSE, auto.reduce.pattern=FALSE)
which.isMatchingStartingAt(pattern, subject, starting.at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
        follow.index=FALSE, auto.reduce.pattern=FALSE)
which.isMatchingEndingAt(pattern, subject, ending.at=1,
        max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
        follow.index=FALSE, auto.reduce.pattern=FALSE)
}

\arguments{
  \item{x}{
    A character vector, or an \link{XString} or \link{XStringSet} object.
  }
  \item{letter}{
    A character string or an \link{XString} object containing the letters
    to check.
  }
  \item{at, starting.at, ending.at}{
    An integer vector specifying the starting (for \code{starting.at}
    and \code{at}) or ending (for \code{ending.at}) positions of the
    pattern relatively to the subject.
    With \code{auto.reduce.pattern} (below), either a single integer or
    a constant vector of length \code{nchar(pattern)} (below), to which
    the former is immediately converted.

    For the \code{hasLetterAt} function, \code{letter} and \code{at}
    must have the same length.
  }
  \item{pattern}{
    The pattern string (but see \code{auto.reduce.pattern}, below).
  }
  \item{subject}{
    A character vector, or an \link{XString} or \link{XStringSet} object
    containing the subject sequence(s).
  }
  \item{max.mismatch, min.mismatch}{
    Integer vectors of length >= 1 recycled to the length of the
    \code{at} (or \code{starting.at}, or \code{ending.at}) argument.
    More details below.
  }
  \item{with.indels}{
    See details below.
  }
  \item{fixed}{
    Only with a \link{DNAString} or \link{RNAString}-based subject can a
    \code{fixed} value other than the default (\code{TRUE}) be used.

    If \code{TRUE} (the default), an IUPAC ambiguity code in the pattern
    can only match the same code in the subject, and vice versa.
    If \code{FALSE}, an IUPAC ambiguity code in the pattern can match
    any letter in the subject that is associated with the code, and
    vice versa.
    See \code{\link{IUPAC_CODE_MAP}} for more information about the
    IUPAC Extended Genetic Alphabet.

    \code{fixed} can also be a character vector, a subset
    of \code{c("pattern", "subject")}.
    \code{fixed=c("pattern", "subject")} is equivalent to \code{fixed=TRUE}
    (the default).
    An empty vector is equivalent to \code{fixed=FALSE}.
    With \code{fixed="subject"}, ambiguities in the pattern only
    are interpreted as wildcards.
    With \code{fixed="pattern"}, ambiguities in the subject only
    are interpreted as wildcards.
  }
  \item{follow.index}{
    Whether the single integer returned by \code{which.isMatchingAt}
    (and related utils) should be the first *value* in \code{at} for
    which a match occurred, or its *index* in \code{at} (the default).
  }
  \item{auto.reduce.pattern}{
    Whether \code{pattern} should be effectively shortened by 1 letter,
    from its beginning for \code{which.isMatchingStartingAt} and from
    its end for \code{which.isMatchingEndingAt}, for each successive
    \code{(at, max.mismatch)} "pair".
  }
}

\details{
  A "match" of pattern P in subject S is a substring S' of S that is considered
  similar enough to P according to some distance (or metric) specified by the
  user. 2 distances are supported by most pattern matching functions in the
  Biostrings package. The first (and simplest) one is the "number of mismatching
  letters". It is defined only when the 2 strings to compare have the same
  length, so when this distance is used, only matches that have the same number
  of letters as P are considered.
  The second one is the "edit distance" (aka Levenshtein distance): it's
  the minimum number of operations needed to transform P into S', where an
  operation is an insertion, deletion, or substitution of a single letter.
  When this metric is used, matches can have a different number of letters
  than P.

  The \code{neditAt} function implements these 2 distances.
  If \code{with.indels} is \code{FALSE} (the default), then the first distance
  is used i.e. \code{neditAt} returns the "number of mismatching letters"
  between the pattern P and the substring S' of S starting at the
  positions specified in \code{at} (note that \code{neditAt} is vectorized
  so a long vector of integers can be passed thru the \code{at} argument).
  If \code{with.indels} is \code{TRUE}, then the "edit distance" is
  used: for each position specified in \code{at}, P is compared to
  all the substrings S' of S starting at this position and the smallest
  distance is returned. Note that this distance is guaranteed to be reached
  for a substring of length < 2*length(P) so, of course, in practice,
  P only needs to be compared to a small number of substrings for every
  starting position.
}

\value{
  \code{hasLetterAt}: A logical matrix with one row per element in \code{x}
  and one column per letter/position to check. When a specified position
  is invalid with respect to an element in \code{x} then the corresponding
  matrix element is set to NA.

  \code{neditAt}: If \code{subject} is an \link{XString} object, then
  return an integer vector of the same length as \code{at}.
  If \code{subject} is an \link{XStringSet} object, then return the
  integer matrix with \code{length(at)} rows and \code{length(subject)}
  columns defined by:
  \preformatted{
    sapply(unname(subject),
           function(x) neditAt(pattern, x, ...))
  }

  \code{neditStartingAt} is identical to \code{neditAt} except
  that the \code{at} argument is now called \code{starting.at}.
  \code{neditEndingAt} is similar to \code{neditAt} except that
  the \code{at} argument is now called \code{ending.at} and must contain
  the ending positions of the pattern relatively to the subject.

  \code{isMatchingAt}: If \code{subject} is an \link{XString} object,
  then return the logical vector defined by:
  \preformatted{
    min.mismatch <= neditAt(...) <= max.mismatch
  }
  If \code{subject} is an \link{XStringSet} object, then return the
  logical matrix with \code{length(at)} rows and \code{length(subject)}
  columns defined by:
  \preformatted{
    sapply(unname(subject),
           function(x) isMatchingAt(pattern, x, ...))
  }

  \code{isMatchingStartingAt} is identical to \code{isMatchingAt} except
  that the \code{at} argument is now called \code{starting.at}.
  \code{isMatchingEndingAt} is similar to \code{isMatchingAt} except that
  the \code{at} argument is now called \code{ending.at} and must contain
  the ending positions of the pattern relatively to the subject.

  \code{which.isMatchingAt}: The default behavior (\code{follow.index=FALSE})
  is as follow. If \code{subject} is an \link{XString} object,
  then return the single integer defined by:
  \preformatted{
    which(isMatchingAt(...))[1]
  }
  If \code{subject} is an \link{XStringSet} object, then return
  the integer vector defined by:
  \preformatted{
    sapply(unname(subject),
           function(x) which.isMatchingAt(pattern, x, ...))
  }
  If \code{follow.index=TRUE}, then the returned value is defined by:
  \preformatted{
    at[which.isMatchingAt(..., follow.index=FALSE)]
  }

  \code{which.isMatchingStartingAt} is identical to \code{which.isMatchingAt}
  except that the \code{at} argument is now called \code{starting.at}.
  \code{which.isMatchingEndingAt} is similar to \code{which.isMatchingAt}
  except that the \code{at} argument is now called \code{ending.at} and must
  contain the ending positions of the pattern relatively to the subject.
}

\seealso{
  \code{\link{nucleotideFrequencyAt}},
  \code{\link{matchPattern}},
  \code{\link{matchPDict}},
  \code{\link{matchLRPatterns}},
  \code{\link{trimLRPatterns}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{XString-class},
  \link{align-utils}
}

\examples{
  ## ---------------------------------------------------------------------
  ## hasLetterAt()
  ## ---------------------------------------------------------------------
  x <- DNAStringSet(c("AAACGT", "AACGT", "ACGT", "TAGGA"))
  hasLetterAt(x, "AAAAAA", 1:6)

  ## hasLetterAt() can be used to answer questions like: "which elements
  ## in 'x' have an A at position 2 and a G at position 4?"
  q1 <- hasLetterAt(x, "AG", c(2, 4))
  which(rowSums(q1) == 2)

  ## or "how many probes in the drosophila2 chip have T, G, T, A at
  ## position 2, 4, 13 and 20, respectively?"
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe)
  q2 <- hasLetterAt(probes, "TGTA", c(2, 4, 13, 20))
  sum(rowSums(q2) == 4)
  ## or "what's the probability to have an A at position 25 if there is
  ## one at position 13?"
  q3 <- hasLetterAt(probes, "AACGT", c(13, 25, 25, 25, 25))
  sum(q3[ , 1] & q3[ , 2]) / sum(q3[ , 1])
  ## Probabilities to have other bases at position 25 if there is an A
  ## at position 13:
  sum(q3[ , 1] & q3[ , 3]) / sum(q3[ , 1])  # C
  sum(q3[ , 1] & q3[ , 4]) / sum(q3[ , 1])  # G
  sum(q3[ , 1] & q3[ , 5]) / sum(q3[ , 1])  # T

  ## See ?nucleotideFrequencyAt for another way to get those results.

  ## ---------------------------------------------------------------------
  ## neditAt() / isMatchingAt() / which.isMatchingAt()
  ## ---------------------------------------------------------------------
  subject <- DNAString("GTATA")

  ## Pattern "AT" matches subject "GTATA" at position 3 (exact match)
  neditAt("AT", subject, at=3)
  isMatchingAt("AT", subject, at=3)

  ## ... but not at position 1
  neditAt("AT", subject)
  isMatchingAt("AT", subject)

  ## ... unless we allow 1 mismatching letter (inexact match)
  isMatchingAt("AT", subject, max.mismatch=1)

  ## Here we look at 6 different starting positions and find 3 matches if
  ## we allow 1 mismatching letter
  isMatchingAt("AT", subject, at=0:5, max.mismatch=1)

  ## No match
  neditAt("NT", subject, at=1:4)
  isMatchingAt("NT", subject, at=1:4)

  ## 2 matches if N is interpreted as an ambiguity (fixed=FALSE)
  neditAt("NT", subject, at=1:4, fixed=FALSE)
  isMatchingAt("NT", subject, at=1:4, fixed=FALSE)

  ## max.mismatch != 0 and fixed=FALSE can be used together
  neditAt("NCA", subject, at=0:5, fixed=FALSE)
  isMatchingAt("NCA", subject, at=0:5, max.mismatch=1, fixed=FALSE)

  some_starts <- c(10:-10, NA, 6)
  subject <- DNAString("ACGTGCA")
  is_matching <- isMatchingAt("CAT", subject, at=some_starts, max.mismatch=1)
  some_starts[is_matching]

  which.isMatchingAt("CAT", subject, at=some_starts, max.mismatch=1)
  which.isMatchingAt("CAT", subject, at=some_starts, max.mismatch=1,
                     follow.index=TRUE)

  ## ---------------------------------------------------------------------
  ## WITH INDELS
  ## ---------------------------------------------------------------------
  subject <- BString("ABCDEFxxxCDEFxxxABBCDE")

  neditAt("ABCDEF", subject, at=9)
  neditAt("ABCDEF", subject, at=9, with.indels=TRUE)
  isMatchingAt("ABCDEF", subject, at=9, max.mismatch=1, with.indels=TRUE)
  isMatchingAt("ABCDEF", subject, at=9, max.mismatch=2, with.indels=TRUE)
  neditAt("ABCDEF", subject, at=17)
  neditAt("ABCDEF", subject, at=17, with.indels=TRUE)
  neditEndingAt("ABCDEF", subject, ending.at=22)
  neditEndingAt("ABCDEF", subject, ending.at=22, with.indels=TRUE)
}

\keyword{methods}
