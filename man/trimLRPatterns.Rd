\name{trimLRPatterns}

\alias{trimLRPatterns}
\alias{trimLRPatterns,XString-method}
\alias{trimLRPatterns,XStringSet-method}
\alias{trimLRPatterns,character-method}


\title{Trim Flanking Patterns from Sequences}

\description{
  The \code{trimLRPatterns} function trims left and/or right flanking patterns
  from sequences.
}

\usage{
  trimLRPatterns(Lpattern = "", Rpattern = "", subject,
                 max.Lmismatch = 0, max.Rmismatch = 0,
                 with.Lindels = FALSE, with.Rindels = FALSE,
                 Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
}

\arguments{
  \item{Lpattern}{
    The left pattern.
  }
  \item{Rpattern}{
    The right pattern.
  }
  \item{subject}{
    An \link{XString} object, \link{XStringSet} object, or character vector
    containing the target sequence(s).
  }
  \item{max.Lmismatch}{
    Either an integer vector of length \code{nLp = nchar(Lpattern)}
    representing an absolute number of mismatches (or edit distance if
    \code{with.Lindels} is \code{TRUE}) or a single numeric value in the
    interval \code{[0, 1)} representing a mismatch rate when aligning
    terminal substrings (suffixes) of \code{Lpattern} with the beginning
    (prefix) of \code{subject} following the conventions set by 
    \code{\link{neditStartingAt}}, \code{\link{isMatchingStartingAt}},
    etc.

    When \code{max.Lmismatch} is \code{0L} or a numeric value in the
    interval \code{[0, 1)}, it is taken as a "rate" and is converted to
    \code{as.integer(1:nLp * max.Lmismatch)}, analogous to \link{agrep}
    (which, however, employs \code{\link{ceiling}}).

    Otherwise, \code{max.Lmismatch} is treated as an integer vector where
    negative numbers are used to prevent trimming at the \code{i}-th
    location. When an input integer vector is shorter than \code{nLp}, it
    is augmented with enough \code{-1}s at the beginning to bring its length
    up to \code{nLp}. Elements of \code{max.Lmismatch} beyond the first
    \code{nLp} are ignored.

    Once the integer vector is constructed using the rules given above, when
    \code{with.Lindels} is \code{FALSE}, \code{max.Lmismatch[i]} is the number
    of acceptable mismatches (errors) between the suffix
    \code{substring(Lpattern, nLp - i + 1, nLp)} of \code{Lpattern} and the
    first \code{i} letters of \code{subject}.
    When \code{with.Lindels} is \code{TRUE}, \code{max.Lmismatch[i]}
    represents the allowed "edit distance" between that suffix of
    \code{Lpattern} and \code{subject}, starting at position \code{1} of
    \code{subject} (as in \code{\link{matchPattern}} and
    \code{\link{isMatchingStartingAt}}).

    For a given element \code{s} of the \code{subject}, the initial segment
    (prefix) \code{substring(s, 1, j)} of \code{s} is trimmed if \code{j} is
    the largest \code{i} for which there is an acceptable match, if any.
  }
  \item{max.Rmismatch}{
    Same as \code{max.Lmismatch} but with \code{Rpattern}, along with
    \code{with.Rindels} (below), and its initial segments (prefixes)
    \code{substring(Rpattern, 1, i)}.

    For a given element \code{s} of the subject, with \code{nS = nchar(s)},
    the terminal segment (suffix) \code{substring(s, nS - j + 1, nS)} of
    \code{s} is trimmed if \code{j} is the largest \code{i} for which there
    is an acceptable match, if any.
  }
  \item{with.Lindels}{
    If \code{TRUE}, indels are allowed in the alignments of the suffixes
    of \code{Lpattern} with the subject, at its beginning.
    See the \code{with.indels} arguments of the \code{\link{matchPattern}}
    and \code{\link{neditStartingAt}} functions for detailed information.
  }
  \item{with.Rindels}{
    Same as \code{with.Lindels} but for alignments of the prefixes of
    \code{Rpattern} with the subject, at its end.
    See the \code{with.indels} arguments of the \code{\link{matchPattern}}
    and \code{\link{neditEndingAt}} functions for detailed information.
  }
  \item{Lfixed, Rfixed}{
    Whether IUPAC extended letters in the left or right pattern should
    be interpreted as ambiguities (see \code{?`\link{lowlevel-matching}`}
    for the details).
  }
  \item{ranges}{
    If \code{TRUE}, then return the ranges to use to trim \code{subject}.
    If \code{FALSE}, then returned the trimmed \code{subject}.
  }
}

\value{
  A new \link{XString} object, \link{XStringSet} object, or character vector
  with the "longest" flanking matches removed, as described above.
}

\author{P. Aboyoun and H. Jaffee}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \link{lowlevel-matching},
  \link{XString-class},
  \link{XStringSet-class}
}

\examples{
  Lpattern <- "TTCTGCTTG"
  Rpattern <- "GATCGGAAG"
  subject <- DNAString("TTCTGCTTGACGTGATCGGA")
  subjectSet <- DNAStringSet(c("TGCTTGACGGCAGATCGG", "TTCTGCTTGGATCGGAAG"))

  ## Only allow for perfect matches on the flanks
  trimLRPatterns(Lpattern = Lpattern, subject = subject)
  trimLRPatterns(Rpattern = Rpattern, subject = subject)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet)

  ## Allow for perfect matches on the flanking overlaps
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 max.Lmismatch = 0, max.Rmismatch = 0)

  ## Allow for mismatches on the flanks
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 0.2, max.Rmismatch = 0.2)
  maxMismatches <- as.integer(0.2 * 1:9)
  maxMismatches
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches)

  ## Produce ranges that can be an input into other functions
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 max.Lmismatch = 0, max.Rmismatch = 0, ranges = TRUE)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 0.2, max.Rmismatch = 0.2, ranges = TRUE)
}

\keyword{methods}
