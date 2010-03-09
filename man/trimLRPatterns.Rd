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
    Either an integer vector of length \code{nLp = nchar(Lpattern)} whose
    elements \code{max.Lmismatch[i]} represent the maximum number of
    acceptable mismatching letters (or edit distance if \code{with.Lindels}
    is \code{TRUE}) when aligning \code{substring(Lpattern, nLp - i + 1, nLp)}
    with \code{substring(subject, 1, i)} or a single numeric value in
    \code{[0, 1)} that represents a constant maximum mismatch rate for each
    of the \code{nL} alignments. Negative numbers in integer vector inputs
    are used to prevent trimming at the \code{i}-th location. If an integer
    vector input has \code{length(max.Lmismatch) < nLp}, then
    \code{max.Lmismatch} will be augmented with enough \code{-1}'s at the
    beginning of the vector to bring it up to length \code{nLp}.

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
