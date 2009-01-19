\name{trimLRPatterns}

\alias{trimLRPatterns}
\alias{trimLRPatterns,XString-method}
\alias{trimLRPatterns,XStringSet-method}


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
    The left part of the pattern.
  }
  \item{Rpattern}{
    The right part of the pattern.
  }
  \item{subject}{
    An \link{XString} or \link{XStringSet} object containing the target
    sequence(s).
  }
  \item{max.Lmismatch}{
    Either an integer vector of length \code{nLp = nchar(Lpattern)} whose
    elements \code{max.Lmismatch[i]} represent the maximum number of acceptable
    mismatching letters when aligning
    \code{substring(Lpattern, nLp - i + 1, nLp)} with
    \code{substring(subject, 1, i)} or a single numeric value in \code{(0, 1)}
    that represents a constant maximum mismatch rate for each of the \code{nL}
    alignments. Negative numbers in integer vector inputs are used to prevent
    trimming at the i-th location. If an integer vector input has
    \code{length(max.Lmismatch) < nLp}, then \code{max.Lmismatch} will be
    augmented with enough -1's at the beginning of the vector to bring it up to
    length \code{nLp}.

    If non-zero, an inexact matching algorithm is used (see the
    \code{\link{matchPattern}} function for more information).
  }
  \item{max.Rmismatch}{
    Either an integer vector of length \code{nRp = nchar(Rpattern)} whose
    elements \code{max.Rmismatch[i]} represent the maximum number of acceptable
    mismatching letters when aligning
    \code{substring(Rpattern, nRp - i + 1, nRp)} with
    \code{substring(subject, 1, i)} or a single numeric value in \code{(0, 1)}
    that represents a constant maximum mismatch rate for each of the \code{nR}
    alignments. Negative numbers in integer vector inputs are used to prevent
    trimming at the i-th location. If an integer vector input has
    \code{length(max.Rmismatch) < nRp}, then \code{max.Rmismatch} will be
    augmented with enough -1's at the beginning of the vector to bring it up to
    length \code{nRp}.

    If non-zero, an inexact matching algorithm is used (see the
    \code{\link{matchPattern}} function for more information).
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

    With \code{Lfixed=FALSE}, ambiguities (i.e. letters from the IUPAC Extended
    Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}}) that are not from the
    base alphabet) in the left pattern \_and\_ in the subject are interpreted
    as wildcards i.e. they match any letter that they stand for.

    See the \code{fixed} argument of the \code{\link{matchPattern}} function
    for more information.
  }
  \item{Rfixed}{
    Same as \code{Lfixed} but for the right part of the pattern.
  }
  \item{ranges}{
    If \code{TRUE}, then return the ranges to use to trim \code{subject}.
    If \code{FALSE}, then returned the trimmed \code{subject}.
  }
}

\value{
  A new \link{XString} or \link{XStringSet} object with the flanking patterns
  within the specified edit distances removed.
}

\author{P. Aboyoun}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \link{match-utils},
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
                 max.Lmismatch = rep(0, 9), max.Rmismatch = rep(0, 9))

  ## Allow for mismatches on the flanks
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 0.2, max.Rmismatch = 0.2)
  maxMismatches <- as.integer(0.2 * 1:9)
  maxMismatches
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches)

  ## Produce ranges that can be an input into other functions
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 max.Lmismatch = rep(0, 9), max.Rmismatch = rep(0, 9),
                 ranges = TRUE)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 0.2, max.Rmismatch = 0.2, ranges = TRUE)
}

\keyword{methods}
