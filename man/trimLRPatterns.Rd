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
                 max.Lmismatch = rep.int(c(-1, 0), c(nchar(Lpattern) - 1, 1)),
                 max.Rmismatch = rep.int(c(-1, 0), c(nchar(Rpattern) - 1, 1)),
                 with.Lindels = FALSE, with.Rindels = FALSE,
                 Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
}

\arguments{
  \item{Lpattern}{
    The left part of the pattern. The \code{Lpattern} adds no restrictions on
    the \code{subject} argument.
  }
  \item{Rpattern}{
    The right part of the pattern. When \code{subject} is an \link{XStringSet},
    the the strings must be of the same length to use \code{Rpattern}.
  }
  \item{subject}{
    An \link{XString} or \link{XStringSet} object containing the target
    sequence(s).
  }
  \item{max.Lmismatch}{
    A non-decreasing integer vector of length \code{nLp = nchar(Lpattern)}
    whose elements \code{max.Lmismatch[i]} represent the maximum number of
    acceptable mismatching letters when aligning
    \code{substring(Lpattern, nLp - i + 1, nLp)} with
    \code{substring(subject, 1, i)}.
    Negative numbers are used to prevent trimming at the i-th location.
    If \code{length(max.Lmismatch) < nLp}, then \code{max.Lmismatch} will be
    augmented with enough -1's at the beginning of the vector to bring it up to
    length \code{nLp}.

    If non-zero, an inexact matching algorithm is used (see the
    \code{\link{matchPattern}} function for more information).
  }
  \item{max.Rmismatch}{
    A non-decreasing integer vector of length \code{nRp = nchar(Rpattern)}
    whose elements \code{max.Rmismatch[i]} represent the maximum number of
    acceptable mismatching letters when aligning
    \code{substring(Rpattern, 1, i)} with
    \code{substring(subject, nchar(subject) - i + 1, nchar(subject))}.
    Negative numbers are used to prevent trimming at the i-th location.
    If \code{length(max.Rmismatch) < nRp}, then \code{max.Rmismatch} will be
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

  ## Allow for mismatches on the flanks
  maxMismatches <- rep(-1:2, c(2,1,3,3))
  maxMismatches
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches)

  ## Produce ranges that can be an input into other functions
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches,
                 ranges = TRUE)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 ranges = TRUE)
}

\keyword{methods}
