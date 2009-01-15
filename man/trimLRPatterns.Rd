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
  trimLRPatterns(Lpattern = NULL, Rpattern = NULL, subject,
                 max.Lmismatch = 0, max.Rmismatch = 0,
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
  trimLRPatterns(Lpattern = Lpattern, subject = subject)
  trimLRPatterns(Rpattern = Rpattern, subject = subject, max.Rmismatch = 2)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 2, max.Rmismatch = 2)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subject,
                 max.Lmismatch = 2, max.Rmismatch = 2, ranges = TRUE)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet)
  trimLRPatterns(Lpattern = Lpattern, Rpattern = Rpattern, subject = subjectSet,
                 ranges = TRUE)
}

\keyword{methods}
