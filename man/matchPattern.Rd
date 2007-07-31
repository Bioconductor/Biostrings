\name{matchPattern}

\alias{matchPattern}
\alias{matchPattern,character-method}
\alias{matchPattern,BString-method}
\alias{matchPattern,BStringViews-method}
\alias{countPattern}
\alias{countPattern,character-method}
\alias{countPattern,BString-method}
\alias{countPattern,BStringViews-method}
\alias{mismatch}
\alias{mismatch,BStringViews-method}

\alias{matchDNAPattern}


\title{String searching functions}

\description{
  Generic that finds all matches of a pattern in a BString.
}

\usage{
  matchPattern(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
  countPattern(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
  mismatch(pattern, x, fixed=TRUE)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    A \link{BString} (or \link{DNAString} or \link{RNAString}) object
    containing the subject string.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-fuzzy"}, \code{"boyer-moore"} or \code{"shift-or"}.
  }
  \item{mismatch}{
    The number of mismatches allowed. If non-zero a fuzzy string searching
    algorithm is used for matching.
  }
  \item{fixed}{
    Only with a \link{DNAString} or \link{RNAString} subject can a \code{fixed}
    value other than the default (\code{TRUE}) be used.

    With \code{fixed=FALSE}, ambiguities (i.e. letters from the IUPAC Extended
    Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}}) that are not from the
    base alphabet) in the pattern \_and\_ in the subject are interpreted as
    wildcards i.e. they match any letter that they stand for.

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
  \item{x}{
    A \link{BStringViews} object (typically, one returned
    by \code{matchPattern(pattern, subject)}).
  }
}

\details{
  Available algorithms are: ``naive exact'', ``naive fuzzy'',
  ``Boyer-Moore-like'' and ``shift-or''. Not all of them can be
  used in all situations: restrictions depend on the length of
  the pattern, the class of the subject and the values of
  \code{mismatch} and \code{fixed}.

  When 2 different algorithms can be used for a given task,
  then choosing one or the other only affects the performance,
  not the result, so there is no "wrong choice" (strictly speaking).
  In short, it is better to just use \code{algorithm="auto"} (the default):
  this way \code{matchPattern} will choose the algo that is best suited
  for the task.
}

\value{
  A \link{BStringViews} object for \code{matchPattern}.

  A single integer for \code{countPattern}.

  A list of integer vectors for \code{mismatch}.
}

\seealso{
  \code{\link{mask}},
  \code{\link{alphabetFrequency}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## A simple fuzzy matching example with a short subject
  x <- DNAString("AAGCGCGATATG")
  m1 <- matchPattern("GCNNNAT", x)
  m1
  m2 <- matchPattern("GCNNNAT", x, fixed=FALSE)
  m2
  as.matrix(m2)

  ## With DNA sequence of yeast chromosome number 1
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)
  PpiI <- "GAACNNNNNCTC" # a restriction enzyme pattern
  match1.PpiI <- matchPattern(PpiI, yeast1, fixed=FALSE)
  match2.PpiI <- matchPattern(PpiI, yeast1, mismatch=1, fixed=FALSE)

  ## With a genome containing isolated Ns
  library(BSgenome.Celegans.UCSC.ce2)
  chrII <- Celegans[["chrII"]]
  alphabetFrequency(chrII)
  matchPattern("N", chrII)
  matchPattern("TGGGTGTCTTT", chrII) # no match
  matchPattern("TGGGTGTCTTT", chrII, fixed=FALSE) # 1 match

  ## Using wildcards ("N") in the pattern on a genome containing N-blocks
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrX <- Dmelanogaster[["X"]]
  noN_chrX <- mask(chrX, "N")
  mask(noN_chrX) # See the N-blocks?
  matchPattern("TTTATGNTTGGTA", noN_chrX, fixed=FALSE)
  ## Can also be achieved with
  matchPattern("TTTATGNTTGGTA", chrX, fixed="subject")
}

\keyword{methods}

