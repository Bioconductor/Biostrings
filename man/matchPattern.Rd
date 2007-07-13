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
    One of the following: \code{"auto"}, \code{"naive"}, \code{"boyer-moore"},
    \code{"forward-search"} or \code{"shift-or"} (note that
    \code{"forward-search"} is not available yet).
  }
  \item{mismatch}{
    The number of mismatches allowed. If non-zero a fuzzy string searching
    algorithm is used for matching.
  }
  \item{fixed}{
    A logical. Only \link{DNAString} or \link{RNAString} support \code{fixed=FALSE}.
    In this case, ambiguities (i.e. letters from the IUPAC Extended Genetic Alphabet
    (see \code{\link{IUPAC_CODE_MAP}}) that are not from the base alphabet)
    match any letter that they stand for.
  }
  \item{x}{
    A \link{BStringViews} object (typically, one returned
    by \code{matchPattern(pattern, subject)}).
  }
}

\details{
  You can select between a ``naive'', a ``Boyer-Moore-like''
  and the ``shift-or'' algorithms for performing the search.
  Note that, depending on the values of the additional arguments
  \code{mismatch} and \code{fixed}, some algorithms can't be used
  so it's better to let \code{matchPattern} choose the right
  algorithm for you.
  When two different algorithms can be used for the same task,
  then choosing one or the other only affect the performance,
  not the semantic. Again it's better to let \code{matchPattern}
  choose because it tries to select the faster for the task.
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
  \code{\link{BStringViews}},
  \code{\link{DNAString}}
}

\examples{
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
}

\keyword{methods}

