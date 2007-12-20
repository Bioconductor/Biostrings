\name{matchPDict}

% Classes
\alias{class:PDict}
\alias{PDict-class}
\alias{class:ULdna_PDict}
\alias{ULdna_PDict-class}

% Generics and methods:
\alias{initialize,ULdna_PDict-method}
\alias{length,ULdna_PDict-method}
\alias{width,ULdna_PDict-method}


\alias{matchPDict}
\alias{matchPDict,BString-method}


\title{Searching a sequence for patterns stored in a dictionary}

\description{
  \code{matchPDict} efficiently finds all occurences in a text (the
  subject) of any pattern from a set of patterns (the dictionary).
  The implementation of \code{matchPDict} is based on the Aho-Corasick
  algorithm and is very fast. However, for now, \code{matchPDict} only
  works with a dictionary of DNA patterns where all the patterns have
  the same number of nucleotides (uniform-length dictionary).
  Also, only exact matching is currently supported.
}

\usage{
  matchPDict(pdict, subject, algorithm="auto", mismatch=0, fixed=TRUE)
}

\arguments{
  \item{pdict}{
    The preprocessed dictionary (see below for the details).
  }
  \item{subject}{
    A \link{BString} (or derived) object containing the subject string.
    For now, only \link{DNAString} subjects are supported.
  }
  \item{algorithm}{
    Not supported yet.
  }
  \item{mismatch}{
    Not supported yet.
  }
  \item{fixed}{
    Not supported yet.
  }
}

\value{
  A list of the same length as the input directory where the i-th
  element is an integer vector containing the \_ending\_ positions
  of every match found for the i-th pattern in the directory.
}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{alphabetFrequency}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## Creating the pattern dictionary
  library(drosophila2probe)
  dict <- drosophila2probe$sequence     # The input dictionary.
  length(dict)                          # Hundreds of thousands of patterns.
  unique(nchar(dict))                   # Patterns are 25-mers.
  dict[1:5]
  pdict <- new("ULdna_PDict", dict)     # Store the input dictionary into a
                                        # ULdna_PDict instance (preprocessing).
  length(pdict)                         # Same as the input dictionary.
  width(pdict)                          # The number of chars per pattern.

  ## Using the pattern dictionary on chromosome 3R
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chr3R <- Dmelanogaster[["3R"]]        # Load chromosome 3R
  chr3R
  mends <- matchPDict(pdict, chr3R)     # Search...

  ## Looking at the results
  length(mends)                         # Same as the input dictionary.
  nmatches <- sapply(mends, length)
  sum(nmatches)                         # Total number of matches.
  table(nmatches)
  i0 <- which(nmatches == max(nmatches))
  dict[i0]                              # The pattern with most occurences.
  views(chr3R, mends[[i0]] - 24, mends[[i0]]) # All the occurences.
}

\keyword{methods}

