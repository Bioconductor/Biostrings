\name{PDict-class}
\docType{class}

% Dups class:
\alias{class:Dups}
\alias{Dups-class}
\alias{Dups}

\alias{length,Dups-method}
\alias{duplicated,Dups-method}
\alias{show,Dups-method}

% PDict class:
\alias{class:PDict}
\alias{PDict-class}

\alias{length,PDict-method}
\alias{width,PDict-method}
\alias{names,PDict-method}
\alias{head,PDict-method}
\alias{tb}
\alias{tb,PDict-method}
\alias{tb.width}
\alias{tb.width,PDict-method}
\alias{tail,PDict-method}

\alias{[[,PDict-method}
\alias{[[<-,PDict-method}

\alias{initialize,PDict-method}
\alias{show,PDict-method}
\alias{duplicated,PDict-method}
\alias{patternFrequency}
\alias{patternFrequency,PDict-method}

% ACtree class:
\alias{class:ACtree}
\alias{ACtree-class}
\alias{ACtree}

\alias{length,ACtree-method}
\alias{show,ACtree-method}
\alias{[,ACtree-method}
\alias{as.matrix,ACtree-method}
\alias{initialize,ACtree-method}

% ACtree_PDict class:
\alias{class:ACtree_PDict}
\alias{ACtree_PDict-class}
\alias{ACtree_PDict}

\alias{initialize,ACtree_PDict-method}

% User level constructor:
\alias{PDict}
\alias{PDict,character-method}
\alias{PDict,DNAStringSet-method}
\alias{PDict,XStringViews-method}
\alias{PDict,AsIs-method}


\title{PDict objects}

\description{
  The PDict class is a container for storing a preprocessed dictionary of DNA
  patterns that can later be passed to the \code{\link{matchPDict}} function
  for fast matching.

  \code{PDict} is the constructor function for creating new PDict objects.
}

\usage{
  PDict(x, tb.start=NA, tb.end=NA, tb.width=NA, type="ACtree", skip.invalid.patterns=FALSE)
}

\arguments{
  \item{x}{
    A character vector, a \link{DNAStringSet} object or
    an \link{XStringViews} object with a \link{DNAString} subject.
  }
  \item{tb.start}{
    [DOCUMENT ME]
  }
  \item{tb.end}{
    [DOCUMENT ME]
  }
  \item{tb.width}{
    [DOCUMENT ME]
  }
  \item{type}{
    [DOCUMENT ME]
  }
  \item{skip.invalid.patterns}{
    [DOCUMENT ME]
  }
}

\details{
  THIS IS STILL WORK IN PROGRESS!

  If the original dictionary \code{x} is a character vector or
  an \link{XStringViews} object with a \link{DNAString} subject,
  then the \code{PDict} constructor will first try to turn it
  into a \link{DNAStringSet} object.

  By default (i.e. if \code{PDict} is called with \code{tb.start=NA},
  \code{tb.end=NA} and \code{tb.width=NA}) the following restrictions
  apply: (1) the original dictionary can only contain As, Cs, Gs and Ts
  (no IUPAC extended letters); (2) all the patterns in the dictionary must
  have the same length ("constant width" dictionary); and (3) later
  \code{matchPdict} can only be used with \code{max.mismatch=0}.

  A Trusted Band can be used in order to relax these restrictions (more
  on this below).

  Two types of preprocessing are currently supported: \code{type="ACtree"}
  and \code{type="Twobit"}. When used with \code{type="ACtree"}, \code{PDict}
  returns an ACtree\_PDict object. When used with \code{type="Twobit"},
  it returns a Twobit\_PDict object. The ACtree\_PDict and Twobit\_PDict
  classes are subclasses of the PDict class.

  Only with ACtree\_PDict objects can \code{matchPdict} be called with
  \code{fixed="pattern"} (instead of \code{fixed=TRUE}, the default)
  so that IUPAC extended letters in the subject are treated as ambiguities.
  Twobit\_PDict objects don't allow this.

  How to specify a Trusted Band?
  The \code{tb.start}, \code{tb.end} and \code{tb.width} arguments can
  be used to specify a Trusted Band i.e. a region in the patterns where
  the restrictions described previously will apply. This region is
  defined by its starting and ending nucleotide position relatively
  to each pattern and must have the same length for all patterns (this
  common length is the width of the Trusted Band).
  You can think of the process of specifying a Trusted Band 
  as drawing 2 vertical lines on the original dictionary
  (note that these 2 lines are not necessarily straight lines but the
  horizontal space between them must be constant).
  So when a Trusted Band is specified, the original dictionary is naturally
  divided in three parts: the head, the Trusted Band, and the tail of the
  dictionary.
  Each part is itself a \link{DNAStringSet} object: the Trusted Band is
  the \link{DNAStringSet} object obtained by extracting the trusted region
  from each input pattern (note that the resulting object has a constant
  width that is >= 1 because Trusted Bands of width 0 are not allowed); the
  head is the \link{DNAStringSet} object obtained by extracting the region
  located before (i.e. to the left of) the trusted region from each input
  pattern; and the tail is the \link{DNAStringSet} object obtained by
  extracting the region located after (i.e. to the right of) the trusted
  region from each input pattern.
  Note that each of these three parts contain the same number of patterns
  than the original dictionary and that the original dictionary can easily
  be reconstructed from those parts.

  Why use a Trusted Band?
  The restrictions described previously will apply only to the Trusted Band
  e.g. the Trusted Band cannot contain IUPAC extended letters but the head
  and the tail can. Also with a Trusted Band, if \code{matchPdict} is called
  with a non-null \code{max.mismatch} value then mismatching letters will
  be allowed in the head and the tail. Or if \code{matchPdict} is called
  with \code{fixed="subject"} then IUPAC extended letters in the head and
  the tail will be treated as ambiguities.

  Finally note that calling \code{PDict} with \code{tb.start=NA},
  \code{tb.end=NA} and \code{tb.width=NA} (the default) is equivalent
  to calling it with \code{tb.start=1}, \code{tb.end=-1} and
  \code{tb.width=NA}, which results in a Trusted Band that covers
  the entire dictionary (no head and no tail).
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a PDict object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of patterns in \code{x}.
    }
    \item{}{
      \code{width(x)}:
      A vector of non-negative integers containing the number
      of letters for each pattern in \code{x}.
    }
    \item{}{
      \code{names(x)}:
      The names of the patterns in \code{x}.
    }
    \item{}{
      \code{head(x)}:
      The head of \code{x} or \code{NULL} if \code{x} has no head.
    }
    \item{}{
      \code{tb(x)}:
      The Trusted Band defined on \code{x}.
    }
    \item{}{
      \code{tb.width(x)}:
      The width of the Trusted Band defined on \code{x}.
      Note that, unlike \code{width(tb(x))}, this is a single integer.
      And because the Trusted Band has a constant width, \code{tb.width(x)}
      is in fact equivalent to \code{unique(width(tb(x)))},
      or to \code{width(tb(x))[1]}.
    }
    \item{}{
      \code{tail(x)}:
      The tail of \code{x} or \code{NULL} if \code{x} has no tail.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is a PDict object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the i-th pattern from \code{x} as a \link{DNAString} object.
    }
  }
}

\section{Other methods}{
  In the code snippet below,
  \code{x} is a PDict object.

  \describe{
    \item{}{
      \code{duplicated(x)}:
      [DOCUMENT ME]
    }
    \item{}{
      \code{patternFrequency(x)}:
      [DOCUMENT ME]
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{matchPDict}},
  \code{\link{DNA_ALPHABET}},
  \link{DNAStringSet-class},
  \link{XStringViews-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. NO HEAD AND NO TAIL (THE DEFAULT)
  ## ---------------------------------------------------------------------
  library(drosophila2probe)
  dict0 <- drosophila2probe$sequence   # The input sequences.
  length(dict0)                        # Hundreds of thousands of patterns.
  unique(nchar(dict0))                 # Patterns are 25-mers.
  dict0[1:5]
  pdict0 <- PDict(dict0)               # Store the input dictionary into a
                                       # PDict object (preprocessing).
  pdict0
  class(pdict0)
  length(pdict0)                       # Same as length(dict0).
  tb.width(pdict0)                     # The width of the (implicit)
                                       # Trusted Band.
  sum(duplicated(pdict0))
  table(patternFrequency(pdict0))      # 9 patterns are repeated 3 times.
  pdict0[[1]]
  pdict0[[5]]

  ## ---------------------------------------------------------------------
  ## B. NO HEAD AND A TAIL
  ## ---------------------------------------------------------------------
  dict1 <- c("ACNG", "GT", "CGT", "AC")
  pdict1 <- PDict(dict1, tb.end=2)
  pdict1
  class(pdict1)
  length(pdict1)
  width(pdict1)
  head(pdict1)
  tb(pdict1)
  tb.width(pdict1)
  width(tb(pdict1))
  tail(pdict1)
  pdict1[[3]]
}

\keyword{methods}
\keyword{classes}
