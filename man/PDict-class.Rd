\name{PDict-class}
\docType{class}

% Classes:
\alias{class:PDict}
\alias{PDict-class}
\alias{class:ULdna_PDict}
\alias{ULdna_PDict-class}
\alias{ULdna_PDict}
\alias{class:TBdna_PDict}
\alias{TBdna_PDict-class}
\alias{TBdna_PDict}

% Low level initialization and accesor methods:
\alias{initialize,ULdna_PDict-method}
\alias{length,ULdna_PDict-method}
\alias{width,ULdna_PDict-method}
\alias{names,ULdna_PDict-method}

\alias{initialize,TBdna_PDict-method}
\alias{head,TBdna_PDict-method}
\alias{tail,TBdna_PDict-method}

% User level constructor:
\alias{PDict}
\alias{PDict,character-method}
\alias{PDict,DNAStringList-method}
\alias{PDict,BStringViews-method}
\alias{PDict,AsIs}

% Other generics and methods:
\alias{show,ULdna_PDict-method}
\alias{show,TBdna_PDict-method}
\alias{duplicated,ULdna_PDict-method}
\alias{duplicated,TBdna_PDict-method}
\alias{patternFrequency}
\alias{patternFrequency,ULdna_PDict-method}
\alias{patternFrequency,TBdna_PDict-method}


\title{PDict objects}

\description{
  The PDict class is a container for storing a preprocessed set of patterns
  (aka the dictionary) that can be used with the \code{\link{matchPDict}}
  function in order to find the occurences of all the input patterns (i.e.
  the patterns stored in the PDict object) in a text in an efficient way.

  Converting a set of input sequences into a PDict object is called
  preprocessing. This operation is done by using the \code{PDict}
  constructor.
}

\usage{
  PDict(dict, tb.start=1, tb.end=NA, tb.width=NA, drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
}

\arguments{
  \item{dict}{
      [DOCUMENT ME]
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
  \item{drop.head}{
      [DOCUMENT ME]
  }
  \item{drop.tail}{
      [DOCUMENT ME]
  }
  \item{skip.invalid.patterns}{
      [DOCUMENT ME]
  }
}

\details{
  This is a work in progress and only 2 types of dictionaries are supported
  at the moment: Uniform Length DNA dictionaries and Trusted Band DNA
  dictionaries.

  A Uniform Length DNA dictionary is a dictionary where all the patterns are
  DNA sequences of the same length (i.e. all the patterns have the same number
  of nucleotides). For now the patterns can only contain As, Cs, Gs and
  Ts (no IUPAC extended letters).
  The container for this particular type of dictionary is the ULdna_PDict
  class (a subclass of the PDict class).

  A Trusted Band DNA dictionary is a dictionary where the patterns are DNA
  sequences with a trusted region i.e. a region that will have to match
  exactly when the dictionary is used with the \code{\link{matchPDict}}
  function.
  This trusted region must have the same length for all patterns. It can be a
  prefix (Trusted Prefix), a suffix (Trusted Suffix) or more generally any
  substring (Trusted Band) of the input patterns.
  The container for this particular type of dictionary is the TBdna_PDict
  class (a subclass of the PDict class). The dictionary stored in a
  TBdna_PDict object is splitted in 3 parts: the Trusted Band, the head
  and the tail of the dictionary. Each of them contains the same number
  of sequences (eventually empty) which is the number of input patterns.
  The Trusted Band is the set of sequences obtained by extracting
  the trusted region from each input pattern.
  The head is the set of sequences obtained by extracting the region located
  before (i.e. to the left of) the trusted region from each input pattern.
  The tail is the set of sequences obtained by extracting the region located
  after (i.e. to the right of) the trusted region from each input pattern.

  Like for a Uniform Length DNA dictionary, the Trusted Band can only contain
  As, Cs, Gs and Ts (no IUPAC extended letters) for now.
  However, the head and the tail of a Trusted Band DNA dictionary can contain
  any valid DNA letter including IUPAC extended letters (see
  \code{\link{DNA_ALPHABET}} for the set of valid DNA letters).

  Note that a Trusted Band DNA dictionary with no head and no tail (i.e. with
  a head and a tail where all the sequences are empty) is in fact a Uniform
  Length DNA dictionary.
  A Trusted Band DNA dictionary with no head is called a Trusted Prefix DNA
  dictionary.
  A Trusted Band DNA dictionary with no tail is called a Trusted Suffix DNA
  dictionary.
  Only Trusted Prefix and Trusted Suffix DNA dictionaries are currently
  supported.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a PDict object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of patterns in the PDict object.
    }
    \item{}{
      \code{width(x)}:
      The number of nucleotides per pattern for a Uniform Length DNA dictionary.
      The width of the Trusted Band (i.e. the number of nucleotides in the
      trusted region of each pattern) for a Trusted Band DNA dictionary.
    }
    \item{}{
      \code{names(x)}:
      The names of the patterns in the PDict object.
    }
    \item{}{
      \code{head(x)}:
      The head of the PDict object (a \code{\link{DNAStringList}} object
      of length \code{length(x)} or \code{NULL}).
    }
    \item{}{
      \code{tail(x)}:
      The tail of the PDict object (a \code{\link{DNAStringList}} object
      of length \code{length(x)} or \code{NULL}).
    }
  }
}

\section{Other functions and generics}{
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
  \link{matchPDict},
  \link{DNA_ALPHABET},
  \link{DNAStringList-class}
}

\examples{
  ## Preprocessing a Uniform Length DNA dictionary
  library(drosophila2probe)
  dict0 <- drosophila2probe$sequence   # The input dictionary.
  length(dict0)                        # Hundreds of thousands of patterns.
  unique(nchar(dict0))                 # Patterns are 25-mers.
  dict0[1:5]
  pdict0 <- PDict(dict0)               # Store the input dictionary into a
                                       # PDict object (preprocessing).
  pdict0
  class(pdict0)
  length(pdict0)                       # Same as length(dict0).
  width(pdict0)                        # The number of chars per pattern.
  sum(duplicated(pdict0))
  table(patternFrequency(pdict0))      # 9 patterns are repeated 3 times.

  ## Creating a Uniform Length DNA dictionary by truncating the input
  ## sequences
  dict1 <- c("ACNG", "GT", "CGT", "AC")
  pdict1a <- PDict(dict1, tb.width=2, drop.tail=TRUE)
  pdict1a
  class(pdict1a)
  length(pdict1a)
  width(pdict1a)
  duplicated(pdict1a)
  patternFrequency(pdict1a)

  ## Preprocessing a Trusted Prefix DNA dictionary
  pdict1b <- PDict(dict1, tb.width=2)
  pdict1b
  class(pdict1b)
  length(pdict1b)
  width(pdict1b)
  head(pdict1b)
  tail(pdict1b)
}

\keyword{methods}
\keyword{classes}
