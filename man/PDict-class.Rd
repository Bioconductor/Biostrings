\name{PDict-class}
\docType{class}

% PreprocessedTB class:
\alias{class:PreprocessedTB}
\alias{PreprocessedTB-class}
\alias{PreprocessedTB}

\alias{length,PreprocessedTB-method}
\alias{width,PreprocessedTB-method}
\alias{tb}
\alias{tb,PreprocessedTB-method}
\alias{tb.width}
\alias{tb.width,PreprocessedTB-method}
\alias{nnodes}
\alias{hasAllFlinks}
\alias{computeAllFlinks}
\alias{initialize,PreprocessedTB-method}
\alias{duplicated,PreprocessedTB-method}

% Twobit class:
\alias{class:Twobit}
\alias{Twobit-class}
\alias{Twobit}

\alias{show,Twobit-method}
\alias{initialize,Twobit-method}

% ACtree2 class:
\alias{class:ACtree2}
\alias{ACtree2-class}
\alias{ACtree2}

\alias{nnodes,ACtree2-method}
\alias{hasAllFlinks,ACtree2-method}
\alias{computeAllFlinks,ACtree2-method}
\alias{show,ACtree2-method}
\alias{initialize,ACtree2-method}

% PDict3Parts class:
\alias{class:PDict3Parts}
\alias{PDict3Parts-class}
\alias{PDict3Parts}

\alias{length,PDict3Parts-method}
\alias{width,PDict3Parts-method}
\alias{head,PDict3Parts-method}
\alias{tb,PDict3Parts-method}
\alias{tb.width,PDict3Parts-method}
\alias{tail,PDict3Parts-method}

% PDict class:
\alias{class:PDict}
\alias{PDict-class}

\alias{length,PDict-method}
\alias{width,PDict-method}
\alias{names,PDict-method}
\alias{names<-,PDict-method}
\alias{[[,PDict-method}
\alias{duplicated,PDict-method}
\alias{patternFrequency}
\alias{patternFrequency,PDict-method}

% TB_PDict class:
\alias{class:TB_PDict}
\alias{TB_PDict-class}
\alias{TB_PDict}

\alias{head,TB_PDict-method}
\alias{tb,TB_PDict-method}
\alias{tb.width,TB_PDict-method}
\alias{tail,TB_PDict-method}
\alias{show,TB_PDict-method}

% MTB_PDict class:
\alias{class:MTB_PDict}
\alias{MTB_PDict-class}
\alias{MTB_PDict}

\alias{as.list,MTB_PDict-method}
\alias{show,MTB_PDict-method}

% Expanded_TB_PDict class:
\alias{class:Expanded_TB_PDict}
\alias{Expanded_TB_PDict-class}
\alias{Expanded_TB_PDict}

% PDict() constructor:
\alias{PDict}
\alias{PDict,character-method}
\alias{PDict,DNAStringSet-method}
\alias{PDict,XStringViews-method}
\alias{PDict,AsIs-method}
\alias{PDict,probetable-method}


\title{PDict objects}

\description{
  The PDict class is a container for storing a preprocessed dictionary of DNA
  patterns that can later be passed to the \code{\link{matchPDict}} function
  for fast matching against a reference sequence (the subject).

  \code{PDict} is the constructor function for creating new PDict objects.
}

\usage{
PDict(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
         algorithm="ACtree2", skip.invalid.patterns=FALSE)
}

\arguments{
  \item{x}{
    A character vector, a \link{DNAStringSet} object or
    an \link{XStringViews} object with a \link{DNAString} subject.
  }
  \item{max.mismatch}{
    A single non-negative integer or \code{NA}. See the "Allowing
    a small number of mismatching letters" section below.
  }
  \item{tb.start,tb.end,tb.width}{
    A single integer or \code{NA}. See the "Trusted Band" section below.
  }
  \item{algorithm}{
    \code{"ACtree2"} (the default) or \code{"Twobit"}.
  }
  \item{skip.invalid.patterns}{
    This argument is not supported yet (and might in fact be replaced
    by the \code{filter} argument very soon).
  }
}

\details{
  THIS IS STILL WORK IN PROGRESS!

  If the original dictionary \code{x} is a character vector or
  an \link{XStringViews} object with a \link{DNAString} subject,
  then the \code{PDict} constructor will first try to turn it
  into a \link{DNAStringSet} object.

  By default (i.e. if \code{PDict} is called with \code{max.mismatch=NA},
  \code{tb.start=NA}, \code{tb.end=NA} and \code{tb.width=NA})
  the following limitations apply: (1) the original dictionary can only
  contain base letters (i.e. only As, Cs, Gs and Ts), therefore IUPAC
  ambiguity codes are not allowed; (2) all the
  patterns in the dictionary must have the same length ("constant width"
  dictionary); and (3) later \code{matchPdict} can only be used with
  \code{max.mismatch=0}.

  A Trusted Band can be used in order to relax these limitations (see
  the "Trusted Band" section below).

  If you are planning to use the resulting \code{PDict} object in order
  to do inexact matching where valid hits are allowed to have a small
  number of mismatching letters, then see the "Allowing a small number
  of mismatching letters" section below.

  Two preprocessing algorithms are currently supported:
  \code{algorithm="ACtree2"} (the default) and \code{algorithm="Twobit"}.
  With the \code{"ACtree2"} algorithm, all the oligonucleotides in the
  Trusted Band are stored in a 4-ary Aho-Corasick tree.
  With the \code{"Twobit"} algorithm, the 2-bit-per-letter
  signatures of all the oligonucleotides in the Trusted Band are computed
  and the mapping from these signatures to the 1-based position of the
  corresponding oligonucleotide in the Trusted Band is stored in a way that
  allows very fast lookup.
  Only PDict objects preprocessed with the \code{"ACtree2"} algo can then
  be used with \code{matchPdict} (and family) and with \code{fixed="pattern"}
  (instead of \code{fixed=TRUE}, the default), so that IUPAC ambiguity codes
  in the subject are treated as ambiguities. PDict objects obtained with the
  \code{"Twobit"} algo don't allow this.
  See \code{?`\link{matchPDict-inexact}`} for more information about support
  of IUPAC ambiguity codes in the subject.
}

\section{Trusted Band}{
  What's a Trusted Band?

  A Trusted Band is a region defined in the original dictionary
  where the limitations described above will apply.

  Why use a Trusted Band?

  Because the limitations described above will apply to the Trusted Band only!
  For example the Trusted Band cannot contain IUPAC ambiguity codes but the
  "head" and the "tail" can (see below for what those are).
  Also with a Trusted Band, if \code{matchPdict} is called with a non-null
  \code{max.mismatch} value then mismatching letters will be allowed in the
  head and the tail. Or, if \code{matchPdict} is called with
  \code{fixed="subject"}, then IUPAC ambiguity codes in the head and the
  tail will be treated as ambiguities.

  How to specify a Trusted Band?

  Use the \code{tb.start}, \code{tb.end} and \code{tb.width} arguments of the
  \code{PDict} constructor in order to specify a Trusted Band.
  This will divide each pattern in the original dictionary into three parts:
  a left part, a middle part and a right part.
  The middle part is defined by its starting and ending nucleotide positions
  given relatively to each pattern thru the \code{tb.start}, \code{tb.end}
  and \code{tb.width} arguments. It must have the same length for all
  patterns (this common length is called the width of the Trusted Band).
  The left and right parts are defined implicitely: they are the
  parts that remain before (prefix) and after (suffix) the middle part,
  respectively.
  Therefore three \link{DNAStringSet} objects result from this division:
  the first one is made of all the left parts and forms the head of the PDict
  object, the second one is made of all the middle parts and forms the Trusted
  Band of the PDict object, and the third one is made of all the right parts
  and forms the tail of the PDict object.

  In other words you can think of the process of specifying a Trusted Band 
  as drawing 2 vertical lines on the original dictionary (note that these
  2 lines are not necessarily straight lines but the horizontal space between
  them must be constant). When doing this, you are dividing the dictionary
  into three regions (from left to right): the head, the Trusted Band and the
  tail. Each of them is a \link{DNAStringSet} object with the same number of
  elements than the original dictionary and the original dictionary could
  easily be reconstructed from those three regions.

  The width of the Trusted Band must be >= 1 because Trusted Bands of
  width 0 are not supported.

  Finally note that calling \code{PDict} with \code{tb.start=NA},
  \code{tb.end=NA} and \code{tb.width=NA} (the default) is equivalent
  to calling it with \code{tb.start=1}, \code{tb.end=-1} and
  \code{tb.width=NA}, which results in a full-width Trusted Band i.e.
  a Trusted Band that covers the entire dictionary (no head and no tail).
}

\section{Allowing a small number of mismatching letters}{
  [TODO]
}

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is a PDict object.

  \describe{
    \item{\code{length(x)}:}{
      The number of patterns in \code{x}.
    }
    \item{\code{width(x)}:}{
      A vector of non-negative integers containing the number
      of letters for each pattern in \code{x}.
    }
    \item{\code{names(x)}:}{
      The names of the patterns in \code{x}.
    }
    \item{\code{head(x)}:}{
      The head of \code{x} or \code{NULL} if \code{x} has no head.
    }
    \item{\code{tb(x)}:}{
      The Trusted Band defined on \code{x}.
    }
    \item{\code{tb.width(x)}:}{
      The width of the Trusted Band defined on \code{x}.
      Note that, unlike \code{width(tb(x))}, this is a single integer.
      And because the Trusted Band has a constant width, \code{tb.width(x)}
      is in fact equivalent to \code{unique(width(tb(x)))},
      or to \code{width(tb(x))[1]}.
    }
    \item{\code{tail(x)}:}{
      The tail of \code{x} or \code{NULL} if \code{x} has no tail.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is a PDict object.

  \describe{
    \item{\code{x[[i]]}:}{
      Extract the i-th pattern from \code{x} as a \link{DNAString} object.
    }
  }
}

\section{Other methods}{
  In the code snippet below,
  \code{x} is a PDict object.

  \describe{
    \item{\code{duplicated(x)}:}{
      [TODO]
    }
    \item{\code{patternFrequency(x)}:}{
      [TODO]
    }
  }
}

\author{H. Pagès}

\references{
  Aho, Alfred V.; Margaret J. Corasick (June 1975). "Efficient string
  matching: An aid to bibliographic search".
  Communications of the ACM 18 (6): 333-340.
}

\seealso{
  \code{\link{matchPDict}},
  \code{\link{DNA_ALPHABET}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{DNAStringSet-class},
  \link{XStringViews-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. NO HEAD AND NO TAIL (THE DEFAULT)
  ## ---------------------------------------------------------------------
  library(drosophila2probe)
  dict0 <- DNAStringSet(drosophila2probe)
  dict0                                # The original dictionary.
  length(dict0)                        # Hundreds of thousands of patterns.
  unique(nchar(dict0))                 # Patterns are 25-mers.

  pdict0 <- PDict(dict0)               # Store the original dictionary in
                                       # a PDict object (preprocessing).
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
