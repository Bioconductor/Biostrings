\name{XString-class}
\docType{class}

% Classes:
\alias{class:XString}
\alias{XString-class}
\alias{class:BString}
\alias{BString-class}

% Constructor-like functions and generics:
\alias{XString}
\alias{BString}

% Accessor methods:
\alias{alphabet}
\alias{alphabet,ANY-method}
\alias{nchar,XString-method}

% Coercion:
\alias{seqtype,BString-method}
\alias{seqtype,DNAString-method}
\alias{seqtype,RNAString-method}
\alias{seqtype,AAString-method}
\alias{seqtype<-,XString-method}
\alias{coerce,XString,BString-method}
\alias{coerce,XString,DNAString-method}
\alias{coerce,XString,RNAString-method}
\alias{coerce,XString,AAString-method}
\alias{coerce,character,BString-method}
\alias{coerce,character,DNAString-method}
\alias{coerce,character,RNAString-method}
\alias{coerce,character,AAString-method}
\alias{coerce,character,XString-method}
\alias{as.character,XString-method}
\alias{toString,XString-method}
\alias{as.vector,XString-method}

% "show" method:
\alias{show,XString-method}
\alias{showAsCell,XString-method}

% Equality:
\alias{==,XString,XString-method}
\alias{==,BString,character-method}
\alias{==,character,BString-method}

\alias{!=,XString,XString-method}
\alias{!=,BString,character-method}
\alias{!=,character,BString-method}

% "substr" and "substring" methods:
\alias{substr,XString-method}
\alias{substring,XString-method}

% "updateObject" method:
\alias{updateObject,XString-method}


\title{BString objects}

\description{
  The BString class is a general container for storing
  a big string (a long sequence of characters) and for making its
  manipulation easy and efficient.

  The \link{DNAString}, \link{RNAString} and \link{AAString} classes are
  similar containers but with the more biology-oriented purpose of storing
  a DNA sequence (\link{DNAString}), an RNA sequence (\link{RNAString}),
  or a sequence of amino acids (\link{AAString}).

  All those containers derive directly (and with no additional slots)
  from the XString virtual class.
}

\details{
  The 2 main differences between an XString object and a standard character
  vector are: (1) the data stored in an XString object are not copied on
  object duplication and (2) an XString object can only store
  a single string (see the \link{XStringSet} container for an
  efficient way to store a big collection of strings in a single object).

  Unlike the \link{DNAString}, \link{RNAString} and \link{AAString}
  containers that accept only a predefined set of letters (the alphabet),
  a BString object can be used for storing any single string based on a
  single-byte character set.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{x} can be a single string (character vector of length 1)
  or an XString object.

  \describe{
    \item{\code{BString(x="", start=1, nchar=NA)}:}{
      Tries to convert \code{x} into a BString object by reading
      \code{nchar} letters starting at position \code{start} in \code{x}.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is an XString object.

  \describe{
    \item{\code{alphabet(x)}:}{
      \code{NULL} for a \code{BString} object.
      See the corresponding man pages when \code{x} is a 
      \link{DNAString}, \link{RNAString} or \link{AAString} object.
    }
    \item{\code{length(x)}:}{ or \code{nchar(x)}:
      Get the length of an XString object, i.e., its number of letters.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is an XString object.

  \describe{
    \item{\code{as.character(x)}:}{
      Converts \code{x} to a character string.
    }
    \item{\code{toString(x)}:}{
      Equivalent to \code{as.character(x)}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is an XString object.

  \describe{
    \item{\code{x[i]}:}{
      Return a new XString object made of the selected letters (subscript
      \code{i} must be an NA-free numeric vector specifying the positions of
      the letters to select).
      The returned object belongs to the same class as \code{x}.

      Note that, unlike \code{subseq}, \code{x[i]} does copy the sequence
      data and therefore will be very inefficient for extracting a big number
      of letters (e.g. when \code{i} contains millions of positions).
    }
  }
}

\section{Equality}{
  In the code snippets below, \code{e1} and \code{e2} are XString objects.

  \describe{
    \item{\code{e1 == e2}:}{
      \code{TRUE} if \code{e1} is equal to \code{e2}.
      \code{FALSE} otherwise.

      Comparison between two XString objects of different base types (e.g.
      a BString object and a \link{DNAString} object) is not supported
      with one exception: a \link{DNAString} object and an \link{RNAString}
      object can be compared (see \link{RNAString-class} for more details
      about this).

      Comparison between a BString object and a character string
      is also supported (see examples below).
    }
    \item{\code{e1 != e2}:}{
      Equivalent to \code{!(e1 == e2)}.
    }
  }
}

\author{H. Pagès}

\seealso{
  \code{\link[XVector]{subseq}},
  \code{\link{letter}},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \code{\link{reverseComplement}},
  \code{\link[XVector]{compact}},
  \link[XVector]{XVector-class}
}

\examples{
b <- BString("I am a BString object")
b
length(b)

## Extracting a linear subsequence:
subseq(b)
subseq(b, start=3)
subseq(b, start=-3)
subseq(b, end=-3)
subseq(b, end=-3, width=5)

## Subsetting:
b2 <- b[length(b):1]       # better done with reverse(b)

as.character(b2)

b2 == b                    # FALSE
b2 == as.character(b2)     # TRUE

## b[1:length(b)] is equal but not identical to b!
b == b[1:length(b)]        # TRUE
identical(b, 1:length(b))  # FALSE
## This is because subsetting an XString object with [ makes a copy
## of part or all its sequence data. Hence, for the resulting object,
## the internal slot containing the memory address of the sequence
## data differs from the original. This is enough for identical() to
## see the 2 objects as different.

## Compacting. As a particular type of XVector objects, XString
## objects can optionally be compacted. Compacting is done typically
## before serialization. See ?compact for more information.
}

\keyword{methods}
\keyword{classes}
