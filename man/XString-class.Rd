\name{XString-class}
\docType{class}

% Classes
\alias{class:XString}
\alias{XString-class}
\alias{class:BString}
\alias{BString-class}

% Constructor-like functions and generics:
\alias{XString}
\alias{XString,character-method}
\alias{XString,AsIs-method}
\alias{XString,XString-method}
\alias{BString}

% Standard generic methods:
\alias{show,XString-method}
\alias{length,XString-method}
\alias{nchar,XString-method}
\alias{[,XString-method}
\alias{[<-,XString-method}

\alias{==,XString,XString-method}
\alias{==,BString,character-method}
\alias{==,character,BString-method}

\alias{!=,XString,XString-method}
\alias{!=,BString,character-method}
\alias{!=,character,BString-method}

\alias{as.character,XString-method}
\alias{toString,XString-method}

% Other functions and generics:
\alias{alphabet}
\alias{alphabet,BString-method}


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
  from the XString virtual class. They are also said to be direct XString
  subtypes with no additional slots.
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
  \code{src} can be a single string (character vector of length 1)
  or an XString object.

  \describe{
    \item{}{
      \code{BString(src, start=1, nchar=NA)}:
      Tries to convert \code{src} into a BString object by reading
      \code{nchar} letters starting at position \code{start} in \code{src}.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2}
  are XString objects,
  and \code{i} is a NA-free numeric vector.

  \describe{
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      Get the length of an XString object, i.e., its number of letters.
    }
    \item{}{
      \code{x[i]}:
      Return a new XString object made of the selected letters.
      The returned object belongs to the same class as \code{x}
      (same XString substype).
    }
    \item{}{
      \code{e1 == e2}:
      \code{TRUE} if \code{e1} is equal to \code{e2}.
      \code{FALSE} otherwise.

      Comparison between two XString objects of different substypes (e.g.
      a BString object and a \link{DNAString} object) is not supported
      with one exception: a \link{DNAString} object and an \link{RNAString}
      object can be compared (see \link{RNAString-class} for more details
      about this).

      Comparison between a BString object and a character string
      is also supported (see examples below).
    }
    \item{}{
      \code{e1 != e2}:
      Equivalent to \code{!(e1 == e2)}.
    }
    \item{}{
      \code{as.character(x)}:
      Converts \code{x} to a character string.
    }
    \item{}{
      \code{toString(x)}:
      Equivalent to \code{as.character(x)}.
    }
  }
}

\section{Other functions and generics}{
  In the code snippet below,
  \code{x} must be an XString object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      \code{NULL} for a \code{BString} object.
      See the corresponding man pages when \code{x} is a 
      \link{DNAString}, an \link{RNAString} or an \link{AAString}
      object.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{subXString}},
  \code{\link{letter}},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{XStringSet-class},
  \link{XStringViews-class}
}

\examples{
  b <- BString("I am a BString object")
  b
  length(b)
  b2 <- b[length(b):1]
  as.character(b2)

  ## b[1:length(b)] is equal but not identical to b!
  b == b[1:length(b)]        # TRUE
  identical(b, 1:length(b))  # FALSE
  ## This is because subsetting an XString object with [ makes a copy
  ## of part or all its sequence data. Hence, for the resulting object,
  ## the internal slot containing the memory address of the sequence
  ## data differs from the original. This is enough for identical() to
  ## see the 2 objects as different. See ?subXString for a much more
  ## efficient substring extraction (does NOT copy the sequence data of
  ## the original object).

  b2 == b                    # FALSE
  b2 == as.character(b2)     # TRUE
}

\keyword{methods}
\keyword{classes}
