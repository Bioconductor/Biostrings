\name{BString-class}
\docType{class}

% Classes
\alias{class:BString}
\alias{BString-class}

% Constructor-like functions and generics:
\alias{initialize,BString-method}
\alias{BString}

% Standard generic methods:
\alias{show,BString-method}
\alias{length,BString-method}
\alias{[,BString-method}
\alias{[<-,BString-method}
\alias{==,BString,BString-method}
\alias{!=,BString,BString-method}
\alias{==,BString,ANY-method}
\alias{!=,BString,ANY-method}
\alias{==,ANY,BString-method}
\alias{!=,ANY,BString-method}
\alias{as.character,BString-method}
\alias{toString,BString-method}
\alias{nchar,BString-method}

% Other functions and generics:
\alias{alphabet,BString-method}


\title{The BString class}

\description{
  The \code{BString} class is a general container for storing
  a big string (a long sequence of characters) and to make its
  manipulation easy and efficient.
  The \code{\link{DNAString}}, \code{\link{RNAString}}
  and \code{\link{AAString}} classes are direct \code{BString}
  derivations with the more specific purpose of storing a DNA sequence
  (\code{\link{DNAString}}), a RNA sequence (\code{\link{RNAString}})
  or a sequence of amino acids (\code{\link{AAString}}).
}

\details{
  The 2 main differences between a \code{BString} object
  and a standard character vector are:
  (1) the data stored in a \code{BString} object are not copied
  on object duplication and (2) a \code{BString} can only store
  a single string (see the \code{\link{BStringViews}} for an
  efficient way to store a collection of big strings in a single object).
  A \code{BString} object can be used to store any non-empty
  string based on a single-byte character set.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a character vector
  or a \code{BString} (or derived) object.

  \describe{
    \item{}{
      \code{BString(src)}:
      Tries to convert \code{src} into a \code{BString} object.
      \code{src} can be a single string (character vector of length 1),
      a \link{DNAString} object, an \link{RNAString} object
      or an \link{AAString} object.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2}
  are BString (or derived) objects,
  and \code{i} is a NA-free numeric vector.

  \describe{
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      Get the length of a BString object, i.e., its number of letters.
    }
    \item{}{
      \code{x[i]}:
      Return a BString (or derived) object made of the selected letters.
      The returned object belongs to the same class as \code{x}.
    }
    \item{}{
      \code{e1 == e2}:
      \code{TRUE} if \code{e1} is equal to \code{e2}.
      \code{FALSE} otherwise.
      Allows comparison between BString objects that are not instances
      of the same class (see \link{RNAString-class}).
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
  \code{x} is a \code{BString} (or derived) object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      \code{NULL} for a \code{BString} object.
      See the corrsponding man pages when \code{x} is of class
      \code{\link{DNAString}}, \code{\link{RNAString}} or code{\link{AAString}}.
    }
  }
}

\seealso{
  \code{\link{letter}},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{BStringViews-class}
}

\examples{
  b <- BString("I am a BString object")
  b
  length(b)
  b2 <- b[length(b):1]
  as.character(b2)
}

\keyword{methods}
\keyword{classes}
