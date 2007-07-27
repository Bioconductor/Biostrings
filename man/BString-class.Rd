\name{BString-class}
\docType{class}

% Classes
\alias{class:BString}
\alias{BString-class}

% Accesor methods:
\alias{letter}
\alias{letter,BString-method}

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

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a \code{BString} (or derived) object
  and \code{i} is a numeric vector.

  \describe{
    \item{}{
      \code{letter(x, i)}:
      [TODO: Document me]
    }
  }
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
      a \code{\link{DNAString}} object, an \code{\link{RNAString}} object
      or an \code{\link{AAString}} object.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2}
  are \code{BString} (or derived) objects,
  and \code{i} is a numeric vector.

  \describe{
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      Used to get the length of a \code{BString}, i.e., its number of letters
    }
    \item{}{
      \code{x[i]}:
      [TODO: Document me]
    }
    \item{}{
      \code{e1 == e2}:
      [TODO: Document me]
    }
    \item{}{
      \code{e1 != e2}:
      [TODO: Document me]
    }
    \item{}{
      \code{as.character(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{toString(x)}:
      [TODO: Document me]
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
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{BStringViews-class}
}

\examples{
  b <- BString("I am a BString object")
  length(b)
}

\keyword{methods}
\keyword{classes}
