\name{BString, DNAString, RNAString, AAString}
\docType{class}

% Classes
\alias{class:BString}
\alias{BString-class}
\alias{class:DNAString}
\alias{DNAString-class}
\alias{class:RNAString}
\alias{RNAString-class}
\alias{class:AAString}
\alias{AAString-class}

% Accesor methods:
\alias{letter}
\alias{letter,BString-method}

% Constructor-like functions and generics:
\alias{initialize,BString-method}
\alias{BString}
\alias{initialize,DNAString-method}
\alias{DNAString}
\alias{initialize,RNAString-method}
\alias{RNAString}
\alias{initialize,AAString-method}
\alias{AAString}

% Standard generic methods:
\alias{show,BString-method}
\alias{length,BString-method}
\alias{[,BString-method}
\alias{[<-,BString-method}
\alias{==,BString,ANY-method}
\alias{==,ANY,BString-method}
\alias{!=,BString,ANY-method}
\alias{!=,ANY,BString-method}
\alias{as.character,BString-method}
\alias{toString,BString-method}
\alias{nchar,BString-method}

% Other functions and generics:
\alias{subBString}
\alias{subBString,BString-method}


\title{The BString class and its 3 direct derivations: the DNAString,
  RNAString and AAString classes}

\description{
  The BString class is a general container for storing a big string
  (a long sequence of characters) and to make its manipulation easy
  and efficient.
  The DNAString, RNAString and AAString classes are direct BString
  derivations with the more specific purpose of storing a DNA sequence
  (DNAString), a RNA sequence (RNAString) or a sequence of amino acids
  (AAString).
}

\details{
  The 2 main differences between a BString and a conventional character
  vector are:
  (1) the data stored in a BString are not copied on object
  duplication and (2) a BString can only store a single string (see
  the \link{BStringViews} class for an efficient way to store a
  collection of big strings in a single object).
  A BString object can be used to store {\it as is} any non-empty
  string based on a single-byte character set.
  The DNAString and RNAString classes derive directly from the BString
  class (they have no additional slot). However, these classes can only
  store strings based on a limited alphabet (the IUPAC extended genetic
  alphabet + the gap letter) and these letters are encoded (before they
  are stored in the object) in a way that optimizes fast search algorithms.
  This makes the DNAString and RNAString classes particularly well suited
  for storing big DNA or RNA sequences.
  All the methods that are defined for a BString object also work with a
  DNAString, RNAString or AAString object.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BString (or derived) object
  and \code{i} is a numeric vector.

  \describe{
    \item{}{
      \code{letter(x, i)}:
      [TODO: Document me]
    }
  }
}

\section{Constructor-like functions and generics}{
  In the code snippets below,
  \code{src} can be a character vector
  or a BString (or derived) object.

  \describe{
    \item{}{
      \code{BString(src)}:
      [TODO: Document me]
    }
    \item{}{
      \code{DNAString(src)}:
      [TODO: Document me]
    }
    \item{}{
      \code{RNAString(src)}:
      [TODO: Document me]
    }
    \item{}{
      \code{AAString(src)}:
      [TODO: Document me]
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2}
  are BString (or derived) objects,
  and \code{i} is a numeric vector.

  \describe{
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      Used to get the length of a BString, i.e., its number of letters
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
  In the code snippets below,
  \code{x} is a BString (or derived) object
  and \code{first} and \code{last} are single numerics.

  \describe{
    \item{}{
      \code{subBString(x, first, last)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BStringViews}
}

\examples{
  b <- BString("I am a BString object")
  length(b)
  d <- DNAString("TTGAAAA-CTC-N")
  length(d)
  r <- RNAString(d)
}

\keyword{methods}
\keyword{classes}
