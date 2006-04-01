\name{BString, DNAString, RNAString}
\docType{class}

% Classes
\alias{class:BString}
\alias{BString-class}
\alias{class:DNAString}
\alias{DNAString-class}
\alias{class:RNAString}
\alias{RNAString-class}

% Constructor-like functions and generics:
\alias{initialize,BString-method}
\alias{BString}
\alias{initialize,DNAString-method}
\alias{DNAString}
\alias{initialize,RNAString-method}
\alias{RNAString}

% Standard generic methods:
\alias{length,BString-method}
\alias{[,BString-method}
\alias{[<-,BString-method}
\alias{show,BString-method}
\alias{==,BString,ANY-method}
\alias{==,ANY,BString-method}
\alias{!=,BString,ANY-method}
\alias{!=,ANY,BString-method}
\alias{as.character,BString-method}
\alias{toString,BString-method}
\alias{nchar,BString-method}

% Other functions and generics:
\alias{letter}
\alias{subBString}
\alias{subBString,BString-method}


\title{The BString, DNAString and RNAString classes}

\description{
  The BString, DNAString and RNAString classes are containers for big
  strings in general, and DNA or RNA sequences in particular.
  BString is the base class designed to store a big string
  and to make its manipulation easy and efficient.
  It's also the common parent of the DNAString and RNAString classes.
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
  When the strings you want to work with are not DNA or RNA sequences,
  then you can use the BString class.
  All the methods that are defined for a BString object also work with a
  DNAString or RNAString object.
}

\section{Constructor-like functions and generics}{
  In the code snippets below,
  \code{src} can be a character vector
  or a BString (or DNAString or RNAString) object.

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
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2}
  are BString (or DNAString or RNAString) objects,
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
      \code{show(object)}:
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
  \code{x} is a BString (or DNAString or RNAString) object,
  \code{i} is a numeric vector,
  and \code{first} and \code{last} are single numerics.

  \describe{
    \item{}{
      \code{letter(x, i)}:
      [TODO: Document me]
    }
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
