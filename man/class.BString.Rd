\name{BString, DNAString, RNAString}
\docType{class}

\alias{class:BString}
\alias{BString-class}

\alias{length,BString-method}
\alias{toString,BString-method}
\alias{initialize,BString-method}
\alias{show,BString-method}
\alias{==,BString,ANY-method}
\alias{==,ANY,BString-method}
\alias{!=,BString,ANY-method}
\alias{!=,ANY,BString-method}
\alias{[,BString-method}
\alias{[<-,BString-method}

\alias{BString}

\alias{class:DNAString}
\alias{DNAString-class}

\alias{initialize,DNAString-method}

\alias{DNAString}

\alias{class:RNAString}
\alias{RNAString-class}

\alias{initialize,RNAString-method}

\alias{RNAString}

\title{Containers for big strings in general, and DNA or RNA sequences
  in particular}

\description{
  BString is the base class designed to store a big string and to make
  its manipulation easy and efficient.
  It's also the common parent of the DNAString and RNAString classes.
}

\details{
  The 2 main differences between a BString and a conventional character
  vector are:
  (1) the data stored in a BString are not copied on object
  duplication and (2) a BString can only store a single string (see
  the \link{BStringViews} class for an efficient way to store a
  collection of big strings in a single object).
  A BString object can be used to store any non-empty string
  based on a single-byte character set.
  The DNAString and RNAString classes derive directly from the BString
  class (they have no additional slot). However, these classes can only
  store strings based on a limited alphabet (the IUPAC extended genetic
  alphabet + the gap letter) and these letters are encoded (before they
  are stored in the object) in a way that optimizes fast search algorithms.
  This makes the DNAString and RNAString classes particularly well suited
  for storing big DNA or RNA sequences.
}

\section{Methods}{
  Standard generic methods:
  \describe{
    \item{\code{length(BString)}:}{Used to get the length of a BString}
    \item{\code{toString(BString)}:}{[TODO: Document me]}
    \item{\code{initialize(BString)}:}{[TODO: Document me]}
    \item{\code{show(BString)}:}{[TODO: Document me]}
    \item{\code{BString1 == BString2}:}{[TODO: Document me]}
    \item{\code{BString1 != BString2}:}{[TODO: Document me]}
    \item{\code{BString[i]}:}{[TODO: Document me]}
    \item{\code{BString[i] <- "attempt to modify the value of a BString"}:}{[TODO: Document me]}
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
