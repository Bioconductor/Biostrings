\name{BStringList-class}
\docType{class}

% Classes:
\alias{class:BStringList}
\alias{BStringList-class}
\alias{BStringList}

% Constructor-like functions and generics:
\alias{initialize,BStringList-method}
\alias{initialize,DNAStringList-method}
\alias{initialize,RNAStringList-method}
\alias{initialize,AAStringList-method}

% Accesor methods:
\alias{as.list,BStringList-method}
\alias{length,BStringList-method}
\alias{nchar,BStringList-method}
\alias{desc}
\alias{desc,BStringList-method}
\alias{desc<-}
\alias{desc<-,BStringList-method}

% Other methods:

\alias{show,BStringList-method}
\alias{[,BStringList-method}
\alias{[<-,BStringList-method}
\alias{[[,BStringList-method}
\alias{[[<-,BStringList-method}


\title{The BStringList, DNAStringList, RNAStringList and AAStringList classes}

\description{
  The BStringList class is a container for storing a list of \code{\link{BString}}
  objects and for making its manipulation easy and efficient.
  The DNAStringList (or RNAStringList, or AAStringList) class is a direct
  BStringList derivation for storing a list of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) objects.
}

\details{
  [TODO]
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringList (or derived) object.

  \describe{
    \item{}{
      \code{as.list(x)}:
      Return the elements in \code{x} as a standard R list.
    }
    \item{}{
      \code{length(x)}:
      Equivalent to \code{length(as.list(x))}.
    }
    \item{}{
      \code{nchar(x)}:
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
    }
    \item{}{
      \code{desc(x)}:
      \code{NULL} or a character vector of the same length as \code{x} containing
      a short user-provided description or comment for each element in \code{x}.
      These are the only data in a BStringList object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is a BStringList (or derived) object,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new BStringList object made of the selected elements.
    }
    \item{}{
      \code{x[[i]]}:
      Equivalent to \code{as.list(x)[[i]]}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{BStringViews-class}
}

\examples{
  seqs <- new("DNAStringList", list(seq1=DNAString("TTGA"), seq2=DNAString("-CTC-N")))
  seqs
  desc(seqs)
  desc(seqs)[2] <- "seqB"
}

\keyword{methods}
\keyword{classes}
