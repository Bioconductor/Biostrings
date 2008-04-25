\name{XStringSet-class}
\docType{class}

% Classes:
\alias{class:XStringSet}
\alias{XStringSet-class}
\alias{XStringSet}

\alias{class:BStringSet}
\alias{BStringSet-class}
\alias{BStringSet}
\alias{class:DNAStringSet}
\alias{DNAStringSet-class}
\alias{DNAStringSet}
\alias{class:RNAStringSet}
\alias{RNAStringSet-class}
\alias{RNAStringSet}
\alias{class:AAStringSet}
\alias{AAStringSet-class}
\alias{AAStringSet}

\alias{class:XStringList}
\alias{XStringList-class}
\alias{XStringList}

\alias{class:BStringList}
\alias{BStringList-class}
\alias{BStringList}
\alias{class:DNAStringList}
\alias{DNAStringList-class}
\alias{DNAStringList}
\alias{class:RNAStringList}
\alias{RNAStringList-class}
\alias{RNAStringList}
\alias{class:AAStringList}
\alias{AAStringList-class}
\alias{AAStringList}

% Accesor methods:
\alias{nchar,XStringSet-method}
\alias{alphabet,XStringSet-method}

\alias{length,XStringList-method}
\alias{nchar,XStringList-method}
\alias{width,XStringList-method}
\alias{names,XStringList-method}
\alias{names<-,XStringList-method}

% Constructor-like functions and generics:
\alias{initialize,XStringSet-method}
\alias{XStringSet,character-method}
\alias{XStringSet,AsIs-method}
\alias{XStringSet,XString-method}
\alias{XStringSet,XStringSet-method}
\alias{XStringSet,XStringViews-method}

\alias{BStringSet}
\alias{DNAStringSet}
\alias{RNAStringSet}
\alias{AAStringSet}

\alias{initialize,BStringList-method}
\alias{initialize,DNAStringList-method}
\alias{initialize,RNAStringList-method}
\alias{initialize,AAStringList-method}
\alias{BStringList,character-method}
\alias{DNAStringList,character-method}
\alias{RNAStringList,character-method}
\alias{AAStringList,character-method}
\alias{BStringList,XStringList-method}
\alias{DNAStringList,XStringList-method}
\alias{RNAStringList,XStringList-method}
\alias{AAStringList,XStringList-method}
\alias{BStringList,vector-method}
\alias{DNAStringList,vector-method}
\alias{RNAStringList,vector-method}
\alias{AAStringList,vector-method}
\alias{BStringList,ANY-method}
\alias{DNAStringList,ANY-method}
\alias{RNAStringList,ANY-method}
\alias{AAStringList,ANY-method}

% Other methods:
\alias{show,XStringSet-method}
\alias{[[,XStringSet-method}
\alias{[[<-,XStringSet-method}
\alias{as.list,XStringSet-method}
\alias{as.character,XStringSet-method}

\alias{show,XStringList-method}
\alias{[,XStringList-method}
\alias{[<-,XStringList-method}
\alias{[[,XStringList-method}
\alias{[[<-,XStringList-method}
\alias{as.list,XStringList-method}


\title{BStringSet, DNAStringSet, RNAStringSet and AAStringSet objects}

\description{
  The BStringSet class is a container for storing a set of \code{\link{BString}}
  objects and for making its manipulation easy and efficient.
  Similarly, the DNAStringSet (or RNAStringSet, or AAStringSet) class is
  a container for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) objects.

  The BStringList class (and family) is equivalent to the BStringSet family.
  It was a first attempt to address the problem of storing a large set of
  sequences but it has been superseded by the BStringSet family.
}

\usage{
  BStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  DNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  RNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  AAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)

  ## The BStringList family has been superseded by the more efficient
  ## BStringSet family
  BStringList(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  DNAStringList(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  RNAStringList(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  AAStringList(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
}

\arguments{
  \item{x}{
    [TODO]
  }
  \item{start}{
    [TODO]
  }
  \item{end}{
    [TODO]
  }
  \item{width}{
    [TODO]
  }
  \item{use.names}{
    [TODO]
  }
  \item{check}{
    [TODO]
  }
}

\details{
  [TODO]
}

\section{Accesor methods}{
  The XStringSet class derives from the \link{IRanges} class hence
  all the accessor methods defined for a \link{IRanges} object can also
  be used on an XStringSet object. In particular, the following methods
  are available (in the code snippets below, \code{x} is an XStringSet
  object:

  \describe{
    \item{}{
      \code{length(x)}:
      The number of sequences in \code{x}.
    }
    \item{}{
      \code{width(x)}:
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
    }
    \item{}{
      \code{nchar(x)}:
      The same as \code{width(x)}.
    }
    \item{}{
      \code{names(x)}:
      \code{NULL} or a character vector of the same length as \code{x} containing
      a short user-provided description or comment for each element in \code{x}.
      These are the only data in an XStringSet object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is an XStringSet object,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new XStringSet object made of the selected elements.
    }
    \item{}{
      \code{x[[i]]}:
      Extract the i-th \code{\link{BString}} object from \code{x}.
    }
  }
}

\section{Other methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{}{
      \code{as.list(x)}:
      Return the elements in \code{x} as a standard R list.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{XStringViews-class},
  \link{IRanges}
}

\examples{
  x1 <- DNAStringSet(c("TTGA", "-CTC-N"))
  x1
  names(x1)
  names(x1)[2] <- "seqB"
  x1

  library(drosophila2probe)
  x2 <- DNAStringSet(drosophila2probe$sequence)
  x2

  RNAStringSet(x2, start=2, end=-5)
}

\keyword{methods}
\keyword{classes}

