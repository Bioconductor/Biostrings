\name{BStringSet-class}
\docType{class}

% Classes:
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
\alias{nchar,BStringSet-method}

\alias{length,BStringList-method}
\alias{nchar,BStringList-method}
\alias{width,BStringList-method}
\alias{names,BStringList-method}
\alias{names<-,BStringList-method}

% Constructor-like functions and generics:
\alias{initialize,BStringSet-method}
\alias{initialize,DNAStringSet-method}
\alias{initialize,RNAStringSet-method}
\alias{initialize,AAStringSet-method}
\alias{BStringSet,character-method}
\alias{DNAStringSet,character-method}
\alias{RNAStringSet,character-method}
\alias{AAStringSet,character-method}
\alias{BStringSet,AsIs-method}
\alias{DNAStringSet,AsIs-method}
\alias{RNAStringSet,AsIs-method}
\alias{AAStringSet,AsIs-method}
\alias{BStringSet,BStringSet-method}
\alias{DNAStringSet,BStringSet-method}
\alias{RNAStringSet,BStringSet-method}
\alias{AAStringSet,BStringSet-method}
\alias{BStringSet,BStringViews-method}
\alias{DNAStringSet,BStringViews-method}
\alias{RNAStringSet,BStringViews-method}
\alias{AAStringSet,BStringViews-method}

\alias{initialize,BStringList-method}
\alias{initialize,DNAStringList-method}
\alias{initialize,RNAStringList-method}
\alias{initialize,AAStringList-method}
\alias{BStringList,character-method}
\alias{DNAStringList,character-method}
\alias{RNAStringList,character-method}
\alias{AAStringList,character-method}
\alias{BStringList,BStringList-method}
\alias{DNAStringList,BStringList-method}
\alias{RNAStringList,BStringList-method}
\alias{AAStringList,BStringList-method}
\alias{BStringList,vector-method}
\alias{DNAStringList,vector-method}
\alias{RNAStringList,vector-method}
\alias{AAStringList,vector-method}
\alias{BStringList,ANY-method}
\alias{DNAStringList,ANY-method}
\alias{RNAStringList,ANY-method}
\alias{AAStringList,ANY-method}

% Other methods:
\alias{show,BStringSet-method}
\alias{[[,BStringSet-method}
\alias{[[<-,BStringSet-method}
\alias{as.list,BStringSet-method}
\alias{as.character,BStringSet-method}

\alias{show,BStringList-method}
\alias{[,BStringList-method}
\alias{[<-,BStringList-method}
\alias{[[,BStringList-method}
\alias{[[<-,BStringList-method}
\alias{as.list,BStringList-method}


\title{BStringSet objects}

\description{
  The BStringSet class is a container for storing a set of \code{\link{BString}}
  instances and for making its manipulation easy and efficient.
  The DNAStringSet (or RNAStringSet, or AAStringSet) class is a direct
  BStringSet derivation for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) instances.

  The BStringList class (and family) is equivalent to the BStringSet family.
  It was a first attempt to address the problem of storing a large set of
  sequences but it has been superseded by the BStringSet family.
}

\usage{
  BStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  DNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  RNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
  AAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)

  ## The BStringList family has been superseded by more efficient BStringSet
  ## family
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
  The BStringSet class derives from the \link{IRanges} class hence
  all the accessor methods defined for a \link{IRanges} object can also
  be used on a BStringSet object. In particular, the following methods
  are available (in the code snippets below, \code{x} is a BStringSet (or
  derived) object):

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
      These are the only data in a BStringSet object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is a BStringSet (or derived) object,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new BStringSet object made of the selected elements.
    }
    \item{}{
      \code{x[[i]]}:
      Extract the i-th \code{\link{BString}} object from \code{x}.
    }
  }
}

\section{Other methods}{
  In the code snippets below,
  \code{x} is a BStringSet (or derived) object.
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
  \link{BStringViews-class},
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

