\name{XStringSetList-class}
\docType{class}

% Classes:

\alias{class:XStringSetList}
\alias{XStringSetList-class}
\alias{XStringSetList}

\alias{class:BStringSetList}
\alias{BStringSetList-class}
\alias{BStringSetList}

\alias{class:DNAStringSetList}
\alias{DNAStringSetList-class}
\alias{DNAStringSetList}

\alias{class:RNAStringSetList}
\alias{RNAStringSetList-class}
\alias{RNAStringSetList}

\alias{class:AAStringSetList}
\alias{AAStringSetList-class}
\alias{AAStringSetList}

% Methods:

\alias{relistToClass,XStringSet-method}

\alias{seqtype,XStringSetList-method}
\alias{seqtype<-,XStringSetList-method}

\alias{show,XStringSetList-method}
\alias{showAsCell,XStringSetList-method}

\alias{coerce,list,XStringSetList-method}
\alias{coerce,List,XStringSetList-method}
\alias{coerce,list,BStringSetList-method}
\alias{coerce,List,BStringSetList-method}
\alias{coerce,list,DNAStringSetList-method}
\alias{coerce,List,DNAStringSetList-method}
\alias{coerce,list,RNAStringSetList-method}
\alias{coerce,List,RNAStringSetList-method}
\alias{coerce,list,AAStringSetList-method}
\alias{coerce,List,AAStringSetList-method}

\alias{nchar,XStringSetList-method}


\title{XStringSetList objects}

\description{
  The XStringSetList class is a virtual container for storing a list
  of \link{XStringSet} objects.
}

\details{
  Concrete flavors of the XStringSetList container are the BStringSetList,
  DNAStringSetList, RNAStringSetList and AAStringSetList containers
  for storing a list of \link{BStringSet}, \link{DNAStringSet},
  \link{RNAStringSet} and \link{AAStringSet} objects, respectively.
  These four containers are direct subclasses of XStringSetList
  with no additional slots. The XStringSetList class itself
  is virtual and has no constructor.
}

\usage{
## Constructors:
BStringSetList(..., use.names=TRUE)
DNAStringSetList(..., use.names=TRUE)
RNAStringSetList(..., use.names=TRUE)
AAStringSetList(..., use.names=TRUE)
}

\arguments{
  \item{\dots}{
    Character vector(s) (with no NAs), or \link{XStringSet} object(s), or
    \link{XStringViews} object(s) to be concatenated into a
    \link{XStringSetList}.
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
}

\section{Methods}{
  The XStringSetList class extends the \link[S4Vectors]{List} class defined
  in the \pkg{IRanges} package. Using a less technical jargon, this just means
  that an XStringSetList object is a list-like object that can be manipulated
  like an ordinary list. Or, said otherwise, most of the operations that work
  on an ordinary list (e.g. \code{length}, \code{names}, \code{[}, \code{[[},
  \code{c}, \code{unlist}, etc...) should work on an XStringSetList object.
  In addition, Bioconductor specific list operations like
  \code{\link[S4Vectors]{elementNROWS}} and
  \code{\link[IRanges]{PartitioningByEnd}} (defined in the \pkg{IRanges}
  package) are supported too.
}

\author{H. Pagès}

\seealso{
  \link{XStringSet-class},
  \link[S4Vectors]{List-class}
}

\examples{
## ------------------------------------------------------------------------
## A. THE XStringSetList CONSTRUCTORS
## ------------------------------------------------------------------------

dna1 <- c("AAA", "AC", "", "T", "GGATA")
dna2 <- c("G", "TT", "C")

x <- DNAStringSetList(dna1, dna2)
x

DNAStringSetList(DNAStringSet(dna1), DNAStringSet(dna2))

DNAStringSetList(dna1, DNAStringSet(dna2))

DNAStringSetList(DNAStringSet(dna1), dna2)

DNAStringSetList(dna1, RNAStringSet(DNAStringSet(dna2)))

DNAStringSetList(list(dna1, dna2))

DNAStringSetList(CharacterList(dna1, dna2))

## Empty object (i.e. zero-length):
DNAStringSetList()

## Not empty (length is 1):
DNAStringSetList(character(0))

## ---------------------------------------------------------------------
## B. UNLISTING AN XStringSetList OBJECT
## ---------------------------------------------------------------------
length(x)
elementNROWS(x)
unlist(x)
x[[1]]
x[[2]]
as.list(x)

names(x) <- LETTERS[1:length(x)]
x[["A"]]
x[["B"]]
as.list(x)  # named list
}

\keyword{methods}
\keyword{classes}
