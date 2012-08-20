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
\alias{partitioning}
\alias{partitioning,XStringSetList-method}

\alias{seqtype,XStringSetList-method}
\alias{seqtype<-,XStringSetList-method}

\alias{[[,XStringSetList-method}

\alias{show,XStringSetList-method}


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
  with no additional slots.

  Currently \code{DNAStringSetList()} is the only XStringSetList constructor. 
  The XStringSetList class itself is virtual and has no constructor.
}

\usage{
## Constructors:
DNAStringSetList(..., use.names=TRUE)
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


\section{Accessor-like methods}{
  In the code snippets below,
  \code{x} is an XStringSetList object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of outer list elements in \code{x}.
    }
    \item{}{
      \code{partitioning(x)}:
      A PartitioningByEnd instance describing the block grouping of the 
      list elements of \code{x}.
    }
    \item{}{
      \code{elementLengths(x)}:
      An integer vector of the length of each list element of \code{x}.
    }
   \item{}{
      \code{names(x)}:
      \code{NULL} or a character vector of the same length as \code{x} 
      containing a short user-provided description or comment for each 
      element in \code{x}. These are the only data in an XStringSetList 
      object that can safely be changed by the user. All the other data 
      are immutable! As a general recommendation, the user should never 
      try to modify an object by accessing its slots directly.
    }
  }
}

\section{Subsetting and appending}{
  In the code snippets below,
  \code{x} and \code{values} are XStringSet objects,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the i-th \code{\link{XStringSet}} object from \code{x}.
    }
    \item{}{
      \code{append(x, values, after=length(x))}:
      Add sequences in \code{values} to \code{x}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{XStringSet-class},
  \link[IRanges]{Grouping-class},
  \link[IRanges]{Vector-class}
}

\examples{
  ## ------------------------------------------------------------------------
  ## A. THE XStringSetList CONSTRUCTORS
  ## ------------------------------------------------------------------------
  ## Currently DNAStringSetList() is the only constructor. Others will 
  ## be developed when the use case arises.

  dna1 <- DNAStringSet(c("AAA", "AC", "", "T", "GGATA"))
  dna2 <- DNAStringSet(c("G", "TT", "C"))
  x <- DNAStringSetList(dna1, dna2)
  x

  y <- DNAStringSetList(c("AAA", "AC", "", "T", "GGATA"), c("G", "TT", "C"))
  stopifnot(identical(x, y))

  ## ---------------------------------------------------------------------
  ## B. UNLISTING AN XStringSetList OBJECT
  ## ---------------------------------------------------------------------
  length(x)
  elementLengths(x)
  unlist(x)
  x[[1]]
  x[[2]]
  as.list(x)

  names(x) <- LETTERS[1:length(x)]
  x[["A"]]
  x[["B"]]
  as.list(x)  # named list

  ## ---------------------------------------------------------------------
  ## B. USING THE GROUPING CORE API ON 'partitioning(x)'
  ## ---------------------------------------------------------------------
  partitioning(x)
  length(partitioning(x))
  nobj(partitioning(x))
  grouplength(partitioning(x))  # same as 'unname(sapply(x, length))'

  ## ---------------------------------------------------------------------
  ## C. USING THE RANGES CORE API ON 'partitioning(x)'
  ## ---------------------------------------------------------------------
  start(partitioning(x))
  end(partitioning(x))
  width(partitioning(x))  # same as 'grouplength(partitioning(x))'
}

\keyword{methods}
\keyword{classes}
