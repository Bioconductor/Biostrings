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
\alias{unlist,XStringSetList-method}
\alias{partitioning}
\alias{partitioning,XStringSetList-method}
\alias{length,XStringSetList-method}
\alias{names,XStringSetList-method}
\alias{names<-,XStringSetList-method}

\alias{xsbasetype,XStringSetList-method}
\alias{xsbasetype<-,XStringSetList-method}

\alias{show,XStringSetList-method}

\alias{[[,XStringSetList-method}
\alias{[[<-,XStringSetList-method}


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
}

\section{Methods}{
  [TODO]
}

\author{H. Pages}

\seealso{
  \link{XStringSet-class},
  \link[IRanges]{Grouping-class},
  \link[IRanges]{Sequence-class}
}

\examples{
  unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
  partitioning <- PartitioningByEnd(c(0, 2, 2, 3))
  x <- new("DNAStringSetList",
           unlisted=unlisted,
           partitioning=partitioning)
  x
  length(x)
  unlist(x)
  x[[1]]
  x[[2]]
  as.list(x)

  names(x) <- LETTERS[1:4]
  x[["A"]]
  x[["B"]]
  as.list(x)  # named list

  ## Using the Grouping core API on 'partitioning(x)':
  partitioning(x)
  length(partitioning(x))
  nobj(partitioning(x))
  grouplength(partitioning(x))  # same as 'unname(sapply(x, length))'

  ## Using the Ranges core API on 'partitioning(x)':
  start(partitioning(x))
  end(partitioning(x))
  width(partitioning(x))  # same as 'grouplength(partitioning(x))'
}

\keyword{methods}
\keyword{classes}
