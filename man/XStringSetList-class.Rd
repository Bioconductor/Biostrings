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
  Concrete flavours of the XStringSetList container are the BStringSetList,
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
  \link[IRanges]{IPartitioning-class},
  \link[IRanges]{ListLike-class}
}

\examples{
  unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
  xsl <- new("DNAStringSetList",
             unlisted=unlisted,
             end=c(0L, 2L, 2L, 3L))
  xsl
  length(xsl)
  width(xsl)  # the element lenths
  end(xsl)  # the element cumulated lengths
  unlist(xsl)
  xsl[[1]]
  xsl[[2]]
  as.list(xsl)

  ## More examples to be added soon...
}

\keyword{methods}
\keyword{classes}
