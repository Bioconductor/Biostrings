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
\alias{length,XStringSetList-method}
\alias{unlist,XStringSetList-method}
\alias{cum_elementLengths}
\alias{cum_elementLengths,XStringSetList-method}
\alias{elementLengths,XStringSetList-method}
\alias{start,XStringSetList-method}
\alias{end,XStringSetList-method}
\alias{width,XStringSetList-method}
\alias{names<-,XStringSetList-method}

\alias{xsbasetype,XStringSetList-method}
\alias{xsbasetype<-,XStringSetList-method}
\alias{coerce,XStringSetList,IRanges-method}
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
  \link{ListLike-class}
}

\examples{
  unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
  xsl <- new("DNAStringSetList",
             unlisted=unlisted,
             cum_eltlength=c(0L, 2L, 2L, 3L))
  xsl
  length(xsl)
  elementLengths(xsl)
  cum_elementLengths(xsl)
  unlist(xsl)
  xsl[[1]]
  xsl[[2]]
  as.list(xsl)

  ## More examples to be added soon...
}

\keyword{methods}
\keyword{classes}
