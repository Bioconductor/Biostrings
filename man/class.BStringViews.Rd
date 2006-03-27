\name{BStringViews}
\docType{class}

\alias{class:BStringViews}
\alias{BStringViews-class}

\alias{length,BStringViews-method}
\alias{toString,BStringViews-method}
\alias{initialize,BStringViews-method}
\alias{show,BStringViews-method}
\alias{==,BStringViews,BString-method}
\alias{==,BString,BStringViews-method}
\alias{==,BStringViews,ANY-method}
\alias{==,ANY,BStringViews-method}
\alias{!=,BStringViews,BString-method}
\alias{!=,BString,BStringViews-method}
\alias{!=,BStringViews,ANY-method}
\alias{!=,ANY,BStringViews-method}
\alias{[,BStringViews-method}
\alias{[<-,BStringViews-method}
\alias{[[,BStringViews-method}
\alias{[[<-,BStringViews-method}

\alias{BStringViews}
\alias{BStringViews,file-method}
\alias{BStringViews,BString-method}
\alias{BStringViews,BStringViews-method}

\alias{subject}
\alias{first}
\alias{last}
\alias{width}
\alias{desc}

\title{A container for a set of views on the same BString
(or DNAString, or RNAString) object}

\description{

}

\details{

}

\section{Methods}{
  Standard generic methods:
  \describe{
    \item{\code{length(BStringViews)}:}{Used to get the length of
      a BStringViews, i.e., its number of views}
    \item{\code{initialize(BStringViews)}:}{[TODO: Document me]}
    \item{\code{show(BStringViews)}:}{[TODO: Document me]}
    \item{\code{BStringViews1 == BStringViews2}:}{[TODO: Document me]}
    \item{\code{BStringViews1 != BStringViews2}:}{[TODO: Document me]}
    \item{\code{BStringViews[i]}:}{[TODO: Document me]}
    \item{\code{BStringViews[i] <- "attempt to modify the value of a BStringViews"}:}{[TODO: Document me]}
    \item{\code{BStringViews[[i]]}:}{[TODO: Document me]}
    \item{\code{BStringViews[[i]] <- "attempt to modify the value of a BStringViews"}:}{[TODO: Document me]}
  }
  Class-specific methods:
  \describe{
    \item{\code{subject(BStringViews)}:}{[TODO: Document me]}
    \item{\code{first(BStringViews)}:}{[TODO: Document me]}
    \item{\code{last(BStringViews)}:}{[TODO: Document me]}
    \item{\code{width(BStringViews)}:}{[TODO: Document me]}
    \item{\code{desc(BStringViews)}:}{[TODO: Document me]}
  }
}

\author{H. Pages}

\seealso{
  \link{BString, DNAString, RNAString, views, adjacentViews}
}

\examples{
  ## A BStringViews object with 4 views
  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
}

\keyword{methods}
\keyword{classes}
