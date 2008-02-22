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

% Accesor methods:
\alias{super}
\alias{super,BStringSet-method}

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

% Other methods:
\alias{show,BStringSet-method}
\alias{[[,BStringSet-method}
\alias{[[<-,BStringSet-method}


\title{BStringSet objects}

\description{
  The BStringSet class is a container for storing a set of \code{\link{BString}}
  instances and for making its manipulation easy and efficient.
  The DNAStringSet (or RNAStringSet, or AAStringSet) class is a direct
  BStringSet derivation for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) instances.
}

\details{
  [TODO]
}

\section{Accesor methods}{
  The BStringSet class derives from the \link{SubstrLocs} class hence
  all the accessor methods defined for a \link{SubstrLocs} object can also
  be used on a BStringSet object.
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

\author{H. Pages}

\seealso{
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{BStringViews-class},
  \link{SubstrLocs}
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
