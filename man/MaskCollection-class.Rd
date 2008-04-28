\name{MaskCollection-class}
\docType{class}

% Classes:
\alias{class:MaskCollection}
\alias{MaskCollection-class}
\alias{MaskCollection}

% Basic accessor methods:
\alias{nirlist}
\alias{nirlist,MaskCollection-method}
\alias{length,MaskCollection-method}
\alias{width,MaskCollection-method}
\alias{names,MaskCollection-method}

% Other methods:
\alias{isEmpty,MaskCollection-method}
\alias{max,MaskCollection-method}
\alias{min,MaskCollection-method}
\alias{maskedwidth}
\alias{maskedwidth,MaskCollection-method}
\alias{maskedratio}
\alias{maskedratio,MaskCollection-method}

% Subsetting:
\alias{[[,MaskCollection-method}
\alias{[[<-,MaskCollection-method}
\alias{[,MaskCollection-method}
\alias{[<-,MaskCollection-method}

% "append" method:
\alias{append,MaskCollection-method}

% Transforming methods:
\alias{narrow,MaskCollection-method}
\alias{reduce,MaskCollection-method}
\alias{gaps,MaskCollection-method}

% Coercion:
\alias{coerce,MaskCollection,NormalIRanges-method}

% "show" method :
\alias{show,MaskCollection-method}


\title{MaskCollection objects}

\description{
  The MaskCollection class is a container for storing a collection of masks
  that can be used to mask regions in a sequence.
}

\details{
  In the context of the Biostrings package, a mask is a set of regions
  in a sequence that need to be excluded from some computation.
  For example, when calling \code{\link{alphabetFrequency}}
  or \code{\link{matchPattern}} on a chromosome sequence, you might want
  to exclude some regions like the centromere or the repeat regions.
  This can be achieved by putting one or several masks on the sequence
  before \code{\link{alphabetFrequency}} is called on it.

  A MaskCollection object is a vector-like object that represents
  such set of masks.
  Like standard R vectors, it has a "length" which is the number of
  masks contained in it. But unlike standard R vectors, it also has
  a "width" which determines the length of the sequences it can be "put on".
  For example, a MaskCollection object of width 20000 can only be put on
  an XString object of 20000 letters.

  Each mask in a MaskCollection object \code{x} is just a finite set of
  integers that are >= 1 and <= \code{width(x)}.
  When "put on" a sequence, these integers indicate the positions of the
  letters to mask.
  Internally, each mask is represented by a \link{NormalIRanges} object.
}

\section{Basic accesor methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of masks in \code{x}.
    }
    \item{}{
      \code{width(x)}:
      The common with of all the masks in \code{x}.
      This determines the length of the sequences that \code{x} can be
      "put on".
    }
    \item{}{
      \code{names(x)}:
      \code{NULL} or a character vector of the same length as \code{x}.
    }
    \item{}{
      \code{nirlist(x)}:
      A list of the same length as \code{x}, where each element is
      a \link{NormalIRanges} object representing a mask in \code{x}.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{isEmpty(x)}:
      Return a logical vector of the same length as \code{x}, indicating,
      for each mask in \code{x}, whether it's empty or not.
    }
    \item{}{
      \code{max(x)}:
      The greatest (or last, or rightmost) masked position for each mask.
      This is a numeric vector of the same length as \code{x}.
    }
    \item{}{
      \code{min(x)}:
      The smallest (or first, or leftmost) masked position for each mask.
      This is a numeric vector of the same length as \code{x}.
    }
    \item{}{
      \code{maskedwidth(x)}:
      The number of masked position for each mask.
      This is an integer vector of the same length as \code{x} where
      all values are >= 0 and <= \code{width(x)}.
    }
    \item{}{
      \code{maskedratio(x)}:
      \code{maskedwidth(x) / width(x)}
    }
  }
}

\section{Subsetting and appending}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new MaskCollection object made of the selected masks.
      \code{i} can be a numeric vector, a logical vector, \code{NULL}
      or missing.
    }
    \item{}{
      \code{append(x, values, after=length(x))}:
      Add masks to \code{x}.
    }
    \item{}{
      \code{x[[i]]}:
      Extract the i-th mask as a \link{NormalIRanges} object.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)}:
      Narrow the masks in \code{x}.
    }
    \item{}{
      \code{reduce(x)}:
      Return a MaskCollection object of length 1 made of the union
      (or merging, or collapsing) of all the masks in \code{x}.
    }
    \item{}{
      \code{gaps(x)}:
      Invert the masks in \code{x}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{MaskedXString-class},
  \code{\link{alphabetFrequency}},
  \code{\link{matchPattern}},
  \link{NormalIRanges-class}
}

\examples{
  ### Making a MaskCollection object:
  nir1 <- new("NormalIRanges", start=c(11, 25, 28), width=c(5, 2, 2))
  nir2 <- new("NormalIRanges", start=c(3, 10, 27), width=c(5, 8, 1))
  nir3 <- new("NormalIRanges", start=c(7, 12), width=c(2, 4))
  mymasks <- new("MaskCollection", nirlist=list(nir1, nir2, nir3), width=29L)
  mymasks
  length(mymasks)
  width(mymasks)

  ### Putting a MaskCollection object on a sequence:
  x <- DNAString("ACTCAAACTAGGATAGNACTTAACCAGGC")
  x
  length(x)  # same as width(mymasks)
  nchar(x)   # same as length(x)
  masks(x) <- mymasks
  x
  length(x)  # has not changed
  nchar(x)   # has changed
  alphabetFrequency(x)

  ### Removing the masks:
  masks(x) <- NULL
  x
  alphabetFrequency(x)

  ### Some operations on the MaskCollection object itself:
  reduce(mymasks)
  mymasks[[2]]
  length(mymasks[[2]])
  mymasks[[2]][-3]
  append(mymasks[-2], gaps(mymasks[2]))
  mymasks2 <- narrow(mymasks, start=8)
  mymasks2
  mymasks2[[2]]
}

\keyword{methods}
\keyword{classes}
