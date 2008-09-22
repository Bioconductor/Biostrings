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

% Accesor methods:
\alias{nchar,XStringSet-method}
\alias{alphabet,XStringSet-method}

% Coercion:
\alias{coerce,XStringSet,BStringSet-method}
\alias{coerce,XStringSet,DNAStringSet-method}
\alias{coerce,XStringSet,RNAStringSet-method}
\alias{coerce,XStringSet,AAStringSet-method}

% Constructor-like functions and generics:
\alias{XStringSet,character-method}
\alias{XStringSet,AsIs-method}
\alias{XStringSet,XString-method}
\alias{XStringSet,XStringSet-method}
\alias{XStringSet,XStringViews-method}

\alias{BStringSet}
\alias{DNAStringSet}
\alias{RNAStringSet}
\alias{AAStringSet}

% "show" method:
\alias{show,XStringSet-method}

% Subsetting and appending:
\alias{[[,XStringSet-method}
\alias{[[<-,XStringSet-method}
\alias{append,XStringSet-method}

% Other methods:
\alias{as.list,XStringSet-method}
\alias{as.character,XStringSet-method}
\alias{as.matrix,XStringSet-method}


\title{BStringSet, DNAStringSet, RNAStringSet and AAStringSet objects}

\description{
  The BStringSet class is a container for storing a set of \code{\link{BString}}
  objects and for making its manipulation easy and efficient.

  Similarly, the DNAStringSet (or RNAStringSet, or AAStringSet) class is
  a container for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) objects.

  All those containers derive directly (and with no additional slots)
  from the XStringSet virtual class. They are also said to be XStringSet
  subtypes.
}

\usage{
  ## Constructors:
  BStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  DNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  RNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  AAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
}

\arguments{
  \item{x}{
    Either a character vector, or an \link{XString}, \link{XStringSet}
    or \link{XStringViews} object.
  }
  \item{start}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges]{narrow}} for the details).
  }
  \item{end}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges]{narrow}} for the details).
  }
  \item{width}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges]{narrow}} for the details).
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
}

\details{
  The \code{BStringSet}, \code{DNAStringSet}, \code{RNAStringSet} and
  \code{AAStringSet} functions are constructors that can be used to
  "naturally" turn \code{x} into an XStringSet object of the desired
  subtype.

  They also allow the user to "narrow" the sequences contained in \code{x}
  via proper use of the \code{start}, \code{end} and/or \code{width}
  arguments. In this context, "narrowing" means dropping unwanted parts
  of \code{x} located at the beginning (prefix) or end (suffix) of each
  sequence in \code{x}.

  The \code{\link[IRanges]{narrow}} function is a generic function
  (defined in the IRanges package) with a method for narrowing
  \link[IRanges]{IRanges} objects. Because XStringSet objects are
  a particular kind of \link[IRanges]{IRanges} objects (the XStringSet
  class is a subclass of the \link[IRanges]{IRanges} class), an XStringSet
  object \code{y} can be narrowed with \code{narrow(y)}.
  Therefore the two following expressions are equivalent:
  \describe{
    \item{}{
      \code{DNAStringSet(x, start=s, end=e, width=w)}
    }
    \item{}{
      \code{narrow(DNAStringSet(x), start=s, end=e, width=w)}
    }
  }
  but, besides being more convenient, the former is also more memory
  efficient on character vectors and would work even if the dropped
  parts contained letters that are not in the DNA alphabet (see
  \code{?\link{DNA_ALPHABET}}).
}

\section{Accesor methods}{
  The XStringSet class derives from the \link[IRanges]{IRanges} class
  hence all the accessor methods defined for a \link[IRanges]{IRanges}
  object can also be used on an XStringSet object. In particular, the
  following methods are available (in the code snippets below, \code{x}
  is an XStringSet object:

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

\section{Subsetting and appending}{
  In the code snippets below,
  \code{x} and \code{values} are XStringSet objects,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new XStringSet object made of the selected elements.
    }
    \item{}{
      \code{x[[i]]}:
      Extract the i-th \code{\link{XString}} object from \code{x}.
    }
    \item{}{
      \code{append(x, values, after=length(x))}:
      Add sequences in \code{values} to \code{x}.
    }
  }
}

\section{Other methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{}{
      \code{as.character(x, use.names)}:
      Convert \code{x} to a character vector of the same length as \code{x}.
      \code{use.names} controls whether or not \code{desc(x)} should be
      used to set the names of the returned vector (default is \code{TRUE}).
    }
    \item{}{
      \code{as.matrix(x, use.names)}:
      Returns a character matrix containing the "exploded" representation of
      the strings. This can only be used on an XStringSet object with
      equal-width strings.
      \code{use.names} controls whether or not \code{desc(x)} should be used
      to set the row names of the returned matrix (default is \code{TRUE}).
    }
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
  \code{\link[IRanges]{narrow}},
  \code{\link{DNA_ALPHABET}}
}

\examples{
  x0 <- c("#TTGA", "#-CTC-N")
  x1 <- DNAStringSet(x0, start=2)
  x1
  names(x1)
  names(x1)[2] <- "seqB"
  x1

  library(drosophila2probe)
  x2 <- DNAStringSet(drosophila2probe$sequence)
  x2

  RNAStringSet(x2, start=2, end=-5)  # does NOT copy the sequence data!
}

\keyword{methods}
\keyword{classes}

