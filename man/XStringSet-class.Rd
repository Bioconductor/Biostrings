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

% Methods:
\alias{length,XStringSet-method}
\alias{width,character-method}
\alias{width,XStringSet-method}
\alias{names,XStringSet-method}
\alias{names<-,XStringSet-method}
\alias{nchar,XStringSet-method}
\alias{narrow,character-method}
\alias{narrow,XStringSet-method}
\alias{threebands,character-method}
\alias{threebands,XStringSet-method}

\alias{XStringSet,character-method}
\alias{XStringSet,AsIs-method}
\alias{XStringSet,XString-method}
\alias{XStringSet,XStringSet-method}
\alias{XStringSet,XStringViews-method}

\alias{coerce,XStringSet,BStringSet-method}
\alias{coerce,XStringSet,DNAStringSet-method}
\alias{coerce,XStringSet,RNAStringSet-method}
\alias{coerce,XStringSet,AAStringSet-method}
\alias{coerce,character,BStringSet-method}
\alias{coerce,character,DNAStringSet-method}
\alias{coerce,character,RNAStringSet-method}
\alias{coerce,character,AAStringSet-method}

\alias{show,XStringSet-method}

\alias{[,XStringSet-method}
\alias{rep,XStringSet-method}
\alias{[[,XStringSet-method}
\alias{[[<-,XStringSet-method}

\alias{append,XStringSet,XStringSet-method}

\alias{union,XStringSet,XStringSet-method}
\alias{intersect,XStringSet,XStringSet-method}
\alias{setdiff,XStringSet,XStringSet-method}
\alias{setequal,XStringSet,XStringSet-method}

\alias{as.character,XStringSet-method}
\alias{toString,XStringSet-method}
\alias{as.matrix,XStringSet-method}

\alias{order,XStringSet-method}
\alias{sort,XStringSet-method}
\alias{rank,XStringSet-method}

\alias{duplicated,XStringSet-method}
\alias{unique,XStringSet-method}


\title{BStringSet, DNAStringSet, RNAStringSet and AAStringSet objects}

\description{
  The BStringSet class is a container for storing a set of \code{\link{BString}}
  objects and for making its manipulation easy and efficient.

  Similarly, the DNAStringSet (or RNAStringSet, or AAStringSet) class is
  a container for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) objects.

  All those containers derive directly (and with no additional slots)
  from the XStringSet virtual class.
}

\usage{
  ## Constructors:
  BStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  DNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  RNAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)
  AAStringSet(x, start=NA, end=NA, width=NA, use.names=TRUE)

  ## Accessor-like methods:
  \S4method{length}{XStringSet}(x)
  \S4method{width}{character}(x)
  \S4method{width}{XStringSet}(x)
  \S4method{names}{XStringSet}(x)
  \S4method{nchar}{XStringSet}(x, type="chars", allowNA=FALSE)

  ## Narrowing:
  \S4method{narrow}{character}(x, start=NA, end=NA, width=NA, use.names=TRUE)
  \S4method{narrow}{XStringSet}(x, start=NA, end=NA, width=NA, use.names=TRUE)

  ## ... and more (see below)
}

\arguments{
  \item{x}{
    Either a character vector (with no NAs), or an \link{XString},
    XStringSet or \link{XStringViews} object.
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
  \item{type,allowNA}{
    Ignored.
  }
}

\details{
  The \code{BStringSet}, \code{DNAStringSet}, \code{RNAStringSet} and
  \code{AAStringSet} functions are constructors that can be used to
  "naturally" turn \code{x} into an XStringSet object of the desired
  base type.

  They also allow the user to "narrow" the sequences contained in \code{x}
  via proper use of the \code{start}, \code{end} and/or \code{width}
  arguments. In this context, "narrowing" means dropping unwanted parts
  of \code{x} located at the beginning (prefix) or end (suffix) of each
  sequence in \code{x}.
  The "narrowing" capabilities of these constructors can be illustrated
  by the following property: if \code{x} is a character vector
  (with no NAs), or an XStringSet or \link{XStringViews} object,
  then the 3 following transformations are equivalent:
  \describe{
    \item{}{
      \code{BStringSet(x, start=mystart, end=myend, width=mywidth)}
    }
    \item{}{
      \code{narrow(BStringSet(x), start=mystart, end=myend, width=mywidth)}
    }
    \item{}{
      \code{BStringSet(narrow(x, start=mystart, end=myend, width=mywidth))}
    }
  }
  Note that, besides being more convenient, the first form is also more
  efficient on character vectors.
}

\section{Accessor-like methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of sequences in \code{x}.
    }
    \item{}{
      \code{width(x)}:
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
      Note that \code{width(x)} is also defined for a character vector
      with no NAs and is equivalent to \code{nchar(x, type="bytes")}.
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
    \item{}{
      \code{alphabet(x)}:
      Return \code{NULL}, \code{\link{DNA_ALPHABET}}, \code{\link{RNA_ALPHABET}} or
      \code{\link{AA_ALPHABET}} depending on whether \code{x} is a BStringSet,
      DNAStringSet, RNAStringSet or AAStringSet object.
    }
    \item{}{
      \code{nchar(x)}:
      The same as \code{width(x)}.
    }
  }
}

\section{Narrowing and related transformations}{
  In the code snippets below,
  \code{x} is a character vector (with no NAs),
  or an XStringSet (or \link{XStringViews}) object.

  \describe{
    \item{}{
      \code{narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)}:
      Narrows the set of sequences in \code{x} by extracting the set of
      subsequences described by the supplied start/end/width values.
      Like with the \code{\link[IRanges]{narrow}} method for
      \link[IRanges]{IRanges} objects, the supplied start/end/width values
      are solved internally by calling
      \code{solveUserSEW(width(x), start=start, end=end, width=width)}
      and therefore must be compliant with the rules of the SEW
      (Start/End/Width) interface (see \code{?solveUserSEW} for the details).

      Note that this narrowing transformation is similar to what
      \code{\link{substr}} does on a character vector. However some
      noticeable differences are: (1) the names of the arguments (\code{start}
      and \code{stop} for \code{\link{substr}}), (2) the richer SEW interface
      for \code{narrow} (e.g. support for negative start or end values),
      and (3) strict checking of valid start/end/width values for \code{narrow}
      (unlike \code{\link{substr}}, \code{narrow} will throw an error if
      the specified start/end/width values define "out of limits" subsequences
      or subsequences with a negative width).
    }
    \item{}{
      \code{threebands(x, start=NA, end=NA, width=NA)}:
      Like the method for \link[IRanges]{IRanges} objects, the
      \code{threebands} methods for character vectors and XStringSet
      objects extend the capability of \code{narrow} by returning the 3
      set of subsequences (the left, middle and right subsequences)
      associated to the narrowing operation.
      See \code{?\link[IRanges]{threebands}} in the IRanges package for
      the details.
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

\section{Ordering and related methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{}{
      \code{order(x)}:
      Return a permutation which rearranges \code{x} into ascending or
      descending order.
    }
    \item{}{
      \code{sort(x)}:
      Sort \code{x} into ascending order (equivalent to \code{x[order(x)]}).
    }
    \item{}{
      \code{rank(x)}:
      Rank \code{x} in ascending order.
    }
  }
}

\section{Duplicated and unique methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{}{
      \code{duplicated(x)}:
      Return a logical vector whose elements denotes duplicates in \code{x}.
    }
    \item{}{
      \code{unique(x)}:
      Return an XStringSet containing the unique values in \code{x}.
    }
  }
}

\section{Set Operations}{
  In the code snippets below,
  \code{x} and \code{y} are XStringSet objects

  \describe{
    \item{}{
      \code{union(x, y)}:
      Union of \code{x} and \code{y}.
    }
    \item{}{
      \code{intersect(x, y)}:
      Intersection of \code{x} and \code{y}.
    }
    \item{}{
      \code{setdiff(x, y)}:
      Asymmetric set difference of \code{x} and \code{y}.
    }
    \item{}{
     \code{setequal(x, y)}:
      Set equality of \code{x} to \code{y}.
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
      \code{use.names} controls whether or not \code{names(x)} should be
      used to set the names of the returned vector (default is \code{TRUE}).
    }
    \item{}{
      \code{as.matrix(x, use.names)}:
      Return a character matrix containing the "exploded" representation of
      the strings. This can only be used on an XStringSet object with
      equal-width strings.
      \code{use.names} controls whether or not \code{names(x)} should be used
      to set the row names of the returned matrix (default is \code{TRUE}).
    }
    \item{}{
      \code{toString(x)}:
      Equivalent to \code{toString(as.character(x))}.
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
  \code{\link{substr}},
  \code{\link[IRanges]{narrow}}
}

\examples{
  x0 <- c("#CTC-NACCAGTAT", "#TTGA", "TACCTAGAG")
  width(x0)
  x1 <- BStringSet(x0)
  x1

  ## 3 equivalent ways to obtain the same BStringSet object:
  BStringSet(x0, start=4, end=-3)
  narrow(x1, start=4, end=-3)
  BStringSet(narrow(x0, start=4, end=-3))

  dna0 <- DNAStringSet(x0, start=4, end=-3)
  dna0
  names(dna0)
  names(dna0)[2] <- "seqB"
  dna0

  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe$sequence)
  probes

  RNAStringSet(probes, start=2, end=-5)  # does NOT copy the sequence data!
}

\keyword{methods}
\keyword{classes}
