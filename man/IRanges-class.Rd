\name{IRanges-class}
\docType{class}

% Classes:
\alias{class:IRanges}
\alias{IRanges-class}
\alias{IRanges}
\alias{class:UnlockedIRanges}
\alias{UnlockedIRanges-class}
\alias{UnlockedIRanges}
\alias{class:Views}
\alias{Views-class}
\alias{Views}
\alias{class:LockedIRanges}
\alias{LockedIRanges-class}
\alias{LockedIRanges}
\alias{class:NormalIRanges}
\alias{NormalIRanges-class}
\alias{NormalIRanges}

% Accesor methods:
\alias{length,IRanges-method}
\alias{start,IRanges-method}
\alias{width}
\alias{width,IRanges-method}
\alias{end,IRanges-method}
\alias{names,IRanges-method}
\alias{desc}

\alias{isNormal}
\alias{isNormal,IRanges-method}
\alias{whichFirstNotNormal}
\alias{whichFirstNotNormal,IRanges-method}
\alias{isEmpty}
\alias{isEmpty,IRanges-method}
\alias{isEmpty,NormalIRanges-method}
\alias{max,NormalIRanges-method}
\alias{min,NormalIRanges-method}

% Initialization:
\alias{initialize,IRanges-method}
\alias{initialize,Views-method}
\alias{initialize,NormalIRanges-method}

% Coercion:
\alias{coerce,IRanges,NormalIRanges-method}
\alias{coerce,UnlockedIRanges,NormalIRanges-method}
\alias{coerce,Views,NormalIRanges-method}
\alias{coerce,LockedIRanges,NormalIRanges-method}

% "show" method:
\alias{as.data.frame,IRanges-method}
\alias{show,IRanges-method}

% Subsetting:
\alias{[,IRanges-method}
\alias{[<-,IRanges-method}

% "duplicated" method:
\alias{duplicated,IRanges-method}

% Other methods:
\alias{as.matrix,IRanges-method}

% Replacement methods:
\alias{start<-}
\alias{start<-,UnlockedIRanges-method}
\alias{width<-}
\alias{width<-,UnlockedIRanges-method}
\alias{end<-}
\alias{end<-,UnlockedIRanges-method}
\alias{names<-,IRanges-method}
\alias{desc<-}
\alias{update,UnlockedIRanges-method}

% Deprecated methods:
\alias{first}
\alias{first,IRanges-method}
\alias{last}
\alias{last,IRanges-method}


\title{IRanges and NormalIRanges objects}

\description{
  The IRanges class is a simple container for storing a set of integer
  ranges.

  A NormalIRanges object is an IRanges object that is "normal".
  See the Normality section below for the definition and properties
  of normal IRanges objects.
}

\details{
  An IRanges object is a data frame-like object where each row describes
  a "range" of integers.

  A "range" of integers is a finite set of consecutive integer values.
  Each range can be fully described with exactly 2 integers which can be
  arbitrarily picked up among the 3 following integers:
  its "start" i.e. its smallest (or first, or leftmost) value;
  its "end" i.e. its greatest (or last, or rightmost) value;
  and its "width" i.e. the number of values in the range.
  For example the set of integers that are greater than or equal to -20
  and less than or equal to 400 is the range of integers that starts
  at -20 and has a width of 421.

  The start can be any integer (see \code{start} below) but the
  width must be a nonnegative integer (see \code{width} below).
  The end of a range is its start plus its width minus one (see
  \code{end} below).
  An "empty" range is a range that contains no value i.e. a range that
  has a null width.
  Note that for an empty range, the end is smaller than the start.

  Two ranges are considered equal iff they share the same start and width.
  Note that with this definition, 2 empty ranges are generally
  not equal (they need to share the same start to be considered equal).

  The length of an IRanges object is the number of ranges in it
  i.e. the number of rows in the object.

  An IRanges object is considered empty iff all its ranges are empty.

  Note that it is unlikely that the user will have to create or manipulate
  directly an IRanges instance when using the Biostrings package. However
  the IRanges class being a superclass of the \link{XStringViews}
  class, any \link{XStringViews} object is also an IRanges object and
  can be manipulated as such. Therefore all the methods described here also
  work with an \link{XStringViews} object.
}

\section{IRanges object vs data frame}{
  An important difference with standard R data frames is that
  IRanges objects only support single subscript subsetting i.e.
  subsetting by row, whereas standard R data frames can be subsetted
  by row and by column.
  As a consequence, the length of an IRanges object is its number of rows,
  whereas the length of a standard R data frame object is its number of
  columns.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is an IRanges object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of ranges in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The start values of the ranges.
      This is a vector of positive integers of the same length as \code{x}.
    }
    \item{}{
      \code{width(x)}:
      The number of integers in each range.
      This is a vector of nonnegative integers of the same length as \code{x}.
    }
    \item{}{
      \code{end(x)}:
      \code{start(x) + width(x) - 1L}
    }
    \item{}{
      \code{names(x)}:
      \code{NULL} or a character vector of the same length as \code{x}.
    }
    \item{}{
      \code{desc(x)}:
      \code{desc} is an alias for \code{names}.
    }
  }
}

\section{Subsetting}{
  In the code snippet below, \code{x} is an IRanges object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new IRanges object (of the same type as \code{x})
      made of the selected ranges.
      \code{i} can be a numeric vector, a logical vector, \code{NULL}
      or missing. If \code{x} is a NormalIRanges object and \code{i}
      a positive numeric subscript (i.e. a numeric vector of positive
      values), then \code{i} must be strictly increasing.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is an IRanges object.

  \describe{
    \item{}{
      \code{isEmpty(x)}:
      Return a logical value indicating whether \code{x} is empty or not.
    }
    \item{}{
      \code{as.data.frame(x, row.names=NULL, optional=FALSE, ...)}:
      Converts \code{x} into a standard R data frame object.
      \code{row.names} must be \code{NULL} or a character vector giving
      the row names for the data frame, and \code{optional} and any
      additional argument (\code{...}) is ignored.
      See \code{?\link{as.data.frame}} for more information about these
      arguments.
    }
    \item{}{
      \code{duplicated(x)}:
      Determines which elements of \code{x} are equal to elements
      with smaller subscripts, and returns a logical vector indicating
      which elements are duplicates.
      It is semantically equivalent to \code{duplicated(as.data.frame(x))}
      (see \code{?\link{duplicated}} for more information).
    }
    \item{}{
      \code{as.matrix(x, ...)}:
      Converts \code{x} into a 2-column integer matrix containing
      \code{start(x)} and \code{width(x)}.
      Extra arguments (\code{...}) are ignored.
    }
  }
}

\section{Normality}{
  A NormalIRanges object is an IRanges object that is "normal".

  An IRanges object is said to be "normal" when its ranges are:
    (a) not empty (i.e. they have a non-null width);
    (b) not overlapping;
    (c) ordered from left to right;
    (d) not even adjacent (i.e. there must be a non empty gap between 2
        consecutive ranges).
  If \code{x} is an IRanges object with more than one element (i.e.
  \code{length >= 2}), then \code{x} is normal iff:
  \preformatted{  start(x)[i] <= end(x)[i] < start(x)[i+1] <= end(x)[i+1]}
  for every 1 <= \code{i} < \code{length(x)}.
  If \code{length(x) == 1}, then \code{x} is normal iff \code{width(x)[1] >= 1}.
  If \code{length(x) == 0}, then \code{x} is normal.

  An IRanges object can be used to represent an arbitrary finite set of
  integers (that are not necessarily consecutive).
  Now the 2 most interesting properties of normal IRanges objects are that:
  (1) they are the "best" (in terms of storage space) IRanges objects for
  representing arbitrary finite sets of integers and (2) the mapping between
  finite sets of integers and normal IRanges objects is one-to-one.
  More precisely, if \code{x} is an IRanges object, then it can be seen
  as representing the set of integers obtained by taking the union of
  all its ranges.
  Inverserly, since any finite set of integers can be obtained by a finite
  union of ranges, then it can be represented by an IRanges object, but this
  representation is clearly not unique.
  However, among all the IRanges objects that represent (or map) the
  same finite set of integers, only one is normal, and this normal
  representation is minimal in terms of length (and therefore in terms
  of storage space).

  Subsetting \code{x} is currently not supported although it could be but
  should only accept strictly increasing subscripts in order to preserve
  normality.

  Use the \code{isNormal} method to check whether an IRanges object is
  normal or not. In the code snippet below, \code{x} is an IRanges object.

  \describe{
    \item{}{
      \code{isNormal(x)}:
      Return a logical value indicating whether \code{x} is normal or not.
    }
    \item{}{
      \code{whichFirstNotNormal(x)}:
      Return \code{NA} if \code{x} is normal, or the smallest valid indice
      \code{i} in \code{x} for which \code{x[1:i]} is not normal.
    }
    \item{}{
      \code{max(x)}:
      (Defined for NormalIRanges objects only.)
      The maximum value in the finite set of integers represented by \code{x}.
    }
    \item{}{
      \code{min(x)}:
      (Defined for NormalIRanges objects only.)
      The minimum value in the finite set of integers represented by \code{x}.
    }
  }
}

\section{Deprecated methods}{
    \item{}{
      \code{first(x)}:
      deprecated. Use \code{start} instead.
    }
    \item{}{
      \code{last(x)}:
      deprecated. Use \code{end} instead.
    }
}

\author{H. Pages}

\seealso{
  \link{IRanges-utils},
  \link{XStringViews-class},
  \code{\link{as.data.frame}},
  \code{\link{duplicated}},
  \code{\link{as.matrix}}
}

\examples{
  x <- new("IRanges", start=c(2:-1, 13:15), width=c(0:3, 2:0))
  x
  length(x)
  start(x)
  width(x)
  end(x)
  isEmpty(x)
  as.data.frame(x)
  as.matrix(x)

  ## Subsetting:
  x[4:2]                  # 3 ranges
  x[-1]                   # 6 ranges
  x[FALSE]                # 0 range
  x0 <- x[width(x) == 0]  # 2 ranges
  isEmpty(x0)

  ## Unlock the IRanges instance and use replacement methods to slide
  ## or resize its elements:
  x <- as(x, "UnlockedIRanges")
  width(x) <- width(x) * 2  + 1  # resize elements
  x
  start(x) <- end(x)             # slide elements
  x
  start(x)[4] <- end(x)[4]       # slide the 4th element
  x
  end(x)[1] <- start(x)[3]       # slide the first element
  x
  width(x) <- c(2, 0)            # resize elements
  x
  duplicated(x)

  ## Name the elements:
  names(x)
  names(x) <- c("range1", "range2")
  x
  x[names(x) == ""]  # 5 ranges
  x[names(x) != ""]  # 2 ranges

  ## Using an IRanges object for storing a big set of ranges is more
  ## efficient than using a standard R data frame:
  N <- 2000000L  # nb of ranges
  W <- 180L      # width of each range
  start <- 1L
  end <- 50000000L
  set.seed(777)
  range_starts <- sort(sample(end-W+1L, N))
  range_widths <- rep.int(W, N)
  ## Instantiation is faster
  system.time(x <- new("IRanges", start=range_starts, width=range_widths))
  system.time(y <- data.frame(start=range_starts, width=range_widths))
  ## Subsetting is faster
  system.time(x16 <- x[c(TRUE, rep.int(FALSE, 15))])
  system.time(y16 <- y[c(TRUE, rep.int(FALSE, 15)), ])
  ## Internal representation is more compact
  object.size(x16)
  object.size(y16)

  ## Normality:
  isNormal(x16)                        # FALSE
  if (interactive())
      x16 <- as(x16, "NormalIRanges")  # Error!
  whichFirstNotNormal(x16)             # 57
  isNormal(x16[1:56])                  # TRUE
  xx <- as(x16[1:56], "NormalIRanges")
  class(xx)
  max(xx)
  min(xx)
}

\keyword{methods}
\keyword{classes}
