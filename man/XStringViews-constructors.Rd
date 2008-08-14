\name{XStringViews-constructors}

% Constructor-like functions and generics:
\alias{XStringViews-constructors}
\alias{views}
\alias{views,ANY-method}
\alias{views,XString-method}
\alias{adjacentViews}
\alias{XStringViews}
\alias{XStringViews,ANY-method}
\alias{XStringViews,file-method}
\alias{XStringViews,XString-method}
\alias{XStringViews,XStringViews-method}

% Other functions:
\alias{restrict,XStringViews-method}
\alias{trim}
\alias{narrow,XStringViews-method}
\alias{subviews}
\alias{gaps,XStringViews-method}

% Deprecated:
\alias{BStringViews}
\alias{BStringViews,ANY-method}
\alias{BStringViews,file-method}
\alias{BStringViews,XString-method}
\alias{BStringViews,XStringViews-method}


\title{Basic functions for creating or modifying XStringViews objects}

\description{
  A set of basic functions for creating or modifying XStringViews objects.
}

\usage{
  views(subject, start=NA, end=NA)
  adjacentViews(subject, width, gapwidth=0)
  XStringViews(x, subjectClass, collapse="")
  trim(x, use.names=TRUE)
  subviews(x, start=NA, end=NA, width=NA, use.names=TRUE)
}

\arguments{
  \item{subject}{
    An \link{XString} object or a single string.
  }
  \item{start}{
    An integer vector containing the start positions of the views.
  }
  \item{end}{
    An integer vector containing the end positions of the views.
  }
  \item{width}{
    An integer vector containing the widths of the views.
  }
  \item{gapwidth}{
    An integer vector containing the widths of the gaps between the views.
  }
  \item{x}{
    An \link{XString} object or a character vector for \code{XStringViews}.
    An XStringViews object for \code{trim} and \code{subviews}.
  }
  \item{subjectClass}{
    The class to be given to the subject of the \link{XStringViews} object
    created and returned by the function.
    Must be the name of one of the direct XString subtypes i.e.
    \code{"BString"}, \code{"DNAString"}, \code{"RNAString"}
    or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string to be inserted between the views of the
    \link{XStringViews} object created and returned by the function.
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
}

\details{
  The \code{views} function returns an XStringViews object
  containing views on \code{subject} with start and end positions
  given in the \code{start} and \code{end} vectors.

  The \code{adjacentViews} function returns an XStringViews object
  containing views on \code{subject} with widths given in the \code{width}
  vector and separated by gaps of width \code{gapwidth}.
  The first view starts at position 1.

  The \code{XStringViews} constructor will try to create an XStringViews
  object from the value passed to its \code{x} argument.
  If \code{x} itself is an XStringViews object, the returned
  object is obtained by coercing its subject to the class specified by
  \code{subjectClass}.
  If \code{x} is an \link{XString} object, the returned object is made of
  a single view that starts at the first letter and ends at the last letter
  of \code{x} (in addition \code{x} itself is coerced to the class
  specified by \code{subjectClass} when specified).
  If \code{x} is a character vector, the returned object has one view
  per character string in \code{x} (and its subject is an instance of
  the class specified by \code{subjectClass}).
}

\value{
  All these functions return an XStringViews object \code{y} with the following
  number of views (\code{length(y)}):
  \code{max(length(start), length(end)} for the \code{views} function
  and \code{length(width)} for the \code{adjacentViews} function.
  For the \code{XStringViews} constructor, \code{length(y)} is \code{1}
  when \code{x} is an \link{XString} object and \code{length(x)} otherwise.
}

\seealso{
  \link{XStringViews-class},
  \link{XString-class},
  \code{\link[IRanges]{gaps}}
}

\examples{
  ## 'start' and 'end' are recycled
  views("abcdefghij", 2:1, 4)
  views("abcdefghij", 5:7, )
  views("abcdefghij", , 5:7)

  ## Views can be "out of limits"
  views("abcdefghij", 2:0, 4)

  ## Views on a DNAString object
  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)

  ## Creating adjacent views
  v3 <- adjacentViews("abcdefghij", 4:2, gapwidth=1)

  XStringViews(v12, subjectClass="RNAString")
  XStringViews(AAString("MARKSLEMSIR*"))
  XStringViews("abcdefghij", subjectClass="BString")

  ## Applying gaps() to an XStringViews object
  v2 <- views("abCDefgHIJK", start=c(8, 3), end=c(14, 4))
  gaps(v2)
}

\keyword{methods}
