\name{XStringViews-constructors}

% Constructor-like functions and generics:
\alias{XStringViews-constructors}
\alias{adjacentViews}
\alias{XStringViews}
\alias{XStringViews,ANY-method}
\alias{XStringViews,XString-method}
\alias{XStringViews,XStringViews-method}

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
  adjacentViews(subject, width, gapwidth=0)
  XStringViews(x, subjectClass, collapse="")
}

\arguments{
  \item{subject}{
    An \link{XString} object or a single string.
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
    Must be the name of one of the direct XString subclasses i.e.
    \code{"BString"}, \code{"DNAString"}, \code{"RNAString"}
    or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string to be inserted between the views of the
    \link{XStringViews} object created and returned by the function.
  }
}

\details{
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
  These functions return an XStringViews object \code{y}.
  \code{length(y)} (the number of views in \code{y})
  is \code{length(width)} for the \code{adjacentViews} function.
  For the \code{XStringViews} constructor, \code{length(y)} is \code{1}
  when \code{x} is an \link{XString} object and \code{length(x)} otherwise.
}

\seealso{
  \link{XStringViews-class},
  \link{XString-class}
}

\examples{
  adjacentViews("abcdefghij", 4:2, gapwidth=1)

  v12 <- Views(DNAString("TAATAATG"), start=-2:9, end=0:11)
  XStringViews(v12, subjectClass="RNAString")
  XStringViews(AAString("MARKSLEMSIR*"))
  XStringViews("abcdefghij", subjectClass="BString")
}

\keyword{methods}
