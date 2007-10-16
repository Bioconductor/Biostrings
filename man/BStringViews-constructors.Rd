\name{BStringViews-constructors}

% Constructor-like functions and generics:
\alias{BStringViews-constructors}
\alias{views}
\alias{adjacentViews}
\alias{BStringViews}
\alias{BStringViews,ANY-method}
\alias{BStringViews,file-method}
\alias{BStringViews,BString-method}
\alias{BStringViews,BStringViews-method}


\title{Basic functions for creating BStringViews objects}

\description{
  A set of basic functions for creating BStringViews objects.
}

\usage{
  views(subject, start=NA, end=NA)
  adjacentViews(subject, width, gapwidth=0)
  BStringViews(src, subjectClass, collapse="")
}

\arguments{
  \item{subject}{
    A \link{BString} (or derived) object or a single string.
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
  \item{src}{
    A \link{BString} (or derived) object or a character vector.
  }
  \item{subjectClass}{
    The class to be given to the subject of the \link{BStringViews} object
    created and returned by the function.
    Must be \code{"BString"} or one of its subclasses: \code{"DNAString"},
    \code{"RNAString"} or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string to be inserted between the views of the
    \link{BStringViews} object created and returned by the function.
  }

}

\details{
  The \code{views} function returns a BStringViews object
  containing views on \code{subject} with start and end positions
  given in the \code{start} and \code{end} vectors.

  The \code{adjacentViews} function returns a BStringViews object
  containing views on \code{subject} with widths given in the \code{width}
  vector and separated by gaps of width \code{gapwidth}.
  The first view starts at position 1.

  The \code{BStringViews} constructor will try to create a BStringViews
  object from the value passed to its \code{src} argument.
  If \code{src} itself is a BStringViews object, the returned
  object is obtained by coercing its subject to the class specified by
  \code{subjectClass}.
  If \code{src} is a \link{BString} object, the returned object is made of
  a single view that starts at the first letter and ends at the last letter
  of \code{src} (in addition \code{src} itself is coerced to the class
  specified by \code{subjectClass} when specified).
  If \code{src} is a character vector, the returned object has one view
  per character string in \code{src} (and its subject is an instance of
  the class specified by \code{subjectClass}).
}

\value{
  All these functions return a BStringViews object \code{y} with the following
  number of views (\code{length(y)}):
  \code{max(length(start), length(end)} for the \code{views} function
  and \code{length(width)} for the \code{adjacentViews} function.
  For the \code{BStringViews} constructor, \code{length(y)} is \code{1}
  when \code{src} is a \link{BString} object and \code{length(src)} otherwise.
}

\seealso{
  \link{BStringViews-class}, \link{BString}, \link{DNAString}, \link{RNAString}
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

  BStringViews(v12, subjectClass="RNAString")
  BStringViews(AAString("MARKSLEMSIR*"))
  BStringViews("abcdefghij", subjectClass="BString")
}

\keyword{methods}
