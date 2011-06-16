\name{XStringViews-constructor}

% Constructor-like functions and generics:
\alias{XStringViews-constructor}
\alias{XStringViews}
\alias{XStringViews,ANY-method}
\alias{XStringViews,XString-method}
\alias{XStringViews,XStringViews-method}


\title{The XStringViews() constructor}

\description{
  WARNING: The \code{XStringViews} constructor function is deprecated!

  The \code{XStringViews} constructor function can be used for creating
  or modifying an XStringViews object.
}

\usage{
XStringViews(x, subjectClass, collapse="")
}

\arguments{
  \item{x}{
    An \link{XString} object or a character vector for \code{XStringViews}.
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
  An XStringViews object \code{y}.
  \code{length(y)} (the number of views in \code{y})
  is \code{1} when \code{x} is an \link{XString} object
  and \code{length(x)} otherwise.
}

\seealso{
  \link{XStringViews-class},
  \link{XStringSet-class},
  \link{XString-class}
}

\examples{
v12 <- Views(DNAString("TAATAATG"), start=-2:9, end=0:11)

## Do not do this anymore:
if (interactive()) {
  XStringViews(v12, subjectClass="RNAString")  # DEPRECATED!
}
## Do this instead:
xsbasetype(v12) <- "RNA"

## Do not do this anymore:
if (interactive()) {
  XStringViews(AAString("MARKSLEMSIR*"))  # DEPRECATED!
}
## Do this instead:
as(AAString("MARKSLEMSIR*"), "Views")

## Do not do this anymore:
x <- c("abcdefghij", "KLMN")
if (interactive()) {
  XStringViews(x, subjectClass="BString")  # DEPRECATED!
}
## Do this instead (if you really want views):
successiveViews(unlist(BStringSet(x)), nchar(x))
## otherwise just do:
BStringSet(x)
}

\keyword{methods}
