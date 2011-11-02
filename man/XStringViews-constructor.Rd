\name{XStringViews-constructor}

% Constructor-like functions and generics:
\alias{XStringViews-constructor}
\alias{XStringViews}
\alias{XStringViews,ANY-method}
\alias{XStringViews,XString-method}
\alias{XStringViews,XStringViews-method}


\title{The XStringViews() constructor}

\description{
  The \code{XStringViews} constructor function is now defunct.
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
  XStringViews(v12, subjectClass="RNAString")  # DEFUNCT!
}
## Do this instead:
xsbasetype(v12) <- "RNA"

## Do not do this anymore:
if (interactive()) {
  XStringViews(AAString("MARKSLEMSIR*"))  # DEFUNCT!
}
## Do this instead:
as(AAString("MARKSLEMSIR*"), "Views")

x <- c("abcdefghij", "KLMN")
## Do not do this anymore:
if (interactive()) {
  XStringViews(x, subjectClass="BString")  # DEFUNCT!
}
## Do this instead (if you really want views):
as(BStringSet(x), "Views")
## otherwise just do:
BStringSet(x)
}

\keyword{methods}
