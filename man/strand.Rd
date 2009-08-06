\name{strand}
\docType{methods}
\alias{strand-methods}

\alias{strand}
\alias{strand,missing-method}
\alias{strand,character-method}
\alias{strand,DataTable-method}

\title{Accessing strand information}

\description{The \code{strand} generic is meant as an accessor for
  strand information. Three methods are defined by the \code{BSgenome}
  package, described below.
}

\usage{strand(x)}

\arguments{
  \item{x}{The object from which to obtain a strand factor, can be
    missing.}
}

\details{
  If \code{x} is missing, returns an empty factor with the
  standard levels that any strand factor should have: \code{+},
  \code{-}, and \code{*} (for either).

  If \code{x} is a \code{character} vector, \code{x} is coerced to
  a factor with the levels listed above.

  If \code{x} inherits from \code{DataTable}, the \code{"strand"}
  column is returned as a factor with the levels listed above. If
  \code{x} has no \code{"strand"} column, this return value is
  populated with \code{NA}s.
}

\author{ Michael Lawrence }

\examples{
strand()
strand(c("+", "-", NA, "*"))
}

\keyword{methods}