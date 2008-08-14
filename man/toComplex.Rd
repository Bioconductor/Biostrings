\name{toComplex}

\alias{toComplex}
\alias{toComplex,DNAString-method}

\alias{as.complex,DNAString-method}


\title{Turning a DNA sequence into a vector of complex numbers}

\description{
  The \code{toComplex} utility function turns a \link{DNAString} object
  into a complex vector.
}

\usage{
  toComplex(x, baseValues)
}

\arguments{
  \item{x}{
    A \link{DNAString} object.
  }
  \item{baseValues}{
    A named complex vector containing the values associated
    to each base e.g.
      \code{c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)}
  }
}

\value{
  A complex vector of the same length as \code{x}.
}

\author{H. Pages}

\seealso{
  \link{DNAString}
}

\examples{
  seq <- DNAString("accacctgaccattgtcct")
  baseValues1 <- c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
  toComplex(seq, baseValues1)

  ## GC content:
  baseValues2 <- c(A=0, C=1, G=1, T=0)
  sum(as.integer(toComplex(seq, baseValues2)))
  ## Note that there are better ways to do this (see ?alphabetFrequency)
}

\keyword{methods}
