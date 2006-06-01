\name{reverseComplement}

\alias{reverse}
\alias{reverse,BString-method}
\alias{complement}
\alias{complement,DNAString-method}
\alias{reverseComplement}
\alias{reverseComplement,DNAString-method}


\title{Sequence reversing and complementing}

\description{
  These functions can reverse a \link{BString}, \link{DNAString}
  or \link{RNAString} object and complement each base
  of a \link{DNAString} object.
}

\usage{
  reverse(x)
  complement(x)
  reverseComplement(x)
}

\arguments{
  \item{x}{
    A \link{BString}, \link{DNAString} or \link{RNAString} object
    for \code{reverse}.
    A \link{DNAString} object only for \code{complement}
    and \code{reverseComplement}.
  }
}

\details{
  Given an object \code{x} of class \link{BString}, \link{DNAString}
  or \link{RNAString}, \code{reverse(x)} returns an object of the same class
  where letters in \code{x} are reordered in the reverse ordered.
  If \code{x} is a \link{DNAString} object, \code{complement(x)} returns
  an object where each base in \code{x} is "complemented" i.e.
  A, C, G, T are replaced by T, G, C, A respectively.
  Letters belonging to the "IUPAC extended genetic alphabet"
  are also replaced by their complement (M <-> K, R <-> Y, S <-> S, V <-> B,
  W <-> W, H <-> D, N <-> N) and the gap symbol (-) is unchanged.
  \code{reverseComplement(x)} is equivalent to \code{reverse(complement(x))}.
}

\value{
  An object of the same class and length as the original object.
}

\author{H. Pages}

\seealso{}

\examples{
  reverseComplement(DNAString("ACGT-YN-"))
}

\keyword{manip}

