\name{BioString-class}
\docType{class}
\alias{BioString-class}
\alias{as.character,BioString-method}
\alias{initialize,BioString-method}
\alias{matchDNAPattern,BioString,BioString-method}
\alias{substr,BioString-method}

\title{Class "BioString", represents a biological sequence}
\description{Class "BioString", contains an encoded string representing a
  biological sequence for a particular alphabet (RNA, DNA or amino acid). }
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("BioString", alphabet, end, start, values, initialized, ...)}.
  However, it is recommended that users should not call this directly.
  For now, use the function \code{\link{NucleotideString}} to create
  objects of class "BioString" that uses a nucleotide alphabet (RNA or DNA).
}
\section{Slots}{
  \describe{
    \item{\code{alphabet}:}{Object of class \code{"BioAlphabet"},
      the alphabet used in the sequence. }
    \item{\code{initialized}:}{Object of class \code{"logical"},
      \code{TRUE} if the sequence initialized with values. Users should
      not modify this slot directly. }
    \item{\code{start}:}{Object of class \code{"integer"}, start of the
      sequence stored in the \code{values} slot - this uses zero-based
      indexing.}
    \item{\code{end}:}{Object of class \code{"integer"}, zero-based
      index of the position after the last value - so for the full
      sequence (for which \code{start} is \code{0}) this is same as the
      length of the sequence. }
    \item{\code{values}:}{Object of class \code{"externalptr"}, this
      internally stores the actual encoded sequence as a vector. As
      objects of class "externalptr" are passed by value in R, this
      saves copying of long sequences. }
  }
}
\section{Methods}{
  \describe{
    \item{initialize(.Object, alphabet,
      end = 0,
      start = 0,
      values=BioStringNewValues(alphabet, end),
      initialized=!missing(values))}{Construct an object of class
      "BioString". Usually not called directly by users. }
    \item{as.character(x)}{Convert a "BioString" object to a character
      string using its native alphabet.}
    \item{substr(x, start, stop)}{Return another BioString object which
      is the substring of \code{x} starting at \code{start} and ending
      at \code{stop}. If \code{start} is less than \code{0}, it is
      increased to \code{1}. If \code{stop} is more than the length of
      \code{x}, it is reduced to the length of \code{x}. If \code{start}
      is more than \code{end}, an empty string is returned.}
    \item{matchDNAPattern(pattern, x, algorithm)}{Match the DNA string
      \code{x} against \code{pattern} using \code{algorithm}. Currently
      only a modified Boyer-Moore algorithm is implemented. The pattern
      can use the letters A,C,G,T,- and also the wildcards N (matching
      A,C,G,T), V (matching A,G,C), R (matching A,G) and Y
      (matching C,T).}
  }
}
\author{Saikat DebRoy}
\section{The structure of the values slot}{
  The \code{values} slot of the "BioString" class is of class
  "externalptr". It always contains an R vector object in its tag
  field. The other fields are not used at present. The vector in the tag
  field is either a \code{CHARSXP} or an \code{INTSXP}. The exact type
  depends on the length of the alphabet. \code{INTSXP} is used if it is
  more than the number of bits in a C \code{char} type and
  \code{CHARSXP} is used otherwise.

  We use the \code{i}-th bit in the \code{char} or \code{int} (depending
  on whether the vector is of type CHARSXP or INTSXP) to represent the
  \code{i}-th letter in the alphabet where \code{i=0} represents the
  first bit. This effectively means that we can have at most \code{32}
  letters (including gap) in our alphabets for almost all standard
  machines.
}
\seealso{
  \code{\link{BioAlphabet-class}} and its subclasses for valid alphabet
  objects.
  \code{\link{DNAString}} for creating objects of class "BioString"
  representing DNA sequences.
  \code{\link{NucleotideString}} for creating objects of class "BioString"
  representing DNA or RNA sequences.
}
\examples{
new("BioString", DNAAlphabet()) # creates an empty DNA string
x <- DNAString("AAGCTANA", gap="N")
as.character(x)
y <- substr(x, 2, 4)
as.character(y)
matchDNAPattern("GC", x)
}
\keyword{classes}
