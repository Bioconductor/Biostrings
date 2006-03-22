\name{NucleotideAlphabet-class}
\docType{class}
\alias{NucleotideAlphabet-class}
\alias{BioStringNewValues,NucleotideAlphabet,numeric-method}
\alias{initialize,NucleotideAlphabet-method}

\title{Class "NucleotideAlphabet" represents alphabets used in DNA or RNA}
\description{ Each object of class "NucleotideAlphabet" represents a
  nucleotide alphabet. Usually this is the RNA or DNA alphabet and
  consists of \sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{T} for DNA and \sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{U} for RNA
  along with \sQuote{-} which represents a gap in the sequence.}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{letters}:}{Object of class \code{"character"}
      representing the letters of the alphabet. Usually roman
      upper case letters are used. The one exception is the gap
      character which is always \sQuote{-}.}
    \item{\code{mapping}:}{Object of class \code{"integer"} representing
      the encoding used to represent these objects internally. }
  }
}
\section{Extends}{
Class \code{"BioAlphabet"}, directly.
}
\section{Methods}{
  \describe{
    \item{BioStringNewValues(alphabet, length.string)}{Given \code{alphabet}
      of class "NucleotideAlphabet" return a new uninitialized
      "BioString" object for that alphabet with length
      \code{length.string}. See \code{\link{BioString-class}} for more
      details.}
    \item{initialize(.Object, letters)}{Initialize an nucleotide alphabet
      object for the given letters. If one of the letters is not the
      gap character \sQuote{-}, then it is added as the first letter.
      With the gap character, the length of the alphabet should be five.
      An encoding also created which is stored in the \code{mapping}
      slot. }
  }
}
\author{ Saikat DebRoy }
\examples{
myDNAAlph <- new("NucleotideAlphabet", letters=c('A', 'G', 'C', 'T'))
myDNAAlph
BioStringNewValues(myDNAAlph, 20)
}
\keyword{classes}
