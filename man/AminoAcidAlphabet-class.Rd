\name{AminoAcidAlphabet-class}
\docType{class}
\alias{AminoAcidAlphabet-class}
\alias{BioStringNewValues,AminoAcidAlphabet,numeric-method}

\title{Class "NucleotideAlphabet" represents alphabets used in DNA or RNA}
\description{ Each object of class "NucleotideAlphabet" represents an
  amino acid alphabet. This usually also contains the gap character
  represented by \sQuote{-}.}
\section{Slots}{
  \describe{
    \item{\code{letters}:}{Object of class \code{"character"}
      representing the letters of the alphabet. Usually roman
      upper case letters are used. The one exception is the gap
      character which is always \sQuote{-}}
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
      of class "AminoAcidAlphabet" return a new uninitialized
      BioString object for that alphabet with length
      \code{length.string}. See \code{\link{BioString-class}} for more
      details.}
  }
}
\author{ Saikat DebRoy }
%\examples{ FIXME: add example
%}
\keyword{classes}
