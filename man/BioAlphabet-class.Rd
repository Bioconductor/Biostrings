\name{BioAlphabet-class}
\docType{class}
\alias{BioAlphabet-class}
\alias{initialize,BioAlphabet-method}

\title{Class "BioAlphabet" represents alphabets used in biology }
\description{ Each object of class "BioAlphabet" represents a particular
  alphabet used to represent some biological sequence.}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{letters}:}{Object of class \code{"character"}
      representing the letters of the alphabet. Usually roman
      upper case letters are used. The one exception is the gap
      character which is always \code{'-'}}
    \item{\code{mapping}:}{Object of class \code{"integer"} representing
      the encoding used to represent these objects internally. }
    \item{\code{gap}:}{Object of class \code{"character"}, the gap
      character in the alphabet. }
  }
}
\section{Methods}{
  \describe{
    \item{initialize(.Object, letters)}{Initialize an alphabet
      object. for the given letters. If one of the letters is not the
      gap character \code{'-'}, then it is added as the first letter.
      An encoding also created which is stored in the \code{mapping}
      slot. }
  }
}
\author{ Saikat DebRoy }
\seealso{
  \code{\link{NucleotideAlphabet-class}} and
  \code{\link{AminoAcidAlphabet-class}}, for two subclasses of the
  \code{"BioAlphabet"} class.
}
\examples{
new("NucleotideAlphabet", letters=c('A', 'G', 'C', 'T'))
new("NucleotideAlphabet", letters=c('-', 'A', 'G', 'C', 'T'))
}
\keyword{classes}
