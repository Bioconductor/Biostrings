\name{BioPatternAlphabet-class}
\docType{class}
\alias{BioPatternAlphabet-class}
\alias{BioStringNewValues,BioPatternAlphabet,numeric-method}
\alias{gapletter<-,BioPatternAlphabet,character-method}
\alias{initialize,BioPatternAlphabet-method}

\title{Class representing alphabets used in patterns for matching
  biological sequences }
\description{ Each object of class "BioPatternAlphabet" represents an
  alphabet for patterns that is used to match a specific type of
  biological sequence.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("BioPatternAlphabet", baseAlphabet, letters)} where
  \code{letters} is a named character vector with the names being single
  letters not in the base alphabet and the value corresponding to
  each name are strings made of letters from the base alphabet. Each
  such name value pair defines a new letter in the pattern alphabet
  which matches all the letters in its value string.
}
\section{Slots}{
  \describe{
    \item{\code{baseAlphabet}:}{Object of class \code{"BioAlphabet"},
      the base alphabet that for the letters matched by members of this
      pattern alphabet. }
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
\section{Extends}{
Class \code{"BioAlphabet"}, directly.
}
\section{Methods}{
  \describe{
    \item{BioStringNewValues(alphabet, length.string)}{Given \code{alphabet}
      of class "BioPatternAlphabet" return a new uninitialized
      "BioString" object for that alphabet with length
      \code{length.string}. See \code{\link{BioString-class}} for more
      details.}
    \item{gapletter(x) <- value}{Change the letter representing gaps in
      \code{x} of class "BioPatternAlphabet" to \code{value} which must
      be a single letter not in the base alphabet of \code{x}.}
    \item{initialize(.Object, ...)}{Method to initialize an object of
      class "BioPatternAlphabet". Usually used indirectly via a call to
      \code{\link{new}}}
  }
}
\author{ Saikat DebRoy }
\seealso{
  \code{\link{BioAlphabet-class}}
}
\examples{
dnaAlph <- new("BioPatternAlphabet",
               DNAAlphabet(), c(N="AGCT",
                                B="CGT",
                                D="AGT",
                                H="ACT",
                                K="GT",
                                M="AC",
                                R="AG",
                                S="CG",
                                V="ACG",
                                W="AT",
                                Y="CT"))
dnaAlph
}
\keyword{classes}
