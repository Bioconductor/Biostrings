\name{NucleotideString}
\alias{NucleotideString}
\alias{DNAString}
\title{Functions to convert a character string to a BioString object }
\description{
  The function \code{NucleotideString} expects a character string
  representing either an RNA or DNA sequence, and converts it to a
  BioString object with the appropriate alphabet.

  The function \code{DNAString} expects a character string
  representing either a DNA sequence, and converts it to a
  BioString object with the standard DNA alphabet.
}
\usage{
NucleotideString(src, type = c("DNA", "RNA"),
                 srctype = c("character", "connection"),
                 alphabet = if (type == "DNA")
                            DNAAlphabet() else RNAAlphabet(),
                 gap = alphabet@gap)
DNAString(src, gap="-")
}
\arguments{
  \item{src}{ A character string }
  \item{type}{ Either "DNA" (the default) or "RNA". }
  \item{srctype}{ Currently the only valid value (and the default) is
    "character". In future, we may allow "connection" as another valid
    value denoting that the \code{src} argument is a file name.}
  \item{alphabet}{Alphabet to be used for the string - usually one of
    DNAAlphabet() or RNAAlphabet(). }
  \item{gap}{The character in the input that represents a gap. All
    occurences of this character in the input is
    converted to the gap character of the alphabet (which is usually
    \sQuote{-}).}
}
\value{
  An object of class "BioString".
}
\author{Saikat DebRoy}
\seealso{\code{\link{BioString-class}} for the class of the return
  value.

  \code{\link{DNAAlphabet}} and \code{\link{RNAAlphabet}} alphabet for
  creating the standard DNA and RNA alphabets.
}
\examples{
NucleotideString("ACTGAACT")
DNAString("ACTGAACT-GC")
DNAString("ACTGAACTNGC", gap="N")
}
\keyword{classes}
