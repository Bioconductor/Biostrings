\name{BioString-class}
\docType{class}
\alias{BioString-class}
\alias{[,BioString-method}
\alias{[[,BioString-method}
\alias{allSameLetter,BioString,character-method}
\alias{allSameLetter,BioString,BioString-method}
\alias{as.character,BioString-method}
\alias{as.matrix,BioString-method}
\alias{initialize,BioString-method}
\alias{length,BioString-method}
\alias{matchDNAPattern,BioString,BioString-method}
\alias{nchar,BioString-method}
\alias{show,BioString-method}
\alias{substr,BioString-method}
\alias{substring,BioString-method}

\title{Class "BioString", represents a biological sequence}
\description{Class "BioString", contains an encoded string representing a
  biological sequence for a particular alphabet (RNA, DNA or amino
  acid). It represents zero or more substrings of the full string.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("BioString", alphabet, end, start, values, initialized, ...)}.
  However, it is recommended that users should not call this directly.
  For now, use the function \code{\link{NucleotideString}} to create
  objects of class "BioString" that uses a nucleotide alphabet (RNA or
  DNA) and the function \code{\link{DNAString}} for objects using DNA
  alphabet.}
}
\section{Slots}{
  \describe{
    \item{\code{alphabet}:}{Object of class \code{"BioAlphabet"},
      the alphabet used in the sequence. }
    \item{\code{initialized}:}{Object of class \code{"logical"},
      \code{TRUE} if the sequence initialized with values. Users should
      not modify this slot directly. }
    \item{\code{offsets}:}{Object of class \code{"matrix"} and storage
      mode "integer", this stores (in two columns) the start and end
      points of the substrings represented in \code{x}. Rows with the
      first value \code{1} and the second value{0} represent empty
      substrings.}
    \item{\code{values}:}{Object of class \code{"externalptr"}, this
      internally stores the actual encoded sequence as a vector. As
      objects of class "externalptr" are passed by value in R, this
      saves copying of long sequences. }
  }
}
\section{Methods}{
  \describe{
    \item{initialize(.Object, alphabet,
      offsets=cbind(1, 0),
      values=BioStringNewValues(alphabet, end),
      initialized=!missing(values))}{Construct an object of class
      "BioString". Usually not called directly by users. }
    \item{length(x)}{Return the number of substrings represented by
      \code{x}.}
    \item{x[i]}{Return the substrings in \code{x} corresponding to index
      \code{i}.}
    \item{x[[i]]}{Return the substring in \code{x} corresponding to the
      index \code{i}. The index \code{i} must be of length \code{1}.}
    \item{nchar(x)}{Return the number of characters in each substring
      represented in \code{x}.}
    \item{show(object)}{Display \code{object} of class "BioString".}
    \item{as.character(x)}{Convert a "BioString" object to a character
      vector using its native alphabet.}
    \item{as.matrix(x)}{Return a two-column matrix of integers, the
      first column representing the start index and the scond column
      representing the end index of the substrings in the full
      string.}
    \item{substr(x, start, stop)}{Return another BioString object with
      value equivalent to \code{substr(as.character(x), start, stop)}.}
    \item{substring(text, first, last)}{Return another BioString object with
      value equivalent to
      \code{substring(as.character(text), first, last)}.}
    \item{matchDNAPattern(pattern, x, algorithm)}{Match the DNA string
      \code{x} against \code{pattern} using \code{algorithm}. Currently
      only a modified Boyer-Moore algorithm is implemented. The pattern
      can use the letters A,C,G,T,- (the last being the gap character)
      and also the wildcards N (matching A,C,G,T), V (matching A,G,C),
      R (matching A,G) and Y (matching C,T).}
    \item{allSameLetter(x, letter)}{Return a logical vetor indicating
      which of the elements of \code{x} are entirely made up of the
      letter \code{letter}.}
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
  letters (including gap) in our alphabets for all standard computer
  architectures.
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
x
as.character(x)
substr(x, 2, 4)
substring(x, 1, seq(length=nchar(x))) # all prefixes of x
substring(x, seq(length=nchar(x)), nchar(x)) # all suffixes of x
matchDNAPattern("GC", x)
x <- substring(x, 1:3, 3:5)
x[1:2]
x[-3] # same as x[1:2]
x[[3]]
}
\keyword{classes}
