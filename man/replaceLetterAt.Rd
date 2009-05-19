\name{replaceLetterAt}

\alias{replaceLetterAt}
\alias{replaceLetterAt,DNAString-method}
\alias{replaceLetterAt,DNAStringSet-method}
\alias{.inplaceReplaceLetterAt}

% Deprecated:
\alias{replaceLetterAtLoc}


\title{Replacing letters in a sequence (or set of sequences)
       at some specified locations}

\description{
  \code{replaceLetterAt} first makes a copy of a sequence (or set of
  sequences) and then replaces some of the original letters by new
  letters at the specified locations.

  \code{.inplaceReplaceLetterAt} is the IN PLACE version of
  \code{replaceLetterAt}: it will modify the original
  sequence in place i.e. without copying it first.
  Note that in place modification of a sequence is fundamentally
  dangerous because it alters all objects defined in your session
  that make reference to the modified sequence.
  NEVER use \code{.inplaceReplaceLetterAt}, unless you know what
  you are doing!
}

\usage{
  replaceLetterAt(x, at, letter, if.not.extending="replace", verbose=FALSE)

  ## NEVER USE THIS FUNCTION!
  .inplaceReplaceLetterAt(x, at, letter)
}

\arguments{
  \item{x}{
    A \link{DNAString} or rectangular \link{DNAStringSet} object.
  }
  \item{at}{
    The locations where the replacements must occur.

    If \code{x} is a \link{DNAString} object, then \code{at} is typically
    an integer vector with no NAs but a logical vector or \link[IRanges]{Rle}
    object is valid too. Locations can be repeated and in this case the
    last replacement to occur at a given location prevails.

    If \code{x} is a rectangular \link{DNAStringSet} object, then
    \code{at} must be a matrix of logicals with the same dimensions
    as \code{x}.
  }
  \item{letter}{
    The new letters.

    If \code{x} is a \link{DNAString} object, then \code{letter} must be
    a \link{DNAString} object or a character vector (with no NAs) with a
    total number of letters (\code{sum(nchar(letter))}) equal to the number
    of locations specified in \code{at}.

    If \code{x} is a rectangular \link{DNAStringSet} object, then
    \code{letter} must be a \link{DNAStringSet} object or a character
    vector of the same length as \code{x}. In addition, the number of
    letters in each element of \code{letter} must match the number of
    locations specified in the corresponding row of \code{at}
    (\code{all(width(letter) == rowSums(at))}).
  }
  \item{if.not.extending}{
    What to do if the new letter is not "extending" the old letter?
    The new letter "extends" the old letter if both are IUPAC letters
    and the new letter is as specific or less specific than the old
    one (e.g. M extends A, Y extends Y, but Y doesn't extend S).
    Possible values are \code{"replace"} (the default) for replacing
    in all cases, \code{"skip"} for not replacing when the new letter
    does not extend the old letter, \code{"merge"} for merging the new
    IUPAC letter with the old one, and \code{"error"} for raising
    an error.

    Note that the gap (\code{"-"}) and hard masking (\code{"+"})
    letters are not extending or extended by any other letter.

    Also note that \code{"merge"} is the only value for the
    \code{if.not.extending} argument that guarantees the final result
    to be independent on the order the replacement is performed (although
    this is only relevant when \code{at} contains duplicated locations,
    otherwise the result is of course always independent on the order,
    whatever the value of \code{if.not.extending} is).
  }
  \item{verbose}{
    When \code{TRUE}, a warning will report the number of skipped or merged
    letters.
  }
}

\details{
  \code{.inplaceReplaceLetterAt} semantic is equivalent to calling
  \code{replaceLetterAt} with \code{if.not.extending="merge"}
  and \code{verbose=FALSE}.

  Never use \code{.inplaceReplaceLetterAt}!
  It is used by the \code{\link[BSgenome]{injectSNPs}} function
  in the BSgenome package, as part of the "lazy sequence loading"
  mechanism, for altering the original sequences of a
  \link[BSgenome]{BSgenome} object at "sequence-load time".
  This alteration consists in injecting the IUPAC ambiguity letters
  representing the SNPs into the just loaded sequence, which is the
  only time where in place modification of the external data of
  an \link{XString} object is safe.
}

\value{
  A \link{DNAString} or \link{DNAStringSet} object of the same shape
  (i.e. length and width) as the orignal object \code{x}
  for \code{replaceLetterAt}.
}

\author{H. Pages}

\seealso{
  \code{\link{IUPAC_CODE_MAP}},
  \code{\link{chartr}},
  \code{\link{injectHardMask}},
  \link{DNAString},
  \link{DNAStringSet},
  \code{\link[BSgenome]{injectSNPs}},
  \link[BSgenome]{BSgenome}
}

\examples{
  ## Replace letters of a DNAString object:
  replaceLetterAt(DNAString("AAMAA"), c(5, 1, 3, 1), "TYNC")
  replaceLetterAt(DNAString("AAMAA"), c(5, 1, 3, 1), "TYNC", if.not.extending="merge")

  ## Replace letters of a DNAStringSet object (sorry for the totally
  ## artificial example with absolutely no biological meaning):
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe$sequence)
  at <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
               nrow=length(probes), ncol=width(probes)[1],
               byrow=TRUE)
  letter_subject <- DNAString(paste(rep.int("-", width(probes)[1]), collapse=""))
  letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
  replaceLetterAt(probes, at, letter)
}

\keyword{utilities}
\keyword{manip}
