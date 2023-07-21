\name{QualityScaledXStringSet-class}
\docType{class}

% Classes:
\alias{class:QualityScaledXStringSet}
\alias{QualityScaledXStringSet-class}
\alias{QualityScaledXStringSet}

\alias{parallel_slot_names,QualityScaledXStringSet-method}

\alias{class:QualityScaledBStringSet}
\alias{QualityScaledBStringSet-class}
\alias{QualityScaledBStringSet}
\alias{class:QualityScaledDNAStringSet}
\alias{QualityScaledDNAStringSet-class}
\alias{QualityScaledDNAStringSet}
\alias{class:QualityScaledRNAStringSet}
\alias{QualityScaledRNAStringSet-class}
\alias{QualityScaledRNAStringSet}
\alias{class:QualityScaledAAStringSet}
\alias{QualityScaledAAStringSet-class}
\alias{QualityScaledAAStringSet}

% Accessor methods:
\alias{quality}
\alias{quality,QualityScaledXStringSet-method}

\alias{windows,QualityScaledXStringSet-method}
\alias{narrow,QualityScaledXStringSet-method}
\alias{reverse,QualityScaledXStringSet-method}
\alias{reverseComplement,QualityScaledDNAStringSet-method}
\alias{reverseComplement,QualityScaledRNAStringSet-method}

\alias{show,QualityScaledXStringSet-method}

\alias{readQualityScaledDNAStringSet}
\alias{writeQualityScaledXStringSet}

\title{QualityScaledBStringSet, QualityScaledDNAStringSet, QualityScaledRNAStringSet and QualityScaledAAStringSet objects}

\description{
  The QualityScaledBStringSet class is a container for storing a
  \code{\link{BStringSet}} object with an \code{\link{XStringQuality}} object.

  Similarly, the QualityScaledDNAStringSet (or QualityScaledRNAStringSet, or
  QualityScaledAAStringSet) class is a container for storing a
  \code{\link{DNAStringSet}} (or \code{\link{RNAStringSet}}, or
  \code{\link{AAStringSet}}) objects with an \code{\link{XStringQuality}}
  object.
}

\usage{
## Constructors:
QualityScaledBStringSet(x, quality)
QualityScaledDNAStringSet(x, quality)
QualityScaledRNAStringSet(x, quality)
QualityScaledAAStringSet(x, quality)

## Read/write a QualityScaledXStringSet object from/to a FASTQ file:
readQualityScaledDNAStringSet(filepath,
                quality.scoring=c("phred", "solexa", "illumina"),
                nrec=-1L, skip=0L, seek.first.rec=FALSE,
                use.names=TRUE)

writeQualityScaledXStringSet(x, filepath, append=FALSE,
                compress=FALSE, compression_level=NA)
}

\arguments{
  \item{x}{
    For the \code{QualityScaled*StringSet} constructors:
    Either a character vector, or an \link{XString}, \link{XStringSet}
    or \link{XStringViews} object.

    For \code{writeQualityScaledXStringSet}:
    A QualityScaledDNAStringSet object or other QualityScaledXStringSet
    derivative.
  }
  \item{quality}{
    An \link{XStringQuality} derivative.
  }
  \item{filepath, nrec, skip, seek.first.rec, use.names,
        append, compress, compression_level}{
    See \code{?`\link{XStringSet-io}`}.
  }
  \item{quality.scoring}{
    Specify the quality scoring used in the FASTQ file. Must be one of
    \code{"phred"} (the default), \code{"solexa"}, or \code{"illumina"}.
    If set to \code{"phred"} (or \code{"solexa"} or \code{"illumina"}),
    the qualities will be stored in a \link{PhredQuality} (or
    \link{SolexaQuality} or \link{IlluminaQuality}, respectively)
    object.
  }
}

\details{
  The \code{QualityScaledBStringSet}, \code{QualityScaledDNAStringSet},
  \code{QualityScaledRNAStringSet} and \code{QualityScaledAAStringSet}
  functions are constructors that can be used to "naturally" turn
  \code{x} into an QualityScaledXStringSet object of the desired base type.
}

\section{Accessor methods}{
  The QualityScaledXStringSet class derives from the \link{XStringSet}
  class hence all the accessor methods defined for an \link{XStringSet}
  object can also be used on an QualityScaledXStringSet object. Common
  methods include (in the code snippets below, \code{x} is an
  QualityScaledXStringSet object):

  \describe{
    \item{\code{length(x)}:}{
      The number of sequences in \code{x}.
    }
    \item{\code{width(x)}:}{
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
    }
    \item{\code{nchar(x)}:}{
      The same as \code{width(x)}.
    }
    \item{\code{names(x)}:}{
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user-provided description or comment for each
      element in \code{x}.
    }
    \item{\code{quality(x)}:}{
      The quality of the strings.
    }
  }
}

\section{Subsetting and appending}{
  In the code snippets below,
  \code{x} and \code{values} are XStringSet objects,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{\code{x[i]}:}{
      Return a new QualityScaledXStringSet object made of the selected elements.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \itemize{
    \item \link{BStringSet}, \link{DNAStringSet}, \link{RNAStringSet},
          and \link{AAStringSet} objects.

    \item \link{XStringQuality} objects.

    \item \code{\link{readDNAStringSet}} and \code{\link{writeXStringSet}}
          for reading/writing a \link{DNAStringSet} object (or other
          \link{XStringSet} derivative) from/to a FASTA or FASTQ file.
  }
}

\examples{
## ---------------------------------------------------------------------
## QualityScaled*StringSet() CONSTRUCTORS
## ---------------------------------------------------------------------

x1 <- DNAStringSet(c("TTGA", "CTCN"))
q1 <- PhredQuality(c("*+,-", "6789"))
qdna1 <- QualityScaledDNAStringSet(x1, q1)
qdna1

## ---------------------------------------------------------------------
## READ/WRITE A QualityScaledDNAStringSet OBJECT FROM/TO A FASTQ FILE
## ---------------------------------------------------------------------

filepath <- system.file("extdata", "s_1_sequence.txt",
                        package="Biostrings")

## By default, readQualityScaledDNAStringSet() assumes that the FASTQ
## file contains "Phred quality scores" (this is the standard Sanger
## variant to assess reliability of a base call):
qdna2 <- readQualityScaledDNAStringSet(filepath)
qdna2

outfile2a <- tempfile()
writeQualityScaledXStringSet(qdna2, outfile2a)

outfile2b <- tempfile()
writeQualityScaledXStringSet(qdna2, outfile2b, compress=TRUE)

## Use 'quality.scoring="solexa"' or 'quality.scoring="illumina"' if the
## quality scores are Solexa quality scores:
qdna3 <- readQualityScaledDNAStringSet(filepath, quality.scoring="solexa")
qdna3

outfile3a <- tempfile()
writeQualityScaledXStringSet(qdna3, outfile3a)

outfile3b <- tempfile()
writeQualityScaledXStringSet(qdna3, outfile3b, compress=TRUE)

## Sanity checks:
stopifnot(identical(readLines(outfile2a), readLines(filepath)))
stopifnot(identical(readLines(outfile2a), readLines(outfile2b)))
stopifnot(identical(readLines(outfile3a), readLines(filepath)))
stopifnot(identical(readLines(outfile3a), readLines(outfile3b)))
}

\keyword{methods}
\keyword{classes}
