\name{MultipleAlignment-class}
\docType{class}

% Classes:
\alias{class:MultipleAlignment}
\alias{MultipleAlignment-class}
\alias{MultipleAlignment}
\alias{class:DNAMultipleAlignment}
\alias{DNAMultipleAlignment-class}
\alias{DNAMultipleAlignment}
\alias{class:RNAMultipleAlignment}
\alias{RNAMultipleAlignment-class}
\alias{RNAMultipleAlignment}
\alias{class:AAMultipleAlignment}
\alias{AAMultipleAlignment-class}
\alias{AAMultipleAlignment}

% Accessor-like methods:
\alias{unmasked,MultipleAlignment-method}
\alias{rownames,MultipleAlignment-method}
\alias{rownames<-,MultipleAlignment-method}
\alias{rowmask}
\alias{rowmask,MultipleAlignment-method}
\alias{rowmask<-}
\alias{rowmask<-,MultipleAlignment,NULL-method}
\alias{rowmask<-,MultipleAlignment,NormalIRanges-method}
\alias{rowmask<-,MultipleAlignment,IRanges-method}
\alias{colmask}
\alias{colmask,MultipleAlignment-method}
\alias{colmask<-}
\alias{colmask<-,MultipleAlignment,NULL-method}
\alias{colmask<-,MultipleAlignment,NormalIRanges-method}
\alias{colmask<-,MultipleAlignment,IRanges-method}
\alias{nrow,MultipleAlignment-method}
\alias{ncol,MultipleAlignment-method}
\alias{dim,MultipleAlignment-method}
\alias{maskednrow}
\alias{maskednrow,MultipleAlignment-method}
\alias{maskedncol}
\alias{maskedncol,MultipleAlignment-method}
\alias{maskeddim}
\alias{maskeddim,MultipleAlignment-method}
\alias{maskedratio,MultipleAlignment-method}
\alias{nchar,MultipleAlignment-method}
\alias{xsbasetype,MultipleAlignment-method}

% Read method:
\alias{read.DNAMultipleAlignment}

% Coercion:
\alias{coerce,MultipleAlignment,DNAStringSet-method}
\alias{coerce,MultipleAlignment,RNAStringSet-method}
\alias{coerce,MultipleAlignment,AAStringSet-method}
\alias{coerce,MultipleAlignment,BStringSet-method}
\alias{coerce,character,DNAMultipleAlignment-method}
\alias{coerce,character,RNAMultipleAlignment-method}
\alias{coerce,character,AAMultipleAlignment-method}
\alias{as.character,MultipleAlignment-method}
\alias{as.matrix,MultipleAlignment-method}

% Utilities:
\alias{consensusMatrix,MultipleAlignment-method}
\alias{consensusString,MultipleAlignment-method}
\alias{alphabetFrequency,MultipleAlignment-method}

% show method:
\alias{show,MultipleAlignment-method}


\title{MultipleAlignment objects}

\description{
   The MultipleAlignment class is a container for storing multiple sequence
   alignments.
}

\usage{
  ## Constructors:
  DNAMultipleAlignment(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
  RNAMultipleAlignment(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
  AAMultipleAlignment(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)

  ## Read method:
  read.DNAMultipleAlignment(filepath, format=c("fasta", "clustalw"))

  ## ... and more (see below)
}

\arguments{
  \item{x}{
    Either a character vector (with no NAs), or an \link{XString},
    XStringSet or \link{XStringViews} object.
  }
  \item{start,end,width}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges:Ranges-utils]{narrow}} for the details).
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
  \item{filepath}{
    A character vector (of arbitrary length when reading, of length 1
    when writing) containing the paths to the files to read or write.
    Note that special values like \code{""} or \code{"|cmd"} (typically
    supported by other I/O functions in R) are not supported here. Also
    \code{filepath} cannot be a connection.
  }
  \item{format}{
    File format.
  }
}

\details{
  The MultipleAlignment class is designed to hold and represent multiple
  sequence alignments. The rows and columns within an alignment can be
  masked for ad hoc analyses.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a MultipleAlignment object.

  \describe{
    \item{}{
      \code{unmasked(x)}:
      The underlying \link{XStringSet} object containing the multiple
      sequence alignment.
    }
    \item{}{
      \code{rownames(x)}:
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user-provided description or comment for each
      sequence in \code{x}.
    }
    \item{}{
      \code{rowmask(x)}, \code{rowmask(x) <- value}:
      Gets and sets the \link{NormalIRanges} object representing the
      masked rows in \code{x}.
    }
    \item{}{
      \code{colmask(x)}, \code{colmask(x) <- value}:
      Gets and sets the \link{NormalIRanges} object representing the
      masked columns in \code{x}.
    }
    \item{}{
      \code{nrow(x)}:
      Returns the number of sequences aligned in \code{x}.
    }
    \item{}{
      \code{ncol(x)}:
      Returns the number of characters for each alignment in \code{x}.
    }
    \item{}{
      \code{dim(x)}:
      Equivalent to \code{c(nrow(x), ncol(x))}.
    }
    \item{}{
      \code{maskednrow(x)}:
      Returns the number of masked aligned sequences in \code{x}.
    }
    \item{}{
      \code{maskedncol(x)}:
      Returns the number of masked aligned characters in \code{x}.
    }
    \item{}{
      \code{maskeddim(x)}:
      Equivalent to \code{c(maskednrow(x), maskedncol(x))}.
    }
    \item{}{
      \code{maskedratio(x)}:
      Equivalent to \code{maskeddim(x) / dim(x)}.
    }
    \item{}{
      \code{nchar(x)}:
      Returns the number of unmasked aligned characters in \code{x},
      i.e. \code{ncol(x) - maskedncol(x)}.
    }
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unmasked(x))}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a MultipleAlignment object.

  \describe{
    \item{}{
      \code{as(from, "DNAStringSet")}, \code{as(from, "RNAStringSet")},
      \code{as(from, "AAStringSet")}, \code{as(from, "BStringSet")}:
      Creates an instance of the specified \link{XStringSet} object subtype
      that contains the unmasked regions of the multiple sequence alignment
      in \code{x}.
    }
    \item{}{
      \code{as.character(x, use.names)}:
      Convert \code{x} to a character vector containing the unmasked
      regions of the multiple sequence alignment. \code{use.names}
      controls whether or not \code{rownames(x)} should be used to set
      the names of the returned vector (default is \code{TRUE}).
    }
    \item{}{
      \code{as.matrix(x, use.names)}:
      Returns a character matrix containing the "exploded" representation
      of the unmasked regions of the multiple sequence alignment.
      \code{use.names} controls whether or not \code{rownames(x)} should
      be used to set the row names of the returned matrix (default is
      \code{TRUE}).
    }
  }
}

\section{Utilities}{
  In the code snippets below, x is a MultipleAlignment object.

  \describe{
    \item{}{
      consensusMatrix(x, as.prob, baseOnly):
      Creates an integer matrix containing the column frequencies of
      the underlying alphabet. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{baseOnly} is \code{TRUE}, then the non-base
      letters are collapsed into an \code{"other"} category.
    }
    \item{}{
      consensusString(x, ...):
      Creates a consensus string for \code{x}.
    }
    \item{}{
      alphabetFrequency(x, as.prob, collapse):
      Creates an integer matrix containing the row frequencies of
      the underlying alphabet. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{collapse} is \code{TRUE}, then returns the
      overall frequency instead of the frequency by row.
    }
  }
}


\author{M. Carlson and P. Aboyoun}

\seealso{
  \link{XStringSet-class},
  \link{MaskedXString-class}
}

\examples{
## Create an object from file
multAlign <-
  read.DNAMultipleAlignment(filepath =
                            system.file("extdata","clustalw2.aln",
                                        package="Biostrings"),
                            format="clustalw")

## mask columns
modified <- multAlign
colmask(modified) <- IRanges(start=c(1,20), end=c(4,30))
modified

## mask rows
modified <- multAlign
rowmask(modified) <- IRanges(start=1, end=2)
modified

## remove row masks
rowmask(modified) <- NULL
modified

## List the names of the sequences in the alignment
rownames(multAlign)
}

\keyword{methods}
\keyword{classes}
