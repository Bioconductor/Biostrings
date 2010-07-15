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
\alias{automask}
\alias{automask,MultipleAlignment}
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

% Read functions:
\alias{read.DNAMultipleAlignment}
\alias{read.RNAMultipleAlignment}
\alias{read.AAMultipleAlignment}

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
\alias{consensusString,DNAMultipleAlignment-method}
\alias{consensusString,RNAMultipleAlignment-method}
\alias{consensusString,AAMultipleAlignment-method}
\alias{consensusViews}
\alias{consensusViews,MultipleAlignment-method}
\alias{consensusViews,DNAMultipleAlignment-method}
\alias{consensusViews,RNAMultipleAlignment-method}
\alias{consensusViews,AAMultipleAlignment-method}
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

  ## Read functions:
  read.DNAMultipleAlignment(filepath, format)
  read.RNAMultipleAlignment(filepath, format)
  read.AAMultipleAlignment(filepath, format)

  ## ... and more (see below)
}

\arguments{
  \item{x}{
    Either a character vector (with no NAs), or an \link{XString},
    \link{XStringSet} or \link{XStringViews} object containing
    strings with the same number of characters.
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
    Either \code{"fasta"} (the default), \code{stockholm}, or
    \code{"clustal"}.
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
      \code{automask(x, max.gappercent, max.gaplength, append)}:
      Returns a MultipleAlignment object with a modified column mask
      based upon gaps in the columns. In particular, this mask is defined
      by \code{max.gaplength} or more consecutive columns that have more
      than \code{max.gappercent} of their non-masked rows containing
      gap codes.
      \describe{
        \item{max.gappercent}{A value from 0 to 1 that indicates the
          tolerable percentage of gaps.}
        \item{max.gaplength}{An non-negative integer that indicates
          the tolerable number of consecutive gaps, as defined by
          \code{max.gappercent}.}
        \item{append}{Overwrite \code{colmask(x)} or append (via
          \code{union}) to it?}
      }
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
      \code{consensusMatrix(x, as.prob, baseOnly)}:
      Creates an integer matrix containing the column frequencies of
      the underlying alphabet with masked columns being represented
      with \code{NA} values. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{baseOnly} is \code{TRUE}, then the non-base
      letters are collapsed into an \code{"other"} category.
    }
    \item{}{
      \code{consensusString(x, ...)}:
      Creates a consensus string for \code{x} with the symbol \code{"#"}
      representing a masked column. See \code{\link{consensusString}}
      for details on the arguments.
    }
    \item{}{
      \code{consensusViews(x, ...)}:
      Similar to the \code{consensusString} method. It returns a 
      \link{XStringViews} on the consensus string containing subsequence
      contigs of non-masked columns. Unlike the \code{consensusString}
      method, the masked columns in the underlying string contain a
      consensus value rather than the \code{"#"} symbol.
    }
    \item{}{
      \code{alphabetFrequency(x, as.prob, collapse)}:
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
## create an object from file
origMAlign <-
  read.DNAMultipleAlignment(filepath =
                            system.file("extdata","clustalw2.aln",
                                        package="Biostrings"),
                            format="clustal")

## list the names of the sequences in the alignment
rownames(origMAlign)

## mask columns
modMAlign <- origMAlign
colmask(modMAlign) <- IRanges(start=c(1,20), end=c(4,30))
modMAlign

## mask rows
modMAlign <- origMAlign
rowmask(modMAlign) <- IRanges(start=1, end=2)
modMAlign

## remove row masks
rowmask(modMAlign) <- NULL
modMAlign

## automatically mask rows based on consecutive gaps
modMAlign <- automask(modMAlign, max.gappercent=.5, max.gaplength=3)
modMAlign

## calculate frequencies
alphabetFrequency(modMAlign)
consensusMatrix(modMAlign, baseOnly=TRUE)[, 1:8]

## get consensus values
consensusString(modMAlign)
consensusViews(modMAlign)
}

\keyword{methods}
\keyword{classes}
