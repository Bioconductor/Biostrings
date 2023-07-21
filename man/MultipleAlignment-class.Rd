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
\alias{rowmask<-,MultipleAlignment,ANY-method}
\alias{colmask}
\alias{colmask,MultipleAlignment-method}
\alias{colmask<-}
\alias{colmask<-,MultipleAlignment,NULL-method}
\alias{colmask<-,MultipleAlignment,ANY-method}
\alias{maskMotif,MultipleAlignment,ANY-method}
\alias{maskGaps}
\alias{maskGaps,MultipleAlignment-method}
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
\alias{seqtype,MultipleAlignment-method}

% Read functions:
\alias{readDNAMultipleAlignment}
\alias{readRNAMultipleAlignment}
\alias{readAAMultipleAlignment}

% Write functions:
\alias{write.phylip}

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

% show style methods:
\alias{show,MultipleAlignment-method}
\alias{detail,MultipleAlignment-method}


\title{MultipleAlignment objects}

\description{
   The MultipleAlignment class is a container for storing multiple sequence
   alignments.
}

\usage{
## Constructors:
DNAMultipleAlignment(x=character(), start=NA, end=NA, width=NA,
    use.names=TRUE, rowmask=NULL, colmask=NULL)
RNAMultipleAlignment(x=character(), start=NA, end=NA, width=NA,
    use.names=TRUE, rowmask=NULL, colmask=NULL)
AAMultipleAlignment(x=character(), start=NA, end=NA, width=NA,
    use.names=TRUE, rowmask=NULL, colmask=NULL)

## Read functions:
readDNAMultipleAlignment(filepath, format)
readRNAMultipleAlignment(filepath, format)
readAAMultipleAlignment(filepath, format)

## Write funtions:
write.phylip(x, filepath)

## ... and more (see below)
}

\arguments{
  \item{x}{
    Either a character vector (with no NAs), or an \link{XString},
    \link{XStringSet} or \link{XStringViews} object containing
    strings with the same number of characters.  If writing out a Phylip
    file, then x would be a \link{MultipleAlignment} object
  }
  \item{start,end,width}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges]{narrow}} in the \pkg{IRanges} package for
    the details).
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
  \item{rowmask}{
    a NormalIRanges object that will set masking for rows
  }
  \item{colmask}{
    a NormalIRanges object that will set masking for columns
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
    \item{\code{unmasked(x)}:}{
      The underlying \link{XStringSet} object containing the multiple
      sequence alignment.
    }
    \item{\code{rownames(x)}:}{
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user-provided description or comment for each
      sequence in \code{x}.
    }
    \item{\code{rowmask(x)}, \code{rowmask(x, append, invert) <- value}:}{
      Gets and sets the \link{NormalIRanges} object representing the
      masked rows in \code{x}. The \code{append} argument takes
      \code{union}, \code{replace} or \code{intersect} to indicate how
      to combine the new \code{value} with \code{rowmask(x)}. The
      \code{invert} argument takes a logical argument to indicate
      whether or not to invert the new mask. The \code{value} argument
      can be of any class that is coercible to a \link{NormalIRanges}
      via the \code{as} function.
    }
    \item{\code{colmask(x)}, \code{colmask(x, append, invert) <- value}:}{
      Gets and sets the \link{NormalIRanges} object representing the
      masked columns in \code{x}. The \code{append} argument takes
      \code{union}, \code{replace} or \code{intersect} to indicate how
      to combine the new \code{value} with \code{colmask(x)}. The
      \code{invert} argument takes a logical argument to indicate
      whether or not to invert the new mask. The \code{value} argument
      can be of any class that is coercible to a \link{NormalIRanges}
      via the \code{as} function.
    }
    \item{\code{maskMotif(x, motif, min.block.width=1, ...)}:}{
      Returns a MultipleAlignment object with a modified column mask
      based upon motifs found in the consensus string where the consensus
      string keeps all the columns but drops the masked rows.
      \describe{
        \item{motif}{The motif to mask.}
        \item{min.block.width}{The minimum width of the blocks to mask.}
        \item{...}{Additional arguments for \code{matchPattern}.}
      }
    }
    \item{\code{maskGaps(x, min.fraction, min.block.width)}:}{
      Returns a MultipleAlignment object with a modified column mask
      based upon gaps in the columns. In particular, this mask is defined
      by \code{min.block.width} or more consecutive columns that have
      \code{min.fraction} or more of their non-masked rows containing
      gap codes.
      \describe{
        \item{min.fraction}{A value in \code{[0, 1]} that indicates
          the minimum fraction needed to call a gap in the consensus string
          (default is \code{0.5}).}
        \item{min.block.width}{A positive integer that indicates the
          minimum number of consecutive gaps to mask, as defined by
          \code{min.fraction} (default is \code{4}).}
      }
    }
    \item{\code{nrow(x)}:}{
      Returns the number of sequences aligned in \code{x}.
    }
    \item{\code{ncol(x)}:}{
      Returns the number of characters for each alignment in \code{x}.
    }
    \item{\code{dim(x)}:}{
      Equivalent to \code{c(nrow(x), ncol(x))}.
    }
    \item{\code{maskednrow(x)}:}{
      Returns the number of masked aligned sequences in \code{x}.
    }
    \item{\code{maskedncol(x)}:}{
      Returns the number of masked aligned characters in \code{x}.
    }
    \item{\code{maskeddim(x)}:}{
      Equivalent to \code{c(maskednrow(x), maskedncol(x))}.
    }
    \item{\code{maskedratio(x)}:}{
      Equivalent to \code{maskeddim(x) / dim(x)}.
    }
    \item{\code{nchar(x)}:}{
      Returns the number of unmasked aligned characters in \code{x},
      i.e. \code{ncol(x) - maskedncol(x)}.
    }
    \item{\code{alphabet(x)}:}{
      Equivalent to \code{alphabet(unmasked(x))}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a MultipleAlignment object.

  \describe{
    \item{\code{as(from, "DNAStringSet")}, \code{as(from, "RNAStringSet")},
      \code{as(from, "AAStringSet")}, \code{as(from, "BStringSet")}:}{
      Creates an instance of the specified \link{XStringSet} object subtype
      that contains the unmasked regions of the multiple sequence alignment
      in \code{x}.
    }
    \item{\code{as.character(x, use.names)}:}{
      Convert \code{x} to a character vector containing the unmasked
      regions of the multiple sequence alignment. \code{use.names}
      controls whether or not \code{rownames(x)} should be used to set
      the names of the returned vector (default is \code{TRUE}).
    }
    \item{\code{as.matrix(x, use.names)}:}{
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
    \item{\code{consensusMatrix(x, as.prob, baseOnly)}:}{
      Creates an integer matrix containing the column frequencies of
      the underlying alphabet with masked columns being represented
      with \code{NA} values. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{baseOnly} is \code{TRUE}, then the non-base
      letters are collapsed into an \code{"other"} category.
    }
    \item{\code{consensusString(x, ...)}:}{
      Creates a consensus string for \code{x} with the symbol \code{"#"}
      representing a masked column. See \code{\link{consensusString}}
      for details on the arguments.
    }
    \item{\code{consensusViews(x, ...)}:}{
      Similar to the \code{consensusString} method. It returns a 
      \link{XStringViews} on the consensus string containing subsequence
      contigs of non-masked columns. Unlike the \code{consensusString}
      method, the masked columns in the underlying string contain a
      consensus value rather than the \code{"#"} symbol.
    }
    \item{\code{alphabetFrequency(x, as.prob, collapse)}:}{
      Creates an integer matrix containing the row frequencies of
      the underlying alphabet. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{collapse} is \code{TRUE}, then returns the
      overall frequency instead of the frequency by row.
    }
    \item{\code{detail(x, invertColMask, hideMaskedCols)}:}{ Allows for a full
      pager driven display of the object so that masked cols and rows
      can be removed and the entire sequence can be visually
      inspected. If \code{hideMaskedCols} is set to it's default value
      of \code{TRUE} then the output will hide all the the masked
      columns in the output.  Otherwise, all columns will be displayed
      along with a row to indicate the masking status.  If
      \code{invertColMask} is \code{TRUE} then any displayed mask will
      be flipped so as to represent things in a way consistent with
      Phylip style files instead of the mask that is actually stored in
      the \code{MultipleAlignment} object.  Please notice that
      \code{invertColMask} will be ignored if \code{hideMaskedCols} is
      set to its default value of \code{TRUE} since in that case it will
      not make sense to show any masking information in the output.
      Masked rows are always hidden in the output.
    }
  }
}

\section{Display}{
  The letters in a DNAMultipleAlignment or RNAMultipleAlignment object
  are colored when displayed by the \code{show()} method. Set global
  option \code{Biostrings.coloring} to FALSE to turn off this coloring.
}

\author{P. Aboyoun and M. Carlson}

\seealso{
  \link{XStringSet-class},
  \link{MaskedXString-class}
}

\examples{
## create an object from file
origMAlign <-
  readDNAMultipleAlignment(filepath =
                           system.file("extdata",
                                       "msx2_mRNA.aln",
                                       package="Biostrings"),
                           format="clustal")

## list the names of the sequences in the alignment
rownames(origMAlign)

## rename the sequences to be the underlying species for MSX2
rownames(origMAlign) <- c("Human","Chimp","Cow","Mouse","Rat",
                          "Dog","Chicken","Salmon")
origMAlign

## See a detailed pager view
if (interactive()) {
detail(origMAlign)
}

## operations to mask rows
## For columns, just use colmask() and do the same kinds of operations
rowMasked <- origMAlign
rowmask(rowMasked) <- IRanges(start=1,end=3)
rowMasked

## remove rowumn masks
rowmask(rowMasked) <- NULL
rowMasked

## "select" rows of interest
rowmask(rowMasked, invert=TRUE) <- IRanges(start=4,end=7)
rowMasked

## or mask the rows that intersect with masked rows
rowmask(rowMasked, append="intersect") <- IRanges(start=1,end=5)
rowMasked

## TATA-masked
tataMasked <- maskMotif(origMAlign, "TATA")
colmask(tataMasked)

## automatically mask rows based on consecutive gaps
autoMasked <- maskGaps(origMAlign, min.fraction=0.5, min.block.width=4)
colmask(autoMasked)
autoMasked

## calculate frequencies
alphabetFrequency(autoMasked)
consensusMatrix(autoMasked, baseOnly=TRUE)[, 84:90]

## get consensus values
consensusString(autoMasked)
consensusViews(autoMasked)

## cluster the masked alignments
sdist <- stringDist(as(autoMasked,"DNAStringSet"), method="hamming")
clust <- hclust(sdist, method = "single")
plot(clust)
fourgroups <- cutree(clust, 4)
fourgroups

## write out the alignement object (with current masks) to Phylip format
write.phylip(x = autoMasked, filepath = tempfile("foo.txt",tempdir()))

}

\keyword{methods}
\keyword{classes}
