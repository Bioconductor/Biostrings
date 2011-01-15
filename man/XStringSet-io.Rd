\name{XStringSet-io}

\alias{read.BStringSet}
\alias{read.DNAStringSet}
\alias{read.RNAStringSet}
\alias{read.AAStringSet}

\alias{fasta.info}
\alias{fastq.geometry}

\alias{write.XStringSet}

\alias{save.XStringSet}

% Old stuff:
\alias{FASTArecordsToCharacter}
\alias{CharacterToFASTArecords}
\alias{FASTArecordsToXStringViews}
\alias{XStringSetToFASTArecords}
\alias{read.XStringViews}
\alias{write.XStringViews}


\title{Read/write an XStringSet object from/to a file}

\description{
  Functions to read/write an \link{XStringSet} object from/to a file.
}

\usage{
## Read FASTA (or FASTQ) files in an XStringSet object:
read.BStringSet(filepath, format="fasta",
                nrec=-1L, skip=0L, use.names=TRUE)
read.DNAStringSet(filepath, format="fasta",
                  nrec=-1L, skip=0L, use.names=TRUE)
read.RNAStringSet(filepath, format="fasta",
                  nrec=-1L, skip=0L, use.names=TRUE)
read.AAStringSet(filepath, format="fasta",
                 nrec=-1L, skip=0L, use.names=TRUE)

## Extract basic information about FASTA (or FASTQ) files
## without loading them:
fasta.info(filepath, nrec=-1L, skip=0L, use.names=TRUE)
fastq.geometry(filepath, nrec=-1L, skip=0L)

## Write an XStringSet object to a FASTA (or FASTQ) file:
write.XStringSet(x, filepath, append=FALSE, format="fasta", ...)

## Serialize an XStringSet object:
save.XStringSet(x, objname, dirpath=".", save.dups=FALSE, verbose=TRUE)
}

\arguments{
  \item{filepath}{
    A character vector (of arbitrary length when reading, of length 1
    when writing) containing the path(s) to the file(s) to read or write.
    Note that special values like \code{""} or \code{"|cmd"} (typically
    supported by other I/O functions in R) are not supported here.
    Also \code{filepath} cannot be a connection.
  }
  \item{format}{
    Either \code{"fasta"} (the default) or \code{"fastq"}.
  }
  \item{nrec}{
    Single integer. The maximum of number of records to read in.
    Negative values are ignored.
  }
  \item{skip}{
    Single non-negative integer. The number of records of the data file(s)
    to skip before beginning to read in records.
  }
  \item{use.names}{
    Should the returned vector be named? For FASTA the names are taken
    from the record description lines. For FASTQ they are taken from
    the record sequence ids.
    Dropping the names can help reducing memory footprint e.g. for
    a FASTQ file containing millions of reads.
  }
  \item{x}{
    For \code{write.XStringSet}, the object to write to \code{file}.

    For \code{save.XStringSet}, the object to serialize.
  }
  \item{append}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} output will be
    appended to \code{file}; otherwise, it will overwrite the contents
    of \code{file}. See \code{?\link[base]{cat}} for the details.
  }
  \item{...}{
    Further format-specific arguments.
    If \code{format="fasta"}, the \code{width} argument (single integer)
    can be used to specify the maximum number of letters per line of
    sequence.
    If \code{format="fastq"}, the \code{qualities} argument (\link{BStringSet}
    object) \emph{must} be used to specify the qualities.
  }
  \item{objname}{
    The name of the serialized object.
  }
  \item{dirpath}{
    The path to the directory where to save the serialized object.
  }
  \item{save.dups}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE} then the \code{\link[IRanges:Grouping-class]{Dups}}
    object describing 
    how duplicated elements in \code{x} are related to each other is
    saved too. For advanced users only.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
}

\details{
  Only FASTA and FASTQ files are supported for now. The qualities
  stored in the FASTQ records are ignored.

  Reading functions \code{read.BStringSet}, \code{read.DNAStringSet},
  \code{read.RNAStringSet} and \code{read.AAStringSet} load sequences
  from an input file (or set of input files) into an \link{XStringSet}
  object.
  When multiple input files are specified, they are read in the
  corresponding order and their data are stored in the returned object
  in that order. Note that when multiple input FASTQ files are specified,
  all must have the same "width" (i.e. all their sequences must have the
  same length).

  The \code{fasta.info} utility returns an integer vector with one element
  per FASTA record in the input files. Each element is the length of the
  sequence found in the corresponding record.

  The \code{fastq.geometry} utility returns an integer vector describing
  the "geometry" of the FASTQ files i.e. a vector of length 2 where the
  first element is the total number of FASTQ records in the files and
  the second element the common "width" of these files (this width is
  \code{NA} if the files contain no FASTQ records or records with
  different widths).

  \code{write.XStringSet} writes an \link{XStringSet} object to a file.
  WARNING: Please be aware that using \code{write.XStringSet} on a
  \link{BStringSet} object that contains the '\n' (LF) or '\r' (CR)
  characters or the FASTA markup characters '>' or ';' is almost
  guaranteed to produce a broken FASTA file!

  Serializing an \link{XStringSet} object with \code{save.XStringSet}
  is equivalent to using the standard \code{save} mechanism. But it will
  try to reduce the size of \code{x} in memory first before calling
  \code{save}. Most of the times this leads to a much reduced size on disk.
}

\seealso{
  \code{\link{readFASTA}},
  \code{\link{writeFASTA}},
  \link{XStringSet-class},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. READ/WRITE FASTA FILES
  ## ---------------------------------------------------------------------
  filepath <- system.file("extdata", "someORF.fa", package="Biostrings")
  fasta.info(filepath)
  x <- read.DNAStringSet(filepath)
  x
  out <- tempfile()
  write.XStringSet(x, out)

  ## ---------------------------------------------------------------------
  ## B. READ/WRITE FASTQ FILES
  ## ---------------------------------------------------------------------
  filepath <- system.file("extdata", "s_1_sequence.txt",
                          package="Biostrings")
  fastq.geometry(filepath)
  read.DNAStringSet(filepath, format="fastq")

  library(BSgenome.Celegans.UCSC.ce2)
  ## Create a "sliding window" on chr I:
  sw_start <- seq.int(1, length(Celegans$chrI)-50, by=50)
  sw <- Views(Celegans$chrI, start=sw_start, width=10)
  my_fake_shortreads <- as(sw, "XStringSet")
  my_fake_ids <- sprintf("ID\%06d",  seq_len(length(my_fake_shortreads)))
  names(my_fake_shortreads) <- my_fake_ids
  my_fake_shortreads
  my_fake_quals <- rep.int(BStringSet("]]]]]]]]]]"),
                           length(my_fake_shortreads))
  my_fake_quals
  out <- tempfile()
  write.XStringSet(my_fake_shortreads, out, format="fastq",
                   qualities=my_fake_quals)

  ## ---------------------------------------------------------------------
  ## C. SERIALIZATION
  ## ---------------------------------------------------------------------
  save.XStringSet(my_fake_shortreads, "my_fake_shortreads", dirpath=tempdir())
}

\keyword{utilities}
\keyword{manip}
