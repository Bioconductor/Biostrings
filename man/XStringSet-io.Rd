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
\alias{FASTArecordsToBStringViews}
\alias{read.BStringViews}
\alias{write.BStringViews}
\alias{read.XStringViews}
\alias{write.XStringViews}


\title{Read/write an XStringSet object from/to a file}

\description{
  Functions to read/write an \link{XStringSet} object from/to a file.
}

\usage{
  ## Read FASTA (or FASTQ) files in an XStringSet object:
  read.BStringSet(filepath, format="fasta")
  read.DNAStringSet(filepath, format="fasta")
  read.RNAStringSet(filepath, format="fasta")
  read.AAStringSet(filepath, format="fasta")

  ## Extract basic information about FASTA (or FASTQ) files
  ## without loading them:
  fasta.info(filepath, use.descs=TRUE)
  fastq.geometry(filepath)

  ## Write an XStringSet object to a FASTA (or FASTQ) file:
  write.XStringSet(x, filepath, append=FALSE, format="fasta", width=80)

  ## Serialize an XStringSet object:
  save.XStringSet(x, objname, dirpath=".", save.dups=FALSE, verbose=TRUE)
}

\arguments{
  \item{filepath}{
    A character vector containing the paths to the files to read or write.
  }
  \item{format}{
    Either \code{"fasta"} (the default) or \code{"fastq"}.
    Note that \code{write.XStringSet} only supports \code{"fasta"} for now.
  }
  \item{use.descs}{
    Should the returned vector be named with the description lines
    found in the FASTA records?
  }
  \item{x}{
    For \code{write.XStringSet}, the object to write to \code{file}.

    For \code{save.XStringSet}, the object to serialize.
  }
  \item{file}{
    A connection, or a character string naming the file to write
    to.  If \code{""} (the default), print to the standard output
    connection (generally the console) unless redirected by \code{sink}.
  }
  \item{append}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} output will be
    appended to \code{file}; otherwise, it will overwrite the contents
    of \code{file}. See \code{?\link[base]{cat}} for the details.
  }
  \item{width}{
    Only relevant if \code{format} is \code{"fasta"}.
    The maximum number of letters per line of sequence.
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
  Only FASTA and FASTQ files are supported for now. The identifiers and
  qualities stored in the FASTQ records are ignored (only the sequences
  are returned).

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
  If \code{use.descs} is \code{TRUE} (the default) then the returned
  vector is named with the description lines found in the FASTA records.

  The \code{fastq.geometry} utility returns an integer vector describing
  the "geometry" of the FASTQ files i.e. a vector of length 2 where the
  first element is the total number of FASTQ records in the files and
  the second element the common "width" of these files (this width is
  \code{NA} if the files contain no FASTQ records or records with
  different "widths").

  \code{write.XStringSet} writes an \link{XStringSet} object to a file.
  Its only supports the FASTA format for now.

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
  ## B. READ FASTQ FILES
  ## ---------------------------------------------------------------------
  filepath <- system.file("extdata", "s_1_sequence.txt", package="Biostrings")
  fastq.geometry(filepath)
  ## Only the FASTQ sequences are returned (identifiers and qualities
  ## are dropped):
  read.DNAStringSet(filepath, format="fastq")

  ## ---------------------------------------------------------------------
  ## C. SERIALIZATION
  ## ---------------------------------------------------------------------
  library(BSgenome.Celegans.UCSC.ce2)
  ## Create a "sliding window" on chr I:
  sw_start <- seq.int(1, length(Celegans$chrI)-50, by=50)
  sw <- Views(Celegans$chrI, start=sw_start, width=10)
  my_fake_shortreads <- as(sw, "XStringSet")
  save.XStringSet(my_fake_shortreads, "my_fake_shortreads", dirpath=tempdir())
}

\keyword{utilities}
\keyword{manip}
