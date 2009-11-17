\name{XStringSet-io}

\alias{read.XStringViews}
\alias{read.BStringSet}
\alias{read.DNAStringSet}
\alias{read.RNAStringSet}
\alias{read.AAStringSet}
\alias{write.XStringSet}
\alias{write.XStringViews}
\alias{save.XStringSet}

\alias{fastq.geometry}
\alias{FASTArecordsToCharacter}
\alias{CharacterToFASTArecords}
\alias{FASTArecordsToXStringViews}
\alias{XStringSetToFASTArecords}

% Deprecated:
\alias{read.BStringViews}
\alias{write.BStringViews}
\alias{FASTArecordsToBStringViews}


\title{Read/write an XStringSet or XStringViews object from/to a file}

\description{
  Functions to read/write an \link{XStringSet} or \link{XStringViews} object
  from/to a file.
}

\usage{
  ## XStringSet object:
  read.BStringSet(filepath, format)
  read.DNAStringSet(filepath, format)
  read.RNAStringSet(filepath, format)
  read.AAStringSet(filepath, format)
  write.XStringSet(x, file="", append=FALSE, format, width=80)

  ## Serializing an XStringSet object:
  save.XStringSet(x, objname, dirpath=".", save.dups=FALSE, verbose=TRUE)

  ## XStringViews object:
  read.XStringViews(filepath, format, subjectClass, collapse="")
  write.XStringViews(x, file="", append=FALSE, format, width=80)

  ## FASTQ utilities:
  fastq.geometry(filepath)

  ## Some related helper functions:
  FASTArecordsToCharacter(FASTArecs, use.names=TRUE)
  CharacterToFASTArecords(x)
  FASTArecordsToXStringViews(FASTArecs, subjectClass, collapse="")
  XStringSetToFASTArecords(x)
}

\arguments{
  \item{filepath}{
    A character vector with no NAs.
  }
  \item{format}{
    Either \code{"fasta"} or \code{"fastq"}.
    Note that \code{write.XStringSet} and \code{write.XStringViews} only
    support \code{"fasta"} for now.
  }
  \item{x}{
    For \code{write.XStringSet} and \code{write.XStringViews}, the object to
    write to \code{file}.
    For \code{CharacterToFASTArecords}, the (possibly named) character
    vector to be converted to a list of FASTA records as one returned
    by \code{\link{readFASTA}}.
    For \code{XStringSetToFASTArecords}, the \link{XStringSet} object
    to be converted to a list of FASTA records as one returned
    by \code{\link{readFASTA}}.
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
  \item{subjectClass}{
    The class to be given to the subject of the \link{XStringViews} object
    created and returned by the function.
    Must be the name of one of the direct XString subclasses i.e.
    \code{"BString"}, \code{"DNAString"}, \code{"RNAString"}
    or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string to be inserted between the views of the
    \link{XStringViews} object created and returned by the function.
  }
  \item{FASTArecs}{
    A list of FASTA records as one returned by \code{\link{readFASTA}}.
  }
  \item{use.names}{
    Whether or not the description line preceding each FASTA records should be
    used to set the names of the returned object.
  }
}

\details{
  Only FASTA and FASTQ files are supported for now. The identifiers and
  qualities stored in the FASTQ records are ignored (only the sequences
  are returned).

  Reading functions \code{read.BStringSet}, \code{read.DNAStringSet},
  \code{read.RNAStringSet}, \code{read.AAStringSet} and \code{read.XStringViews}
  load sequences from a file into an \link{XStringSet} or \link{XStringViews}
  object. Only one FASTA file, but more than one FASTQ file, can be read at a
  time (by passing a character vector of length > 1). In that case, all the
  FASTQ files must have the same "width" (i.e. all their sequences must have
  the same length) and the sequences from all the files are stored in the
  returned object in the order they were read.

  The \code{fastq.geometry} utility returns an integer vector describing
  the "geometry" of the FASTQ files i.e. a vector of length 2 where the
  first element is the total number of FASTQ records in the files and
  the second element the common "width" of these files (this width is
  \code{NA} if the files contain no FASTQ records or records with
  different "widths").

  Writing functions \code{write.XStringSet} and \code{write.XStringViews}
  write an \link{XStringSet} or \link{XStringViews} object to a file or
  connection. They only support the FASTA format for now.

  Serializing an \link{XStringSet} object with \code{save.XStringSet}
  is equivalent to using the standard \code{save} mechanism. But it will
  try to reduce the size of \code{x} in memory first before calling
  \code{save}. Most of the times this leads to a much reduced size on disk.

  \code{FASTArecordsToCharacter}, \code{CharacterToFASTArecords},
  \code{FASTArecordsToXStringViews} and \code{XStringSetToFASTArecords}
  are helper functions used internally by \code{write.XStringSet} and
  \code{read.XStringViews} for switching between different
  representations of the same object.
}

\seealso{
  \code{\link{fasta.info}},
  \code{\link{readFASTA}},
  \code{\link{writeFASTA}},
  \link{XStringSet-class},
  \link{XStringViews-class},
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
  x <- read.DNAStringSet(filepath, "fasta")
  x
  write.XStringSet(x, format="fasta")  # writes to the console

  ## ---------------------------------------------------------------------
  ## B. READ FASTQ FILES
  ## ---------------------------------------------------------------------
  filepath <- system.file("extdata", "s_1_sequence.txt", package="Biostrings")
  fastq.geometry(filepath)
  ## Only the FASTQ sequences are returned (identifiers and qualities
  ## are dropped):
  read.DNAStringSet(filepath, "fastq")

  ## ---------------------------------------------------------------------
  ## C. SERIALIZATION
  ## ---------------------------------------------------------------------
  library(BSgenome.Celegans.UCSC.ce2)
  ## Create a "sliding window" on chr I:
  sw_start <- seq.int(1, length(Celegans$chrI)-50, by=50)
  sw <- Views(Celegans$chrI, start=sw_start, width=10)
  my_fake_shortreads <- as(sw, "XStringSet")
  save.XStringSet(my_fake_shortreads, "my_fake_shortreads", dirpath=tempdir())

  ## ---------------------------------------------------------------------
  ## D. SOME RELATED HELPER FUNCTIONS
  ## ---------------------------------------------------------------------
  ## Converting 'x'...
  ## ... to a list of FASTA records (as one returned by the "readFASTA" function)
  x1 <- XStringSetToFASTArecords(x)
  ## ... to a named character vector
  x2 <- FASTArecordsToCharacter(x1) # same as 'as.character(x)'
}

\keyword{utilities}
\keyword{manip}
