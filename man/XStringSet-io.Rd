\name{XStringSet-io}

\alias{read.XStringViews}
\alias{read.BStringSet}
\alias{read.DNAStringSet}
\alias{read.RNAStringSet}
\alias{read.AAStringSet}
\alias{write.XStringSet}
\alias{write.XStringViews}

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
  ## XStringSet object
  read.BStringSet(file, format)
  read.DNAStringSet(file, format)
  read.RNAStringSet(file, format)
  read.AAStringSet(file, format)
  write.XStringSet(x, file="", format, width=80)

  ## XStringViews object
  read.XStringViews(file, format, subjectClass, collapse="")
  write.XStringViews(x, file="", format, width=80)

  ## Some related helper functions
  FASTArecordsToCharacter(FASTArecs, use.names=TRUE)
  CharacterToFASTArecords(x)
  FASTArecordsToXStringViews(FASTArecs, subjectClass, collapse="")
  XStringSetToFASTArecords(x)
}

\arguments{
  \item{file}{
    A character vector with no NAs.
    If \code{""} (the default for \code{write.XStringSet} and
    \code{write.XStringViews}), then the functions write to the standard
    output connection (the console) unless redirected by \code{sink}.
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
  \item{width}{
    Only relevant if \code{format} is \code{"fasta"}.
    The maximum number of letters per line of sequence.
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
  qualities stored in the FASTQ records are ignored: only the sequences
  are returned.

  Reading functions \code{read.BStringSet}, \code{read.DNAStringSet},
  \code{read.RNAStringSet}, \code{read.AAStringSet} and \code{read.XStringViews}
  load sequences from a file into an \link{XStringSet} or \link{XStringViews}
  object. Only one FASTA file, but more than one FASTQ file, can be read at a
  time (by passing a character vector of length > 1). In that case, all the
  FASTQ files must have the same "width" (i.e. all their sequences must have
  the same length) and the sequences from all the files are stored in the
  returned object in the order they were read.

  Writing functions \code{write.XStringSet} and \code{write.XStringViews}
  write an \link{XStringSet} or \link{XStringViews} object to a file or
  connection. They only support the FASTA format for now.

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
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  x <- read.DNAStringSet(file, "fasta")
  x
  write.XStringSet(x, format="fasta")  # writes to the console

  ## ---------------------------------------------------------------------
  ## B. READ FASTQ FILES
  ## ---------------------------------------------------------------------
  file <- system.file("extdata", "s_1_sequence.txt", package="Biostrings")
  read.DNAStringSet(file, "fastq")  # FASTQ qualities are dropped

  ## ---------------------------------------------------------------------
  ## C. SOME RELATED HELPER FUNCTIONS
  ## ---------------------------------------------------------------------
  ## Converting 'x'...
  ## ... to a list of FASTA records (as one returned by the "readFASTA" function)
  x1 <- XStringSetToFASTArecords(x)
  ## ... to a named character vector
  x2 <- FASTArecordsToCharacter(x1) # same as 'as.character(x)'
}

\keyword{utilities}
\keyword{manip}
