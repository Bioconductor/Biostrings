\name{XStringSet-io}

\alias{read.BStringViews}
\alias{read.BStringSet}
\alias{read.DNAStringSet}
\alias{read.RNAStringSet}
\alias{read.AAStringSet}
\alias{write.XStringSet}
\alias{write.BStringViews}

\alias{FASTArecordsToCharacter}
\alias{CharacterToFASTArecords}
\alias{FASTArecordsToBStringViews}
\alias{XStringSetToFASTArecords}

\title{Read/write an XStringSet or BStringViews object from/to a file}

\description{
  Functions to read/write an \link{XStringSet} or \link{BStringViews} object
  from/to a file.
}

\usage{
  ## XStringSet object
  read.BStringSet(file, format)
  read.DNAStringSet(file, format)
  read.RNAStringSet(file, format)
  read.AAStringSet(file, format)
  write.XStringSet(x, file="", format, width=80)

  ## BStringViews object
  read.BStringViews(file, format, subjectClass, collapse="")
  write.BStringViews(x, file="", format, width=80)

  ## Some related helper functions
  FASTArecordsToCharacter(FASTArecs, use.names=TRUE)
  CharacterToFASTArecords(x)
  FASTArecordsToBStringViews(FASTArecs, subjectClass, collapse="")
  XStringSetToFASTArecords(x)
}

\arguments{
  \item{file}{
    Either a character string naming a file or a connection open
    for reading or writing.
    If \code{""} (the default for \code{write.XStringSet} and
    \code{write.BStringViews}), then the functions write to the standard
    output connection (the console) unless redirected by \code{sink}.
  }
  \item{format}{
    Only \code{"fasta"} is supported for now.
  }
  \item{x}{
    For \code{write.XStringSet} and \code{write.BStringViews}, the object to
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
    The class to be given to the subject of the \link{BStringViews} object
    created and returned by the function.
    Must be \code{"BString"} or one of its subclasses: \code{"DNAString"},
    \code{"RNAString"} or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string to be inserted between the views of the
    \link{BStringViews} object created and returned by the function.
  }
  \item{FASTArecs}{
    A list of FASTA records as one returned by \code{\link{readFASTA}}.
  }
  \item{use.names}{
    Whether or not the description line preceding each FASTA records should be
    used to set the names of the returned vector.
  }
}

\details{
  Only FASTA files are supported for now.

  Reading functions \code{read.BStringSet}, \code{read.DNAStringSet},
  \code{read.RNAStringSet}, \code{read.AAStringSet} and \code{read.BStringViews}
  load sequences from a file into an \link{XStringSet} or \link{BStringViews}
  object.

  Writing functions \code{write.XStringSet} and \code{write.BStringViews}
  write an \link{XStringSet} or \link{BStringViews} object to a file or
  connection.

  \code{FASTArecordsToCharacter}, \code{CharacterToFASTArecords},
  \code{FASTArecordsToBStringViews} and \code{XStringSetToFASTArecords}
  are helper functions used internally by \code{write.XStringSet} and
  \code{read.BStringViews} for switching between different
  representations of the same object.
}

\seealso{
  \code{\link{readFASTA}},
  \code{\link{writeFASTA}},
  \link{XStringSet-class},
  \link{BStringViews-class},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class}
}

\examples{
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  x <- read.DNAStringSet(file, "fasta")
  x
  write.XStringSet(x, format="fasta") # writes to the console

  ## Converting 'x'...
  ## ... to a list of FASTA records (as one returned by the "readFASTA" function)
  x1 <- XStringSetToFASTArecords(x)
  ## ... to a named character vector
  x2 <- FASTArecordsToCharacter(x1) # same as 'as.character(x)'
}

\keyword{manip}
