\name{BStringViews-io}

\alias{read.BStringViews}
\alias{write.BStringViews}

\alias{FASTArecordsToCharacter}
\alias{CharacterToFASTArecords}
\alias{FASTArecordsToBStringViews}
\alias{BStringViewsToFASTArecords}

\title{Input/output of BStringViews objects}

\description{
  Functions to create a \link{BStringViews} object by reading a file
  or to write a \link{BStringViews} object to a file.
}

\usage{
  read.BStringViews(file, format, subjectClass, collapse="")
  write.BStringViews(x, file="", format, width=80)

  FASTArecordsToCharacter(FASTArecs, use.names=TRUE)
  CharacterToFASTArecords(x)
  FASTArecordsToBStringViews(FASTArecs, subjectClass, collapse="")
  BStringViewsToFASTArecords(x)
}

\arguments{
  \item{file}{
    Either a character string naming a file or a connection open
    for reading or writing.
    If \code{""} (the default for \code{write.BStringViews}),
    \code{write.BStringViews} writes to the standard
    output connection (the console) unless redirected by \code{sink}.
  }
  \item{format}{
    Only \code{"fasta"} is supported for now.
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
  \item{x}{
    For \code{write.BStringViews}, the \link{BStringViews} object to be
    written to \code{file}.
    For \code{CharacterToFASTArecords}, the (possibly named) character
    vector to be converted to a list of FASTA records as one returned
    by \code{\link{readFASTA}}.
    For \code{BStringViewsToFASTArecords}, the \link{BStringViews} object
    to be converted to a list of FASTA records as one returned
    by \code{\link{readFASTA}}.
  }
  \item{width}{
    Only relevant if \code{format} is \code{"fasta"}.
    The maximum number of letters per line of sequence.
  }
  \item{FASTArecs}{
    A list of FASTA records as one returned by \code{\link{readFASTA}}.
  }
  \item{use.names}{
    Whether or not the \code{desc} component of the FASTA records should be
    used to set the names of the returned vector.
  }
}

\details{
  Only FASTA files are supported for now.

  \code{read.BStringViews} creates a \link{BStringViews} object
  by reading a file.

  \code{write.BStringViews} writes a \link{BStringViews} object
  to a file or connection.

  \code{FASTArecordsToCharacter}, \code{CharacterToFASTArecords},
  \code{FASTArecordsToBStringViews} and \code{BStringViewsToFASTArecords}
  are convenience functions used internally by \code{read.BStringViews}
  and \code{write.BStringViews} for switching between different
  representations of the same object.
}

\seealso{
  \code{\link{readFASTA}},
  \code{\link{writeFASTA}},
  \link{BStringViews-class},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class}
}

\examples{
  file <- system.file("extdata", "someORF.fsa", package="Biostrings")
  x <- read.BStringViews(file, "fasta", "DNAString")
  desc(x)
  write.BStringViews(x, format="fasta") # writes to the console

  ## Converting 'x'...
  ## ... to a list of FASTA records (as one returned by the "readFASTA" function)
  x1 <- BStringViewsToFASTArecords(x)
  ## ... to a named character vector
  x2 <- FASTArecordsToCharacter(x1) # same as 'as.character(x)'
}

\keyword{manip}
