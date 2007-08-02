\name{BStringViews-io}

\alias{read.BStringViews}
\alias{write.BStringViews}

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
    The class to be given to the subject of the returned \link{BStringViews}
    object.
    Must be \code{"BString"} or one of its derivations: \code{"DNAString"},
    \code{"RNAString"} or \code{"AAString"}.
  }
  \item{collapse}{
    An optional character string used to separate the views of the
    returned \link{BStringViews} object.
  }
  \item{x}{
    The \link{BStringViews} object to be written.
  }
  \item{width}{
    Only relevant if \code{format} is \code{"fasta"}.
    The maximum number of letters per line of sequence.
  }
  \item{FASTArecs}{
    A list as one returned by \code{\link{readFASTA}}.
  }
}

\details{
  Only FASTA files are supported for now.

  \code{read.BStringViews} creates a \link{BStringViews} object
  by reading a file.

  \code{write.BStringViews} writes a \link{BStringViews} object
  to a file or connection.

  \code{FASTArecordsToBStringViews} and \code{BStringViewsToFASTArecords}
  are lower-level functions that are used by \code{read.BStringViews}
  and \code{write.BStringViews} for the conversion between a "FASTA
  records object" (a list as one returned by \code{\link{readFASTA}})
  and a \link{BStringViews} object.
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
  file <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
  x <- read.BStringViews(file, "fasta", "DNAString")
  write.BStringViews(x, format="fasta") # writes to the console
}

\keyword{manip}
