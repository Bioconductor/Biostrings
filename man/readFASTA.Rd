\name{readFASTA}

\alias{readFASTA}
\alias{writeFASTA}

\title{Functions to read/write FASTA formatted files}
\description{
\code{readFASTA} and \code{writeFASTA} read from and write to a FASTA file.
Note that the object returned by \code{readFASTA} or passed to
\code{writeFASTA} is a standard list. For faster and more memory efficient
alternatives that return/accept an \link{XStringSet} object, see the
\code{\link{read.DNAStringSet}} function and family.
}
\usage{
  readFASTA(file, checkComments=TRUE, strip.descs=TRUE)
  writeFASTA(x, file="", desc=NULL, append=FALSE, width=80)
}
\arguments{
  \item{file}{
    Either a character string naming a file or a connection.
    If \code{""} (the default for \code{writeFASTA}),
    then the function writes to the standard output connection (the console)
    unless redirected by \code{sink}.
  }
  \item{checkComments}{
    Whether or not comments, lines beginning with a semi-colon
    should be found and removed.
  }
  \item{strip.descs}{
    Whether or not the ">" marking the beginning of the description
    lines should be removed. Note that this argument is new
    in Biostrings >= 2.8. In previous versions \code{readFASTA}
    was keeping the ">".
  }
  \item{x}{
    A list as one returned by \code{readFASTA} if \code{desc} is
    not specified (i.e. \code{NULL}). If \code{desc} is specified
    (see below) then \code{x} can also be a list-like object with
    \link{XString} elements (for example it can be an
    \link{XStringSet}, \link{XStringViews} or
    \link[BSgenome:BSgenome-class]{BSgenome}
    object) or just a character vector.
  }
  \item{desc}{
    If \code{NULL} (the default) then \code{x} must be a list as
    one returned by \code{readFASTA} and all the sequences in \code{x}
    are written to the file.
    Otherwise \code{desc} must be a character vector no longer than
    the number of sequences in \code{x} containing the descriptions
    of the sequences in \code{x} that must be written to the file.
  }
  \item{append}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} output will be
    appended to \code{file}; otherwise, it will overwrite the contents
    of \code{file}. See \code{?\link[base]{cat}} for the details.
  }
  \item{width}{
    The maximum number of letters per line of sequence.
  }
}
\details{
  FASTA is a simple file format for biological sequence data. A file may
  contain one or more sequences, for each sequence there is a description 
  line which begins with a \code{>}.

  FASTA is a widely used format in biology. It is a relatively simple markup.
  I am not aware of a standard. It might be nice to check to see if the 
  data that were parsed are sequences of some appropriate type, but without
  a standard that does not seem possible.

  There are many other packages that provide similar, but different 
  capabilities.  The one in the package seqinr seems most similar but they
  separate the biological sequence into single character strings, which
  is too inefficient for large problems.
}

\value{
  For \code{readFASTA}: A list with one element per FASTA record in the file.
  Each element is in two parts, one is the description of the record
  and the second a character string of the biological sequence.
}

\author{R. Gentleman, H. Pages. Improvements to writeFASTA by Kasper D. Hansen}

\seealso{
  \code{\link{read.DNAStringSet}},
  \code{\link{fasta.info}},
  \code{\link{write.XStringSet}},
  \code{\link{read.table}},
  \code{\link{scan}},
  \code{\link{write.table}},
  \code{\link[BSgenome:BSgenome-class]{BSgenome-class}}
}

\examples{
  f1 <- system.file("extdata", "someORF.fa", package="Biostrings")
  ff <- readFASTA(f1, strip.descs=TRUE)
  desc <- sapply(ff, function(x) x$desc)
  desc

  ## Keep the "reverse complement" sequences only:
  ff2 <- ff[grep("reverse complement", desc, fixed=TRUE)]

  ## Write them to a FASTA file:
  temp_file <- file.path(tempdir(), "temp.fa")
  writeFASTA(ff2, file=temp_file)

  ## Write the first 2 to a FASTA file with a modified description:
  writeFASTA(ff2, file=temp_file, desc=c("a", "b"))

  ## Write a genome to a FASTA file:
  library(BSgenome.Celegans.UCSC.ce2)
  writeFASTA(Celegans, file=temp_file, desc=seqnames(Celegans))
}

\keyword{utilities}
\keyword{manip}
