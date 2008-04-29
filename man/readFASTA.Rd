\name{readFASTA}
\alias{readFASTA}
\alias{writeFASTA}

\title{Functions to read/write FASTA formatted files}
\description{
 FASTA is a simple file format for biological sequence data. A file may
contain one or more sequences, for each sequence there is a description 
line which begins with a \code{>}.
}
\usage{
  readFASTA(file, checkComments=TRUE, strip.desc=TRUE)
  writeFASTA(x, file="", width=80)
}
\arguments{
  \item{file}{
    Either a character string naming a file or a connection open
    for reading or writing.
    If \code{""} (the default for \code{writeFASTA}),
    then the function writes to the standard output connection (the console)
    unless redirected by \code{sink}.
  }
  \item{checkComments}{
    Whether or not comments, lines beginning with a semi-colon
    should be found and removed.
  }
  \item{strip.desc}{
    Whether or not the ">" marking the beginning of the description
    lines should be removed. Note that this argument is new
    in Biostrings >= 2.8. In previous versions \code{readFASTA}
    was keeping the ">".
  }
  \item{x}{
    A list as one returned by \code{readFASTA}.
  }
  \item{width}{
    The maximum number of letters per line of sequence.
  }
}
\details{
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
 A list with one element for each sequence in the file. The elements 
are in two parts, one the description and the second a character string
of the biological sequence.
}
\author{ R. Gentleman }

\seealso{
  \code{\link{read.BStringSet}},
  \code{\link{read.DNAStringSet}},
  \code{\link{read.RNAStringSet}},
  \code{\link{read.AAStringSet}},
  \code{\link{write.XStringSet}},
  \code{\link{read.table}},
  \code{\link{scan}},
  \code{\link{write.table}}
}

\examples{
f1 <- system.file("extdata/someORF.fa", package="Biostrings")
ff <- readFASTA(f1, strip.desc=TRUE)
desc <- sapply(ff, function(x) x$desc)
## Keep the "reverse complement" sequences only
ff2 <- ff[grep("reverse complement", desc, fixed=TRUE)]
writeFASTA(ff2, file.path(tempdir(), "someORF2.fa"))
}
\keyword{ manip}

