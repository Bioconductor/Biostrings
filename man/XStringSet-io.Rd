\name{XStringSet-io}

\alias{XStringSet-io}

\alias{readBStringSet}
\alias{readDNAStringSet}
\alias{readRNAStringSet}
\alias{readAAStringSet}

\alias{fasta.info}
\alias{fastq.geometry}

\alias{writeXStringSet}

\alias{saveXStringSet}


\title{Read/write an XStringSet object from/to a file}

\description{
  Functions to read/write an \link{XStringSet} object from/to a file.
}

\usage{
## Read FASTA (or FASTQ) files in an XStringSet object:
readBStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readDNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readRNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readAAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

## Extract basic information about FASTA (or FASTQ) files
## without actually loading the sequence data:
fasta.info(filepath,
           nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE,
           seqtype="B")
fastq.geometry(filepath, nrec=-1L, skip=0L, seek.first.rec=FALSE)

## Write an XStringSet object to a FASTA (or FASTQ) file:
writeXStringSet(x, filepath, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta", ...)

## Serialize an XStringSet object:
saveXStringSet(x, objname, dirpath=".", save.dups=FALSE, verbose=TRUE)
}

\arguments{
  \item{filepath}{
    A character vector (of arbitrary length when reading, of length 1
    when writing) containing the path(s) to the file(s) to read or write.
    Reading files in gzip format (which usually have the '.gz' extension) is
    supported.

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
  \item{seek.first.rec}{
    \code{TRUE} or \code{FALSE} (the default). If \code{TRUE}, then the
    reading function starts by setting the file position indicator at the
    beginning of the first line in the file that looks like the beginning of
    a FASTA (if \code{format} is \code{"fasta"}) or FASTQ (if \code{format}
    is \code{"fastq"}) record. More precisely this is the first line in the
    file that starts with a '>' (for FASTA) or a '@' (for FASTQ). An error
    is raised if no such line is found.

    Normal parsing then starts from there, and everything happens like if
    the file actually started there. In particular it will be an error if
    this first record is not a valid FASTA or FASTQ record.

    Using \code{seek.first.rec=TRUE} is useful for example to parse GFF3
    files with embedded FASTA data.
  }
  \item{use.names}{
    Should the returned vector be named? For FASTA the names are taken
    from the record description lines. For FASTQ they are taken from
    the record sequence ids.
    Dropping the names can help reducing memory footprint e.g. for
    a FASTQ file containing millions of reads.
  }
  \item{seqtype}{
    A single string specifying the type of sequences contained in the
    FASTA file(s). Supported sequence types:
    \itemize{
      \item \code{"B"} for anything i.e. any letter is a valid one-letter
            sequence code.
      \item \code{"DNA"} for DNA sequences i.e. only letters in
            \code{\link{DNA_ALPHABET}} (case ignored) are valid
            one-letter sequence codes.
      \item \code{"RNA"} for RNA sequences i.e. only letters in
            \code{\link{RNA_ALPHABET}} (case ignored) are valid
            one-letter sequence codes.
      \item \code{"AA"} for Amino Acid sequences. Currently treated as
            \code{"B"} but this will change in the near future i.e. only
            letters in \code{\link{AA_ALPHABET}} (case ignored) will be
            valid one-letter sequence codes.
    }
    Invalid one-letter sequence codes are ignored with a warning.
  }
  \item{x}{
    For \code{writeXStringSet}, the object to write to \code{file}.

    For \code{saveXStringSet}, the object to serialize.
  }
  \item{append}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} output will be
    appended to \code{file}; otherwise, it will overwrite the contents
    of \code{file}. See \code{?\link[base]{cat}} for the details.
  }
  \item{compress}{
    Like for the \code{save} function in base R, must be \code{TRUE} or
    \code{FALSE} (the default), or a single string specifying whether writing
    to the file is to use compression.
    The only type of compression supported at the moment is \code{"gzip"}.
    Except on Windows where they're not supported at all, bzip2-compressed
    files are supported by the reading functions but not by
    \code{writeXStringSet}.

    Passing \code{TRUE} is equivalent to passing \code{"gzip"}.
  }
  \item{compression_level}{
    Not implemented yet.
  }
  \item{...}{
    Further format-specific arguments.
    If \code{format="fasta"}, the \code{width} argument (single integer)
    can be used to specify the maximum number of letters per line of
    sequence.
    If \code{format="fastq"}, the \code{qualities} argument (\link{BStringSet}
    object) can be used to specify the qualities. If the qualities are
    omitted, then the fake quality ';' is assigned to each letter in
    \code{x} and written to the file.
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
  gzip compression is supported by reading and writing functions on all
  platforms. bzip2 compression is supported by reading functions only,
  and on Unix/Linux/MacOSX only (granted the bzip2 library and header
  files are installed in the standard locations). bzip2 compression is
  not supported on Windows.

  \code{readDNAStringSet} and family (i.e. \code{readBStringSet},
  \code{readDNAStringSet}, \code{readRNAStringSet} and \code{readAAStringSet})
  load sequences from an input file (or multiple input files) into an
  \link{XStringSet} object. When multiple input files are specified,
  all must have the same format (i.e. FASTA or FASTQ) and files with
  different compression types can be mixed with non-compressed files. 
  The files are read in the order they were specified and the sequences
  are stored in the returned object in the order they were read.

  Only FASTA and FASTQ files are supported for now. The read qualities
  stored in FASTQ files are ignored by \code{readDNAStringSet} and family.
  When multiple input FASTQ files are specified, all must have the same
  "width" (i.e. all their sequences must have the same length).

  The \code{fasta.info} utility returns an integer vector with one element
  per FASTA record in the input files. Each element is the length of the
  sequence found in the corresponding record, that is, the number of
  valid one-letter sequence codes in the record. See description of the
  \code{seqtype} argument above for how to control the set of valid
  one-letter sequence codes.

  The \code{fastq.geometry} utility returns an integer vector describing
  the "geometry" of the FASTQ files i.e. a vector of length 2 where the
  first element is the total number of FASTQ records in the files and
  the second element the common "width" of these files (this width is
  \code{NA} if the files contain no FASTQ records or records with
  different widths).

  \code{writeXStringSet} writes an \link{XStringSet} object to a file.
  Like with \code{readDNAStringSet} and family, only FASTA and FASTQ
  files are supported for now.
  WARNING: Please be aware that using \code{writeXStringSet} on a
  \link{BStringSet} object that contains the '\\n' (LF) or '\\r' (CR)
  characters or the FASTA markup characters '>' or ';' is almost
  guaranteed to produce a broken FASTA file!

  Serializing an \link{XStringSet} object with \code{saveXStringSet}
  is equivalent to using the standard \code{save} mechanism. But it will
  try to reduce the size of \code{x} in memory first before calling
  \code{save}. Most of the times this leads to a much reduced size on disk.
}

\references{
  \url{http://en.wikipedia.org/wiki/FASTA_format}
}

\seealso{
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

## Read a non-compressed FASTA files:
filepath1 <- system.file("extdata", "someORF.fa", package="Biostrings")
fasta.info(filepath1, seqtype="DNA")
x1 <- readDNAStringSet(filepath1)
x1

## Read a gzip-compressed FASTA file:
filepath2 <- system.file("extdata", "someORF.fa.gz", package="Biostrings")
fasta.info(filepath2, seqtype="DNA")
x2 <- readDNAStringSet(filepath2)
x2

## Sanity check:
stopifnot(identical(as.character(x1), as.character(x2)))

## Read 2 FASTA files at once:
filepath3 <- system.file("extdata", "fastaEx.fa", package="Biostrings")
fasta.info(c(filepath2, filepath3), seqtype="DNA")
x23 <- readDNAStringSet(c(filepath2, filepath3))
x23

## Sanity check:
x3 <- readDNAStringSet(filepath3)
stopifnot(identical(as.character(x23), as.character(c(x2, x3))))

## Write FASTA files:
out23a <- tempfile()
writeXStringSet(x23, out23a)
out23b <- tempfile()
writeXStringSet(x23, out23b, compress=TRUE)
file.info(c(out23a, out23b))$size

## Sanity checks:
stopifnot(identical(as.character(readDNAStringSet(out23a)),
                    as.character(x23)))
stopifnot(identical(readLines(out23a), readLines(out23b)))

## ---------------------------------------------------------------------
## B. READ/WRITE FASTQ FILES
## ---------------------------------------------------------------------
filepath <- system.file("extdata", "s_1_sequence.txt",
                        package="Biostrings")
fastq.geometry(filepath)
readDNAStringSet(filepath, format="fastq")

library(BSgenome.Celegans.UCSC.ce2)
## Create a "sliding window" on chr I:
sw_start <- seq.int(1, length(Celegans$chrI)-50, by=50)
sw <- Views(Celegans$chrI, start=sw_start, width=10)
my_fake_shortreads <- as(sw, "XStringSet")
my_fake_ids <- sprintf("ID\%06d",  seq_len(length(my_fake_shortreads)))
names(my_fake_shortreads) <- my_fake_ids
my_fake_shortreads

## Fake quality ';' will be assigned to each base in 'my_fake_shortreads':
out2 <- tempfile()
writeXStringSet(my_fake_shortreads, out2, format="fastq")

## Passing qualities thru the 'qualities' argument:
my_fake_quals <- rep.int(BStringSet("DCBA@?>=<;"),
                         length(my_fake_shortreads))
my_fake_quals
out3 <- tempfile()
writeXStringSet(my_fake_shortreads, out3, format="fastq",
                qualities=my_fake_quals)

## ---------------------------------------------------------------------
## C. SERIALIZATION
## ---------------------------------------------------------------------
saveXStringSet(my_fake_shortreads, "my_fake_shortreads", dirpath=tempdir())
}

\keyword{utilities}
\keyword{manip}
