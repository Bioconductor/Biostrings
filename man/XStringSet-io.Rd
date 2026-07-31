\name{XStringSet-io}

\alias{XStringSet-io}

\alias{readBStringSet}
\alias{readDNAStringSet}
\alias{readRNAStringSet}
\alias{readAAStringSet}

\alias{fasta.seqlengths}
\alias{fasta.index}
\alias{fastq.seqlengths}
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
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               use.names=TRUE, with.qualities=FALSE)
readDNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               use.names=TRUE, with.qualities=FALSE)
readRNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               use.names=TRUE, with.qualities=FALSE)
readAAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               use.names=TRUE, with.qualities=FALSE)

## Extract basic information about FASTA (or FASTQ) files
## without actually loading the sequence data:
fasta.seqlengths(filepath,
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               seqtype="B", use.names=TRUE)
fasta.index(filepath,
               nrec=-1L, skip=0L, seek.first.rec=FALSE,
               seqtype="B")

fastq.seqlengths(filepath,
               nrec=-1L, skip=0L, seek.first.rec=FALSE)
fastq.geometry(filepath,
               nrec=-1L, skip=0L, seek.first.rec=FALSE)

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

    Also \code{filepath} cannot be a standard \link[base]{connection}.
    However \code{filepath} can be an object as returned by
    \code{open_input_files}. This object can be used to read files
    by chunks. See "READ FILES BY CHUNK" in the examples section for
    the details.
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
    \code{TRUE} (the default) or \code{FALSE}. If \code{TRUE}, then the
    returned vector is named. For FASTA the names are taken from the record
    description lines. For FASTQ they are taken from the record sequence ids.
    Dropping the names with \code{use.names=FALSE} can help reduce memory
    footprint e.g. for a FASTQ file containing millions of reads.
  }
  \item{with.qualities}{
    \code{TRUE} or \code{FALSE} (the default).
    This argument is only supported when reading a FASTQ file.
    If \code{TRUE}, then the quality strings are also read and returned
    in the \code{qualities} metadata column of the returned \link{DNAStringSet}
    object. Note that by default the quality strings are ignored. This
    helps reduce memory footprint if the FASTQ file contains millions of reads.

    Note that using \code{\link{readQualityScaledDNAStringSet}()} is the
    preferred way to load a set of DNA sequences and their qualities from
    a FASTQ file into Bioconductor. Its main advantage is that it will return
    a \link{QualityScaledDNAStringSet} object instead of a \link{DNAStringSet}
    object, which makes handling of the qualities more convenient and less
    error prone. See "READ A FASTQ FILE AS A QualityScaledDNAStringSet OBJECT"
    in the Examples section below for more information.
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
      \item \code{"AA"} for Amino Acid sequences i.e. only
            letters in \code{\link{AA_ALPHABET}} (case ignored) are
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

    Passing \code{TRUE} is equivalent to passing \code{"gzip"}.
  }
  \item{compression_level}{
    Not implemented yet.
  }
  \item{...}{
    Further format-specific arguments.

    If \code{format="fasta"}, the \code{width} argument can be used to
    specify the maximum number of letters per line of sequence.
    \code{width} must be a single integer, and uses an unbounded width if set to a value less than 1.

    If \code{format="fastq"}, the \code{qualities} argument can be used to
    specify the quality strings. \code{qualities} must be a \link{BStringSet}
    object. If the argument is omitted, then the quality strings are
    taken from the \code{qualities} metadata column of \code{x} (i.e. from
    \code{mcols(x)$qualities}). If \code{x} has no \code{qualities} metadata
    column and the \code{qualities} argument is omitted, then the fake
    quality ';' is assigned to each letter in \code{x} and written to
    the FASTQ file. If \code{x} is a \code{\link{QualityScaledXStringSet}} and
    \code{qualities} is not defined, the \code{qualities} contained in \code{x}
    are used automatically.
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
  platforms.

  \code{readDNAStringSet} and family (i.e. \code{readBStringSet},
  \code{readDNAStringSet}, \code{readRNAStringSet} and \code{readAAStringSet})
  load sequences from an input file (or multiple input files) into an
  \link{XStringSet} object. When multiple input files are specified,
  all must have the same format (i.e. FASTA or FASTQ) and files with
  different compression types can be mixed with non-compressed files.
  The files are read in the order they were specified and the sequences
  are stored in the returned object in the order they were read.

  Only FASTA and FASTQ files are supported for now.

  The \code{fasta.seqlengths} utility returns an integer vector with one
  element per FASTA record in the input files. Each element is the length
  of the sequence found in the corresponding record, that is, the number of
  valid one-letter sequence codes in the record. See description of the
  \code{seqtype} argument above for how to control the set of valid
  one-letter sequence codes.

  The \code{fasta.index} utility returns a data frame with 1 row per
  FASTA record in the input files and the following columns:
  \itemize{
    \item \code{recno}: The rank of the record in the (virtually) concatenated
          input files.
    \item \code{fileno}: The rank of the file where the record is located.
    \item \code{offset}: The offset of the record relative to the start of the
          file where it's located. Measured in bytes.
    \item \code{desc}: The description line (a.k.a. header) of the record.
    \item \code{seqlength}: The length of the sequence in the record (not
          counting invalid letters).
    \item \code{filepath}: The path to the file where the record is located.
          Always a local file, so if the user specified a remote file, this
          column will contain the path to the downloaded file.
  }
  A subset of this data frame can be passed to \code{readDNAStringSet}
  and family for direct access to an arbitrary subset of sequences. More
  precisely, if \code{fai} is a FASTA index that was obtained with
  \code{fasta.index(filepath, ..., seqtype="DNA")}, then
  \code{readDNAStringSet(fai[i, ])} is equivalent to
  \code{readDNAStringSet(filepath, ...)[i]} for any valid subscript \code{i},
  except that the former only loads the requested sequences in memory
  and thus will be more memory efficient if only a small subset of sequences
  is requested.

  The \code{fastq.seqlengths} utility returns the read lengths in an integer
  vector with one element per FASTQ record in the input files.

  The \code{fastq.geometry} utility is a convenience wrapper around
  \code{fastq.seqlengths} that returns an integer vector of length 2
  describing the \emph{geometry} of the FASTQ files. The first integer gives
  the total number of FASTQ records in the files and the second element the
  common length of the reads (this common length is set to \code{NA} in case
  of variable length reads or if no FASTQ record was found). This compact
  representation of the geometry can be useful if the FASTQ files are known
  to contain fixed length reads.

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
  \itemize{
    \item \link{BStringSet}, \link{DNAStringSet}, \link{RNAStringSet},
          and \link{AAStringSet} objects.

    \item \code{\link{readQualityScaledDNAStringSet}} and
          \code{\link{writeQualityScaledXStringSet}} for reading/writing
          a \link{QualityScaledDNAStringSet} object (or other
          \link{QualityScaledXStringSet} derivative) from/to a FASTQ file.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. READ/WRITE FASTA FILES
## ---------------------------------------------------------------------

## Read a non-compressed FASTA files:
filepath1 <- system.file("extdata", "someORF.fa", package="Biostrings")
fasta.seqlengths(filepath1, seqtype="DNA")
x1 <- readDNAStringSet(filepath1)
x1

## Read a gzip-compressed FASTA file:
filepath2 <- system.file("extdata", "someORF.fa.gz", package="Biostrings")
fasta.seqlengths(filepath2, seqtype="DNA")
x2 <- readDNAStringSet(filepath2)
x2

## Sanity check:
stopifnot(identical(as.character(x1), as.character(x2)))

## Read 2 FASTA files at once:
filepath3 <- system.file("extdata", "fastaEx.fa", package="Biostrings")
fasta.seqlengths(c(filepath2, filepath3), seqtype="DNA")
x23 <- readDNAStringSet(c(filepath2, filepath3))
x23

## Sanity check:
x3 <- readDNAStringSet(filepath3)
stopifnot(identical(as.character(x23), as.character(c(x2, x3))))

## Use a FASTA index to load only an arbitrary subset of sequences:
filepath4 <- system.file("extdata", "dm3_upstream2000.fa.gz",
                         package="Biostrings")
fai <- fasta.index(filepath4, seqtype="DNA")
head(fai)
head(fai$desc)
i <- sample(nrow(fai), 10)  # randomly pick up 10 sequences
x4 <- readDNAStringSet(fai[i, ])

## Sanity check:
stopifnot(identical(as.character(readDNAStringSet(filepath4)[i]),
                    as.character(x4)))

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

filepath5 <- system.file("extdata", "s_1_sequence.txt",
                         package="Biostrings")

fastq.geometry(filepath5)

## The quality strings are ignored by default:
reads <- readDNAStringSet(filepath5, format="fastq")
reads
mcols(reads)

## Use 'with.qualities=TRUE' to load them:
reads <- readDNAStringSet(filepath5, format="fastq", with.qualities=TRUE)
reads
mcols(reads)
mcols(reads)$qualities

## Each quality string contains one letter per nucleotide in the
## corresponding read:
stopifnot(identical(width(mcols(reads)$qualities), width(reads)))

## Write the reads to a FASTQ file:
outfile <- tempfile()
writeXStringSet(reads, outfile, format="fastq")
outfile2 <- tempfile()
writeXStringSet(reads, outfile2, compress=TRUE, format="fastq")

## Sanity checks:
stopifnot(identical(readLines(outfile), readLines(filepath5)))
stopifnot(identical(readLines(outfile), readLines(outfile2)))

## ---------------------------------------------------------------------
## C. READ FILES BY CHUNK
## ---------------------------------------------------------------------
## readDNAStringSet() supports reading an arbitrary number of FASTA or
## FASTQ records at a time in a loop. This can be useful to process
## big FASTA or FASTQ files by chunk and thus avoids loading the entire
## file in memory. To achieve this the files to read from need to be
## opened with open_input_files() first. Note that open_input_files()
## accepts a vector of file paths and/or URLs.

## With FASTA files:
files <- open_input_files(filepath4)
i <- 0
while (TRUE) {
    i <- i + 1
    ## Load 4000 records at a time. Each new call to readDNAStringSet()
    ## picks up where the previous call left.
    dna <- readDNAStringSet(files, nrec=4000)
    if (length(dna) == 0L)
        break
    cat("processing chunk", i, "...\n")
    ## do something with 'dna' ...
}

## With FASTQ files:
files <- open_input_files(filepath5)
i <- 0
while (TRUE) {
    i <- i + 1
    ## Load 75 records at a time.
    reads <- readDNAStringSet(files, format="fastq", nrec=75)
    if (length(reads) == 0L)
        break
    cat("processing chunk", i, "...\n")
    ## do something with 'reads' ...
}

## IMPORTANT NOTE: Like connections, the object returned by
## open_input_files() can NOT be shared across workers in the
## context of parallelization!

## ---------------------------------------------------------------------
## D. READ A FASTQ FILE AS A QualityScaledDNAStringSet OBJECT
## ---------------------------------------------------------------------

## Use readQualityScaledDNAStringSet() if you want the object to be
## returned as a QualityScaledDNAStringSet instead of a DNAStringSet
## object. See ?readQualityScaledDNAStringSet for more information.

## Note that readQualityScaledDNAStringSet() is a simple wrapper around
## readDNAStringSet() that does the following if the file contains
## "Phred quality scores" (which is the standard Sanger variant to
## assess reliability of a base call):
reads <- readDNAStringSet(filepath5, format="fastq", with.qualities=TRUE)
quals <- PhredQuality(mcols(reads)$qualities)
QualityScaledDNAStringSet(reads, quals)

## The call to PhredQuality() is replaced with a call to SolexaQuality()
## or IlluminaQuality() if the quality scores are Solexa quality scores.

## ---------------------------------------------------------------------
## E. GENERATE FAKE READS AND WRITE THEM TO A FASTQ FILE
## ---------------------------------------------------------------------

library(BSgenome.Celegans.UCSC.ce2)

## Create a "sliding window" on chr I:
sw_start <- seq.int(1, length(Celegans$chrI)-50, by=50)
sw <- Views(Celegans$chrI, start=sw_start, width=10)
my_fake_reads <- as(sw, "XStringSet")
my_fake_ids <- sprintf("ID\%06d",  seq_len(length(my_fake_reads)))
names(my_fake_reads) <- my_fake_ids
my_fake_reads

## Fake quality ';' will be assigned to each base in 'my_fake_reads':
out2 <- tempfile()
writeXStringSet(my_fake_reads, out2, format="fastq")

## Passing qualities thru the 'qualities' argument:
my_fake_quals <- rep.int(BStringSet("DCBA@?>=<;"), length(my_fake_reads))
my_fake_quals
out3 <- tempfile()
writeXStringSet(my_fake_reads, out3, format="fastq",
                qualities=my_fake_quals)

## ---------------------------------------------------------------------
## F. SERIALIZATION
## ---------------------------------------------------------------------
saveXStringSet(my_fake_reads, "my_fake_reads", dirpath=tempdir())
}

\keyword{utilities}
\keyword{manip}
