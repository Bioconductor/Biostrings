\name{XStringSet-class}
\docType{class}

% Classes:
\alias{class:XStringSet}
\alias{XStringSet-class}
\alias{XStringSet}

\alias{class:BStringSet}
\alias{BStringSet-class}
\alias{BStringSet}
\alias{class:DNAStringSet}
\alias{DNAStringSet-class}
\alias{DNAStringSet}
\alias{class:RNAStringSet}
\alias{RNAStringSet-class}
\alias{RNAStringSet}
\alias{class:AAStringSet}
\alias{AAStringSet-class}
\alias{AAStringSet}

% Methods:
\alias{width,character-method}
\alias{nchar,XStringSet-method}
\alias{windows,character-method}
\alias{narrow,character-method}
\alias{subseq,character-method}
\alias{threebands,character-method}
\alias{subseq<-,character-method}
\alias{subseq<-,XStringSet-method}

\alias{seqtype,XStringSet-method}
\alias{seqtype<-,XStringSet-method}

\alias{coerce,ANY,BStringSet-method}
\alias{coerce,ANY,DNAStringSet-method}
\alias{coerce,ANY,RNAStringSet-method}
\alias{coerce,ANY,AAStringSet-method}
\alias{coerce,ANY,XStringSet-method}

\alias{show,XStringSet-method}
\alias{showAsCell,XStringSet-method}

\alias{unlist,XStringSet-method}
\alias{as.character,XStringSet-method}
\alias{as.factor,XStringSet-method}
\alias{as.data.frame,XStringSet-method}
\alias{as.vector,XStringSet-method}
\alias{toString,XStringSet-method}
\alias{as.matrix,XStringSet-method}

\alias{updateObject,XStringSet-method}


\title{XStringSet objects}

\description{
  The BStringSet class is a container for storing a set of
  \code{\link{BString}} objects and for making its manipulation
  easy and efficient.

  Similarly, the DNAStringSet (or RNAStringSet, or AAStringSet) class is
  a container for storing a set of \code{\link{DNAString}}
  (or \code{\link{RNAString}}, or \code{\link{AAString}}) objects.

  All those containers derive directly (and with no additional slots)
  from the XStringSet virtual class.
}

\usage{
## Constructors:
BStringSet(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
DNAStringSet(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
RNAStringSet(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
AAStringSet(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)

## Accessor-like methods:
\S4method{width}{character}(x)
\S4method{nchar}{XStringSet}(x, type="chars", allowNA=FALSE)

## ... and more (see below)
}

\arguments{
  \item{x}{
    Either a character vector (with no NAs), or an \link{XString},
    XStringSet or \link{XStringViews} object.
  }
  \item{start,end,width}{
    Either \code{NA}, a single integer, or an integer vector of the same
    length as \code{x} specifying how \code{x} should be "narrowed"
    (see \code{?\link[IRanges]{narrow}} for the details).
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
  \item{type,allowNA}{
    Ignored.
  }
}

\details{
  The \code{BStringSet}, \code{DNAStringSet}, \code{RNAStringSet} and
  \code{AAStringSet} functions are constructors that can be used to
  turn input \code{x} into an XStringSet object of the desired base type.

  They also allow the user to "narrow" the sequences contained in \code{x}
  via proper use of the \code{start}, \code{end} and/or \code{width}
  arguments. In this context, "narrowing" means dropping a prefix or/and
  a suffix of each sequence in \code{x}.
  The "narrowing" capabilities of these constructors can be illustrated
  by the following property: if \code{x} is a character vector
  (with no NAs), or an XStringSet (or \link{XStringViews}) object,
  then the 3 following transformations are equivalent:
  \describe{
    \item{\code{BStringSet(x, start=mystart, end=myend, width=mywidth)}:}{
    }
    \item{\code{subseq(BStringSet(x), start=mystart, end=myend, width=mywidth)}:}{
    }
    \item{\code{BStringSet(subseq(x, start=mystart, end=myend, width=mywidth))}:}{
    }
  }
  Note that, besides being more convenient, the first form is also more
  efficient on character vectors.
}

\section{Accessor-like methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{\code{length(x)}:}{
      The number of sequences in \code{x}.
    }
    \item{\code{width(x)}:}{
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
      Note that \code{width(x)} is also defined for a character vector
      with no NAs and is equivalent to \code{nchar(x, type="bytes")}.
    }
    \item{\code{names(x)}:}{
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user-provided description or comment for each
      element in \code{x}.
      These are the only data in an XStringSet object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
    \item{\code{alphabet(x)}:}{
      Return \code{NULL}, \code{\link{DNA_ALPHABET}},
      \code{\link{RNA_ALPHABET}} or \code{\link{AA_ALPHABET}} depending on
      whether \code{x} is a BStringSet, DNAStringSet, RNAStringSet or
      AAStringSet object.
    }
    \item{\code{nchar(x)}:}{
      The same as \code{width(x)}.
    }
  }
}

\section{Subsequence extraction and related transformations}{
  In the code snippets below,
  \code{x} is a character vector (with no NAs),
  or an XStringSet (or \link{XStringViews}) object.

  \describe{
    \item{\code{subseq(x, start=NA, end=NA, width=NA)}:}{
      Applies \code{subseq} on each element in \code{x}.
      See \code{?\link[XVector]{subseq}} for the details.

      Note that this is similar to what \code{\link{substr}} does on a
      character vector. However there are some noticeable differences:

      (1) the arguments are \code{start} and \code{stop} for
      \code{\link{substr}};

      (2) the SEW interface (start/end/width) interface of \code{subseq}
      is richer (e.g. support for negative start or end values);
      and (3) \code{subseq} checks that the specified start/end/width values
      are valid i.e., unlike \code{\link{substr}}, it throws an error if
      they define "out of limits" subsequences or subsequences with a
      negative width.
    }
    \item{\code{narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)}:}{
      Same as \code{subseq}. The only differences are: (1) \code{narrow}
      has a \code{use.names} argument; and (2) all the things \code{narrow}
      and \code{subseq} work on
      (\link[IRanges]{IRanges}, XStringSet or
      \link{XStringViews} objects for \code{narrow},
      \link[XVector]{XVector} or XStringSet objects for
      \code{subseq}). But they both work and do the same thing on an
      XStringSet object.
    }
    \item{\code{threebands(x, start=NA, end=NA, width=NA)}:}{
      Like the method for \link[IRanges]{IRanges}
      objects, the
      \code{threebands} methods for character vectors and XStringSet
      objects extend the capability of \code{narrow} by returning the 3
      set of subsequences (the left, middle and right subsequences)
      associated to the narrowing operation.
      See \code{?\link[IRanges]{threebands}} in the
      \pkg{IRanges} package for the details.
    }
    \item{\code{subseq(x, start=NA, end=NA, width=NA) <- value}:}{
      A vectorized version of the \code{\link[XVector]{subseq<-}}
      method for \link[XVector]{XVector} objects.
      See \code{?`\link[XVector]{subseq<-}`} for the details.
    }
  }
}

\section{Subsetting and appending}{
  In the code snippets below,
  \code{x} and \code{values} are XStringSet objects,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{\code{x[i]}:}{
      Return a new XStringSet object made of the selected elements.
    }
    \item{\code{x[[i]]}:}{
      Extract the i-th \code{\link{XString}} object from \code{x}.
    }
    \item{\code{append(x, values, after=length(x))}:}{
      Add sequences in \code{values} to \code{x}.
    }
  }
}

\section{Set operations}{
  In the code snippets below,
  \code{x} and \code{y} are XStringSet objects.

  \describe{
    \item{\code{union(x, y)}:}{
      Union of \code{x} and \code{y}.
    }
    \item{\code{intersect(x, y)}:}{
      Intersection of \code{x} and \code{y}.
    }
    \item{\code{setdiff(x, y)}:}{
      Asymmetric set difference of \code{x} and \code{y}.
    }
    \item{\code{setequal(x, y)}:}{
      Set equality of \code{x} to \code{y}.
    }
  }
}

\section{Other methods}{
  In the code snippets below,
  \code{x} is an XStringSet object.

  \describe{
    \item{\code{unlist(x)}:}{
      Turns \code{x} into an \link{XString} object by combining the
      sequences in \code{x} together.
      Fast equivalent to \code{do.call(c, as.list(x))}.
    }
    \item{\code{as.character(x, use.names=TRUE)}:}{
      Converts \code{x} to a character vector of the same length as \code{x}.
      The \code{use.names} argument controls whether or not \code{names(x)}
      should be propagated to the names of the returned vector.
    }
    \item{\code{as.factor(x)}:}{
      Converts \code{x} to a factor, via \code{as.character(x)}.
    }
    \item{\code{as.matrix(x, use.names=TRUE)}:}{
      Returns a character matrix containing the "exploded" representation of
      the strings. Can only be used on an XStringSet object with
      equal-width strings.
      The \code{use.names} argument controls whether or not \code{names(x)}
      should be propagated to the row names of the returned matrix.
    }
    \item{\code{toString(x)}:}{
      Equivalent to \code{toString(as.character(x))}.
    }
    \item{\code{show(x)}:}{
      By default the \code{show} method displays 5 head and 5 tail
      lines. The number of lines can be altered by setting the global
      options \code{showHeadLines} and \code{showTailLines}. If the
      object length is less than the sum of the options, the full object
      is displayed. These options affect GRanges, GAlignments, IRanges,
      and XStringSet objects.
    }
  }
}

\section{Display}{
  The letters in a DNAStringSet, RNAStringSet, or AAStringSet
  object are colored when displayed by the \code{show()} method.
  Set global option \code{Biostrings.coloring} to FALSE to turn off
  this coloring.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link{readDNAStringSet}} and \code{\link{writeXStringSet}}
          for reading/writing a \link{DNAStringSet} object (or other
          \link{XStringSet} derivative) from/to a FASTA or FASTQ file.

    \item \link{XStringSet-comparison}

    \item \link{XString} objects.

    \item \link{XStringViews} objects.

    \item \link{XStringSetList} objects.

    \item \code{\link[XVector]{subseq}}, \code{\link[IRanges]{narrow}},
          and \code{\link{substr}}.

    \item \code{\link[XVector]{compact}}

    \item \link[XVector]{XVectorList} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. USING THE XStringSet CONSTRUCTORS ON A CHARACTER VECTOR OR FACTOR
## ---------------------------------------------------------------------
## Note that there is no XStringSet() constructor, but an XStringSet
## family of constructors: BStringSet(), DNAStringSet(), RNAStringSet(),
## etc...
x0 <- c("#CTC-NACCAGTAT", "#TTGA", "TACCTAGAG")
width(x0)
x1 <- BStringSet(x0)
x1

## 3 equivalent ways to obtain the same BStringSet object:
BStringSet(x0, start=4, end=-3)
subseq(x1, start=4, end=-3)
BStringSet(subseq(x0, start=4, end=-3))

dna0 <- DNAStringSet(x0, start=4, end=-3)
dna0  # 'options(Biostrings.coloring=FALSE)' to turn off coloring

names(dna0)
names(dna0)[2] <- "seqB"
dna0

## When the input vector contains a lot of duplicates, turning it into
## a factor first before passing it to the constructor will produce an
## XStringSet object that is more compact in memory:
library(hgu95av2probe)
x2 <- sample(hgu95av2probe$sequence, 999000, replace=TRUE)
dna2a <- DNAStringSet(x2)
dna2b <- DNAStringSet(factor(x2))  # slower but result is more compact
object.size(dna2a)
object.size(dna2b)

## ---------------------------------------------------------------------
## B. USING THE XStringSet CONSTRUCTORS ON A SINGLE SEQUENCE (XString
##    OBJECT OR CHARACTER STRING)
## ---------------------------------------------------------------------
x3 <- "abcdefghij"
BStringSet(x3, start=2, end=6:2)  # behaves like 'substring(x3, 2, 6:2)'
BStringSet(x3, start=-(1:6))
x4 <- BString(x3)
BStringSet(x4, end=-(1:6), width=3)

## Randomly extract 1 million 40-mers from C. elegans chrI:
extractRandomReads <- function(subject, nread, readlength)
{
    if (!is.integer(readlength))
        readlength <- as.integer(readlength)
    start <- sample(length(subject) - readlength + 1L, nread,
                    replace=TRUE)
    DNAStringSet(subject, start=start, width=readlength)
}
library(BSgenome.Celegans.UCSC.ce2)
rndreads <- extractRandomReads(Celegans$chrI, 1000000, 40)
## Notes:
## - This takes only 2 or 3 seconds versus several hours for a solution
##   using substring() on a standard character string.
## - The short sequences in 'rndreads' can be seen as the result of a
##   simulated high-throughput sequencing experiment. A non-realistic
##   one though because:
##     (a) It assumes that the underlying technology is perfect (the
##         generated reads have no technology induced errors).
##     (b) It assumes that the sequenced genome is exactly the same as the
##         reference genome.
##     (c) The simulated reads can contain IUPAC ambiguity letters only
##         because the reference genome contains them. In a real
##         high-throughput sequencing experiment, the sequenced genome
##         of course doesn't contain those letters, but the sequencer
##         can introduce them in the generated reads to indicate ambiguous
##         base-calling.
##     (d) The simulated reads come from the plus strand only of a single
##         chromosome.
## - See the getSeq() function in the BSgenome package for how to
##   circumvent (d) i.e. how to generate reads that come from the whole
##   genome (plus and minus strands of all chromosomes).

## ---------------------------------------------------------------------
## C. USING THE XStringSet CONSTRUCTORS ON AN XStringSet OBJECT
## ---------------------------------------------------------------------
library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)
probes

RNAStringSet(probes, start=2, end=-5)  # does NOT copy the sequence data!

## ---------------------------------------------------------------------
## D. USING THE XStringSet CONSTRUCTORS ON AN ORDINARY list OF XString
##    OBJECTS
## ---------------------------------------------------------------------
probes10 <- head(probes, n=10)
set.seed(33)
shuffled_nucleotides <- lapply(probes10, sample)
shuffled_nucleotides

DNAStringSet(shuffled_nucleotides)  # does NOT copy the sequence data!

## Note that the same result can be obtained in a more compact way with
## just:
set.seed(33)
endoapply(probes10, sample)

## ---------------------------------------------------------------------
## E. USING subseq() ON AN XStringSet OBJECT
## ---------------------------------------------------------------------
subseq(probes, start=2, end=-5)

subseq(probes, start=13, end=13) <- "N"
probes

## Add/remove a prefix:
subseq(probes, start=1, end=0) <- "--"
probes
subseq(probes, end=2) <- ""
probes

## Do more complicated things:
subseq(probes, start=4:7, end=7) <- c("YYYY", "YYY", "YY", "Y")
subseq(probes, start=4, end=6) <- subseq(probes, start=-2:-5)
probes

## ---------------------------------------------------------------------
## F. UNLISTING AN XStringSet OBJECT
## ---------------------------------------------------------------------
library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)
unlist(probes)

## ---------------------------------------------------------------------
## G. COMPACTING AN XStringSet OBJECT
## ---------------------------------------------------------------------
## As a particular type of XVectorList objects, XStringSet objects can
## optionally be compacted. Compacting is done typically before
## serialization. See ?compact for more information.
library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)

y <- subseq(probes[1:12], start=5)
probes@pool
y@pool
object.size(probes)
object.size(y)

y0 <- compact(y)
y0@pool
object.size(y0)
}

\keyword{methods}
\keyword{classes}
