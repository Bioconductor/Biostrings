\name{extractTranscripts}

\alias{transcriptWidths}
\alias{extractTranscripts}
\alias{transcriptLocs2refLocs}

\title{Extract a set of transcripts}

\description{
  \code{extractTranscripts} allows the user to extract a set of
  transcripts specified by the starts and ends of their exons
  as well as the strand from which the transcript is coming.

  \code{transcriptWidths} only returns the lengths of the
  transcripts (called the "widths" in this context) specified
  by the starts and ends of their exons.

  \code{transcriptLocs2refLocs} converts transcript-based
  locations into reference-based locations.
}

\usage{
  extractTranscripts(x, exonStarts=list(), exonEnds=list(),
                     strand=character(0), reorder.exons.on.minus.strand=FALSE)

  transcriptWidths(exonStarts=list(), exonEnds=list())

  transcriptLocs2refLocs(tlocs, exonStarts=list(), exonEnds=list(),
                         strand=character(0),
                         reorder.exons.on.minus.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link{DNAString} or \link{MaskedDNAString} object.
  }
  \item{exonStarts, exonEnds}{
    The starts and ends of the exons, respectively.

    Each argument can be a list of integer vectors,
    an \link[IRanges]{IntegerList} object,
    or a character vector where each element is a
    comma-separated list of integers.
    In addition, the lists represented by \code{exonStarts}
    and \code{exonEnds} must have the same shape i.e.
    have the same lengths and have elements of the same lengths.
    The length of \code{exonStarts} and \code{exonEnds}
    is the number of transcripts.
  }
  \item{strand}{
    A character vector of the same length as \code{exonStarts} and
    \code{exonEnds} specifying the strand (\code{"+"} or \code{"-"})
    from which the transcript is coming.
  }
  \item{reorder.exons.on.minus.strand}{
    \code{TRUE} or \code{FALSE}. Should the order of exons
    for transcripts coming from the minus strand be reversed?
  }
  \item{tlocs}{
    A list of integer vectors of the same length as \code{exonStarts}
    and \code{exonEnds}. Each element in \code{tlocs} must contain
    transcript-based locations.
  }
}

\details{
  \code{extractTranscripts} allows the user to extract a set of
  transcripts specified by the starts and ends of their exons
  as well as the strand from which the transcript is coming.
  See \code{\link[GenomicFeatures:extractTranscriptsFromGenome]{extractTranscriptsFromGenome}}
  in the GenomicFeatures package for extracting transcripts from
  a genome.
}

\value{
  A \link{DNAStringSet} object for \code{extractTranscripts}.

  An integer vector for \code{transcriptWidths}.

  A list of integer vectors of the same shape as \code{tlocs}
  for \code{transcriptLocs2refLocs}.
}

\seealso{
  \code{\link[GenomicFeatures:extractTranscriptsFromGenome]{extractTranscriptsFromGenome}},
  \code{\link{reverseComplement}},
  \link{DNAString-class},
  \link{DNAStringSet-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. EXTRACTING WORM TRANSCRIPTS ZC101.3 AND F37B1.1
  ## ---------------------------------------------------------------------

  ## Transcript ZC101.3 (is on + strand):
  ##   Exons starts/ends relative to transcript:
  rstarts1 <- c(1, 488, 654, 996, 1365, 1712, 2163, 2453)
  rends1 <- c(137, 578, 889, 1277, 1662, 1870, 2410, 2561)
  ##   Exons starts/ends relative to chromosome:
  starts1 <- 14678410 + rstarts1
  ends1 <- 14678410 + rends1

  ## Transcript F37B1.1 (is on - strand):
  ##   Exons starts/ends relative to transcript:
  rstarts2 <- c(1, 325)
  rends2 <- c(139, 815)
  ##   Exons starts/ends relative to chromosome:
  starts2 <- 13611188 - rends2
  ends2 <- 13611188 - rstarts2

  exon_starts <- list(as.integer(starts1), as.integer(starts2))
  exon_ends <- list(as.integer(ends1), as.integer(ends2))

  library(BSgenome.Celegans.UCSC.ce2)
  ## Both transcripts are on chrII:
  chrII <- Celegans$chrII
  transcripts <- extractTranscripts(chrII,
                   exonStarts=exon_starts,
                   exonEnds=exon_ends,
                   strand=c("+","-"))

  ## Same as 'width(transcripts)':
  transcriptWidths(exonStarts=exon_starts, exonEnds=exon_ends)

  transcriptLocs2refLocs(list(c(1:6, 135:140, 1555:1560), c(1:6, 137:142, 625:630)),
                   exonStarts=exon_starts,
                   exonEnds=exon_ends,
                   strand=c("+","-"))

  ## A sanity check:
  ref_locs <- transcriptLocs2refLocs(list(1:1560, 1:630),
                   exonStarts=exon_starts,
                   exonEnds=exon_ends,
                   strand=c("+","-"))
  stopifnot(chrII[ref_locs[[1]]] == transcripts[[1]])
  stopifnot(complement(chrII)[ref_locs[[2]]] == transcripts[[2]])
}

\keyword{methods}
\keyword{manip}
