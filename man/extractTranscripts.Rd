\name{extractTranscripts}

\alias{transcriptWidths}
\alias{extractTranscripts}
\alias{transcriptLocs2refLocs}

\title{DNA/RNA transcription and translation}

\description{
  Functions for transcription and/or translation of DNA or RNA
  sequences, and related utilities.
}

\usage{
  transcriptWidths(exonStarts=list(), exonEnds=list())

  extractTranscripts(x, exonStarts=list(), exonEnds=list(),
                     strand=character(0), reorder.exons.on.minus.strand=FALSE)

  transcriptLocs2refLocs(tlocs, exonStarts=list(), exonEnds=list(),
                         strand=character(0),
                         reorder.exons.on.minus.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link{DNAString} object.
  }
}

\details{
  \code{extractTranscripts} allows the user to extract a set of
  transcripts specified by the starts and ends of their exons
  as well as the strand the transcript is coming from.
  See \code{\link[GenomicFeatures:extractTranscriptsFromGenome]{extractTranscriptsFromGenome}}
  in the GenomicFeatures package for extracting transcripts from
  a genome.
}

\value{
  A \link{DNAStringSet} object for \code{extractTranscripts}.
}

\seealso{
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

  library(BSgenome.Celegans.UCSC.ce2)
  ## Both transcripts are on chrII:
  chrII <- Celegans$chrII
  transcripts <- extractTranscripts(chrII,
                   exonStarts=list(as.integer(starts1), as.integer(starts2)),
                   exonEnds=list(as.integer(ends1), as.integer(ends2)),
                   strand=c("+","-")
                 )
}

\keyword{methods}
\keyword{manip}
