\name{PairwiseAlignments-io}

\alias{PairwiseAlignments-io}

\alias{writePairwiseAlignments}


\title{Write a PairwiseAlignments object to a file}

\description{
  The \code{writePairwiseAlignments} function writes a
  \link{PairwiseAlignments} object to a file.
  Only the "pair" format is supported at the moment.
}

\usage{
writePairwiseAlignments(x, file="", Matrix=NA, block.width=50)
}

\arguments{
  \item{x}{
    A \link{PairwiseAlignments} object, typically returned by the
    \code{\link{pairwiseAlignment}} function.
  }
  \item{file}{
    A connection, or a character string naming the file to print
    to. If \code{""} (the default), \code{writePairwiseAlignments}
    prints to the standard output connection (aka the console)
    unless redirected by \code{sink}.
    If it is \code{"|cmd"}, the output is piped to the command given
    by \code{cmd}, by opening a pipe connection.
  }
  \item{Matrix}{
    A single string containing the name of the substitution matrix
    (e.g. \code{"BLOSUM50"}) used for the alignment.
    See the \code{substitutionMatrix} argument of the
    \code{\link{pairwiseAlignment}} function for the details.
    See \code{?\link{substitution.matrices}} for a list of
    predefined substitution matrices available in the Biostrings package.
  }
  \item{block.width}{
    A single integer specifying the maximum number of sequence
    letters (including the "-" letter, which represents gaps)
    per line.
  }
}

\details{
  The "pair" format is one of the numerous pairwise sequence alignment
  formats supported by the EMBOSS software. See
  \url{http://emboss.sourceforge.net/docs/themes/AlignFormats.html}
  for a brief (and rather informal) description of this format.
}

\note{
  This brief description of the "pair" format suggests that
  it is best suited for \emph{global} pairwise alignments,
  because, in that case, the original pattern and subject sequences
  can be inferred (by just removing the gaps).

  However, even though the "pair" format can also be used for non global
  pairwise alignments (i.e. for \emph{global-local}, \emph{local-global},
  and \emph{local} pairwise alignments), in that case the original
  pattern and subject sequences \emph{cannot} be inferred. This is because
  the alignment written to the file doesn't necessarily span the entire
  pattern (if \code{type(x)} is \code{local-global} or \code{local})
  or the entire subject (if \code{type(x)} is \code{global-local}
  or \code{local}).

  As a consequence, the \code{writePairwiseAlignments} function can be
  used on a \link{PairwiseAlignments} containing non global alignments
  (i.e. with \code{type(x) != "global"}) but with the 2 following caveats:
  \enumerate{
    \item The type of the alignments (\code{type(x)}) is not written
          to the file.
    \item The original pattern and subject sequences cannot be inferred.
          Furthermore, there is no way to infer their lengths (because
          we don't know if they've been trimmed or not).
  }

  Also note that the \code{\link{pairwiseAlignment}} function
  interprets the \code{gapOpening} and \code{gapExtension} arguments
  differently than most other alignment tools. As a consequence
  the values of the Gap_penalty and Extend_penalty fields written to
  the file are not the same as the values that were passed to the
  \code{gapOpening} and \code{gapExtension} arguments. With the
  following relationship:
  \itemize{
    \item Gap_penalty = - gapOpening - gapExtension
    \item Extend_penalty = - gapExtension
  }
}

\author{
  H. Pages
}

\references{
  \url{http://emboss.sourceforge.net/docs/themes/AlignFormats.html}
}

\seealso{
  \itemize{
    \item \code{\link{pairwiseAlignment}}
    \item \link{PairwiseAlignments-class}
    \item \link{substitution.matrices}
  }
}

\examples{
## With 1 pair:
pattern <- DNAString("CGTACGTAACGTTCGT")
subject <- DNAString("CGTCGTCGTCCGTAA")
x1 <- pairwiseAlignment(pattern, subject)
x1
writePairwiseAlignments(x1)
writePairwiseAlignments(x1, block.width=10)
## The 2 bottom-right numbers (16 and 15) are the lengths of
## the original pattern and subject, respectively.

x2 <- pairwiseAlignment(pattern, subject, type="global-local")
x2 # score is different!
writePairwiseAlignments(x2)
## By just looking at the file, we can't tell the length of the
## original subject! Could be 13, could be more...

## With more than 1 pair (and named patterns):
pattern <- DNAStringSet(c(myp1="ACCA", myp2="ACGCA", myp3="ACGGCA"))
x3 <- pairwiseAlignment(pattern, subject)
x3
writePairwiseAlignments(x3)
}

\keyword{utilities}
\keyword{manip}
