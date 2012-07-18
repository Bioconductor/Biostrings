\name{pid}
\alias{pid}

\alias{pid,PairwiseAlignments-method}

\title{Percent Sequence Identity}
\description{
Calculates the percent sequence identity for a pairwise sequence alignment.
}
\usage{
pid(x, type="PID1")
}
\arguments{
  \item{x}{a \code{\link{PairwiseAlignments}} object.}
  \item{type}{one of percent sequence identity. One of \code{"PID1"}, \code{"PID2"},
    \code{"PID3"}, and \code{"PID4"}. See Details for more information.}
}
\details{
Since there is no universal definition of percent sequence identity, the \code{pid} function
calculates this statistic in the following types:
\describe{
  \item{\code{"PID1"}:}{
    100 * (identical positions) / (aligned positions + internal gap positions)
  }
  \item{\code{"PID2"}:}{
    100 * (identical positions) / (aligned positions)
  }
  \item{\code{"PID3"}:}{
    100 * (identical positions) / (length shorter sequence)
  }
  \item{\code{"PID4"}:}{
    100 * (identical positions) / (average length of the two sequences)
  }
}
}
\value{
A numeric vector containing the specified sequence identity measures.
}
\references{
A. May, Percent Sequence Identity: The Need to Be Explicit, Structure 2004, 12(5):737.

G. Raghava and G. Barton, Quantification of the variation in percentage
identity for protein sequence alignments, BMC Bioinformatics 2006, 7:415.
}
\author{P. Aboyoun}
\seealso{
  \link{pairwiseAlignment},
  \link{PairwiseAlignments-class},
  \link{match-utils}
}
\examples{
  s1 <- DNAString("AGTATAGATGATAGAT")
  s2 <- DNAString("AGTAGATAGATGGATGATAGATA")

  palign1 <- pairwiseAlignment(s1, s2)
  palign1
  pid(palign1)

  palign2 <-
    pairwiseAlignment(s1, s2,
      substitutionMatrix =
      nucleotideSubstitutionMatrix(match = 2, mismatch = 10, baseOnly = TRUE))
  palign2
  pid(palign2, type = "PID4")
}
\keyword{methods}
