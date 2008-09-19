\name{stringDist}
\alias{stringDist}

\alias{stringDist,character-method}
\alias{stringDist,XStringSet-method}
\alias{stringDist,QualityScaledXStringSet-method}

\title{String Distance/Alignment Score Matrix}
\description{
Computes the Levenshtein edit distance or pairwise alignment score matrix for a set of strings.
}
\usage{
stringDist(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE, \dots)
\S4method{stringDist}{XStringSet}(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", quality = PhredQuality(22L), substitutionMatrix = NULL,
                   gapOpening = 0, gapExtension = -1)
\S4method{stringDist}{QualityScaledXStringSet}(x, method = "quality", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = -1)
}
\arguments{
  \item{x}{a character vector or an \code{\link{XStringSet}} object.}
  \item{method}{calculation method. One of \code{"levenshtein"}, \code{"quality"}, or
    \code{"substitutionMatrix"}.}
  \item{ignoreCase}{logical value indicating whether to ignore case during scoring.}
  \item{diag}{logical value indicating whether the diagonal of the matrix should be printed
    by \code{print.dist}.}
  \item{upper}{logical value indicating whether the diagonal of the matrix should be printed
    by \code{print.dist}.}
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"}, and \code{"overlap"},
    where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties.
     This argument is ignored if \code{method == "levenshtein"}.}
  \item{quality}{object of class \code{\link{XStringQuality}} representing the quality scores
    for \code{x} that are used in a quality-based method for generating a substitution matrix.
    This argument is ignored if \code{method != "quality"}.}
  \item{substitutionMatrix}{symmetric substitution matrix for the alignment. This argument is
    ignored if \code{method != "substitutionMatrix"}.}
  \item{gapOpening}{penalty for opening a gap in the alignment. This argument is ignored if
    \code{method == "levenshtein"}.}
  \item{gapExtension}{penalty for extending a gap in the alignment. This argument is ignored if
    \code{method == "levenshtein"}.}
  \item{\dots}{optional arguments to generic function to support additional methods.}
}
\details{
Uses the underlying pairwiseAlignment code to compute the distance/alignment score matrix.
}
\value{
Returns an object of class \code{"dist"}.
}
\author{P. Aboyoun}
\seealso{
  \link[stats]{dist},
  \link[base]{agrep},
  \link{pairwiseAlignment},
  \link{substitution.matrices}
}
\examples{
  stringDist(c("lazy", "HaZy", "crAzY"))
  stringDist(c("lazy", "HaZy", "crAzY"), ignoreCase = TRUE)

  data(srPhiX174)
  stringDist(srPhiX174[1:4])
  stringDist(srPhiX174[1:4], method = "quality",
             quality = SolexaQuality(quPhiX174[1:4]),
             gapOpening = -10, gapExtension = -4)
}
\keyword{character}
\keyword{multivariate}
\keyword{cluster}
