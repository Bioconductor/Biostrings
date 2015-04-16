\name{stringDist}
\alias{stringDist}

\alias{stringDist,character-method}
\alias{stringDist,XStringSet-method}
\alias{stringDist,QualityScaledXStringSet-method}

\title{String Distance/Alignment Score Matrix}
\description{
Computes the Levenshtein edit distance or pairwise alignment score matrix for a
set of strings.
}
\usage{
stringDist(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE, \dots)
\S4method{stringDist}{XStringSet}(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE,
                   upper = FALSE, type = "global", quality = PhredQuality(22L),
                   substitutionMatrix = NULL, fuzzyMatrix = NULL, gapOpening = 0,
                   gapExtension = 1)
\S4method{stringDist}{QualityScaledXStringSet}(x, method = "quality", ignoreCase = FALSE,
                   diag = FALSE, upper = FALSE, type = "global", substitutionMatrix = NULL,
                   fuzzyMatrix = NULL, gapOpening = 0, gapExtension = 1)
}
\arguments{
  \item{x}{a character vector or an \code{\link{XStringSet}} object.}
  \item{method}{calculation method. One of \code{"levenshtein"},
    \code{"hamming"}, \code{"quality"}, or \code{"substitutionMatrix"}.}
  \item{ignoreCase}{logical value indicating whether to ignore case during
    scoring.}
  \item{diag}{logical value indicating whether the diagonal of the matrix
    should be printed by \code{print.dist}.}
  \item{upper}{logical value indicating whether the upper triangle of the matrix
    should be printed by \code{print.dist}.}
  \item{type}{(applicable when \code{method = "quality"} or
    \code{method = "substitutionMatrix"}).
    type of alignment. One of \code{"global"}, \code{"local"}, and
    \code{"overlap"}, where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties.}
  \item{quality}{(applicable when \code{method = "quality"}).
    object of class \code{\link{XStringQuality}} representing the quality scores
    for \code{x} that are used in a quality-based method for generating a
    substitution matrix.}
  \item{substitutionMatrix}{(applicable when
    \code{method = "substitutionMatrix"}).
    symmetric matrix representing the fixed substitution scores in the
    alignment.}
  \item{fuzzyMatrix}{(applicable when \code{method = "quality"}).
    fuzzy match matrix for quality-based alignments. It takes values between 0
    and 1; where 0 is an unambiguous mismatch, 1 is an unambiguous match, and
    values in between represent a fraction of "matchiness".}
  \item{gapOpening}{(applicable when \code{method = "quality"} or
    \code{method = "substitutionMatrix"}).
    penalty for opening a gap in the alignment.}
  \item{gapExtension}{(applicable when \code{method = "quality"} or
    \code{method = "substitutionMatrix"}).
    penalty for extending a gap in the alignment}
  \item{\dots}{optional arguments to generic function to support additional
    methods.}
}
\details{
When \code{method = "hamming"}, uses the underlying \code{neditStartingAt} code
to calculate the distances, where the Hamming distance is defined as the number
of substitutions between two strings of equal length. Otherwise, uses the
underlying \code{pairwiseAlignment} code to compute the distance/alignment
score matrix.
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

  data(phiX174Phage)
  plot(hclust(stringDist(phiX174Phage), method = "single"))

  data(srPhiX174)
  stringDist(srPhiX174[1:4])
  stringDist(srPhiX174[1:4], method = "quality",
             quality = SolexaQuality(quPhiX174[1:4]),
             gapOpening = 10, gapExtension = 4)
}
\keyword{character}
\keyword{multivariate}
\keyword{cluster}
