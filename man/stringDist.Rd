\name{stringDist}
\alias{stringDist}

\alias{stringDist,character-method}
\alias{stringDist,XStringSet-method}

\title{String Distance/Alignment Score Matrix}
\description{
Computes the Levenshtein edit distance or pairwise alignment score matrix for a set of strings.
}
\usage{
stringDist(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE, \dots)
\S4method{stringDist}{XStringSet}(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", quality = 22L, qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = 0, gapExtension = -1)
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
      \code{"overlap"} = align whole strings without end gap penalties.}
  \item{quality}{quality scores for \code{x} that are used in a quality-based method for
    generating a substitution matrix.
    When \code{qualityType = "Phred"}, these scores must either be represented by [0 - 99]
    integer vectors, character vectors, \code{\link{BString}}, or, in the case of
    \code{patternQuality}, \code{\link{BStringSet}} objects. When \code{qualityType = "Phred"},
    characters are interpreted as [0 - 99] quality measures by subtracting 33 from their ASCII
    decimal representation (e.g. ! = 0, " = 1, \# = 2, ...).
    When \code{qualityType = "Solexa"}, these scores must either be represented by [-5 - 99]
    integer vectors, character vectors, \code{\link{BString}}, or, in the case of
    \code{patternQuality}, \code{\link{BStringSet}} objects. When \code{qualityType = "Solexa"},
    characters are interpreted as [-5 - 99] quality measures by subtracting 64 from their ASCII
    decimal representation (e.g. ; = -5, < = -4, = = -3, ...).
    These two arguments are ignored if \code{!is.null(substitutionMatrix)}.}
  \item{qualityType}{type of quality measure. Either of \code{"Phred"} or \code{"Solexa"}.
    See \code{substitution.matrices} help page for more details.}
  \item{substitutionMatrix}{constant substitution matrix for the alignment. Do not
    use \code{substitutionMatrix} in conjunction with \code{patternQuality} and
    \code{subjectQuality} arguments.}
  \item{gapOpening}{penalty for opening a gap in the alignment.}
  \item{gapExtension}{penalty for extending a gap in the alignment.}
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
  stringDist(srPhiX174[1:4], method = "quality", quality = quPhiX174[1:4],
             gapOpening = -10, gapExtension = -4)
}
\keyword{character}
\keyword{multivariate}
\keyword{cluster}
