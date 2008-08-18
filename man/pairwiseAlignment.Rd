\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,character,character-method}
\alias{pairwiseAlignment,character,XString-method}
\alias{pairwiseAlignment,XString,character-method}
\alias{pairwiseAlignment,XString,XString-method}
\alias{pairwiseAlignment,XStringSet,character-method}
\alias{pairwiseAlignment,XStringSet,XString-method}

\title{Optimal Pairwise Alignment}
\description{
Solves (Needleman-Wunsch) global alignment, (Smith-Waterman) local alignment, and (ends-free) overlap alignment problems.
}
\usage{
pairwiseAlignment(pattern, subject, \dots)
\S4method{pairwiseAlignment}{XStringSet,XString}(pattern, subject,
                  patternQuality = 22L, subjectQuality = 22L, type = "global",
                  qualityType = "Phred", substitutionMatrix = NULL, gapOpening = -10,
                  gapExtension = -4, scoreOnly = FALSE)
}
\arguments{
  \item{pattern}{a character vector of any length, an \code{\link{XString}}, or an
    \code{\link{XStringSet}} object.}
  \item{subject}{a character vector of length 1 or an \code{\link{XString}} object.}
  \item{patternQuality, subjectQuality}{respective quality scores for \code{pattern} and
    \code{subject} that are used in a quality-based method for generating a substitution
    matrix.
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
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"},
    \code{"overlap"}, \code{"patternOverlap"}, and \code{"subjectOverlap"} where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties,
      \code{"patternOverlap"} = align whole strings without end gap penalties on
        \code{pattern} and with end gap penalties on \code{subject},
      \code{"subjectOverlap"} = align whole strings with end gap penalties on
        \code{pattern} and without end gap penalties on \code{subject}.}
  \item{qualityType}{type of quality measure. Either of \code{"Phred"} or \code{"Solexa"}.
    See \code{substitution.matrices} help page for more details.}
  \item{substitutionMatrix}{constant substitution matrix for the alignment. Do not
    use \code{substitutionMatrix} in conjunction with \code{patternQuality} and
    \code{subjectQuality} arguments.}
  \item{gapOpening}{the cost for opening a gap in the alignment.}
  \item{gapExtension}{the incremental cost incurred along the length of the gap in the alignment.}
  \item{scoreOnly}{logical to denote whether or not to return just the scores of the
    optimal pairwise alignment.}
  \item{\dots}{optional arguments to generic function to support additional methods.}
}
\details{
General implementation based on Chapter 2 of Haubold and Wiehe (2006).
Quality-based method for generating a substitution matrix based on the Bioinformatics
article by Ketil Malde given below.
}
\value{
If \code{scoreOnly == FALSE}, an instance of class \code{\link{PairwiseAlignment}} is returned.
If \code{scoreOnly == TRUE}, a numeric vector containing the scores for the optimal
pairwise alignments is returned.
}
\references{
R. Durbin, S. Eddy, A. Krogh, G. Mitchison, Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.
B. Haubold, T. Wiehe, Introduction to Computational Biology, Birkhauser Verlag 2006, Chapter 2.
K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics 2008 24(7):897-900.
}
\author{P. Aboyoun and H. Pages}
\seealso{
  \link{stringDist},
  \link{PairwiseAlignment-class},
  \link{substitution.matrices}
}
\examples{
  ## Nucleotide global, local, and overlap alignments
  s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
  s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")

  # First use a constant substitution matrix
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  globalAlign <-
    pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  localAlign <-
    pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  overlapAlign <-
    pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)

  # Then use quality-based method for generating a substitution matrix
  pairwiseAlignment(s1, s2,
                    patternQuality = rep(c(22L, 12L), times = c(36, 18)),
                    subjectQuality = rep(c(22L, 12L), times = c(40, 20)),
                    scoreOnly = TRUE)


  ## Amino acid global alignment
  pairwiseAlignment(AAString("PAWHEAE"), AAString("HEAGAWGHEE"), substitutionMatrix = "BLOSUM50",
                    gapOpening = 0, gapExtension = -8)
}
\keyword{models}
\keyword{methods}
