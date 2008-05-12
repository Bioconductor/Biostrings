\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,character,character-method}
\alias{pairwiseAlignment,character,XString-method}
\alias{pairwiseAlignment,XString,character-method}
\alias{pairwiseAlignment,XString,XString-method}

\title{Optimal Pairwise Alignment}
\description{
Solves the (Needleman-Wunsch) global alignment and (Smith-Waterman) local alignment problems.
}
\usage{
pairwiseAlignment(string1, string2, quality1 = 22L, quality2 = 22L, type = "global",
                  substitutionMatrix = NULL, gapOpening = -10, gapExtension = -2,
                  scoreOnly = FALSE)
}
\arguments{
  \item{string1, string2}{a character vector of length 1 or an \code{\link{XString}} object.}
  \item{quality1, quality2}{respective quality scores for \code{string1} and \code{string2}
    that are used in a quality-based method for generating a substitution matrix.
    These scores must either be represented by [0 - 99] integer vectors, character
    vectors, or \code{\link{BString}} objects. Characters are interpreted as [0 - 99]
    quality measures by subtracting 33 from their ASCII decimal representation
    (e.g. ! = 0, " = 1, \# = 2, ...). These two arguments are ignored if
    \code{!is.null(substitutionMatrix)}.}
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"},
    \code{"overlap"}, \code{"overlap1"}, and \code{"overlap2"} where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties,
      \code{"overlap1"} = align whole strings without end gap penalties on
        \code{string1} and with end gap penalties on \code{string2},
      \code{"overlap2"} = align whole strings with end gap penalties on
        \code{string1} and without end gap penalties on \code{string2}.}
  \item{substitutionMatrix}{constant substitution matrix for the alignment. Do not
    use \code{substitutionMatrix} in conjunction with \code{quality1} and \code{quality2}
    arguments.}
  \item{gapOpening}{penalty for opening a gap in the alignment.}
  \item{gapExtension}{penalty for extending a gap in the alignment.}
  \item{scoreOnly}{logical to denote whether or not to only return the scores of the
    optimal pairwise alignment. (See Value section below.)}
}
\details{
General implementation based on Chapter 2 of Haubold and Wiehe (2006).
Quality-based method for generating a substitution matrix based on the Bioinformatics
article by Ketil Malde given below.
}
\value{
If \code{scoreOnly == FALSE}, an instance of class \code{\link{XStringAlign}} is returned.
If \code{scoreOnly == TRUE}, a numeric vector containing the scores for the optimal
pairwise alignments is returned.
}
\references{
R. Durbin, S. Eddy, A. Krogh, G. Mitchison, Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.
B. Haubold, T. Wiehe, Introduction to Computational Biology, Birkhauser Verlag 2006, Chapter 2.
K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics, Feb 23, 2008.
}
\author{Patrick Aboyoun and Herve Pages.}
\seealso{
  \link{XStringAlign-class},
  \link{substitution.matrices}
}
\examples{
  ## Nucleotide global, local, and overlap alignments
  s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
  s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")

  # First use a constant substitution matrix
  mat <- matrix(-3, nrow = 4, ncol = 4)
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
  globalAlign <-
    pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  localAlign <-
    pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  overlapAlign <-
    pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)

  # Then use quality-based method for generating a substitution matrix
  pairwiseAlignment(s1, s2,
                    quality1 = rep(c(22L, 12L), times = c(36, 18)),
                    quality2 = rep(c(22L, 12L), times = c(40, 20)),
                    scoreOnly = TRUE)


  ## Amino acid global alignment
  pairwiseAlignment(AAString("PAWHEAE"), AAString("HEAGAWGHEE"), substitutionMatrix = "BLOSUM50",
                    gapOpening = 0, gapExtension = -8)
}
\keyword{models}
\keyword{methods}
