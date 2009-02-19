\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,character,character-method}
\alias{pairwiseAlignment,character,XString-method}
\alias{pairwiseAlignment,character,XStringSet-method}
\alias{pairwiseAlignment,character,QualityScaledXStringSet-method}
\alias{pairwiseAlignment,XString,character-method}
\alias{pairwiseAlignment,XString,XString-method}
\alias{pairwiseAlignment,XString,XStringSet-method}
\alias{pairwiseAlignment,XString,QualityScaledXStringSet-method}
\alias{pairwiseAlignment,XStringSet,character-method}
\alias{pairwiseAlignment,XStringSet,XString-method}
\alias{pairwiseAlignment,XStringSet,XStringSet-method}
\alias{pairwiseAlignment,XStringSet,QualityScaledXStringSet-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,character-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,XString-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,XStringSet-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,QualityScaledXStringSet-method}

\title{Optimal Pairwise Alignment}
\description{
Solves (Needleman-Wunsch) global alignment, (Smith-Waterman) local alignment, and (ends-free) overlap alignment problems.
}
\usage{
pairwiseAlignment(pattern, subject, \dots)
\S4method{pairwiseAlignment}{XStringSet,XStringSet}(pattern, subject,
                  patternQuality = PhredQuality(22L), subjectQuality = PhredQuality(22L),
                  type = "global", substitutionMatrix = NULL, fuzzyMatrix = NULL,
                  gapOpening = -10, gapExtension = -4, scoreOnly = FALSE)
\S4method{pairwiseAlignment}{QualityScaledXStringSet,QualityScaledXStringSet}(pattern, subject,
                  type = "global", substitutionMatrix = NULL, fuzzyMatrix = NULL, 
                  gapOpening = -10, gapExtension = -4, scoreOnly = FALSE)
}
\arguments{
  \item{pattern}{a character vector of any length, an \code{\link{XString}}, or an
    \code{\link{XStringSet}} object.}
  \item{subject}{a character vector of length 1 or an \code{\link{XString}} object.}
  \item{patternQuality, subjectQuality}{objects of class \code{\link{XStringQuality}}
    representing the respective quality scores for \code{pattern} and \code{subject}
    that are used in a quality-based method for generating a substitution matrix.
    These two arguments are ignored if \code{!is.null(substitutionMatrix)} or if its
    respective string set (\code{pattern}, \code{subject}) is of class
    \code{\link{QualityScaledXStringSet}}.}
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"},
    \code{"overlap"}, \code{"global-local"}, and \code{"local-global"} where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties,
      \code{"global-local"} = align whole strings with end gap penalties on
        \code{pattern} and without end gap penalties on \code{subject}
      \code{"local-global"} = align whole strings without end gap penalties on
        \code{pattern} and with end gap penalties on \code{subject}.}
  \item{substitutionMatrix}{substitution matrix for a non-quality based alignment.
    It cannot be used in conjunction with \code{patternQuality} and
    \code{subjectQuality} arguments.}
  \item{fuzzyMatrix}{fuzzy match matrix for quality-based alignments. It takes
    values between 0 and 1; where 0 is an unambiguous mismatch, 1 is an unambiguous
    match, and values in between represent a fraction of "matchiness".}
  \item{gapOpening}{the cost for opening a gap in the alignment.}
  \item{gapExtension}{the incremental cost incurred along the length of the gap in the alignment.}
  \item{scoreOnly}{logical to denote whether or not to return just the scores of the
    optimal pairwise alignment.}
  \item{\dots}{optional arguments to generic function to support additional methods.}
}
\details{
If \code{scoreOnly == FALSE}, the pairwise alignment with the maximum alignment
score is returned. If more than one pairwise alignment has the maximum alignment
score exists, the first alignment along the subject is returned. If there are
multiple pairwise alignments with the maximum alignment score at the chosen subject
location, then at each location along the alignment mismatches are given preference
to insertions/deletions. For example, \code{pattern: [1] ATTA; subject: [1] AT-A}
is chosen above \code{pattern: [1] ATTA; subject: [1] A-TA} if they both have the
maximum alignment score.

General implementation based on Chapter 2 of Haubold and Wiehe (2006).
Quality-based method for generating a substitution matrix based on the Bioinformatics
article by Ketil Malde given below.
}
\value{
If \code{scoreOnly == FALSE}, an instance of class \code{\link{PairwiseAlignedXStringSet}}
or \code{\link{PairwiseAlignedFixedSubject}} is returned.
If \code{scoreOnly == TRUE}, a numeric vector containing the scores for the optimal
pairwise alignments is returned.
}
\references{
R. Durbin, S. Eddy, A. Krogh, G. Mitchison, Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.

B. Haubold, T. Wiehe, Introduction to Computational Biology, Birkhauser Verlag 2006, Chapter 2.

K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics 2008 24(7):897-900.
}
\note{
Use \code{\link{matchPattern}} or \code{\link{vmatchPattern}} if you need to
find all the occurences (eventually with indels) of a given pattern in a
reference sequence or set of sequences.

Use \code{\link{matchPDict}} if you need to match a (big) set of patterns
against a reference sequence.
}
\author{P. Aboyoun and H. Pages}
\seealso{
  \code{\link{stringDist}},
  \link{PairwiseAlignedXStringSet-class},
  \link{XStringQuality-class},
  \link{substitution.matrices},
  \code{\link{matchPattern}}
}
\examples{
  ## Nucleotide global, local, and overlap alignments
  s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
  s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")

  # First use a fixed substitution matrix
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  globalAlign <-
    pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  localAlign <-
    pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
  overlapAlign <-
    pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)

  # Then use quality-based method for generating a substitution matrix
  pairwiseAlignment(s1, s2,
                    patternQuality = SolexaQuality(rep(c(22L, 12L), times = c(36, 18))),
                    subjectQuality = SolexaQuality(rep(c(22L, 12L), times = c(40, 20))),
                    scoreOnly = TRUE)


  ## Amino acid global alignment
  pairwiseAlignment(AAString("PAWHEAE"), AAString("HEAGAWGHEE"), substitutionMatrix = "BLOSUM50",
                    gapOpening = 0, gapExtension = -8)
}
\keyword{models}
\keyword{methods}
