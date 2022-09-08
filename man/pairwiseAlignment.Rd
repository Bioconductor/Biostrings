\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,ANY,ANY-method}
\alias{pairwiseAlignment,ANY,QualityScaledXStringSet-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,ANY-method}
\alias{pairwiseAlignment,QualityScaledXStringSet,QualityScaledXStringSet-method}

\title{Optimal Pairwise Alignment}
\description{
Solves (Needleman-Wunsch) global alignment, (Smith-Waterman) local alignment,
and (ends-free) overlap alignment problems.
}
\usage{
pairwiseAlignment(pattern, subject, \dots)

\S4method{pairwiseAlignment}{ANY,ANY}(pattern, subject,
                  patternQuality=PhredQuality(22L),
                  subjectQuality=PhredQuality(22L),
                  type="global",
                  substitutionMatrix=NULL, fuzzyMatrix=NULL,
                  gapOpening=10, gapExtension=4,
                  scoreOnly=FALSE)

\S4method{pairwiseAlignment}{QualityScaledXStringSet,QualityScaledXStringSet}(pattern, subject,
                  type="global",
                  substitutionMatrix=NULL, fuzzyMatrix=NULL, 
                  gapOpening=10, gapExtension=4,
                  scoreOnly=FALSE)
}

\arguments{
  \item{pattern}{a character vector or \code{\link{XStringSet}} derivative
    of any length, or an \code{\link{XString}} derivative.}
  \item{subject}{a character vector or \code{\link{XStringSet}} derivative
    of length 1 or \code{length(pattern)}, or an \code{\link{XString}}
    derivative.}
  \item{patternQuality, subjectQuality}{objects of class
    \code{\link{XStringQuality}} representing the respective quality scores for
    \code{pattern} and \code{subject} that are used in a quality-based method
    for generating a substitution matrix. These two arguments are ignored if
    \code{!is.null(substitutionMatrix)} or if its respective string set
    (\code{pattern}, \code{subject}) is of class
    \code{\link{QualityScaledXStringSet}}.}
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"},
    \code{"overlap"}, \code{"global-local"}, and \code{"local-global"} where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties,
      \code{"global-local"} = align whole strings in \code{pattern} with
        consecutive subsequence of \code{subject},
      \code{"local-global"} = align consecutive subsequence of \code{pattern}
        with whole strings in \code{subject}.}
  \item{substitutionMatrix}{substitution matrix representing the fixed
    substitution scores for an alignment. It cannot be used in conjunction with
    \code{patternQuality} and \code{subjectQuality} arguments.}
  \item{fuzzyMatrix}{fuzzy match matrix for quality-based alignments. It takes
    values between 0 and 1; where 0 is an unambiguous mismatch, 1 is an
    unambiguous match, and values in between represent a fraction of
    "matchiness". (See details section below.)}
  \item{gapOpening}{the cost for opening a gap in the alignment.}
  \item{gapExtension}{the incremental cost incurred along the length of the gap
    in the alignment.}
  \item{scoreOnly}{logical to denote whether or not to return just the scores of
    the optimal pairwise alignment.}
  \item{\dots}{optional arguments to generic function to support additional
    methods.}
}
\details{
Quality-based alignments are based on the paper the Bioinformatics article by
Ketil Malde listed in the Reference section below. Let \eqn{\epsilon_i} be the
probability of an error in the base read. For \code{"Phred"} quality measures
\eqn{Q} in \eqn{[0, 99]}, these error probabilities are given by
\eqn{\epsilon_i = 10^{-Q/10}}. For \code{"Solexa"} quality measures \eqn{Q} in
\eqn{[-5, 99]}, they are given by \eqn{\epsilon_i = 1 - 1/(1 + 10^{-Q/10})}.
Assuming independence within and between base reads, the combined error
probability of a mismatch when the underlying bases do match is
\eqn{\epsilon_c = \epsilon_1 + \epsilon_2 - (n/(n-1)) * \epsilon_1 * \epsilon_2},
where \eqn{n} is the number of letters in the underlying alphabet (i.e.
\eqn{n = 4} for DNA input, \eqn{n = 20} for amino acid input, otherwise
\eqn{n} is the number of distinct letters in the input).
Using \eqn{\epsilon_c}, the substitution score is given by
\eqn{b * \log_2(\gamma_{x,y} * (1 - \epsilon_c) * n + (1 - \gamma_{x,y}) * \epsilon_c * (n/(n-1)))},
where \eqn{b} is the bit-scaling for the scoring and \eqn{\gamma_{x,y}} is the
probability that characters \eqn{x} and \eqn{y} represents the same underlying
information (e.g. using IUPAC, \eqn{\gamma_{A,A} = 1} and \eqn{\gamma_{A,N} = 1/4}.
In the arguments listed above \code{fuzzyMatch} represents \eqn{\gamma_{x,y}}
and \code{patternQuality} and \code{subjectQuality} represents \eqn{\epsilon_1}
and \eqn{\epsilon_2} respectively.

If \code{scoreOnly == FALSE}, a pairwise alignment with the maximum alignment
score is returned. If more than one pairwise alignment produces the maximum
alignment score, then the alignment with the smallest initial deletion whose
mismatches occur before its insertions and deletions is chosen. For example,
if \code{pattern = "AGTA"} and \code{subject = "AACTAACTA"}, then the alignment
\code{pattern: [1] AG-TA; subject: [1] AACTA} is chosen over
\code{pattern: [1] A-GTA; subject: [1] AACTA} or
\code{pattern: [1] AG-TA; subject: [5] AACTA} if they all achieve the maximum
alignment score.
}
\value{
If \code{scoreOnly == FALSE} (the default), the function returns a
\code{\link{PairwiseAlignmentsSingleSubject}} object (if a single subject
was supplied) or a \code{\link{PairwiseAlignments}} object (if more
than one subject was supplied). In both cases, the returned object
contains N \emph{optimal pairwise alignments} where N is the number of
supplied patterns, that is, N = \code{length(pattern)} if \code{pattern}
is a character vector or \code{\link{XStringSet}} derivative, or N = 1 if
it's an \code{\link{XString}} derivative. If more than one subject
was supplied, the alignments in the returned \code{\link{PairwiseAlignments}}
object are obtained by aligning \code{pattern[[1]]} to \code{subject[[1]]},
\code{pattern[[2]]} to \code{subject[[2]]}, \code{pattern[[3]]} to
\code{subject[[3]]}, etc...

If \code{scoreOnly == TRUE}, a numeric vector containing the scores for
the N \emph{optimal pairwise alignments} is returned.
}
\references{
R. Durbin, S. Eddy, A. Krogh, G. Mitchison, Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.

B. Haubold, T. Wiehe, Introduction to Computational Biology, Birkhauser Verlag 2006, Chapter 2.

K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics 2008 24(7):897-900.
}
\note{
Use \code{\link{matchPattern}} or \code{\link{vmatchPattern}} if you need to
find all the occurrences (eventually with indels) of a given pattern in a
reference sequence or set of sequences.

Use \code{\link{matchPDict}} if you need to match a (big) set of patterns
against a reference sequence.
}
\author{P. Aboyoun and H. Pagès}
\seealso{
  \code{\link{writePairwiseAlignments}},
  \code{\link{stringDist}},
  \link{PairwiseAlignments-class},
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
    pairwiseAlignment(s1, s2, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
  localAlign <-
    pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
  overlapAlign <-
    pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)

  # Then use quality-based method for generating a substitution matrix
  pairwiseAlignment(s1, s2,
                    patternQuality = SolexaQuality(rep(c(22L, 12L), times = c(36, 18))),
                    subjectQuality = SolexaQuality(rep(c(22L, 12L), times = c(40, 20))),
                    scoreOnly = TRUE)

  # Now assume can't distinguish between C/T and G/A
  pairwiseAlignment(s1, s2,
                    patternQuality = SolexaQuality(rep(c(22L, 12L), times = c(36, 18))),
                    subjectQuality = SolexaQuality(rep(c(22L, 12L), times = c(40, 20))),
                    type = "local")
  mapping <- diag(4)
  dimnames(mapping) <- list(DNA_BASES, DNA_BASES)
  mapping["C", "T"] <- mapping["T", "C"] <- 1
  mapping["G", "A"] <- mapping["A", "G"] <- 1
  pairwiseAlignment(s1, s2,
                    patternQuality = SolexaQuality(rep(c(22L, 12L), times = c(36, 18))),
                    subjectQuality = SolexaQuality(rep(c(22L, 12L), times = c(40, 20))),
                    fuzzyMatrix = mapping,
                    type = "local")

  ## Amino acid global alignment
  pairwiseAlignment(AAString("PAWHEAE"), AAString("HEAGAWGHEE"),
                    substitutionMatrix = "BLOSUM50",
                    gapOpening = 0, gapExtension = 8)
}
\keyword{models}
\keyword{methods}
