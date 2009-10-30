\name{matchPWM}

\alias{maxWeights}
\alias{maxWeights,matrix-method}
\alias{minWeights}
\alias{minWeights,matrix-method}
\alias{maxScore}
\alias{maxScore,ANY-method}
\alias{minScore}
\alias{minScore,ANY-method}
\alias{unitScale}

\alias{PWM}
\alias{PWM,character-method}
\alias{PWM,DNAStringSet-method}

\alias{PWMscoreStartingAt}

\alias{matchPWM}
\alias{matchPWM,character-method}
\alias{matchPWM,DNAString-method}
\alias{matchPWM,XStringViews-method}
\alias{matchPWM,MaskedDNAString-method}
\alias{countPWM}
\alias{countPWM,character-method}
\alias{countPWM,DNAString-method}
\alias{countPWM,XStringViews-method}
\alias{countPWM,MaskedDNAString-method}

\alias{reverseComplement,matrix-method}

% Deprecated:
\alias{PWMscore}


\title{PWM creating, matching, and related utilities}

\description{
  Position Weight Matrix (PWM) creating, matching, and related utilities
  for DNA data. (PWM for amino acid sequences are not supported.)
}

\usage{
  PWM(x, type = c("log2probratio", "prob"),
      prior.params = c("A"=0.25, "C"=0.25, "G"=0.25, "T"=0.25))

  matchPWM(pwm, subject, min.score="80\%", ...)
  countPWM(pwm, subject, min.score="80\%", ...)
  PWMscoreStartingAt(pwm, subject, starting.at=1)

  ## Utility functions for basic manipulation of the Position Weight Matrix
  maxWeights(x)
  minWeights(x)
  maxScore(x)
  minScore(x)
  unitScale(x)
  \S4method{reverseComplement}{matrix}(x, ...)
}

\arguments{
  \item{x}{
    For \code{PWM} a character string or DNAStringSet whose elements all have
    the same number of characters.

    For \code{maxWeights}, \code{minWeights}, \code{maxScore}, \code{minScore},
    \code{unitScale} , and \code{reverseComplement} a numeric matrix with row
    names A, C, G and T representing a Position Weight Matrix.
  }
  \item{type}{
    The type of position weight matrix, either "log2probratio" or "prob".
    See Details section for more information.
  }
  \item{prior.params}{
    A positive numeric vector, which represents the parameters of the
    Dirichlet conjugate prior, with names A, C, G, and T.
    See Details section for more information.
  }
  \item{pwm}{
    A numeric matrix with row names A, C, G and T representing a Position
    Weight Matrix.
  }
  \item{subject}{
    An \link{DNAString}, \link{XStringViews} or \link{MaskedDNAString}
    object for \code{matchPWM} and \code{countPWM}.

    A \link{DNAString} object containing the subject sequence.
  }
  \item{min.score}{
    The minimum score for counting a match.
    Can be given as a character string containing a percentage (e.g.
    \code{"85\%"}) of the highest possible score or as a single number.
  }
  \item{starting.at}{
    An integer vector specifying the starting positions of the
    Position Weight Matrix relatively to the subject.
  }
  \item{...}{
    Additional arguments for methods.
  }
}

\details{
  The \code{PWM} function uses a multinomial model with a Dirichlet conjugate
  prior to calculate the estimated probability of base b at position i. As
  mentioned in the Arguments section, \code{prior.params} supplies the
  parameters for the DNA bases A, C, G, and T in the Dirichlet prior. These
  values result in a position independent initial estimate of the probabilities
  for the bases to be
  \code{priorProbs = prior.params/sum(prior.params)} and the
  posterior (data infused) estimate for the probabilities for the bases in each
  of the positions to be
  \code{postProbs = (consensusMatrix(x) + prior.params)/(length(x) + sum(prior.params))}.
  When \code{type = "log2probratio"}, the PWM = \code{unitScale(log2(postProbs/priorProbs))}.
  When \code{type = "prob"}, the PWM = \code{unitScale(postProbs)}.  
}

\value{
  A numeric matrix representing the Position Weight Matrix for \code{PWM}.

  A numeric vector containing the Position Weight Matrix-based scores
  for \code{PWMscoreStartingAt}.

  An \link{XStringViews} object for \code{matchPWM}.

  A single integer for \code{countPWM}.

  A vector containing the max weight for each position in \code{pwm}
  for \code{maxWeights}.

  A vector containing the min weight for each position in \code{pwm}
  for \code{minWeights}.

  The highest possible score for a given Position Weight Matrix for
  \code{maxScore}.

  The lowest possible score for a given Position Weight Matrix for
  \code{maxScore}.

  The modified numeric matrix given by
  \code{(x - minScore(x)/ncol(x))/(maxScore(x) - minScore(x))} for
  \code{unitScale}.

  A PWM obtained by reverting the column order in PWM \code{x} and by
  reassigning each row to its complementary nucleotide
  for \code{reverseComplement}.
}

\references{
Wasserman, WW, Sandelin, A., (2004) Applied bioinformatics for the identification of
regulatory elements, Nat Rev Genet., 5(4):276-87.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{reverseComplement}},
  \link{DNAString-class},
  \link{XStringViews-class}
}

\examples{
  ## Data setup
  data(HNF4alpha)
  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chr3R <- Dmelanogaster$chr3R
  chr3R

  ## Create a PWM and perform some general routines
  pwm <- PWM(HNF4alpha)
  round(pwm, 2)
  maxWeights(pwm)
  maxScore(pwm)
  reverseComplement(pwm)

  ## Score the first 5 positions
  PWMscoreStartingAt(pwm, unmasked(chr3R), starting.at=1:5)

  ## Match the plus strand
  matchPWM(pwm, chr3R)
  countPWM(pwm, chr3R)

  ## Match the minus strand
  matchPWM(reverseComplement(pwm), chr3R)
}

\keyword{methods}
\keyword{manip}
\keyword{utilities}
