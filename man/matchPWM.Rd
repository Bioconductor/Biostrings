\name{matchPWM}

\alias{maxWeights}
\alias{maxScore}
\alias{PWMscoreStartingAt}
\alias{matchPWM}
\alias{countPWM}

\alias{reverseComplement,matrix-method}

% Deprecated:
\alias{PWMscore}


\title{A simple PWM matching function and related utilities}

\description{
  A function implementing a simple algorithm for matching a set of
  patterns represented by a Position Weight Matrix (PWM) to a DNA sequence.
  PWM for amino acid sequences are not supported.
}

\usage{
  matchPWM(pwm, subject, min.score="80\%")
  countPWM(pwm, subject, min.score="80\%")
  PWMscoreStartingAt(pwm, subject, starting.at=1)

  ## Utility functions for basic manipulation of the Position Weight Matrix
  maxWeights(pwm)
  maxScore(pwm)
  \S4method{reverseComplement}{matrix}(x, ...)
}

\arguments{
  \item{pwm, x}{
    A Position Weight Matrix (numeric matrix with row names A, C, G and T).
  }
  \item{subject}{
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
    Additional arguments are currently ignored by the \code{reverseComplement}
    method for matrix objects.
  }
}

\value{
  An \link{XStringViews} object for \code{matchPWM}.

  A single integer for \code{countPWM}.

  A numeric vector containing the Position Weight Matrix-based scores
  for \code{PWMscoreStartingAt}.

  A vector containing the max weight for each position in \code{pwm}
  for \code{maxWeights}.

  The highest possible score for a given Position Weight Matrix for
  \code{maxScore}.

  A PWM obtained by reverting the column order in PWM \code{x} and by
  reassigning each row to its complementary nucleotide
  for \code{reverseComplement}.
}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{reverseComplement}},
  \link{DNAString-class},
  \link{XStringViews-class}
}

\examples{
  pwm <- rbind(A=c( 1,  0, 19, 20, 18,  1, 20,  7),
               C=c( 1,  0,  1,  0,  1, 18,  0,  2),
               G=c(17,  0,  0,  0,  1,  0,  0,  3),
               T=c( 1, 20,  0,  0,  0,  1,  0,  8))
  maxWeights(pwm)
  maxScore(pwm)
  reverseComplement(pwm)

  subject <- DNAString("AGTAAACAA")
  PWMscoreStartingAt(pwm, subject, starting.at=c(2:1, NA))

  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chr3R <- unmasked(Dmelanogaster$chr3R)
  chr3R

  ## Match the plus strand
  matchPWM(pwm, chr3R)
  countPWM(pwm, chr3R)

  ## Match the minus strand
  matchPWM(reverseComplement(pwm), chr3R)
}

\keyword{methods}
\keyword{manip}
\keyword{utilities}
