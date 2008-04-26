\name{matchPWM}

\alias{maxWeights}
\alias{maxScore}
\alias{PWMscore}
\alias{matchPWM}
\alias{countPWM}

\alias{reverseComplement,matrix-method}


\title{A simple PWM matching function and related utilities}

\description{
  A function implementing a simple algorithm for matching a set of
  patterns represented by a Position Weight Matrix (PWM) to a DNA sequence.
  PWM for amino acid sequences are not supported.
}

\usage{
  matchPWM(pwm, subject, min.score="80\%")
  countPWM(pwm, subject, min.score="80\%")

  ## Utility functions for basic manipulation of the Position Weight Matrix
  maxWeights(pwm)
  maxScore(pwm)
  reverseComplement(x, ...) # S4 method for matrix objects
}

\arguments{
  \item{pwm}{
    A Position Weight Matrix (integer matrix with row names A, C, G and T).
  }
  \item{subject}{
    A \link{DNAString} object containing the subject sequence.
  }
  \item{min.score}{
    The minimum score for counting a match.
    Can be given as a percentage (e.g. \code{"85\%"}) of the highest possible
    score or as an integer.
  }
  \item{x}{
    A Position Weight Matrix.
  }
  \item{...}{
    Additional arguments are currently ignored by the \code{reverseComplement}
    method for matrix objects.
  }
}

\value{
  An \link{XStringViews} object for \code{matchPWM}.

  A single integer for \code{countPWM}.

  An integer vector containing the max weight for each position in \code{pwm}
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
  PWMscore(pwm, subject, c(2:1, NA))

  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chr3R <- Dmelanogaster[["3R"]]
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

