\name{findPalindromes}

\alias{findPalindromes}
\alias{findPalindromes,XString-method}
\alias{findPalindromes,DNAString-method}
\alias{findPalindromes,RNAString-method}
\alias{findPalindromes,XStringViews-method}
\alias{findPalindromes,MaskedXString-method}

\alias{palindromeArmLength}
\alias{palindromeArmLength,XString-method}
\alias{palindromeArmLength,DNAString-method}
\alias{palindromeArmLength,RNAString-method}
\alias{palindromeArmLength,XStringViews-method}

\alias{palindromeLeftArm}
\alias{palindromeRightArm}
\alias{palindromeLeftArm,XString-method}
\alias{palindromeRightArm,XString-method}
\alias{palindromeLeftArm,XStringViews-method}
\alias{palindromeRightArm,XStringViews-method}

% Old stuff (Deprecated or Defunct):
\alias{findComplementedPalindromes}
\alias{findComplementedPalindromes,DNAString-method}
\alias{findComplementedPalindromes,RNAString-method}
\alias{findComplementedPalindromes,XStringViews-method}
\alias{findComplementedPalindromes,MaskedXString-method}
\alias{complementedPalindromeArmLength}
\alias{complementedPalindromeArmLength,DNAString-method}
\alias{complementedPalindromeArmLength,RNAString-method}
\alias{complementedPalindromeArmLength,XStringViews-method}
\alias{complementedPalindromeLeftArm}
\alias{complementedPalindromeLeftArm,XString-method}
\alias{complementedPalindromeLeftArm,XStringViews-method}
\alias{complementedPalindromeRightArm}
\alias{complementedPalindromeRightArm,XString-method}
\alias{complementedPalindromeRightArm,XStringViews-method}


\title{Searching a sequence for palindromes}

\description{
  The \code{findPalindromes} function can be used to find palindromic
  regions in a sequence.

  \code{palindromeArmLength}, \code{palindromeLeftArm}, and
  \code{palindromeRightArm} are utility functions for operating on
  palindromic sequences.
}

\usage{
findPalindromes(subject, min.armlength=4,
                max.looplength=1, min.looplength=0, max.mismatch=0)
palindromeArmLength(x, max.mismatch=0, ...)
palindromeLeftArm(x, max.mismatch=0, ...)
palindromeRightArm(x, max.mismatch=0, ...)
}

\arguments{
  \item{subject}{
    An \link{XString} object containing the subject string,
    or an \link{XStringViews} object.
  }
  \item{min.armlength}{
    An integer giving the minimum length of the arms of the palindromes
    to search for.
  }
  \item{max.looplength}{
    An integer giving the maximum length of "the loop" (i.e the sequence
    separating the 2 arms) of the palindromes to search for.
    Note that by default (\code{max.looplength=1}), \code{findPalindromes}
    will search for strict palindromes only.
  }
  \item{min.looplength}{
    An integer giving the minimum length of "the loop" of the palindromes
    to search for.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed between the 2 arms of
    the palindromes to search for.
  }
  \item{x}{
    An \link{XString} object containing a 2-arm palindrome, or an
    \link{XStringViews} object containing a set of 2-arm palindromes.
  }
  \item{...}{
    Additional arguments to be passed to or from methods.
  }
}

\details{
  The \code{findPalindromes} function finds palindromic substrings in a subject
  string. The palindromes that can be searched for are either strict
  palindromes or 2-arm palindromes (the former being a particular case of
  the latter) i.e. palindromes where the 2 arms are separated by an arbitrary
  sequence called "the loop".

  If the subject string is a nucleotide sequence (i.e. DNA or RNA), the 2 arms
  must contain sequences that are reverse complement from each other.
  Otherwise, they must contain sequences that are the same.
}

\value{
  \code{findPalindromes} returns an \link{XStringViews} object containing all
  palindromes found in \code{subject} (one view per palindromic substring
  found).

  \code{palindromeArmLength} returns the arm length (integer) of the 2-arm
  palindrome \code{x}. It will raise an error if \code{x} has no arms. Note
  that any sequence could be considered a 2-arm palindrome if we were OK with
  arms of length 0 but we are not: \code{x} must have arms of length greater
  or equal to 1 in order to be considered a 2-arm palindrome.
  When applied to an \link{XStringViews} object \code{x},
  \code{palindromeArmLength} behaves in a vectorized fashion by returning
  an integer vector of the same length as \code{x}.

  \code{palindromeLeftArm} returns an object of the same class as the original
  object \code{x} and containing the left arm of \code{x}.

  \code{palindromeRightArm} does the same as \code{palindromeLeftArm} but on
  the right arm of \code{x}.

  Like \code{palindromeArmLength}, both \code{palindromeLeftArm} and
  \code{palindromeRightArm} will raise an error if \code{x} has no arms.
  Also, when applied to an \link{XStringViews} object \code{x}, both behave
  in a vectorized fashion by returning an \link{XStringViews} object of the
  same length as \code{x}.
}

\author{H. Pages}

\seealso{
  \code{\link{maskMotif}},
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \code{\link{matchProbePair}},
  \link{XStringViews-class},
  \link{DNAString-class}
}

\examples{
x0 <- BString("abbbaabbcbbaccacabbbccbcaabbabacca")

pals0a <- findPalindromes(x0, min.armlength=3, max.looplength=5)
pals0a
palindromeArmLength(pals0a)
palindromeLeftArm(pals0a)
palindromeRightArm(pals0a)

pals0b <- findPalindromes(x0, min.armlength=9, max.looplength=5,
                          max.mismatch=3)
pals0b
palindromeArmLength(pals0b, max.mismatch=3)
palindromeLeftArm(pals0b, max.mismatch=3)
palindromeRightArm(pals0b, max.mismatch=3)

## Whitespaces matter:
x1 <- BString("Delia saw I was aileD")
palindromeArmLength(x1)
palindromeLeftArm(x1)
palindromeRightArm(x1)

x2 <- BString("was it a car or a cat I saw")
palindromeArmLength(x2)
palindromeLeftArm(x2)
palindromeRightArm(x2)

## On a DNA or RNA sequence:
x3 <- DNAString("CCGAAAACCATGATGGTTGCCAG")
findPalindromes(x3)
findPalindromes(RNAString(x3))

## Note that palindromes can be nested:
x4 <- DNAString("ACGTTNAACGTCCAAAATTTTCCACGTTNAACGT")
findPalindromes(x4, max.looplength=19)

## A real use case:
library(BSgenome.Dmelanogaster.UCSC.dm3)
chrX <- Dmelanogaster$chrX
chrX_pals0 <- findPalindromes(chrX, min.armlength=40, max.looplength=80)
chrX_pals0
palindromeArmLength(chrX_pals0)  # 251 70 262

## Allowing up to 2 mismatches between the 2 arms:
chrX_pals2 <- findPalindromes(chrX, min.armlength=40, max.looplength=80,
                              max.mismatch=2)
chrX_pals2
palindromeArmLength(chrX_pals2, max.mismatch=2)  # 254 77 44 48 40 264
}

\keyword{methods}
