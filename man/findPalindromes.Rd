\name{findPalindromes}

\alias{findPalindromes}
\alias{findPalindromes,XString-method}
\alias{findPalindromes,BStringViews-method}

\alias{palindromeArmLength}
\alias{palindromeArmLength,XString-method}
\alias{palindromeArmLength,BStringViews-method}

\alias{palindromeLeftArm}
\alias{palindromeLeftArm,XString-method}
\alias{palindromeLeftArm,BStringViews-method}

\alias{palindromeRightArm}
\alias{palindromeRightArm,XString-method}
\alias{palindromeRightArm,BStringViews-method}

\alias{findComplementedPalindromes}
\alias{findComplementedPalindromes,DNAString-method}
\alias{findComplementedPalindromes,BStringViews-method}

\alias{complementedPalindromeArmLength}
\alias{complementedPalindromeArmLength,DNAString-method}
\alias{complementedPalindromeArmLength,BStringViews-method}

\alias{complementedPalindromeLeftArm}
\alias{complementedPalindromeLeftArm,DNAString-method}
\alias{complementedPalindromeLeftArm,BStringViews-method}

\alias{complementedPalindromeRightArm}
\alias{complementedPalindromeRightArm,DNAString-method}
\alias{complementedPalindromeRightArm,BStringViews-method}


\title{Searching a sequence for palindromes or complemented palindromes}

\description{
  The \code{findPalindromes} function finds palindromic substrings in a subject
  string. The palindromes that can be searched for are either strict palindromes
  or 2-arm palindromes (the former being a particular case of the latter) i.e.
  palindromes where the 2 arms are separated by an arbitrary sequence called
  "the loop".
  Use the \code{findComplementedPalindromes} function to find complemented
  palindromic substrings in a \link{DNAString} subject (in a complemented
  palindrome the 2 arms are reverse-complementary sequences).
}

\usage{
  findPalindromes(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
  palindromeArmLength(x, max.mismatch=0, ...)
  palindromeLeftArm(x, max.mismatch=0, ...)
  palindromeRightArm(x, max.mismatch=0, ...)

  findComplementedPalindromes(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
  complementedPalindromeArmLength(x, max.mismatch=0, ...)
  complementedPalindromeLeftArm(x, max.mismatch=0, ...)
  complementedPalindromeRightArm(x, max.mismatch=0, ...)
}

\arguments{
  \item{subject}{
    An \link{XString} object containing the subject string,
    or a \link{BStringViews} object.
  }
  \item{min.armlength}{
    An integer giving the minimum length of the arms of the palindromes
    (or complemented palindromes) to search for.
  }
  \item{max.looplength}{
    An integer giving the maximum length of "the loop" (i.e the sequence
    separating the 2 arms) of the palindromes (or complemented palindromes)
    to search for.
    Note that by default (\code{max.looplength=1}), \code{findPalindromes}
    will search for strict palindromes (or complemented palindromes) only.
  }
  \item{min.looplength}{
    An integer giving the minimum length of "the loop" of the palindromes
    (or complemented palindromes) to search for.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed between the 2 arms of
    the palindromes (or complemented palindromes) to search for.
  }
  \item{x}{
    An \link{XString} object containing a 2-arm palindrome
    or complemented palindrome, or a \link{BStringViews} object containing a set
    of 2-arm palindromes or complemented palindromes.
  }
  \item{...}{
    Additional arguments to be passed to or from methods.
  }
}

\details{
}

\value{
  \code{findPalindromes} and \code{findComplementedPalindromes} return a
  \link{BStringViews} object containing all palindromes (or complemented
  palindromes) found in \code{subject} (one view per palindromic substring
  found).

  \code{palindromeArmLength} and \code{complementedPalindromeArmLength} return
  the arm length (integer) of the 2-arm palindrome (or complemented palindrome)
  \code{x}.
  It will raise an error if \code{x} has no arms. Note that any sequence could
  be considered a 2-arm palindrome if we were OK with arms of length 0 but we
  are not: \code{x} must have arms of length greater or equal to 1 in order to
  be considered a 2-arm palindrome. The same apply to 2-arm complemented
  palindromes.
  When applied to a \link{BStringViews} object \code{x},
  \code{palindromeArmLength} and \code{complementedPalindromeArmLength} behave
  in a vectorized fashion by returning an integer vector of the same length
  as \code{x}.

  \code{palindromeLeftArm} and \code{complementedPalindromeLeftArm} return an
  object of the same class as the original object \code{x} and containing the
  left arm of \code{x}.

  \code{palindromeRightArm} does the same as \code{palindromeLeftArm} but on
  the right arm of \code{x}.

  Like \code{palindromeArmLength}, both \code{palindromeLeftArm} and
  \code{palindromeRightArm} will raise an error if \code{x} has no arms.
  Also, when applied to a \link{BStringViews} object \code{x}, both behave
  in a vectorized fashion by returning a \link{BStringViews} object of the
  same length as \code{x}.
}

\author{H. Pages}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \code{\link{matchProbePair}},
  \code{\link{masks}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## Note that complemented palindromes (like palindromes) can be nested
  findComplementedPalindromes(DNAString("ACGTTNAACGT-ACGTTNAACGT"))

  ## A real use case
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrX <- Dmelanogaster[["X"]]
  masks(chrX) <- "N"
  chrX_pals <- findComplementedPalindromes(chrX, min.armlength=50, max.looplength=20)
  complementedPalindromeArmLength(chrX_pals) # 251

  ## Of course, whitespaces matter
  palindromeArmLength(BString("was it a car or a cat I saw"))

  ## Note that the 2 arms of a strict palindrome (or strict complemented
  ## palindrome) are equal to the full sequence.
  palindromeLeftArm(BString("Delia saw I was aileD"))
  complementedPalindromeLeftArm(DNAString("N-ACGTT-AACGT-N"))
  palindromeLeftArm(DNAString("N-AAA-N-N-TTT-N"))
}

\keyword{methods}

