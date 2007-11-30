\name{findPalindromes}

\alias{findPalindromes}
\alias{findPalindromes,BString-method}
\alias{findPalindromes,BStringViews-method}

\alias{palindromeArmLength}
\alias{palindromeArmLength,BString-method}
\alias{palindromeArmLength,BStringViews-method}

\alias{palindromeLeftArm}
\alias{palindromeLeftArm,BString-method}
\alias{palindromeLeftArm,BStringViews-method}

\alias{palindromeRightArm}
\alias{palindromeRightArm,BString-method}
\alias{palindromeRightArm,BStringViews-method}


\title{Searching a sequence for palindromes}

\description{
  The \code{findPalindromes} function finds palindromic substrings in a subject
  string. The palindromes that can be searched for are either strict palindromes
  or 2-arm palindromes (the former being a particular case of the latter).
  With a \link{DNAString} or \link{RNAString} subject, they can also be
  antipalindromes i.e. the 2 arms are reverse-complementary sequences.
}

\usage{
  findPalindromes(subject, min.armlength=4, max.ngaps=1, anti=NULL)
  palindromeArmLength(x, anti=NULL, ...)
  palindromeLeftArm(x, anti=NULL, ...)
  palindromeRightArm(x, anti=NULL, ...)
}

\arguments{
  \item{subject}{
    A \link{BString} (or derived) object containing the subject string,
    or a \link{BStringViews} object.
  }
  \item{min.armlength}{
    An integer giving the minimum length of the arms of the palindromes
    (or antipalindromes) to search for.
  }
  \item{max.ngaps}{
    An integer giving the maximum length of the sequence separating the
    2 arms of the palindromes (or antipalindromes) to search for.
    Note that by default (\code{max.ngaps=1}), \code{findPalindromes} will
    search for strict palindromes (or antipalindromes) only.
  }
  \item{anti}{
    Either \code{NULL}, \code{TRUE} or \code{FALSE}.
    Values \code{anti=TRUE} or \code{anti=FALSE} can only be used with
    a \link{DNAString} or \link{RNAString} subject to specify wether or not
    the substrings to search for must be antipalindromes.
    By default (\code{anti=NULL}), antipalindromes are searched for.
  }
  \item{x}{
    A \link{BString} (or derived) object containing a 2-arm palindrome
    or antipalindrome, or a \link{BStringViews} object containing a set
    of 2-arm palindromes or antipalindromes.
  }
  \item{...}{
    Additional arguments to be passed to or from methods.
  }
}

\details{
}

\value{
  \code{findPalindromes} returns a \link{BStringViews} object containing all
  palindromes (or antipalindromes) found in \code{subject} (one view per
  palindromic substring found).

  \code{palindromeArmLength} returns the arm length (integer) of the 2-arm
  palindrome or antipalindrome \code{x}.
  It will raise an error if \code{x} has no arms. Note that any sequence could
  be considered a 2-arm palindrome if we were OK with arms of length 0 but we
  are not: \code{x} must have arms of length greater or equal to 1 in order to
  be considered a 2-arm palindrome. The same apply to 2-arm antipalindromes.
  When applied to a \link{BStringViews} object \code{x},
  \code{palindromeArmLength} behaves in a vectorized fashion by returning an
  integer vector of the same length as \code{x}.

  \code{palindromeLeftArm} returns an object of the same class as the original
  object \code{x} and containing the left arm of \code{x}.

  \code{palindromeRightArm} does the same as \code{palindromeLeftArm} but on
  the right arm of \code{x}.

  Like \code{palindromeArmLength}, both \code{palindromeLeftArm} and
  \code{palindromeRightArm} will raise an error if \code{x} has no arms.
  Also, when applied to a \link{BStringViews} object \code{x}, both behave
  in a vectorized fashion by returning a \link{BStringViews} object of the
  same length as \code{x}.
}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \code{\link{matchProbePair}},
  \code{\link{mask}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## Note that palindromes can be nested
  findPalindromes(DNAString("ACGTTNAACGT-ACGTTNAACGT"))

  ## A real use case
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrX <- Dmelanogaster[["X"]]
  noN_chrX <- mask(chrX, "N")
  chrX_pals <- findPalindromes(noN_chrX, min.armlength=50, max.ngaps=20)
  palindromeArmLength(chrX_pals) # 251

  ## Of course, spaces count
  palindromeArmLength(BString("was it a car or a cat I saw"))
  ## Note that the 2 arms of a strict palindrome (or antipalindrome) are equal
  ## to the full sequence.
  palindromeLeftArm(BString("Delia saw I was aileD"))
  palindromeLeftArm(DNAString("N-ACGTT-AACGT-N"))
  palindromeLeftArm(DNAString("N-AAA-N-N-TTT-N"), anti=FALSE)
}

\keyword{methods}

