\name{XStringAlign-class}
\docType{class}

% Classes
\alias{class:XStringAlign}
\alias{XStringAlign-class}
\alias{XStringAlign}

% Constructor-like functions and generics:
\alias{initialize,XStringAlign-method}

% Accesor methods:
\alias{align1}
\alias{align1,XStringAlign-method}
\alias{align2}
\alias{align2,XStringAlign-method}
\alias{type}
\alias{type,XStringAlign-method}
\alias{score}
\alias{score,XStringAlign-method}
\alias{profile1}
\alias{profile1,XStringAlign-method}
\alias{profile2}
\alias{profile2,XStringAlign-method}
\alias{length,XStringAlign-method}
\alias{nchar,XStringAlign-method}
\alias{alphabet,XStringAlign-method}

% Standard generic methods:
\alias{show,XStringAlign-method}
\alias{as.character,XStringAlign-method}


\title{XStringAlign objects}

\description{
  The \code{XStringAlign} class is a container for storing
  an alignment between 2 \link{XString} objects of the same subtype.
}

\details{
  Before we define the notion of alignment, we introduce the notion of
  "filled-with-gaps subsequence". A "filled-with-gaps subsequence" of
  a string string1 is obtained by inserting 0 or any number of gaps in
  a subsequence of s1. For example L-A--ND and A--N-D are
  "filled-with-gaps subsequences" of LAND. An alignment between two strings
  string1 and string2 results in two strings (align1 and align2) that have
  the same length and are "filled-with-gaps subsequences" of string1 and string2.

  For example, this is an alignment between LAND and LEAVES:
  \preformatted{
    L-A
    LEA
  }

  An alignment can be seen as a compact representation of one set of basic
  operations that transforms string1 into align1. There are 3 different kinds
  of basic operations: "insertions" (gaps in align1), "deletions" (gaps in
  align2),  "replacements".
  The above alignment represents the following basic operations:
  \preformatted{
    insert E at pos 2
    insert V at pos 4
    insert E at pos 5
    replace by S at pos 6 (N is replaced by S)
    delete at pos 7 (D is deleted)
  }
  Note that "insert X at pos i" means that all letters at a position >= i
  are moved 1 place to the right before X is actually inserted.

  There are many possible alignments between two given strings string1 and
  string2 and a common problem is to find the one (or those ones) with the
  highest score, i.e. with the lower total cost in terms of basic operations.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a \code{XStringAlign} object.

  \describe{
    \item{}{
      \code{align1(x)} and \code{align2(x)}:
      The "filled-with-gaps subsequences" of the original strings to align.
      Note that \code{align1(x)} and \code{align2(x)} are \link{XString}
      objects of the same subtype and length.
    }
    \item{}{
      \code{type(x)}:
      The type of the alignment (\code{"global"}, \code{"local"},
      \code{"overlap"}, \code{"overlap1"}, or \code{"overlap2"}).
    }
    \item{}{
      \code{score(x)}:
      The score of the alignment (integer).
    }
    \item{}{
      \code{profile1(x)} and \code{profile2(x)}:
      A numeric vector whose values are the maximum score that would result
      by starting the alignment at the corresponding location of the original
      string. Note that \code{length(profile1(x)) = nchar(string1)} and
      \code{length(profile1(x)) = nchar(string1)}.
    }
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      The length of the alignment i.e. the common length of \code{align1(x)}
      and \code{align2(x)}.
    }
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(align1(x))} (or \code{alphabet(align2(x))}).
    }
    \item{}{
      \code{as.character(x)}:
      Converts \code{x} to a named character vector of length 2.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{pairwiseAlignment}},
  \link{XString-class}
}

\examples{
  s1 <- AAString("LAND")
  s2 <- AAString("LEAVES")
  nw1 <- pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM50", gapOpening = -3, gapExtension = -1)
  nw1
  length(nw1)
  nw0 <- pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM50", gapOpening = 0, gapExtension = 0)
  nw0
  length(nw0)
  ## Low gap penalties tend to produce longer alignments!

  as.character(nw0)
}

\keyword{methods}
\keyword{classes}
