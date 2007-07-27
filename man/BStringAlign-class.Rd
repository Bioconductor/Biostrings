\name{BStringAlign-class}
\docType{class}

% Classes
\alias{class:BStringAlign}
\alias{BStringAlign-class}
\alias{BStringAlign}

% Constructor-like functions and generics:
\alias{initialize,BStringAlign-method}

% Accesor methods:
\alias{align1}
\alias{align1,BStringAlign-method}
\alias{align2}
\alias{align2,BStringAlign-method}
\alias{score}
\alias{score,BStringAlign-method}
\alias{length,BStringAlign-method}
\alias{nchar,BStringAlign-method}
\alias{alphabet,BStringAlign-method}

% Standard generic methods:
\alias{show,BStringAlign-method}


\title{The BStringAlign class}

\description{
  The \code{BStringAlign class} is a container for storing
  an alignment between 2 \code{\link{BString}} (or derived) objects.
}

\details{
  Before we define the notion of alignment, we introduce
  the notion of "filled-with-gaps supersequence".
  A "filled-with-gaps supersequence" of a string s1 is a string S1 that is
  obtained by inserting 0 or any number of gaps in s1.
  For example L-A--ND is a "filled-with-gaps supersequence" of LAND.
  An alignment between 2 strings s1 and s2 is made of 2 strings align1 and
  align2 that are "filled-with-gaps supersequences" of s1 and s2, and that
  have the same length.
  Note that this common length must be greater or equal to the lengths
  of s1 and s2: nchar(align1) = nchar(align2) >= max(nchar(s1), nchar(s2))

  For example, this is an alignment between LAND and LEAVES:
  \preformatted{
    L-A--ND
    LEAVES-
  }

  An alignment can be seen as a compact representation of one set of basic
  operations that transforms s1 into s2. There are 3 different kinds of
  basic operations: "insertions" (gaps in align1), "deletions" (gaps in
  align2),  "replacements".
  The above alignement represents the following basic operations:
  \preformatted{
    insert E at pos 2
    insert V at pos 4
    insert E at pos 5
    replace by S at pos 6 (N is replaced by S)
    delete at pos 7 (D is deleted)
  }
  Note that "insert X at pos i" means that all letters at a position >= i
  are moved 1 place to the right before X is actually inserted.

  There are many possible alignments between 2 given strings s1 and s2
  and a common problem is to find the one (or those ones) with the highest
  score i.e. with the lower total cost in terms of basic operations.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a \code{BStringAlign} object.

  \describe{
    \item{}{
      \code{align1(x)} and \code{align2(x)}:
      The "filled-with-gaps supersequences" of the original strings to align.
      Note that \code{align1(x)} and \code{align2(x)} are \code{\link{BString}}
      (or derived) objects of the same class and of the same length.
    }
    \item{}{
      \code{score(x)}:
      The score of the alignment (integer).
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
  }
}

\author{H. Pages}

\seealso{
  \code{\link{needwunsQS}},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class}
}

\examples{
  s1 <- AAString("LAND")
  s2 <- AAString("LEAVES")
  ## With the needwunsQS function, the cost of an insertion or deletion
  ## is controlled by the gappen (gap penalty) arg, the cost of a replacement
  ## is controlled by the "substitution scoring matrix" passed thru the substmat
  ## arg
  nw1 <- needwunsQS(s1, s2, substmat="BLOSUM50", gappen=1)
  nw1
  length(nw1)
  nw0 <- needwunsQS(s1, s2, substmat="BLOSUM50", gappen=0)
  nw0
  length(nw0)
  ## Low gap penalties tend to produce longer alignments!
}

\keyword{methods}
\keyword{classes}
