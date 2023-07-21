\name{XStringSet-comparison}

\alias{XStringSet-comparison}

\alias{pcompare,XStringSet,XStringSet-method}
\alias{pcompare,XStringSet,ANY-method}
\alias{pcompare,ANY,XStringSet-method}

\alias{match,XStringSet,XStringSet-method}
\alias{match,XStringSet,ANY-method}
\alias{match,ANY,XStringSet-method}

\alias{is.na,XStringSet-method}
\alias{anyNA,XStringSet-method}

\title{Comparing and ordering the elements in one or more XStringSet objects}

\description{
  Methods for comparing and ordering the elements in one or more
  \link{XStringSet} objects.
}

\details{
  Element-wise (aka "parallel") comparison of 2 \link{XStringSet} objects
  is based on the lexicographic order between 2 \link{BString},
  \link{DNAString}, \link{RNAString}, or \link{AAString} objects.

  For \link{DNAStringSet} and \link{RNAStringSet} objects, the letters in
  the respective alphabets (i.e. \link{DNA_ALPHABET} and \link{RNA_ALPHABET})
  are ordered based on a predefined code assigned to each letter. The code
  assigned to each letter can be retrieved with:
  \preformatted{
  dna_codes <- as.integer(DNAString(paste(DNA_ALPHABET, collapse="")))
  names(dna_codes) <- DNA_ALPHABET

  rna_codes <- as.integer(RNAString(paste(RNA_ALPHABET, collapse="")))
  names(rna_codes) <- RNA_ALPHABET
  }
  Note that this order does NOT depend on the locale in use. Also note that
  comparing DNA sequences with RNA sequences is supported and in that case T
  and U are considered to be the same letter.

  For \link{BStringSet} and \link{AAStringSet} objects, the alphabetical
  order is defined by the C collation. Note that, at the moment,
  \link{AAStringSet} objects are treated like \link{BStringSet} objects
  i.e. the alphabetical order is NOT defined by the order of the letters
  in \link{AA_ALPHABET}. This might change at some point.
}

\section{\code{pcompare()} and related methods}{
  In the code snippets below,
  \code{x} and \code{y} are \link{XStringSet} objects.

  \describe{
    \item{\code{pcompare(x, y)}:}{
      Performs element-wise (aka "parallel") comparison of \code{x} and
      \code{y}, that is, returns an integer vector where the i-th element
      is less than, equal to, or greater than zero if the i-th element in
      \code{x} is considered to be respectively less than, equal to, or
      greater than the i-th element in \code{y}.
      If \code{x} and \code{y} don't have the same length, then the shortest
      is recycled to the length of the longest (the standard recycling rules
      apply).
    }
    \item{\code{x == y}, \code{x != y}, \code{x <= y}, \code{x >= y},
      \code{x < y}, \code{x > y}:}{
      Equivalent to \code{pcompare(x, y) == 0}, \code{pcompare(x, y) != 0},
      \code{pcompare(x, y) <= 0}, \code{pcompare(x, y) >= 0},
      \code{pcompare(x, y) < 0}, and \code{pcompare(x, y) > 0}, respectively.
    }
  }
}

\section{\code{order()} and related methods}{
  In the code snippets below, \code{x} is an \link{XStringSet} object.

  \describe{
    \item{\code{is.unsorted(x, strictly=FALSE)}:}{
      Return a logical values specifying if \code{x} is unsorted. The
      \code{strictly} argument takes logical value indicating if the check
      should be for _strictly_ increasing values.
    }
    \item{\code{order(x, decreasing=FALSE)}:}{
      Return a permutation which rearranges \code{x} into ascending or
      descending order.
    }
    \item{\code{rank(x, ties.method=c("first", "min"))}:}{
      Rank \code{x} in ascending order.
    }
    \item{\code{sort(x, decreasing=FALSE)}:}{
      Sort \code{x} into ascending or descending order.
    }
  }
}

\section{\code{duplicated()} and \code{unique()}}{
  In the code snippets below, \code{x} is an \link{XStringSet} object.

  \describe{
    \item{\code{duplicated(x)}:}{
      Return a logical vector whose elements denotes duplicates in \code{x}.
    }
    \item{\code{unique(x)}:}{
      Return the subset of \code{x} made of its unique elements.
    }
  }
}

\section{\code{match()} and \code{\%in\%}}{
  In the code snippets below,
  \code{x} and \code{table} are \link{XStringSet} objects.

  \describe{
    \item{\code{match(x, table, nomatch=NA_integer_)}:}{
      Returns an integer vector containing the first positions of an identical
      match in \code{table} for the elements in \code{x}.
    }
    \item{\code{x \%in\% table}:}{
      Returns a logical vector indicating which elements in \code{x} match
      identically with an element in \code{table}.
    }
  }
}

\section{\code{is.na()} and related methods}{
  In the code snippets below, \code{x} is an \link{XStringSet}
  object. An \code{XStringSet} object never contains missing values
  (these methods exist for compatibility).

  \describe{
    \item{\code{is.na(x)}:}{ Returns \code{FALSE} for every element.
    }
    \item{\code{anyNA(x)}:}{ Returns \code{FALSE}.
    }
  }
}

\author{H. Pagès}

\seealso{
  \link{XStringSet-class},
  \code{\link{==}},
  \code{\link{is.unsorted}},
  \code{\link{order}},
  \code{\link{rank}},
  \code{\link{sort}},
  \code{\link{duplicated}},
  \code{\link{unique}},
  \code{\link{match}},
  \code{\link{\%in\%}}
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLES
## ---------------------------------------------------------------------

dna <- DNAStringSet(c("AAA", "TC", "", "TC", "AAA", "CAAC", "G"))
match(c("", "G", "AA", "TC"), dna)

library(drosophila2probe)
fly_probes <- DNAStringSet(drosophila2probe)
sum(duplicated(fly_probes))  # 481 duplicated probes

is.unsorted(fly_probes)  # TRUE
fly_probes <- sort(fly_probes)
is.unsorted(fly_probes)  # FALSE
is.unsorted(fly_probes, strictly=TRUE)  # TRUE, because of duplicates
is.unsorted(unique(fly_probes), strictly=TRUE)  # FALSE

## Nb of probes that are the reverse complement of another probe:
nb1 <- sum(reverseComplement(fly_probes) \%in\% fly_probes)
stopifnot(identical(nb1, 455L))  # 455 probes

## Probes shared between drosophila2probe and hgu95av2probe:
library(hgu95av2probe)
human_probes <- DNAStringSet(hgu95av2probe)
m <- match(fly_probes, human_probes)
stopifnot(identical(sum(!is.na(m)), 493L))  # 493 shared probes

## ---------------------------------------------------------------------
## B. AN ADVANCED EXAMPLE
## ---------------------------------------------------------------------
## We want to compare the first 5 bases with the 5 last bases of each
## probe in drosophila2probe. More precisely, we want to compute the
## percentage of probes for which the first 5 bases are the reverse
## complement of the 5 last bases.

library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)

first5 <- narrow(probes, end=5)
last5 <- narrow(probes, start=-5)
nb2 <- sum(first5 == reverseComplement(last5))
stopifnot(identical(nb2, 17L))

## Percentage:
100 * nb2 / length(probes)  # 0.0064 %

## If the probes were random DNA sequences, a probe would have 1 chance
## out of 4^5 to have this property so the percentage would be:
100 / 4^5  # 0.098 %

## With randomly generated probes:
set.seed(33)
random_dna <- sample(DNAString(paste(DNA_BASES, collapse="")),
                     sum(width(probes)), replace=TRUE)
random_probes <- successiveViews(random_dna, width(probes))
random_probes
random_probes <- as(random_probes, "XStringSet")
random_probes

random_first5 <- narrow(random_probes, end=5)
random_last5 <- narrow(random_probes, start=-5)

nb3 <- sum(random_first5 == reverseComplement(random_last5))
100 * nb3 / length(random_probes)  # 0.099 %
}

\keyword{methods}
