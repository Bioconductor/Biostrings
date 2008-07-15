\name{matchPDict}

\alias{matchPDict-exact}

% Classes
\alias{MIndex}
\alias{class:MIndex}
\alias{MIndex-class}
\alias{class:ByPos_MIndex}
\alias{ByPos_MIndex-class}
\alias{class:ByName_MIndex}
\alias{ByName_MIndex-class}

% Generics and methods:
\alias{startIndex}
\alias{startIndex,MIndex-method}
\alias{endIndex}
\alias{countIndex}
\alias{countIndex,MIndex-method}
\alias{[[,MIndex-method}
\alias{$,MIndex-method}
\alias{unlist}
\alias{unlist,MIndex-method}

\alias{length,ByPos_MIndex-method}
\alias{names,ByPos_MIndex-method}
\alias{show,ByPos_MIndex-method}
\alias{[[,ByPos_MIndex-method}
\alias{startIndex,ByPos_MIndex-method}
\alias{endIndex,ByPos_MIndex-method}
\alias{countIndex,ByPos_MIndex-method}

\alias{length,ByName_MIndex-method}
\alias{names,ByName_MIndex-method}
\alias{show,ByName_MIndex-method}
\alias{[[,ByName_MIndex-method}
\alias{startIndex,ByName_MIndex-method}
\alias{endIndex,ByName_MIndex-method}
\alias{countIndex,ByName_MIndex-method}

\alias{matchPDict}
\alias{matchPDict,XString-method}
\alias{matchPDict,XStringViews-method}
\alias{matchPDict,MaskedXString-method}
\alias{countPDict}
\alias{countPDict,XString-method}
\alias{countPDict,XStringViews-method}
\alias{countPDict,MaskedXString-method}
\alias{whichPDict}
\alias{whichPDict,XString-method}

% Functions:
\alias{extractAllMatches}


\title{Searching a sequence for patterns stored in a preprocessed dictionary}

\description{
  The \code{matchPDict}, \code{countPDict} and \code{whichPDict} functions
  efficiently find the occurrences in a text (the subject) of all patterns
  stored in a preprocessed dictionary.

  The three functions differ in what they return: \code{matchPDict} returns
  the "where" information i.e. the positions in the subject of all the
  occurrences of every pattern; \code{countPDict} returns the "how many
  times" information i.e. the number of occurrences for each pattern;
  and \code{whichPDict} returns the "who" information i.e. which patterns
  in the preprocessed dictionary have at least one match.

  This man page shows how to use \code{matchPDict}/\code{countPDict}/\code{whichPDict}
  for exact matching of a constant width dictionary i.e. a dictionary where all the
  patterns have the same length (same number of nucleotides).

  See \code{?`\link{matchPDict-inexact}`} for how to use these functions
  for inexact matching or when the original dictionary has a variable width.
}

\usage{
  matchPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
  countPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
  whichPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)

  ## Manipulation of the MIndex object returned by matchPDict()
  startIndex(x, all.names=FALSE)
  endIndex(x, all.names=FALSE)
  countIndex(x, all.names=FALSE)
  unlist(x, recursive=TRUE, use.names=TRUE)
  extractAllMatches(subject, mindex)
}

\arguments{
  \item{pdict}{
    A \link{PDict} object containing the preprocessed dictionary.
  }
  \item{subject}{
    An \link{XString} object containing the subject string.
    For now, only \link{XString} subjects of subtype \link{DNAString}
    are supported.
  }
  \item{algorithm}{
    Not supported yet.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed (see
    \code{?\link{isMatching}} for the details).
    When used on a constant width dictionary,
    \code{matchPDict}/\code{countPDict}/\code{whichPDict}
    only support exact matching so \code{max.mismatch} must be zero.
    See \code{?`\link{matchPDict-inexact}`} for inexact matching.
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{?\link{isMatching}} for the details).
    When used on a constant width dictionary,
    \code{matchPDict}/\code{countPDict}/\code{whichPDict}
    only support exact matching so \code{fixed} must be \code{TRUE}.
    See \code{?`\link{matchPDict-inexact}`} for inexact matching.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{x}{
    A PDict or MIndex object for \code{names}.
    An MIndex object for the other methods.
  }
  \item{all.names}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{recursive}{
    ignored.
  }
  \item{use.names}{
    ignored.
  }
  \item{mindex}{
    An MIndex object returned by a previous call to \code{matchPDict}.
  }
}

\details{
  In this man page, we assume that you know how to preprocess a dictionary
  of DNA patterns that can then be used with \code{matchPDict}, \code{countPDict}
  or \code{whichPDict}.
  Please see \code{?\link{PDict}} if you don't.

  When using \code{matchPDict}, \code{countPDict} or \code{whichPDict}
  for exact matching of a constant width dictionary, the standard way to
  preprocess the original dictionary is by calling the \code{\link{PDict}}
  constructor on it with no extra arguments. This returns the preprocessed
  dictionary in a \link{PDict} object that can be used with
  \code{matchPDict}/\code{countPDict}/\code{whichPDict}.
}

\value{
  \code{matchPDict} returns an MIndex object with \code{length} equal to the
  number of patterns in the \code{pdict} argument.

  \code{countPDict} returns an integer vector with \code{length} equal to the
  number of patterns in the \code{pdict} argument.

  \code{whichPDict} returns an integer vector made of the indices of the
  patterns in the \code{pdict} argument that have at least one match.

  \code{startIndex}, \code{endIndex} and \code{countIndex} return vectors of
  the same length as the original dictionary: \code{startIndex} and
  \code{endIndex} return a list of integer vectors, and \code{countIndex}
  returns an integer vector.

  \code{extractAllMatches} returns an \link{XStringViews} object with names.
}

\author{H. Pages}

\seealso{
  \link{PDict},
  \link{matchPDict-inexact},
  \code{\link{isMatching}},
  \code{\link{coverage}},
  \code{\link{matchPattern}},
  \code{\link{alphabetFrequency}},
  \link{XStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. WITH UNNAMED PATTERNS
  ## ---------------------------------------------------------------------

  ## Creating the pattern dictionary
  library(drosophila2probe)
  dict0 <- DNAStringSet(drosophila2probe$sequence)
  dict0                                # The original dictionary.
  length(dict0)                        # Hundreds of thousands of patterns.
  pdict0 <- PDict(dict0)               # Store the original dictionary in
                                       # a PDict object (preprocessing).

  ## Using the pattern dictionary on chromosome 3R
  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chr3R <- Dmelanogaster$chr3R         # Load chromosome 3R
  chr3R
  mindex <- matchPDict(pdict0, chr3R)  # Search...

  ## Looking at the matches
  start_index <- startIndex(mindex)    # Get the start index.
  length(start_index)                  # Same as the original dictionary.
  start_index[[8220]]                  # Starts of the 8220th pattern.
  end_index <- endIndex(mindex)        # Get the end index.
  end_index[[8220]]                    # Ends of the 8220th pattern.
  count_index <- countIndex(mindex)    # Get the number of matches per pattern.
  count_index[[8220]]
  mindex[[8220]]                       # Get the matches for the 8220th pattern.
  start(mindex[[8220]])                # Equivalent to startIndex(mindex)[[8220]].
  sum(count_index)                     # Total number of matches.
  table(count_index)
  i0 <- which(count_index == max(count_index))
  pdict0[[i0]]                         # The pattern with most occurrences.
  mindex[[i0]]                         # Its matches as a Views object.
  views(chr3R, start_index[[i0]], end_index[[i0]]) # And as an XStringViews object.

  ## Get the coverage of the original subject
  cov3R <- coverage(mindex, 1, length(chr3R))
  max(cov3R)
  mean(cov3R)
  sum(cov3R != 0) / length(cov3R)      # Only 2.44% of chr3R is covered.
  if (interactive()) {
    plotCoverage <- function(coverage, start, end)
    {
      plot.new()
      plot.window(c(start, end), c(0, 20))
      axis(1)
      axis(2)
      axis(4)
      lines(start:end, coverage[start:end], type="l")
    }
    plotCoverage(cov3R, 27600000, 27900000)
  }

  ## ---------------------------------------------------------------------
  ## B. WITH NAMED PATTERNS
  ## ---------------------------------------------------------------------

  dict1 <- dict0[8211:8236]
  names(dict1) <- LETTERS
  dict1                                # The original dictionary.
  pdict1 <- PDict(dict1)
  
  length(pdict1)                       # Same as the original dictionary.
  names(pdict1)                        # Same as names(dict1).
  mindex <- matchPDict(pdict1, chr3R)  # Search...

  ## Looking at the matches
  names(mindex)                        # Same as names(dict1).
  start_index <- startIndex(mindex)
  start_index
  length(start_index)                  # NOT the same as the original
                                       # dictionary.
  start_index <- startIndex(mindex, all.names=TRUE)
  length(start_index)                  # Same as the original dictionary.
  countIndex(mindex)
  unlist(mindex)
  length(unlist(mindex))               # Total number of matches.
  all_matches <- extractAllMatches(chr3R, mindex)
  all_matches
  names(all_matches)

  ## ---------------------------------------------------------------------
  ## C. PERFORMANCE
  ## ---------------------------------------------------------------------

  ## If getting the number of matches is what matters only (without
  ## regarding their positions), then countPDict() will be faster,
  ## especially when there is a high number of matches
  count_index1 <- countPDict(pdict1, chr3R)
  identical(count_index1, count_index) # TRUE

  if (interactive()) {
    ## What's the impact of the dictionary width on performance?
    ## Below is some code that can be used to figure out (will take a long
    ## time to run). For different widths of the original dictionary, we
    ## look at:
    ##   o pptime: preprocessing time (in sec.) i.e. time needed for
    ##             building the PDict object from the truncated input
    ##             sequences;
    ##   o nnodes: nb of nodes in the resulting Aho-Corasick tree;
    ##   o nupatt: nb of unique truncated input sequences;
    ##   o matchtime: time (in sec.) needed to find all the matches;
    ##   o totalcount: total number of matches.
    getPDictStats <- function(dict0, subject)
    {
      list(
        width=width(dict0[1]),
        pptime=system.time(pdict <- PDict(dict0))[["elapsed"]],
        nnodes=length(pdict@actree),
        nupatt=sum(!duplicated(pdict)),
        matchtime=system.time(mindex <- matchPDict(pdict, subject))[["elapsed"]],
        totalcount=sum(countIndex(mindex))
      )
    }
    stats <- lapply(6:25,
                 function(width)
                     getPDictStats(DNAStringSet(dict0, end=width), chr3R))
    stats <- data.frame(do.call(rbind, stats))
    stats
  }
}

\keyword{methods}
