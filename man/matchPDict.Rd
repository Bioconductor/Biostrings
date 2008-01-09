\name{matchPDict}

% Classes
\alias{class:P2Views}
\alias{P2Views-class}
\alias{class:P2ViewsWithoutIDs}
\alias{P2ViewsWithoutIDs-class}
\alias{class:P2ViewsWithIDs}
\alias{P2ViewsWithIDs-class}
\alias{class:PDict}
\alias{PDict-class}
\alias{class:ULdna_PDict}
\alias{ULdna_PDict-class}

% Generics and methods:
\alias{pids}
\alias{p2starts}
\alias{p2starts,P2Views-method}
\alias{p2ends}
\alias{p2nmatch}
\alias{p2nmatch,P2Views-method}
\alias{[[,P2Views-method}
\alias{$,P2Views-method}

\alias{length,P2ViewsWithoutIDs-method}
\alias{pids,P2ViewsWithoutIDs-method}
\alias{show,P2ViewsWithoutIDs-method}
\alias{[[,P2ViewsWithoutIDs-method}
\alias{p2starts,P2ViewsWithoutIDs-method}
\alias{p2ends,P2ViewsWithoutIDs-method}
\alias{p2nmatch,P2ViewsWithoutIDs-method}

\alias{length,P2ViewsWithIDs-method}
\alias{pids,P2ViewsWithIDs-method}
\alias{show,P2ViewsWithIDs-method}
\alias{[[,P2ViewsWithIDs-method}
\alias{p2starts,P2ViewsWithIDs-method}
\alias{p2ends,P2ViewsWithIDs-method}
\alias{p2nmatch,P2ViewsWithIDs-method}

\alias{length,ACtree-method}
\alias{show,ACtree-method}
\alias{[,ACtree-method}
\alias{as.matrix,ACtree-method}
\alias{initialize,ACtree-method}

\alias{initialize,ULdna_PDict-method}
\alias{length,ULdna_PDict-method}
\alias{width,ULdna_PDict-method}
\alias{pids,ULdna_PDict-method}
\alias{show,ULdna_PDict-method}

\alias{matchPDict}
\alias{matchPDict,BString-method}

% Functions:
\alias{extractAllMatches}


\title{Searching a sequence for patterns stored in a dictionary}

\description{
  \code{matchPDict} efficiently finds all occurences in a text (the
  subject) of any pattern from a set of patterns (the dictionary).
  The implementation of \code{matchPDict} is based on the Aho-Corasick
  algorithm and is very fast. However, for now, \code{matchPDict} only
  works with a dictionary of DNA patterns where all the patterns have
  the same number of nucleotides (uniform-length dictionary).
  Also, only exact matching is currently supported.
}

\usage{
  matchPDict(pdict, subject, algorithm="auto", mismatch=0, fixed=TRUE)

  ## Manipulation of the P2Views object returned by matchPDict()
  p2starts(x, all.pids=FALSE)
  p2ends(x, all.pids=FALSE)
  p2nmatch(x, all.pids=FALSE)
  extractAllMatches(subject, pdictmatches)
}

\arguments{
  \item{pdict}{
    The preprocessed dictionary (see below for the details).
  }
  \item{subject}{
    A \link{BString} (or derived) object containing the subject string.
    For now, only \link{DNAString} subjects are supported.
  }
  \item{algorithm}{
    Not supported yet.
  }
  \item{mismatch}{
    Not supported yet.
  }
  \item{fixed}{
    Not supported yet.
  }
  \item{x}{
    A P2Views object returned by a previous call to \code{matchPDict}.
  }
  \item{all.pids}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{pdictmatches}{
    A P2Views object returned by a previous call to \code{matchPDict}.
  }
}

\value{
  \code{matchPDict} returns a P2Views object.

  \code{p2starts}, \code{p2ends} and \code{p2nmatch} return vectors of the same
  length as the input dictionary: \code{p2starts} and \code{p2ends} return a
  list of integer vectors, and \code{p2nmatch} returns an integer vector.

  \code{extractAllMatches} returns a \link{BStringViews} object.
}

\author{H. Pages}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{alphabetFrequency}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## Creating the pattern dictionary
  library(drosophila2probe)
  dict0 <- drosophila2probe$sequence  # The input dictionary.
  length(dict0)                       # Hundreds of thousands of patterns.
  unique(nchar(dict0))                # Patterns are 25-mers.
  dict0[1:5]
  pdict <- new("ULdna_PDict", dict0)  # Store the input dictionary into a
                                      # ULdna_PDict instance (preprocessing).
  length(pdict)                       # Same as the input dictionary.
  width(pdict)                        # The number of chars per pattern.

  ## Using the pattern dictionary on chromosome 3R
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chr3R <- Dmelanogaster[["3R"]]      # Load chromosome 3R
  chr3R
  matches <- matchPDict(pdict, chr3R) # Search...

  ## Looking at the matches
  starts <- p2starts(matches)         # Get the starts of the matches
  length(starts)                      # Same as the input dictionary.
  starts[[8220]]                      # Starts of the 8220th pattern.
  ends <- p2ends(matches)             # Get the ends of the matches.
  ends[[8220]]                        # Ends of the 8220th pattern.
  nmatch <- p2nmatch(matches)         # Get the number of matches per pattern.
  nmatch[[8220]]
  matches[[8220]]                     # Get the matches for the 8220th pattern.
  sum(nmatch)                         # Total number of matches.
  table(nmatch)
  i0 <- which(nmatch == max(nmatch))
  dict0[i0]                           # The pattern with most occurences.
  matches[[i0]]                       # Its matches as a Views object.
  views(chr3R, starts[[i0]], ends[[i0]]) # And as a BStringViews object.

  ## NOT RUN
  if (FALSE) {
    ## What's the impact of the dictionary width on performance?
    ## Below is some code that can be used to figure out (will take a long
    ## time to run). For different widths of the input dictionary, we look at:
    ##   o pptime: preprocessing time (in sec.) i.e. time needed for building
    ##             the PDict object from the truncated input sequences
    ##   o nnodes: nb of nodes in the resulting Aho-Corasick tree
    ##   o nupatt: nb of unique truncated input sequences
    ##   o matchtime: time (in sec.) needed to find all the matches
    ##   o sum_nmatch: total number of matches
    getPDictStats <- function(dict0, width, subject)
    {
      dict <- substr(dict0, 1, width)
      list(
        width=width,
        pptime=system.time(pdict <- new("ULdna_PDict", dict))[["elapsed"]],
        nnodes=length(pdict@actree),
        nupatt=sum(pdict@dups == 0),
        matchtime=system.time(matches <- matchPDict(pdict, subject))[["elapsed"]],
        sum_nmatch=sum(p2nmatch(matches))
      )
    }
    stats <- lapply(6:25, function(width) getPDictStats(dict0, width, chr3R))
    stats <- data.frame(do.call(rbind, stats))
    stats
  }
}

\keyword{methods}

