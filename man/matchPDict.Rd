\name{matchPDict}

% Classes
\alias{class:PDict}
\alias{PDict-class}
\alias{class:ULdna_PDict}
\alias{ULdna_PDict-class}
\alias{class:P2Views}
\alias{P2Views-class}
\alias{class:P2ViewsWithoutIDs}
\alias{P2ViewsWithoutIDs-class}
\alias{class:P2ViewsWithIDs}
\alias{P2ViewsWithIDs-class}

% Generics and methods:
\alias{pids}

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

\alias{p2start}
\alias{p2start,P2Views-method}
\alias{p2end}
\alias{p2nview}
\alias{p2nview,P2Views-method}
\alias{[[,P2Views-method}
\alias{$,P2Views-method}
\alias{unlist,P2Views-method}

\alias{length,P2ViewsWithoutIDs-method}
\alias{pids,P2ViewsWithoutIDs-method}
\alias{show,P2ViewsWithoutIDs-method}
\alias{[[,P2ViewsWithoutIDs-method}
\alias{p2start,P2ViewsWithoutIDs-method}
\alias{p2end,P2ViewsWithoutIDs-method}
\alias{p2nview,P2ViewsWithoutIDs-method}

\alias{length,P2ViewsWithIDs-method}
\alias{pids,P2ViewsWithIDs-method}
\alias{show,P2ViewsWithIDs-method}
\alias{[[,P2ViewsWithIDs-method}
\alias{p2start,P2ViewsWithIDs-method}
\alias{p2end,P2ViewsWithIDs-method}
\alias{p2nview,P2ViewsWithIDs-method}

\alias{matchPDict}
\alias{matchPDict,BString-method}
\alias{countPDict}
\alias{countPDict,BString-method}

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
  matchPDict(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
  countPDict(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)

  ## Manipulation of the PDict object passed to matchPDict()
  pids(x)

  ## Manipulation of the P2Views object returned by matchPDict()
  p2start(x, all.pids=FALSE)
  p2end(x, all.pids=FALSE)
  p2nview(x, all.pids=FALSE)
  unlist(x, recursive=TRUE, use.names=TRUE)
  extractAllMatches(subject, p2v)
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
  \item{max.mismatch}{
    Not supported yet.
  }
  \item{fixed}{
    Not supported yet.
  }
  \item{x}{
    A PDict or P2Views object for \code{pids}.
    A P2Views object for the other methods.
  }
  \item{all.pids}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{recursive}{
    ignored.
  }
  \item{use.names}{
    ignored.
  }
  \item{p2v}{
    A P2Views object returned by a previous call to \code{matchPDict}.
  }
}

\value{
  \code{matchPDict} returns a P2Views object.

  \code{countPDict} returns an integer vector.

  \code{p2start}, \code{p2end} and \code{p2nview} return vectors of the same
  length as the input dictionary: \code{p2start} and \code{p2end} return a
  list of integer vectors, and \code{p2nview} returns an integer vector.

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
  ## -------------------------------------------------------------------------
  ## A. WITHOUT USING PATTERN IDS
  ## -------------------------------------------------------------------------

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
  p2v <- matchPDict(pdict, chr3R)     # Search...

  ## Looking at the matches
  p2start <- p2start(p2v)             # Get the starts of the matches
  length(p2start)                     # Same as the input dictionary.
  p2start[[8220]]                     # Starts of the 8220th pattern.
  p2end <- p2end(p2v)                 # Get the ends of the matches.
  p2end[[8220]]                       # Ends of the 8220th pattern.
  nmatch <- p2nview(p2v)              # Get the number of matches per pattern.
  nmatch[[8220]]
  p2v[[8220]]                         # Get the matches for the 8220th pattern.
  start(p2v[[8220]])                  # Equivalent to p2start(p2v)[[8220]].
  sum(nmatch)                         # Total number of matches.
  table(nmatch)
  i0 <- which(nmatch == max(nmatch))
  dict0[i0]                           # The pattern with most occurences.
  p2v[[i0]]                           # Its matches as a Views object.
  views(chr3R, p2start[[i0]], p2end[[i0]]) # And as a BStringViews object.

  ## -------------------------------------------------------------------------
  ## B. USING PATTERN IDS
  ## -------------------------------------------------------------------------

  dict1 <- dict0[8211:8236]           # The input dictionary.
  names(dict1) <- LETTERS
  dict1[1:5]
  pdict <- new("ULdna_PDict", dict1)
  length(pdict)                       # Same as the input dictionary.
  pids(pdict)                         # Same as names(dict1).
  p2v <- matchPDict(pdict, chr3R)     # Search...

  ## Looking at the matches
  pids(p2v)                           # Same as names(dict1).
  p2start <- p2start(p2v)
  p2start
  length(p2start)                     # NOT the same as the input dictionary.
  p2start <- p2start(p2v, all.pids=TRUE)
  length(p2start)                     # Same as the input dictionary.
  p2nview(p2v)
  unlist(p2v)
  nview(unlist(p2v))                  # Total number of matches.
  extractAllMatches(chr3R, p2v)

  ## -------------------------------------------------------------------------
  ## C. PERFORMANCE
  ## -------------------------------------------------------------------------

  ## If getting the number of matches is what matters only (without regarding
  ## their position), then countPDict() will be faster, especially when there
  ## is a high number of matches
  nmatch1 <- countPDict(pdict, chr3R)
  identical(nmatch1, nmatch)          # TRUE

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
      list(
        width=width,
        pptime=system.time(pdict <- new("ULdna_PDict", dict0, width))[["elapsed"]],
        nnodes=length(pdict@actree),
        nupatt=sum(pdict@dups == 0),
        matchtime=system.time(p2v <- matchPDict(pdict, subject))[["elapsed"]],
        sum_nmatch=sum(p2nview(p2v))
      )
    }
    stats <- lapply(6:25, function(width) getPDictStats(dict0, width, chr3R))
    stats <- data.frame(do.call(rbind, stats))
    stats
  }
}

\keyword{methods}

