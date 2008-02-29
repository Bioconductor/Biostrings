\name{matchPDict}

\alias{matchPDict-exact}

% Classes
\alias{class:ViewsIndex}
\alias{ViewsIndex-class}
\alias{class:ByPos_ViewsIndex}
\alias{ByPos_ViewsIndex-class}
\alias{class:ByName_ViewsIndex}
\alias{ByName_ViewsIndex-class}

% Generics and methods:
\alias{length,ACtree-method}
\alias{show,ACtree-method}
\alias{[,ACtree-method}
\alias{as.matrix,ACtree-method}
\alias{initialize,ACtree-method}

\alias{startIndex}
\alias{startIndex,ViewsIndex-method}
\alias{endIndex}
\alias{countIndex}
\alias{countIndex,ViewsIndex-method}
\alias{[[,ViewsIndex-method}
\alias{$,ViewsIndex-method}
\alias{unlist,ViewsIndex-method}

\alias{length,ByPos_ViewsIndex-method}
\alias{names,ByPos_ViewsIndex-method}
\alias{show,ByPos_ViewsIndex-method}
\alias{[[,ByPos_ViewsIndex-method}
\alias{startIndex,ByPos_ViewsIndex-method}
\alias{endIndex,ByPos_ViewsIndex-method}
\alias{countIndex,ByPos_ViewsIndex-method}

\alias{length,ByName_ViewsIndex-method}
\alias{names,ByName_ViewsIndex-method}
\alias{show,ByName_ViewsIndex-method}
\alias{[[,ByName_ViewsIndex-method}
\alias{startIndex,ByName_ViewsIndex-method}
\alias{endIndex,ByName_ViewsIndex-method}
\alias{countIndex,ByName_ViewsIndex-method}

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
  algorithm and is very fast.

  Note that, for now, \code{matchPDict} only works with a dictionary
  of DNA patterns.

  This man page shows how to use \code{matchPDict} for exact matching
  of a constant width dictionary i.e. a dictionary where all the
  patterns have the same length (same number of nucleotides).

  See \link{matchPDict-inexact} for inexact matching (or for using a
  variable width dictionary).
}

\usage{
  matchPDict(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
  countPDict(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)

  ## Manipulation of the ViewsIndex object returned by matchPDict()
  startIndex(x, all.names=FALSE)
  endIndex(x, all.names=FALSE)
  countIndex(x, all.names=FALSE)
  unlist(x, recursive=TRUE, use.names=TRUE)
  extractAllMatches(subject, vindex)
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
    The maximum number of mismatching letters allowed (see
    \code{\link{isMatching}} for the details).
    When used on a constant width dictionary, \code{matchPDict} (and
    \code{countPDict}) only support exact matching so \code{max.mismatch}
    must be zero. See \link{matchPDict-inexact} for inexact matching.
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{\link{isMatching}} for the details).
    When used on a constant width dictionary, \code{matchPDict} (and
    \code{countPDict}) only support exact matching so \code{fixed}
    must be \code{TRUE}. See \link{matchPDict-inexact} for inexact matching.
  }
  \item{x}{
    A PDict or ViewsIndex object for \code{names}.
    A ViewsIndex object for the other methods.
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
  \item{vindex}{
    A ViewsIndex object returned by a previous call to \code{matchPDict}.
  }
}

\value{
  \code{matchPDict} returns a ViewsIndex object.

  \code{countPDict} returns an integer vector.

  \code{startIndex}, \code{endIndex} and \code{countIndex} return vectors of the same
  length as the input dictionary: \code{startIndex} and \code{endIndex} return a
  list of integer vectors, and \code{countIndex} returns an integer vector.

  \code{extractAllMatches} returns a \link{BStringViews} object.
}

\author{H. Pages}

\seealso{
  \link{PDict},
  \link{matchPDict-inexact},
  \code{\link{isMatching}},
  \code{\link{matchPattern}},
  \code{\link{alphabetFrequency}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  ## -------------------------------------------------------------------------
  ## A. WITH UNNAMED PATTERNS
  ## -------------------------------------------------------------------------

  ## Creating the pattern dictionary
  library(drosophila2probe)
  dict0 <- drosophila2probe$sequence   # The input dictionary.
  length(dict0)                        # Hundreds of thousands of patterns.
  pdict0 <- PDict(dict0)               # Store the input dictionary into a
                                       # PDict object (preprocessing).

  ## Using the pattern dictionary on chromosome 3R
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chr3R <- Dmelanogaster[["3R"]]       # Load chromosome 3R
  chr3R
  vindex <- matchPDict(pdict0, chr3R)  # Search...

  ## Looking at the matches
  start_index <- startIndex(vindex)    # Get the start index.
  length(start_index)                  # Same as the input dictionary.
  start_index[[8220]]                  # Starts of the 8220th pattern.
  end_index <- endIndex(vindex)        # Get the end index.
  end_index[[8220]]                    # Ends of the 8220th pattern.
  count_index <- countIndex(vindex)    # Get the number of matches per pattern.
  count_index[[8220]]
  vindex[[8220]]                       # Get the matches for the 8220th pattern.
  start(vindex[[8220]])                # Equivalent to startIndex(vindex)[[8220]].
  sum(count_index)                     # Total number of matches.
  table(count_index)
  i0 <- which(count_index == max(count_index))
  dict0[i0]                            # The pattern with most occurences.
  vindex[[i0]]                         # Its matches as a Views object.
  views(chr3R, start_index[[i0]], end_index[[i0]]) # And as a BStringViews object.

  ## -------------------------------------------------------------------------
  ## B. WITH NAMED PATTERNS
  ## -------------------------------------------------------------------------

  dict1 <- dict0[8211:8236]            # The input dictionary.
  names(dict1) <- LETTERS
  dict1[1:5]
  pdict1 <- PDict(dict1)
  length(pdict1)                       # Same as the input dictionary.
  names(pdict1)                        # Same as names(dict1).
  vindex <- matchPDict(pdict1, chr3R)  # Search...

  ## Looking at the matches
  names(vindex)                        # Same as names(dict1).
  start_index <- startIndex(vindex)
  start_index
  length(start_index)                  # NOT the same as the input dictionary.
  start_index <- startIndex(vindex, all.names=TRUE)
  length(start_index)                  # Same as the input dictionary.
  countIndex(vindex)
  unlist(vindex)
  length(unlist(vindex))               # Total number of matches.
  extractAllMatches(chr3R, vindex)

  ## -------------------------------------------------------------------------
  ## C. PERFORMANCE
  ## -------------------------------------------------------------------------

  ## If getting the number of matches is what matters only (without regarding
  ## their positions), then countPDict() will be faster, especially when there
  ## is a high number of matches
  count_index1 <- countPDict(pdict1, chr3R)
  identical(count_index1, count_index) # TRUE

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
    ##   o totalcount: total number of matches
    getPDictStats <- function(dict0, width, subject)
    {
      list(
        width=width,
        pptime=system.time(pdict <- PDict(dict0, tb.end=width, drop.tail=TRUE))[["elapsed"]],
        nnodes=length(pdict@actree),
        nupatt=sum(pdict@dups == 0),
        matchtime=system.time(vindex <- matchPDict(pdict, subject))[["elapsed"]],
        totalcount=sum(countIndex(vindex))
      )
    }
    stats <- lapply(6:25, function(width) getPDictStats(dict0, width, chr3R))
    stats <- data.frame(do.call(rbind, stats))
    stats
  }
}

\keyword{methods}

