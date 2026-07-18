\name{MIndexList-class}
\docType{class}

% Classes
\alias{MIndexList}
\alias{class:MIndexList}
\alias{MIndexList-class}

% Generics and methods:
\alias{coerce,list,MIndexList-method}
\alias{coerce,List,MIndexList-method}
\alias{show,MIndexList-method}
\alias{showAsCell,MIndexList-method}

\title{MIndexList objects}

\description{
  The MIndexList class is a basic container for storing groups of
  \link{MIndex} objects.
}

\details{
  An MIndexList object contains a set of \link{MIndex} objects.
  This holds the result of matching multiple input sequences to multiple target
  sequences.

  The \code{\link{vmatchPDict}} function returns an MIndexList object.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is an MIndexList object.

  \describe{
    \item{\code{length(x)}:}{
      The number of \link{MIndex} objects in \code{x} (corresponding to the
      number of subject sequences matched).
    }
    \item{\code{names(x)}:}{
      The names of the \link{MIndex} objects in \code{x} (corresponding to the
      names of the subject sequences matched).
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below, \code{x} is an MIndexList object.

  \describe{
    \item{\code{x[[i]]}:}{
      Extract the matches for the i-th subject sequence as an
      \link{MIndex} object.
    }
    \item{\code{x[i]}:}{
      Extract a subset of \code{x} as another MIndexList object.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a list or List object.

  \describe{
    \item{\code{as(x, "MIndexList")}:}{
      Coerce a list or List of \link{MIndex} objects to an MIndexList object.
    }
  }
}

\author{A. Lakshman}

\seealso{
  \code{\link{vmatchPDict}},
  \link{MIndex-class},
  \link{PDict-class}
}

\examples{
  ## 1. Create a set of patterns (dictionary) and subject sequences
  patterns <- DNAStringSet(c("AAC", "TTA", "CGT"))
  pdict <- PDict(patterns)
  subjects <- DNAStringSet(c("AACTA", "TTAC", "CGTTAA", "NNN"))
  names(subjects) <- paste0("seq", 1:4)

  ## 2. Match patterns to subjects
  mindex_list <- vmatchPDict(pdict, subjects)
  mindex_list

  ## 3. Basic properties
  length(mindex_list)
  names(mindex_list)

  ## 4. Accessing and subsetting
  mindex_list[[1]]  # extract matches for seq1 as an MIndex object
  mindex_list[2:3]  # extract subset of seqs as another MIndexList object

  ## 5. Coercion
  list_of_mindexes <- as.list(mindex_list)
  as(list_of_mindexes, "MIndexList")
}

\keyword{methods}
\keyword{classes}
