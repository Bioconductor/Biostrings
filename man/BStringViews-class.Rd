\name{BStringViews-class}
\docType{class}

% Classes:
\alias{class:BStringViews}
\alias{BStringViews-class}

% Accesor methods:
\alias{subject}
\alias{subject,BStringViews-method}
\alias{length,BStringViews-method}
\alias{start}
\alias{start,BStringViews-method}
\alias{first}
\alias{first,BStringViews-method}
\alias{end}
\alias{end,BStringViews-method}
\alias{last}
\alias{last,BStringViews-method}
\alias{width}
\alias{width,BStringViews-method}
\alias{nchar,BStringViews-method}
\alias{desc}
\alias{desc,BStringViews-method}
\alias{desc<-}
\alias{desc<-,BStringViews-method}

% Standard generic methods:

\alias{show,BStringViews-method}
\alias{[,BStringViews-method}
\alias{[<-,BStringViews-method}
\alias{[[,BStringViews-method}
\alias{[[<-,BStringViews-method}

\alias{==,BStringViews,BStringViews-method}
\alias{==,BStringViews,BString-method}
\alias{==,BStringViews,character-method}
\alias{==,BString,BStringViews-method}
\alias{==,character,BStringViews-method}

\alias{!=,BStringViews,BStringViews-method}
\alias{!=,BStringViews,BString-method}
\alias{!=,BStringViews,character-method}
\alias{!=,BString,BStringViews-method}
\alias{!=,character,BStringViews-method}

\alias{as.character,BStringViews-method}
\alias{as.matrix,BStringViews-method}
\alias{as.list,BStringViews-method}
\alias{toString,BStringViews-method}


\title{The BStringViews class}

\description{
  A container for a set of views on the same \link{BString} (or
  derived) object.
}

\details{
  A BStringViews object contains a set of views on the same
  \link{BString} (or derived) object called "the subject string"
  or "the subject sequence" or simply "the subject".
  Each view is defined by its start and end locations: both are
  integers such that start <= end.
  Note that two views can overlap and that a view can be "out of limits"
  i.e. it starts before the first letter of the subject or/and ends
  after its last letter. 
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringViews object.

  \describe{
    \item{}{
      \code{subject(x)}:
      The subject of \code{x}. This is always a \link{BString}
      (or derived) object.
    }
    \item{}{
      \code{length(x)}:
      The number of views in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The starting positions of the views.
      An integer vector of the same length as \code{x}.
    }
    \item{}{
      \code{first(x)}:
      deprecated. Use 'start' instead.
    }
    \item{}{
      \code{end(x)}:
      The ending positions of the views.
      An integer vector of the same length as \code{x}.
    }
    \item{}{
      \code{last(x)}:
      deprecated. Use 'end' instead.
    }
    \item{}{
      \code{width(x)}:
      The width of the views.
      Equivalent to \code{end(x) - start(x) + 1L}.
      A vector of positive integers of the same length as \code{x}.
    }
    \item{}{
      \code{nchar(x)}:
      A vector of non-negative integers containing the number
      of letters in each view.
      Values in \code{nchar(x)} coincide with values in \code{width(x)}
      except for "out of limits" views where they are lower.
    }
    \item{}{
      \code{desc(x)}:
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user provided description or comment for each view.
      These are the only data in a BStringViews object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2} are BStringViews objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new BStringViews object made of the selected views.
      \code{i} can be a numeric vector, a logical vector, \code{NULL} or missing.
      The returned object has the same subject as \code{x}.
    }
    \item{}{
      \code{x[[i]]}:
      Extract a view as a \link{BString} (or derived) object.
      \code{i} must be a single numeric value (a numeric vector of length 1).
      Can't be used for extracting a view that is "out of limits" (raise an error).
      The returned object belongs to the same class as \code{subject(x)}.
    }
    \item{}{
      \code{e1 == e2}:
      A vector of logicals indicating the result of the view by
      view comparison. The views in the shorter of the two BStringViews
      object being compared are recycled as necessary.

      Like for comparison between \link{BString} objects, comparison between two
      BStringViews objects with subjects of different classes is not supported
      with one exception: when the subjects are \link{DNAString} and \link{RNAString}
      instances.

      Also, like with \link{BString} objects, comparison between a BStringViews
      object with a subject of class "BString" and a character vector
      is supported (see examples below).
    }
    \item{}{
      \code{e1 != e2}:
      Equivalent to \code{!(e1 == e2)}.
    }
    \item{}{
      \code{as.character(x, use.names)}:
      Converts \code{x} to a character vector of the same length as \code{x}.
      \code{use.names} controls whether or not \code{desc(x)} should be
      used to set the names of the returned vector.
    }
    \item{}{
      \code{as.matrix(x)}:
      Return a 2-col matrix containing \code{start(x)} and \code{end(x)}.
    }
    \item{}{
      \code{as.list(x)}:
      Convert \code{x} to a list of \link{BString} objects.
      Can't be used if \code{x} has "out of limits" views (raise an error).
    }
    \item{}{
      \code{toString(x)}:
      Equivalent to \code{toString(as.character(x))}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BStringViews-constructors},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \code{\link{letter}}
}

\examples{
  s <- DNAString("-CTC-N")

  ## One standard way to create a BStringViews object is to use
  ## the "views" constructor:
  v4 <- views(s, 3:0, 5:8)
  v4
  subject(v4)
  length(v4)
  start(v4)
  end(v4)
  width(v4)

  ## Attach a comment to views #3 and #4:
  desc(v4)[3:4] <- "out of limits"
  desc(v4)

  ## A more programatical way to "tag" the "out of limits" views:
  desc(v4)[start(v4) < 1 | nchar(subject(v4)) < end(v4)] <- "out of limits"
  ## or just:
  desc(v4)[nchar(v4) < width(v4)] <- "out of limits"

  ## Make a new BStringViews object by selecting some views only:
  v4[4:2]               # 3 views
  v4[-1]                # 3 views
  v4[c(2,4)]            # 2 views
  v4[width(v4) >= 5]    # 3 views
  v4[is.na(desc(v4))]   # 2 views
  v4[!is.na(desc(v4))]  # 2 views
  v4[3]                 # 1 view
  v4[FALSE]             # 0 view

  ## Two equivalent ways to extract a view as a BString (or derived) object:
  s2a <- v4[[2]]
  s2b <- subBString(subject(v4), start(v4)[2], end(v4)[2])
  identical(s2a, s2b) # TRUE

  ## It is an error to try to extract an "out of limits" view:
  #v4[[3]] # Error!

  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
  v12 == DNAString("TAA")
  v12[v12 == v12[4]]
  v12[v12 == v12[1]]
  v12[3] == views(RNAString("AU"), 0, 2)
}

\keyword{methods}
\keyword{classes}
