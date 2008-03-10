\name{BStringViews-class}
\docType{class}

% Classes:
\alias{class:BStringViews}
\alias{BStringViews-class}

% Constructor-like functions and generics:
\alias{initialize,BStringViews-method}

% Accesor methods:
\alias{subject}
\alias{subject,BStringViews-method}
\alias{width<-,BStringViews-method}
\alias{nchar,BStringViews-method}

% Other methods:

\alias{show,BStringViews-method}
\alias{[[,BStringViews-method}
\alias{[[<-,BStringViews-method}

\alias{==,BStringViews,BStringViews-method}
\alias{==,BStringViews,XString-method}
\alias{==,BStringViews,character-method}
\alias{==,XString,BStringViews-method}
\alias{==,character,BStringViews-method}

\alias{!=,BStringViews,BStringViews-method}
\alias{!=,BStringViews,XString-method}
\alias{!=,BStringViews,character-method}
\alias{!=,XString,BStringViews-method}
\alias{!=,character,BStringViews-method}

\alias{as.character,BStringViews-method}
\alias{as.matrix,BStringViews-method}
\alias{as.list,BStringViews-method}
\alias{toString,BStringViews-method}


\title{The BStringViews class}

\description{
  The BStringViews class is the basic container for storing a set of views
  (start/end locations) on the same sequence (an \link{XString} object).
}

\details{
  A BStringViews object contains a set of views (start/end locations) on the
  same \link{XString} object called "the subject string"
  or "the subject sequence" or simply "the subject".
  Each view is defined by its start and end locations: both are
  integers such that start <= end.
  A BStringViews object is in fact a particular case of a \link{IRanges} object
  (the BStringViews class contains the \link{IRanges} class) so it can be
  manipulated in the same way: see \link{IRanges} for more information.
  Note that two views can overlap and that a view can be "out of limits"
  i.e. it can start before the first letter of the subject or/and end
  after its last letter.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringViews object.

  \describe{
    \item{}{
      \code{subject(x)}:
      The subject of \code{x}. This is always an \link{XString}
      object.
    }
    \item{}{
      \code{nchar(x)}:
      A vector of non-negative integers containing the number
      of letters in each view.
      Values in \code{nchar(x)} coincide with values in \code{width(x)}
      except for "out of limits" views where they are lower.
    }
  }
}

\section{Other methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2} are BStringViews objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract a view as an \link{XString} object.
      \code{i} must be a single numeric value (a numeric vector of length 1).
      Can't be used for extracting a view that is "out of limits" (raise an
      error). The returned object has the same \link{XString} subtype as
      \code{subject(x)}.
    }
    \item{}{
      \code{e1 == e2}:
      A vector of logicals indicating the result of the view by
      view comparison. The views in the shorter of the two BStringViews
      object being compared are recycled as necessary.

      Like for comparison between \link{XString} objects, comparison between two
      BStringViews objects with subjects of different classes is not supported
      with one exception: when the subjects are \link{DNAString} and
      \link{RNAString} instances.

      Also, like with \link{XString} objects, comparison between a BStringViews
      object with a BString subject and a character vector is supported (see
      examples below).
    }
    \item{}{
      \code{e1 != e2}:
      Equivalent to \code{!(e1 == e2)}.
    }
    \item{}{
      \code{as.character(x, use.names, check.limits)}:
      Convert \code{x} to a character vector of the same length as \code{x}.
      \code{use.names} controls whether or not \code{desc(x)} should be
      used to set the names of the returned vector (default is \code{TRUE}).
      \code{check.limits} controls whether or not an error should be raised
      if \code{x} contains "out of limit" views (default is \code{TRUE}).
      With \code{check.limits=FALSE} then "out of limit" views are padded
      with spaces.
    }
    \item{}{
      \code{as.matrix(x, mode, use.names, check.limits)}:
      Depending on what \code{mode} is choosen (\code{"integer"} or
      \code{"character"}), return either a 2-column integer matrix containing
      \code{start(x)} and \code{end(x)} or a character matrix containing the
      "exploded" representation of the views.
      \code{mode="character"} can only be used on a BStringViews
      object with equal-width views.
      Arguments \code{use.names} and \code{check.limits} are ignored with
      \code{mode="integer"}.
      With \code{mode="character"}, \code{use.names} controls
      whether or not \code{desc(x)} should be used to set the row names
      of the returned matrix (default is \code{TRUE}),
      and \code{check.limits} controls whether or not an error should be raised
      if \code{x} contains "out of limit" views (default is \code{TRUE}).
      With \code{check.limits=FALSE} then "out of limit" views are padded
      with spaces.
    }
    \item{}{
      \code{as.list(x)}:
      Convert \code{x} to a list of \link{XString} objects.
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
  \link{IRanges-class},
  \link{BStringViews-constructors},
  \link{XString-class},
  \link{XStringSet-class},
  \code{\link{letter}}
}

\examples{

  ## One standard way to create a BStringViews object is to use
  ## the "views" constructor:
  s <- DNAString("-CTC-N")
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

  ## Two equivalent ways to extract a view as an XString object:
  s2a <- v4[[2]]
  s2b <- subXString(subject(v4), start(v4)[2], end(v4)[2])
  identical(s2a, s2b) # TRUE

  ## It is an error to try to extract an "out of limits" view:
  #v4[[3]] # Error!

  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
  v12 == DNAString("TAA")
  v12[v12 == v12[4]]
  v12[v12 == v12[1]]
  v12[3] == views(RNAString("AU"), 0, 2)

  ## Here the first view doesn't even overlap with the subject:
  views(BString("aaa--b"), -3:4, -3:4 + c(3:6, 6:3))
}

\keyword{methods}
\keyword{classes}
