\name{XStringViews-class}
\docType{class}

% Classes:
\alias{class:XStringViews}
\alias{XStringViews-class}

% Methods:
\alias{Views,XString-method}
\alias{Views,character-method}

\alias{nchar,XStringViews-method}

\alias{xsbasetype,XStringViews-method}
\alias{xsbasetype<-,XStringViews-method}

\alias{coerce,XStringViews,XStringSet-method}
\alias{coerce,XStringViews,BStringSet-method}
\alias{coerce,XStringViews,DNAStringSet-method}
\alias{coerce,XStringViews,RNAStringSet-method}
\alias{coerce,XStringViews,AAStringSet-method}

\alias{show,XStringViews-method}

\alias{==,XStringViews,XStringViews-method}
\alias{==,XStringViews,XString-method}
\alias{==,XStringViews,character-method}
\alias{==,XString,XStringViews-method}
\alias{==,character,XStringViews-method}

\alias{!=,XStringViews,XStringViews-method}
\alias{!=,XStringViews,XString-method}
\alias{!=,XStringViews,character-method}
\alias{!=,XString,XStringViews-method}
\alias{!=,character,XStringViews-method}

\alias{as.character,XStringViews-method}
\alias{as.matrix,XStringViews-method}
\alias{toString,XStringViews-method}


\title{The XStringViews class}

\description{
  The XStringViews class is the basic container for storing a set of views
  (start/end locations) on the same sequence (an \link{XString} object).
}

\usage{
  ## Constructors:

  \S4method{Views}{character}(subject, start=NA, end=NA, names=NULL)
  \S4method{Views}{XString}(subject, start=NA, end=NA, names=NULL)
}

\arguments{
  \item{subject}{
    The subject sequence.
  }
  \item{start, end}{
    Integer vectors specifying the starting and ending positions of each view.
  }
  \item{names}{
    If not \code{NULL}, the names to assign to each view.
  }
}

\details{
  An XStringViews object contains a set of views (start/end locations) on the
  same \link{XString} object called "the subject string"
  or "the subject sequence" or simply "the subject".
  Each view is defined by its start and end locations: both are
  integers such that start <= end.
  An XStringViews object is in fact a particular case of an \link[IRanges]{Views}
  object (the XStringViews class contains the \link[IRanges]{Views} class) so it
  can be manipulated in a similar manner: see \code{?\link[IRanges]{Views}} for
  more information.
  Note that two views can overlap and that a view can be "out of limits"
  i.e. it can start before the first letter of the subject or/and end
  after its last letter.
}

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is an XStringViews object.

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
  \code{x}, \code{object}, \code{e1} and \code{e2} are XStringViews objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{e1 == e2}:
      A vector of logicals indicating the result of the view by
      view comparison. The views in the shorter of the two XStringViews
      object being compared are recycled as necessary.

      Like for comparison between \link{XString} objects, comparison between two
      XStringViews objects with subjects of different classes is not supported
      with one exception: when the subjects are \link{DNAString} and
      \link{RNAString} instances.

      Also, like with \link{XString} objects, comparison between an XStringViews
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
      \code{use.names} controls whether or not \code{names(x)} should be
      used to set the names of the returned vector (default is \code{TRUE}).
      \code{check.limits} controls whether or not an error should be raised
      if \code{x} contains "out of limit" views (default is \code{TRUE}).
      With \code{check.limits=FALSE} then "out of limit" views are padded
      with spaces.
    }
    \item{}{
      \code{as.matrix(x, mode, use.names, check.limits)}:
      Depending on what \code{mode} is chosen (\code{"integer"} or
      \code{"character"}), return either a 2-column integer matrix containing
      \code{start(x)} and \code{end(x)} or a character matrix containing the
      "exploded" representation of the views.
      \code{mode="character"} can only be used on an XStringViews
      object with equal-width views.
      Arguments \code{use.names} and \code{check.limits} are ignored with
      \code{mode="integer"}.
      With \code{mode="character"}, \code{use.names} controls
      whether or not \code{names(x)} should be used to set the row names
      of the returned matrix (default is \code{TRUE}),
      and \code{check.limits} controls whether or not an error should be raised
      if \code{x} contains "out of limit" views (default is \code{TRUE}).
      With \code{check.limits=FALSE} then "out of limit" views are padded
      with spaces.
    }
    \item{}{
      \code{toString(x)}:
      Equivalent to \code{toString(as.character(x))}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link[IRanges]{Views-class},
  \code{\link[IRanges]{gaps}},
  \link{XStringViews-constructors},
  \link{XString-class},
  \link{XStringSet-class},
  \code{\link{letter}},
  \link{MIndex-class}
}

\examples{
  ## One standard way to create an XStringViews object is to use
  ## the Views() constructor.

  ## Views on a DNAString object:
  s <- DNAString("-CTC-N")
  v4 <- Views(s, start=3:0, end=5:8)
  v4
  subject(v4)
  length(v4)
  start(v4)
  end(v4)
  width(v4)

  ## Attach a comment to views #3 and #4:
  names(v4)[3:4] <- "out of limits"
  names(v4)

  ## A more programatical way to "tag" the "out of limits" views:
  names(v4)[start(v4) < 1 | nchar(subject(v4)) < end(v4)] <- "out of limits"
  ## or just:
  names(v4)[nchar(v4) < width(v4)] <- "out of limits"

  ## Two equivalent ways to extract a view as an XString object:
  s2a <- v4[[2]]
  s2b <- subseq(subject(v4), start=start(v4)[2], end=end(v4)[2])
  identical(s2a, s2b) # TRUE

  ## It is an error to try to extract an "out of limits" view:
  #v4[[3]] # Error!

  v12 <- Views(DNAString("TAATAATG"), start=-2:9, end=0:11)
  v12 == DNAString("TAA")
  v12[v12 == v12[4]]
  v12[v12 == v12[1]]
  v12[3] == Views(RNAString("AU"), start=0, end=2)

  ## Here the first view doesn't even overlap with the subject:
  Views(BString("aaa--b"), start=-3:4, end=-3:4 + c(3:6, 6:3))

  ## 'start' and 'end' are recycled:
  subject <- "abcdefghij"
  Views(subject, start=2:1, end=4)
  Views(subject, start=5:7, end=length(subject))
  Views(subject, start=1, end=5:7)

  ## Applying gaps() to an XStringViews object:
  v2 <- Views("abCDefgHIJK", start=c(8, 3), end=c(14, 4))
  gaps(v2)

  ## Coercion:
  as(v12, "XStringSet")  # same as 'as(v12, "DNAStringSet")'
  as(v12, "RNAStringSet")
}

\keyword{methods}
\keyword{classes}
