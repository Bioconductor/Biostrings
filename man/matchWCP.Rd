\name{matchWCP}

\alias{WCPscoreStartingAt}

\alias{matchWCP}
\alias{matchWCP,character-method}
\alias{matchWCP,XString-method}
\alias{matchWCP,XStringViews-method}
\alias{matchWCP,MaskedXString-method}
\alias{countWCP}
\alias{countWCP,character-method}
\alias{countWCP,XString-method}
\alias{countWCP,XStringViews-method}
\alias{countWCP,MaskedXString-method}


\title{A simple WCP matching function and related utilities}

\description{
  The functions described in this man page are now defunct.
}

\usage{
matchWCP(wcp, subject, min.score="80\%")
countWCP(wcp, subject, min.score="80\%")
WCPscoreStartingAt(wcp, subject, starting.at=1)
}

\arguments{
  \item{wcp}{
    A \link{WCP} object.
  }
  \item{subject}{
    An \link{XString}, \link{XStringViews} or \link{MaskedXString}
    object for \code{matchWCP} and \code{countWCP}.

    A \link{XString} object for \code{WCPscoreStartingAt}.
  }
  \item{min.score}{
    The minimum score for counting a match.
    Can be given as a character string containing a percentage (e.g.
    \code{"85\%"}) of the highest possible score or as a single number.
  }
  \item{starting.at}{
    An integer vector specifying the starting positions of the
    Weighted Clustered Positions relatively to the subject.
  }
}

\value{
  An \link{XStringViews} object for \code{matchWCP}.

  A single integer for \code{countWCP}.

  A numeric vector containing the Weighted Clustered Positions-based scores
  for \code{WCPscoreStartingAt}.
}

\author{P. Aboyoun}

\seealso{
  \code{\link{matchPWM}},
  \code{\link{matchPattern}},
  \link{WCP-class},
  \link{XString-class},
  \link{XStringViews-class}
}

\keyword{methods}
\keyword{manip}
\keyword{utilities}
