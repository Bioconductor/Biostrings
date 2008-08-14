\name{pmatchPattern}

\alias{lcprefix}
\alias{lcprefix,character,character-method}
\alias{lcprefix,character,XString-method}
\alias{lcprefix,XString,character-method}
\alias{lcprefix,XString,XString-method}

\alias{lcsuffix}
\alias{lcsuffix,character,character-method}
\alias{lcsuffix,character,XString-method}
\alias{lcsuffix,XString,character-method}
\alias{lcsuffix,XString,XString-method}

\alias{lcsubstr}
\alias{lcsubstr,character,character-method}
\alias{lcsubstr,character,XString-method}
\alias{lcsubstr,XString,character-method}
\alias{lcsubstr,XString,XString-method}

\alias{pmatchPattern}
\alias{pmatchPattern,character-method}
\alias{pmatchPattern,XString-method}
\alias{pmatchPattern,XStringViews-method}


\title{Longest Common Prefix/Suffix/Substring searching functions}

\description{
  Functions for searching the Longest Common Prefix/Suffix/Substring
  of two strings.

  WARNING: These functions are experimental and might
  not work properly! Full documentation will come later.

  Please send questions/comments to hpages@fhcrc.org

  Thanks for your comprehension!
}

\usage{
  lcprefix(s1, s2)
  lcsuffix(s1, s2)
  lcsubstr(s1, s2)
  pmatchPattern(pattern, subject, maxlength.out=1L)
}

\arguments{
  \item{s1}{
    1st string, a character string or an \link{XString} object.
  }
  \item{s2}{
    2nd string, a character string or an \link{XString} object.
  }
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    An \link{XString} object containing the subject string.
  }
  \item{maxlength.out}{
    The maximum length of the output i.e. the maximum number of views
    in the returned object.
  }
}

\seealso{
  \link{matchPattern},
  \link{XStringViews-class},
  \link{XString-class}
}

\keyword{methods}
