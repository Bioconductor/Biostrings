\name{lcsuffix}

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

\title{Find Longest Common Prefix/Suffix}

\description{
  Functions for searching the Longest Common Prefix/Suffix of two strings.

  WARNING: These functions are experimental and might
  not work properly! Full documentation will come later.

  Thanks for your comprehension!
}

\usage{
lcprefix(s1, s2)
lcsuffix(s1, s2)
}

\arguments{
  \item{s1}{
    1st string, a character string or an \link{XString} object.
  }
  \item{s2}{
    2nd string, a character string or an \link{XString} object.
  }
}

\seealso{
  \link{matchPattern},
  \link{XStringViews-class},
  \link{XString-class}
}

\keyword{methods}
