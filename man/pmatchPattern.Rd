\name{pmatchPattern}

\alias{lcprefix}
\alias{lcprefix,character,character-method}
\alias{lcprefix,character,BString-method}
\alias{lcprefix,BString,character-method}
\alias{lcprefix,BString,BString-method}

\alias{lcsuffix}
\alias{lcsuffix,character,character-method}
\alias{lcsuffix,character,BString-method}
\alias{lcsuffix,BString,character-method}
\alias{lcsuffix,BString,BString-method}

\alias{lcsubstr}
\alias{lcsubstr,character,character-method}
\alias{lcsubstr,character,BString-method}
\alias{lcsubstr,BString,character-method}
\alias{lcsubstr,BString,BString-method}

\alias{pmatchPattern}
\alias{pmatchPattern,character-method}
\alias{pmatchPattern,BString-method}
\alias{pmatchPattern,BStringViews-method}


\title{Longest Common Prefix/Suffix/Substring searching functions}

\description{
  Functions for searching the Longest Common Prefix/Suffix/Substring
  of two strings.

  WARNING: These functions are currently under development and might
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
    1st string, a character string or a \link{BString} object.
  }
  \item{s2}{
    2nd string, a character string or a \link{BString} object.
  }
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    A \link{BString} (or derived) object containing the subject string.
  }
  \item{maxlength.out}{
    The maximum length of the output i.e. the maximum number of views
    in the returned object.
  }
}

\details{
  [TODO]
}

\value{
  [TODO]
}

\seealso{
  \link{matchPattern},
  \link{BStringViews-class},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class}
}

\examples{
  [TODO]
}

\keyword{methods}

