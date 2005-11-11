\name{simpleBioString-class}
\docType{class}
\alias{simpleBioString-class}
\alias{[,simpleBioString-method}
\alias{as.character,simpleBioString-method}
\alias{coerce,simpleBioString,character-method}
\alias{show,simpleBioString-method}
\alias{substr,simpleBioString-method}
\alias{expl}
\alias{extOrChar-class}

\title{Class "simpleBioString"}
\description{Support the BioString class API for R character data}
\section{Objects from the Class}{
Objects can be created by calls of the 
form \code{new("simpleBioString", alphabet, offsets, values, initialized, ...)}.
At this time the offsets component is not meaningful; values
slot holds the character vector of strings.
}
\section{Slots}{
	 \describe{
    \item{\code{values}:}{Object of class \code{"character"}: the biostrings. }
    \item{\code{alphabet}:}{Object of class \code{"BioAlphabet"}: the alphabet from
which the strings are built. }
    \item{\code{offsets}:}{Object of class \code{"matrix"}: currently ignored }
    \item{\code{initialized}:}{Object of class \code{"logical"}: currently ignored }
  }
}
\section{Extends}{
Class \code{"BioString"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "simpleBioString")}: extract from the \code{values} vector;
permits single subscripting only. }
    \item{as.character}{\code{signature(x = "simpleBioString")}: evaluates to the
values slot. }
    \item{coerce}{\code{signature(from = "simpleBioString", to = "character")}: formal
coercion to a character vector. }
    \item{show}{\code{signature(object = "simpleBioString")}: reporting method. }
    \item{substr}{\code{signature(x = "simpleBioString")}: applies base \code{substr}
and returns a \code{simpleBioString}. }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
\note{The \code{extOrChar} class union allows \code{\link{BioString-class}}
instances to have values that are either external pointer (references)
or strings.}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\link{CLASSNAME-class}} for links to other classes
%}
\examples{
lit <- c("AABABCB","BBABCB")
slit <- new("simpleBioString",values=lit, alphabet=StandardPeptideAlphabet())
slit
sslit <- substr(slit,2,3)
sslit
as.character(sslit)
}
\keyword{classes}
