\name{MultipleAlignment-class}
\docType{class}

% Classes:
\alias{class:MultipleAlignment}
\alias{MultipleAlignment-class}
\alias{MultipleAlignment}
\alias{class:DNAMultipleAlignment}
\alias{DNAMultipleAlignment-class}
\alias{DNAMultipleAlignment}
\alias{class:RNAMultipleAlignment}
\alias{RNAMultipleAlignment-class}
\alias{RNAMultipleAlignment}
\alias{class:AAMultipleAlignment}
\alias{AAMultipleAlignment-class}
\alias{AAMultipleAlignment}

% Accessor-like methods:
\alias{unmasked,MultipleAlignment-method}
\alias{rownames,MultipleAlignment-method}
\alias{rownames<-,MultipleAlignment-method}
\alias{rowmask}
\alias{rowmask,MultipleAlignment-method}
\alias{rowmask<-}
\alias{rowmask<-,MultipleAlignment,NULL-method}
\alias{rowmask<-,MultipleAlignment,ANY-method}
\alias{colmask}
\alias{colmask,MultipleAlignment-method}
\alias{colmask<-}
\alias{colmask<-,MultipleAlignment,NULL-method}
\alias{colmask<-,MultipleAlignment,ANY-method}
\alias{maskMotif,MultipleAlignment,ANY-method}
\alias{maskGaps}
\alias{maskGaps,MultipleAlignment-method}
\alias{nrow,MultipleAlignment-method}
\alias{ncol,MultipleAlignment-method}
\alias{dim,MultipleAlignment-method}
\alias{maskednrow}
\alias{maskednrow,MultipleAlignment-method}
\alias{maskedncol}
\alias{maskedncol,MultipleAlignment-method}
\alias{maskeddim}
\alias{maskeddim,MultipleAlignment-method}
\alias{maskedratio,MultipleAlignment-method}
\alias{nchar,MultipleAlignment-method}
\alias{seqtype,MultipleAlignment-method}

% Read functions:
\alias{readDNAMultipleAlignment}
\alias{readRNAMultipleAlignment}
\alias{readAAMultipleAlignment}

% Write functions:
\alias{write.phylip}

% Coercion:
\alias{coerce,MultipleAlignment,DNAStringSet-method}
\alias{coerce,MultipleAlignment,RNAStringSet-method}
\alias{coerce,MultipleAlignment,AAStringSet-method}
\alias{coerce,MultipleAlignment,BStringSet-method}
\alias{coerce,character,DNAMultipleAlignment-method}
\alias{coerce,character,RNAMultipleAlignment-method}
\alias{coerce,character,AAMultipleAlignment-method}
\alias{as.character,MultipleAlignment-method}
\alias{as.matrix,MultipleAlignment-method}

% Utilities:
\alias{consensusMatrix,MultipleAlignment-method}
\alias{consensusString,MultipleAlignment-method}
\alias{consensusString,DNAMultipleAlignment-method}
\alias{consensusString,RNAMultipleAlignment-method}
\alias{consensusString,AAMultipleAlignment-method}
\alias{consensusViews}
\alias{consensusViews,MultipleAlignment-method}
\alias{consensusViews,DNAMultipleAlignment-method}
\alias{consensusViews,RNAMultipleAlignment-method}
\alias{consensusViews,AAMultipleAlignment-method}
\alias{alphabetFrequency,MultipleAlignment-method}

% show style methods:
\alias{show,MultipleAlignment-method}
\alias{detail,MultipleAlignment-method}

\title{MultipleAlignment objects}

\description{
  IMPORTANT NOTE:

  Starting with Biostrings 2.81.3 (Bioconductor 3.24),
  \link[MultipleAlignment]{MultipleAlignment} objects and associated
  functionality are implemented and documented in the \pkg{MultipleAlignment}
  package.

  To access the man page for these objects, load the \pkg{MultipleAlignment}
  package with \code{library(MultipleAlignment)} and use
  \code{?MultipleAlignment::MultipleAlignment}.

  If the \pkg{MultipleAlignment} package is not installed, you will first
  need to install it with \code{BiocManager::install("MultipleAlignment")}.

  If you were trying to access the man page for
  \code{\link[MultipleAlignment]{detail}()}, use
  \code{?MultipleAlignment::detail}.
}

\keyword{methods}
\keyword{classes}
