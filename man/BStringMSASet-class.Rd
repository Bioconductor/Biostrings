\name{BStringMSASet-class}
\docType{class}

% Classes:
\alias{class:BStringMSASet}
\alias{BStringMSASet-class}
\alias{BStringMSASet}

% Accessor-like methods:
\alias{names,BStringMSASet-method}
\alias{names<-,BStringMSASet-method}
\alias{nchar,BStringMSASet-method}
\alias{[,BStringMSASet-method}
\alias{[[,BStringMSASet-method}


% Other methods:
\alias{narrow,BStringMSASet-method}
\alias{xscat,BStringMSASet-method}
\alias{consensusMatrix,BStringMSASet-method}
\alias{consensusString,BStringMSASet-method}
\alias{alphabetFrequency,BStringMSASet-method}
\alias{length,BStringMSASet-method}
\alias{show,BStringMSASet-method}

% Coercion:
\alias{as.list,BStringMSASet-method}
\alias{as.character,BStringMSASet-method}

% Convenience Function:
\alias{subsetCollumns}
\alias{clustReader}

\title{BStringMSASet objects}

\description{
   The BStringMSASet class is a containers for storing multiple sequence
   alignments.  
}

\details{
  The BStringMSASet is meant to hold and represent multiple sequence
  alignments.  Internally, these a series of equivalent length strings,
  each of which has had gaps inserted etc. such that they represent the
  output of an alignment algorithm.  There are also a series of simple
  methods provided so that these alignments can have rows and columns
  removed, and so that the consensus sequence can be obtained etc.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a BStringMSASet object.

  \describe{
    \item{}{
      \code{names(x)}:
      Retrieves the sequence names for the \code{x} object.
    }
    \item{}{
      \code{nchar(x)}:
      Returns the number of characters that each of the strings in the
      \code{x} object contains.
    }
    
  }
}



\section{Other methods}{
  In the code snippets below, \code{x} is a BStringMSASet object.

  \describe{
    \item{}{
      \code{narrow(x,start,end,width,use.names)}:
      Narrows the \code{x} object to the range specified by
      \code{start}, \code{end}, and \code{width}. 
    }
    \item{}{
      \code{xscat(...)}:
      Combines the columns of the BStringMSASet objects passed in.
    }
    \item{}{
      \code{length(x)}:
      Returns the number of sequences aligned in \code{x}.
    }
    \item{}{
      \code{show(x)}:
      Shows an alignment of the 1st section of \code{x}.
    }    
    \item{}{
      \code{consensusMatrix(x)}:
      Calculates that Consensus matrix for \code{x}.
    }    
    \item{}{
      \code{consensusString(x)}:
      Calculates a Consensus string for \code{x}.
    }
    \item{}{
      \code{alphabetFrequency(x)}:
      Calculates the alphabet frequency for \code{x}.
    }        
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a BStringMSASet object.

  \describe{
    \item{}{
      \code{as.list(x, "XStringViews")}:
      Turns \code{x} into an \code{list}.
    }
    \item{}{
      \code{as.character(x, "XStringViews")}:
      Turns \code{x} into an \code{character vector}.
    }
  }
}


\section{Convenience Function}{
  In the code snippets below, \code{x} is a BStringMSASet object.

  \describe{
    \item{}{
      \code{subsetColumns(x, start, end)}:
      Subsets \code{x} by using the vectors contained in \code{start}
      and \code{end} such that the returned BStringMSASet object
      contains only the columns specified.
    }
  }
}



\author{M. Carlson}


\seealso{
  \link{XString-class},
  \link{XStringViews-class}
}


\examples{
## Create an object from file
cl = clustReader(file =
                  system.file("extdata","clustalw2.aln",
                               package="Biostrings"))

## Remove unwanted rows
cl[3]

## Remove all but a range of columns
narrow(cl, start= 9200, end = 10000)

## Select out a series of ranges of columns
subsetCollumns(cl, start=c(1,11), end=c(10,20))

## List the names of the sequences in the alignment
names(cl)
}

\keyword{methods}
\keyword{classes}
