\name{MultipleAlignment-class}
\docType{class}

% Classes:
\alias{class:MultipleAlignment}
\alias{MultipleAlignment-class}
\alias{MultipleAlignment}

% Accessor-like methods:
\alias{names,MultipleAlignment-method}
\alias{names<-,MultipleAlignment-method}
\alias{nchar,MultipleAlignment-method}
\alias{[,MultipleAlignment-method}
\alias{[[,MultipleAlignment-method}

% Other methods:
\alias{narrow,MultipleAlignment-method}
\alias{xscat,MultipleAlignment-method}
\alias{consensusMatrix,MultipleAlignment-method}
\alias{consensusString,MultipleAlignment-method}
\alias{alphabetFrequency,MultipleAlignment-method}
\alias{length,MultipleAlignment-method}
\alias{show,MultipleAlignment-method}

% Coercion:
\alias{as.list,MultipleAlignment-method}
\alias{as.character,MultipleAlignment-method}

% Convenience Function:
\alias{clustReader}


\title{MultipleAlignment objects}

\description{
   The MultipleAlignment class is a containers for storing multiple sequence
   alignments.  
}

\details{
  The MultipleAlignment is meant to hold and represent multiple sequence
  alignments.  Internally, these are a series of equivalent length strings,
  each of which has had gaps inserted etc. such that they represent the
  output of an alignment algorithm.  There are also a series of simple
  methods provided so that these alignments can have rows and columns
  removed, and so that the consensus sequence can be obtained etc.
}

\section{Accessor methods}{
  In the code snippets below, x is a MultipleAlignment object.

  \describe{
    \item{}{
      names(x):
      Retrieves the sequence names for the x object.
    }
    \item{}{
      nchar(x):
      Returns the number of characters that each of the strings in the
      x object contains.
    }
    
  }
}

\section{Other methods}{
  In the code snippets below, x is a MultipleAlignment object.

  \describe{
    \item{}{
      narrow(x,start,end,width,use.names):
      Narrows the x object to the range specified by
      start, end, and width. 
    }
    \item{}{
      xscat(...):
      Combines the columns of the MultipleAlignment objects passed in.
    }
    \item{}{
      length(x):
      Returns the number of sequences aligned in x.
    }
    \item{}{
      show(x):
      Shows an alignment of the 1st section of x.
    }    
    \item{}{
      consensusMatrix(x):
      Calculates that Consensus matrix for x.
    }    
    \item{}{
      consensusString(x):
      Calculates a Consensus string for x.
    }
    \item{}{
      alphabetFrequency(x):
      Calculates the alphabet frequency for x.
    }        
  }
}

\section{Coercion}{
  In the code snippets below, x is a MultipleAlignment object.

  \describe{
    \item{}{
      as.list(x):
      Turns x into a list.
    }
    \item{}{
      as.character(x):
      Turns x into a character vector.
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
cl <- clustReader(file =
                  system.file("extdata","clustalw2.aln",
                               package="Biostrings"))

subset <- cl
alignmentMask(subset) <-  as(IRanges(start=c(1,20),end=c(4,30)),"NormalIRanges")
subset


## List the names of the sequences in the alignment
rownames(cl)
}

\keyword{methods}
\keyword{classes}
