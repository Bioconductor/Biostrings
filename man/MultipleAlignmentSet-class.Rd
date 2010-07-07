\name{MultipleAlignmentSet-class}
\docType{class}

% Classes:
\alias{class:MultipleAlignmentSet}
\alias{MultipleAlignmentSet-class}
\alias{MultipleAlignmentSet}

% Accessor-like methods:
\alias{names,MultipleAlignmentSet-method}
\alias{names<-,MultipleAlignmentSet-method}
\alias{nchar,MultipleAlignmentSet-method}
\alias{[,MultipleAlignmentSet-method}
\alias{[[,MultipleAlignmentSet-method}

% Other methods:
\alias{narrow,MultipleAlignmentSet-method}
\alias{xscat,MultipleAlignmentSet-method}
\alias{consensusMatrix,MultipleAlignmentSet-method}
\alias{consensusString,MultipleAlignmentSet-method}
\alias{alphabetFrequency,MultipleAlignmentSet-method}
\alias{length,MultipleAlignmentSet-method}
\alias{show,MultipleAlignmentSet-method}

% Coercion:
\alias{as.list,MultipleAlignmentSet-method}
\alias{as.character,MultipleAlignmentSet-method}

% Convenience Function:
\alias{subsetColumns}
\alias{clustReader}


\title{MultipleAlignmentSet objects}

\description{
   The MultipleAlignmentSet class is a containers for storing multiple sequence
   alignments.  
}

\details{
  The MultipleAlignmentSet is meant to hold and represent multiple sequence
  alignments.  Internally, these are a series of equivalent length strings,
  each of which has had gaps inserted etc. such that they represent the
  output of an alignment algorithm.  There are also a series of simple
  methods provided so that these alignments can have rows and columns
  removed, and so that the consensus sequence can be obtained etc.
}

\section{Accessor methods}{
  In the code snippets below, x is a MultipleAlignmentSet object.

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
  In the code snippets below, x is a MultipleAlignmentSet object.

  \describe{
    \item{}{
      narrow(x,start,end,width,use.names):
      Narrows the x object to the range specified by
      start, end, and width. 
    }
    \item{}{
      xscat(...):
      Combines the columns of the MultipleAlignmentSet objects passed in.
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
  In the code snippets below, x is a MultipleAlignmentSet object.

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


\section{Convenience Function}{
  In the code snippets below, x is a MultipleAlignmentSet object.

  \describe{
    \item{}{
      subsetColumns(x, start, end):
      Subsets x by using the vectors contained in start
      and end such that the returned MultipleAlignmentSet object
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
subsetColumns(cl, start=c(1,11), end=c(10,20))

## List the names of the sequences in the alignment
names(cl)
}

\keyword{methods}
\keyword{classes}
