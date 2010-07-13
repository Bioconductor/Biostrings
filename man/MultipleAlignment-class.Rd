\name{MultipleAlignment-class}
\docType{class}

% Classes:
\alias{class:MultipleAlignment}
\alias{MultipleAlignment-class}
\alias{MultipleAlignment}
\alias{class:DNAMultipleAlignment}
\alias{DNAMultipleAlignment-class}
\alias{DNAMultipleAlignment}
% \alias{class:AAMultipleAlignment}
% \alias{AAMultipleAlignment-class}
% \alias{AAMultipleAlignment}
% \alias{class:RNAMultipleAlignment}
% \alias{RNAMultipleAlignment-class}
% \alias{RNAMultipleAlignment}

% Accessor-like methods:
\alias{rownames,MultipleAlignment-method}
\alias{rownames<-,MultipleAlignment-method}
\alias{alignmentMask}
\alias{alignmentMask,MultipleAlignment-method}
\alias{alignmentMask<-}
\alias{alignmentMask<-,MultipleAlignment,NormalIRanges-method}
\alias{alignmentMask<-,MultipleAlignment,NULL-method}
\alias{sampleMask}
\alias{sampleMask,MultipleAlignment-method}
\alias{sampleMask<-}
\alias{sampleMask<-,MultipleAlignment,NormalIRanges-method}
\alias{sampleMask<-,MultipleAlignment,NULL-method}
\alias{unmask}
\alias{unmask,MultipleAlignment-method}


% Other methods:
\alias{consensusMatrix,MultipleAlignment-method}
\alias{consensusString,MultipleAlignment-method}
\alias{alphabetFrequency,MultipleAlignment-method}
\alias{ncol,MultipleAlignment-method}
\alias{nrow,MultipleAlignment-method}
\alias{show,MultipleAlignment-method}

% Coercion:
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
      rownames(x):
      Retrieves the sequence names for the x object.
    }
    \item{}{
      alignmentMask(x):
      Allows access to the alignmentMask slot.
    }
    \item{}{
      sampleMask(x):
      Allows access to the sampleMask slot.
    }
    \item{}{
      unmask(x):
      Allows removal of all mask slots.
    }
    \item{}{
      ncol(x):
      Returns the number of characters that each of the strings in the
      x object contains.
    }
    
  }
}

\section{Other methods}{
  In the code snippets below, x is a MultipleAlignment object.

  \describe{
    \item{}{
      nrow(x):
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
      as.character(x):
      Turns x into a character vector.
    }
  }
}

\section{Convenience}{

  \describe{
    \item{}{
      clustReader():
      creates a MultipleAlignment object from a clustalW .aln file.
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
## mask columns
subset <- cl
alignmentMask(subset) <-  as(IRanges(start=c(1,20),end=c(4,30)),"NormalIRanges")
subset

## mask rows
subset <- cl
sampleMask(subset) <-  as(IRanges(start=c(1,2),end=c(1,2)),"NormalIRanges")
subset

## remove masks
sampleMask(subset) <- NULL
subset

## remove all masks
alignmentMask(subset) <-  as(IRanges(start=c(1,20),end=c(4,30)),"NormalIRanges")
subset
subset <- unmask(subset)
subset

## List the names of the sequences in the alignment
rownames(cl)
}

\keyword{methods}
\keyword{classes}
