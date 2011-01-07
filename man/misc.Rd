\name{misc}

\alias{misc}

\alias{N50}


\title{Some miscellaneous stuff}

\description{
  Some miscellaneous stuff.
}

\usage{
N50(csizes)
}

\arguments{
  \item{csizes}{
    A vector containing the contig sizes.
  }
}

\section{The N50 contig size}{
  \bold{Definition}
  The N50 contig size of an assembly (aka the N50 value) is the
  size of the largest contig such that the contigs larger than
  that have at least 50\% the bases of the assembly.

  \bold{How is it calculated?}
  It is calculated by adding the sizes of the biggest contigs until
  you reach half the total size of the contigs.
  The N50 value is then the size of the contig that was added last
  (i.e. the smallest of the big contigs covering 50\% of the genome).

  \bold{What for?}
  The N50 value is a standard measure of the quality of a de novo
  assembly.
}

\value{
  \code{N50}: The N50 value as an integer.
}

\author{Nicolas Delhomme <delhomme@embl.de>}

\seealso{
  \link{XStringSet-class}
}

\examples{
  # Generate 10 random contigs of sizes comprised between 100 and 10000:
  my.contig <- DNAStringSet(
                 sapply(
                   sample(c(100:10000), 10),
                   function(size)
                       paste(sample(DNA_BASES, size, replace=TRUE), collapse="")
                 )
               )

  # Get their sizes:
  my.size <- width(my.contig)

  # Calculate the N50 value of this set of contigs:
  my.contig.N50 <- N50(my.size)
}

\keyword{methods}
\keyword{manip}
