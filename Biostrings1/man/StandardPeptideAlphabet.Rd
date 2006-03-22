\name{StandardPeptideAlphabet}
\alias{StandardPeptideAlphabet}
\alias{AAString}
\alias{rpeptide}
\title{alphabet for proteins  }
\description{alphabet for proteins; generate random peptide strings
}
\usage{
StandardPeptideAlphabet()
rpeptide(nchar, base)
}
\arguments{
  \item{nchar}{ length of string }
  \item{base}{ character string giving codes that can be used }
}

\details{Has codes for amino acids; prepared by
Saikay Debroy.
}
\value{
a \code{\link{BioAlphabet-class}} instance
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
StandardPeptideAlphabet()
AAString("AABB")
rpeptide(12, StandardPeptideAlphabet()@letters)
}
\keyword{ models }
