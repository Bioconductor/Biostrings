\name{phiX174}
\alias{phiX174}
\alias{srPhiX174}
\alias{quPhiX174}
\alias{wtPhiX174}
\title{Coliphage phiX174 complete genome and sample short reads}
\description{
The complete genome for Coliphage phiX174 as well as a small number of Solexa
short reads, qualities associated with those short reads, and counts for the
number times those short reads occurred.
}
\details{
The \code{phiX174} object is a \code{DNAString} containing version
NC\_001422.1, GI:9626372 of the complete genome for Coliphage phiX174.

The \code{srPhiX174} object is a \code{DNAStringSet} containing short reads
from a Solexa machine.

The \code{quPhiX174} object is a \code{BStringSet} containing Solexa quality
scores associated with \code{srPhiX174}.

The \code{wtPhiX174} object is an integer vector containing counts associated
with \code{srPhiX174}.
}
\references{
\url{http://www.genome.jp/dbget-bin/www_bget?refseq+NC_001422}
}
\examples{
data(phiX174)
nchar(phiX174)
phiX174Substring <- substring(phiX174, 2793-34, 2811+34)

data(srPhiX174)
srPhiX174
quPhiX174
summary(wtPhiX174)

alignPhiX174 <-
  pairwiseAlignment(srPhiX174, phiX174Substring,
                    patternQuality = SolexaQuality(quPhiX174),
                    subjectQuality = SolexaQuality(99L),
                    type = "subjectOverlap")
summary(alignPhiX174, weight = wtPhiX174)
}
\keyword{datasets}
