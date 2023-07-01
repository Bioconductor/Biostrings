\name{phiX174Phage}
\alias{phiX174Phage}
\alias{srPhiX174}
\alias{quPhiX174}
\alias{wtPhiX174}
\title{Versions of bacteriophage phiX174 complete genome and sample short reads}
\description{
Six versions of the complete genome for bacteriophage \eqn{\phi} X174 as well
as a small number of Solexa short reads, qualities associated with those short
reads, and counts for the number times those short reads occurred.
}
\details{
The \code{phiX174Phage} object is a \code{DNAStringSet} containing the
following six naturally occurring versions of the bacteriophage \eqn{\phi} X174
genome cited in Smith et al.:

\describe{
  \item{Genbank:}{The version of the genome from GenBank (NC_001422.1, GI:9626372).}
  \item{RF70s:}{A preparation of \eqn{\phi} X double-stranded replicative form (RF) of DNA
    by Clyde A. Hutchison III from the late 1970s.}
  \item{SS78:}{A preparation of \eqn{\phi} X virion single-stranded DNA from 1978.}
  \item{Bull:}{The sequence of wild-type \eqn{\phi} X used by Bull et al.}
  \item{G'97:}{The \eqn{\phi} X replicative form (RF) of DNA from Bull et al.}
  \item{NEB'03:}{A \eqn{\phi} X replicative form (RF) of DNA from New England BioLabs (NEB).}  
}

The \code{srPhiX174} object is a \code{DNAStringSet} containing short reads
from a Solexa machine.

The \code{quPhiX174} object is a \code{BStringSet} containing Solexa quality
scores associated with \code{srPhiX174}.

The \code{wtPhiX174} object is an integer vector containing counts associated
with \code{srPhiX174}.
}
\references{
\url{http://www.genome.jp/dbget-bin/www_bget?refseq+NC_001422}

Bull, J. J., Badgett, M. R., Wichman, H. A., Huelsenbeck, Hillis, D. M., Gulati, A., Ho, C. & Molineux, J.
(1997) Genetics 147, 1497-1507.

Smith, Hamilton O.; Clyde A. Hutchison, Cynthia Pfannkoch, J. Craig Venter (2003-12-23).
"Generating a synthetic genome by whole genome assembly: \{phi\}X174 bacteriophage from synthetic oligonucleotides".
Proceedings of the National Academy of Sciences 100 (26): 15440-15445. doi:10.1073/pnas.2237126100.
}
\examples{
data(phiX174Phage)
nchar(phiX174Phage)
genBankPhage <- phiX174Phage[[1]]
genBankSubstring <- substring(genBankPhage, 2793-34, 2811+34)

data(srPhiX174)
srPhiX174
quPhiX174
summary(wtPhiX174)

alignPhiX174 <-
  pairwiseAlignment(srPhiX174, genBankSubstring,
                    patternQuality = SolexaQuality(quPhiX174),
                    subjectQuality = SolexaQuality(99L),
                    type = "global-local")
summary(alignPhiX174, weight = wtPhiX174)
}
\keyword{datasets}
