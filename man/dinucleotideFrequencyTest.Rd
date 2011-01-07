\name{dinucleotideFrequencyTest}

\alias{dinucleotideFrequencyTest}
\alias{dinucleotideFrequencyTest,DNAStringSet-method}
\alias{dinucleotideFrequencyTest,RNAStringSet-method}


\title{Pearson's chi-squared Test and G-tests for String Position Dependence}

\description{
  Performs Person's chi-squared test, G-test, or William's corrected G-test to
  determine dependence between two nucleotide positions.
}

\usage{
dinucleotideFrequencyTest(x, i, j, test = c("chisq", "G", "adjG"),
                          simulate.p.value = FALSE, B = 2000)
}

\arguments{
  \item{x}{
    A \link{DNAStringSet} or \link{RNAStringSet} object.
  }
  \item{i, j}{
    Single integer values for positions to test for dependence.
  }
  \item{test}{
    One of \code{"chisq"} (Person's chi-squared test), \code{"G"} (G-test), or
    \code{"adjG"} (William's corrected G-test). See Details section.
  }
  \item{simulate.p.value}{a logical indicating whether to compute
    p-values by Monte Carlo simulation.
  }
  \item{B}{an integer specifying the number of replicates used in the
    Monte Carlo test.
  }
}

\details{
  The null and alternative hypotheses for this function are:
  \describe{
    \item{H0: }{positions \code{i} and \code{j} are independent}
    \item{H1: }{otherwise}
  }

  Let O and E be the observed and expected probabilities for base pair
  combinations at positions \code{i} and \code{j} respectively. Then the
  test statistics are calculated as:
  \describe{
    \item{\code{test="chisq"}: }{stat = sum(abs(O - E)^2/E)}
    \item{\code{test="G"}: }{stat = 2 * sum(O * log(O/E))}
    \item{\code{test="adjG"}: }{stat = 2 * sum(O * log(O/E))/q, where
      q = 1 + ((df - 1)^2 - 1)/(6*length(x)*(df - 2))}
  }
  Under the null hypothesis, these test statistics are approximately distributed
  chi-squared(df = ((distinct bases at i) - 1) * ((distinct bases at j) - 1)).
}

\value{
  An htest object. See help(chisq.test) for more details.
}

\references{
  Ellrott, K., Yang, C., Sladek, F.M., Jiang, T. (2002)
  "Identifying transcription factor binding sites through Markov chain optimations",
  Bioinformatics, 18 (Suppl. 2), S100-S109.

  Sokal, R.R., Rohlf, F.J. (2003)
  "Biometry: The Principle and Practice of Statistics in Biological Research",
   W.H. Freeman and Company, New York.

  Tomovic, A., Oakeley, E. (2007)
  "Position dependencies in transcription factor binding sites",
  Bioinformatics, 23, 933-941.

  Williams, D.A. (1976)
  "Improved Likelihood ratio tests for complete contingency tables",
  Biometrika, 63, 33-37.
}

\author{P. Aboyoun}

\seealso{
  \code{\link{nucleotideFrequencyAt}},
  \link{XStringSet-class},
  \code{\link[stats]{chisq.test}}
}

\examples{
  data(HNF4alpha)
  dinucleotideFrequencyTest(HNF4alpha, 1, 2)
  dinucleotideFrequencyTest(HNF4alpha, 1, 2, test = "G")
  dinucleotideFrequencyTest(HNF4alpha, 1, 2, test = "adjG")
}

\keyword{htest}
\keyword{distribution}
