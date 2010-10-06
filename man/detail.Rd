\name{detail}

\alias{detail}

\title{Show (display) detailed object content}

\description{

  This is a variant of \code{\link{show}}, offering a more detailed
  display of object content.

}
\usage{

detail(x, ...)

}

\arguments{

  \item{x}{An object. The default simply invokes \code{\link{show}}.}

  \item{...}{Additional arguments. The default definition makes no use
    of these arguments.}

}
\value{
  
  None; the function is invoked for its side effect (detailed display of
  object content).

}
\author{Martin Morgan}
\examples{
origMAlign <-
  read.DNAMultipleAlignment(filepath =
                            system.file("extdata",
                                        "msx2_mRNA.aln",
                                        package="Biostrings"),
                            format="clustal")
detail(origMAlign)
}
\keyword{manip}
