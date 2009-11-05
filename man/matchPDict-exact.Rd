\name{matchPDict}

\alias{matchPDict-exact}

\alias{matchPDict}
\alias{matchPDict,XString-method}
\alias{matchPDict,XStringSet-method}
\alias{matchPDict,XStringViews-method}
\alias{matchPDict,MaskedXString-method}

\alias{countPDict}
\alias{countPDict,XString-method}
\alias{countPDict,XStringSet-method}
\alias{countPDict,XStringViews-method}
\alias{countPDict,MaskedXString-method}

\alias{whichPDict}
\alias{whichPDict,XString-method}

\alias{vmatchPDict}
\alias{vmatchPDict,ANY-method}
\alias{vmatchPDict,XString-method}
\alias{vmatchPDict,MaskedXString-method}

\alias{vcountPDict}
\alias{vcountPDict,XString-method}
\alias{vcountPDict,XStringSet-method}
\alias{vcountPDict,XStringViews-method}
\alias{vcountPDict,MaskedXString-method}

\alias{vwhichPDict}
\alias{vwhichPDict,XString-method}
\alias{vwhichPDict,XStringSet-method}
\alias{vwhichPDict,XStringViews-method}
\alias{vwhichPDict,MaskedXString-method}

% Functions:
\alias{extractAllMatches}


\title{Searching a sequence for patterns stored in a preprocessed dictionary}

\description{
  A set of functions for finding all the occurrences (aka "matches" or "hits")
  of a set of patterns (aka the dictionary) in a reference
  sequence or set of reference sequences (aka the subject)

  The following functions differ in what they return: \code{matchPDict}
  returns the "where" information i.e. the positions in the subject of all the
  occurrences of every pattern; \code{countPDict} returns the "how many
  times" information i.e. the number of occurrences for each pattern;
  and \code{whichPDict} returns the "who" information i.e. which patterns
  in the preprocessed dictionary have at least one match.

  \code{vcountPDict} and \code{vwhichPDict} are vectorized versions
  of \code{countPDict} and \code{whichPDict}, respectively, that is,
  they work on a set of reference sequences in a vectorized fashion.

  This man page shows how to use these functions (aka the \code{*PDict}
  functions) for exact matching of a constant width dictionary i.e.
  a dictionary where all the patterns have the same length (same number
  of nucleotides).

  See \code{?`\link{matchPDict-inexact}`} for how to use these functions
  for inexact matching or when the original dictionary has a variable width.
}

\usage{
  matchPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, fixed=TRUE, verbose=FALSE)
  countPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, fixed=TRUE, verbose=FALSE)
  whichPDict(pdict, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, fixed=TRUE, verbose=FALSE)

  vcountPDict(pdict, subject, algorithm="auto",
              max.mismatch=0, min.mismatch=0, fixed=TRUE,
              collapse=FALSE, weight=1L, verbose=FALSE)
  vwhichPDict(pdict, subject, algorithm="auto",
              max.mismatch=0, min.mismatch=0, fixed=TRUE, verbose=FALSE)
}

\arguments{
  \item{pdict}{
    A \link{PDict} object containing the preprocessed dictionary.
  }
  \item{subject}{
    An \link{XString} or \link{MaskedXString} object containing the
    subject sequence for \code{matchPDict}, \code{countPDict} and \code{whichPDict}.

    An \link{XStringSet} object containing the subject sequences
    for \code{vcountPDict} and \code{vwhichPDict}.

    For now, only subjects of base class \link{DNAString} are supported.
  }
  \item{algorithm}{
    Not supported yet.
  }
  \item{max.mismatch, min.mismatch}{
    The maximum and minimum number of mismatching letters allowed (see
    \code{?\link{isMatching}} for the details).
    This man page focuses on exact matching of a constant width
    dictionary so \code{max.mismatch=0} in the examples below.
    See \code{?`\link{matchPDict-inexact}`} for inexact matching.
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{?\link{isMatching}} for the details).
    This man page focuses on exact matching of a constant width
    dictionary so \code{fixed=TRUE} in the examples below.
    See \code{?`\link{matchPDict-inexact}`} for inexact matching.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{collapse,weight}{
    \code{collapse} must be \code{FALSE}, \code{1}, or \code{2}.

    If \code{collapse=FALSE} (the default), then \code{weight} is ignored
    and \code{vcountPDict} returns the full matrix of counts (\code{M0}).
    If \code{collapse=1}, then \code{M0} is collapsed "horizontally"
    i.e. it is turned into a vector with \code{length} equal to
    \code{length(pdict)}.
    If \code{weight=1L} (the default), then this vector is defined by
    \code{rowSums(M0)}.
    If \code{collapse=2}, then \code{M0} is collapsed "vertically"
    i.e. it is turned into a vector with \code{length} equal to
    \code{length(subject)}.
    If \code{weight=1L} (the default), then this vector is defined by
    \code{colSums(M0)}.

    If \code{collapse=1} or \code{collapse=2}, then the elements in
    \code{subject} (\code{collapse=1}) or in \code{pdict} (\code{collapse=2})
    can be weighted thru the \code{weight} argument.
    In that case, the returned vector is defined by
    \code{M0 \%*\% rep(weight, length.out=length(subject))}
    and \code{rep(weight, length.out=length(pdict)) \%*\% M0},
    respectively.
  }
}

\details{
  In this man page, we assume that you know how to preprocess a dictionary
  of DNA patterns that can then be used with any of the \code{*PDict}
  functions described here. Please see \code{?\link{PDict}} if you don't.

  When using the \code{*PDict} functions for exact matching of a constant
  width dictionary, the standard way to preprocess the original dictionary
  is by calling the \code{\link{PDict}} constructor on it with no extra
  arguments. This returns the preprocessed dictionary in a \link{PDict}
  object that can be used with any of the \code{*PDict} functions.
}

\value{
  If \code{M} denotes the number of patterns in the \code{pdict}
  argument (\code{M <- length(pdict)}), then \code{matchPDict} returns
  an \link{MIndex} object of length \code{M},
  and \code{countPDict} an integer vector of length \code{M}.

  \code{whichPDict} returns an integer vector made of the indices of the
  patterns in the \code{pdict} argument that have at least one match.

  If \code{N} denotes the number of sequences in the \code{subject}
  argument (\code{N <- length(subject)}), then \code{vcountPDict}
  returns an integer matrix with \code{M} rows and \code{N} columns,
  unless the \code{collapse} argument is used. In that case, depending
  on the type of \code{weight}, an integer or numeric vector is returned
  (see above for the details).

  \code{vwhichPDict} returns a list of \code{N} integer vectors.
}

\author{H. Pages}

\references{
  Aho, Alfred V.; Margaret J. Corasick (June 1975). "Efficient string
  matching: An aid to bibliographic search".
  Communications of the ACM 18 (6): 333-340.
}

\seealso{
  \link{PDict-class},
  \link{MIndex-class},
  \link{matchPDict-inexact},
  \code{\link{isMatching}},
  \code{\link{coverage,MIndex-method}},
  \code{\link{matchPattern}},
  \code{\link{alphabetFrequency}},
  \link{DNAString-class},
  \link{DNAStringSet-class},
  \link{XStringViews-class},
  \link{MaskedDNAString-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. A SIMPLE EXAMPLE OF EXACT MATCHING
  ## ---------------------------------------------------------------------

  ## Creating the pattern dictionary:
  library(drosophila2probe)
  dict0 <- DNAStringSet(drosophila2probe)
  dict0                                # The original dictionary.
  length(dict0)                        # Hundreds of thousands of patterns.
  pdict0 <- PDict(dict0)               # Store the original dictionary in
                                       # a PDict object (preprocessing).

  ## Using the pattern dictionary on chromosome 3R:
  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chr3R <- Dmelanogaster$chr3R         # Load chromosome 3R
  chr3R
  mi0 <- matchPDict(pdict0, chr3R)     # Search...

  ## Looking at the matches:
  start_index <- startIndex(mi0)       # Get the start index.
  length(start_index)                  # Same as the original dictionary.
  start_index[[8220]]                  # Starts of the 8220th pattern.
  end_index <- endIndex(mi0)           # Get the end index.
  end_index[[8220]]                    # Ends of the 8220th pattern.
  count_index <- countIndex(mi0)       # Get the number of matches per pattern.
  count_index[[8220]]
  mi0[[8220]]                          # Get the matches for the 8220th pattern.
  start(mi0[[8220]])                   # Equivalent to startIndex(mi0)[[8220]].
  sum(count_index)                     # Total number of matches.
  table(count_index)
  i0 <- which(count_index == max(count_index))
  pdict0[[i0]]                         # The pattern with most occurrences.
  mi0[[i0]]                            # Its matches as an IRanges object.
  Views(chr3R, mi0[[i0]])              # And as an XStringViews object.

  ## Get the coverage of the original subject:
  cov3R <- as.integer(coverage(mi0, width=length(chr3R)))
  max(cov3R)
  mean(cov3R)
  sum(cov3R != 0) / length(cov3R)      # Only 2.44\% of chr3R is covered.
  if (interactive()) {
    plotCoverage <- function(cx, start, end)
    {
      plot.new()
      plot.window(c(start, end), c(0, 20))
      axis(1)
      axis(2)
      axis(4)
      lines(start:end, cx[start:end], type="l")
    }
    plotCoverage(cov3R, 27600000, 27900000)
  }

  ## ---------------------------------------------------------------------
  ## B. NAMING THE PATTERNS
  ## ---------------------------------------------------------------------

  ## The names of the original patterns, if any, are propagated to the
  ## PDict and MIndex objects:
  names(dict0) <- mkAllStrings(letters, 4)[seq_len(length(dict0))]
  dict0
  dict0[["abcd"]]
  pdict0n <- PDict(dict0)
  names(pdict0n)[1:30]
  pdict0n[["abcd"]]
  mi0n <- matchPDict(pdict0n, chr3R)
  names(mi0n)[1:30]
  mi0n[["abcd"]]

  ## This is particularly useful when unlisting an MIndex object:
  unlist(mi0)[1:10]
  unlist(mi0n)[1:10]  # keep track of where the matches are coming from

  ## ---------------------------------------------------------------------
  ## C. PERFORMANCE
  ## ---------------------------------------------------------------------

  ## If getting the number of matches is what matters only (without
  ## regarding their positions), then countPDict() will be faster,
  ## especially when there is a high number of matches:

  count_index0 <- countPDict(pdict0, chr3R)
  stopifnot(identical(count_index0, count_index))

  if (interactive()) {
    ## What's the impact of the dictionary width on performance?
    ## Below is some code that can be used to figure out (will take a long
    ## time to run). For different widths of the original dictionary, we
    ## look at:
    ##   o pptime: preprocessing time (in sec.) i.e. time needed for
    ##             building the PDict object from the truncated input
    ##             sequences;
    ##   o nnodes: nb of nodes in the resulting Aho-Corasick tree;
    ##   o nupatt: nb of unique truncated input sequences;
    ##   o matchtime: time (in sec.) needed to find all the matches;
    ##   o totalcount: total number of matches.
    getPDictStats <- function(dict, subject)
    {
      ans_width <- width(dict[1])
      ans_pptime <- system.time(pdict <- PDict(dict))[["elapsed"]]
      pptb <- pdict@threeparts@pptb
      ans_nnodes <- nnodes(pptb)
      ans_nupatt <- sum(!duplicated(pdict))
      ans_matchtime <- system.time(
                         mi0 <- matchPDict(pdict, subject)
                       )[["elapsed"]]
      ans_totalcount <- sum(countIndex(mi0))
      list(
        width=ans_width,
        pptime=ans_pptime,
        nnodes=ans_nnodes,
        nupatt=ans_nupatt,
        matchtime=ans_matchtime,
        totalcount=ans_totalcount
      )
    }
    stats <- lapply(8:25,
                 function(width)
                     getPDictStats(DNAStringSet(dict0, end=width), chr3R))
    stats <- data.frame(do.call(rbind, stats))
    stats
  }

  ## ---------------------------------------------------------------------
  ## D. vcountPDict()
  ## ---------------------------------------------------------------------
  subject <- Dmelanogaster$upstream1000[1:200]
  subject
  mat1 <- vcountPDict(pdict0, subject)
  dim(mat1)  # length(pdict0) x length(subject)
  nhit_per_probe <- rowSums(mat1)
  table(nhit_per_probe)

  ## Without vcountPDict(), 'mat1' could have been computed with:
  mat2 <- sapply(unname(subject), function(x) countPDict(pdict0, x))
  stopifnot(identical(mat1, mat2))
  ## but using vcountPDict() is faster (10x or more, depending of the
  ## average length of the sequences in 'subject').

  if (interactive()) {
    ## This will fail (with message "allocMatrix: too many elements
    ## specified") because, on most platforms, vectors and matrices in R
    ## are limited to 2^31 elements:
    subject <- Dmelanogaster$upstream1000
    vcountPDict(pdict0, subject)
    length(pdict0) * length(Dmelanogaster$upstream1000)
    1 * length(pdict0) * length(Dmelanogaster$upstream1000)  # > 2^31
    ## But this will work:
    nhit_per_seq <- vcountPDict(pdict0, subject, collapse=2)
    sum(nhit_per_seq >= 1)  # nb of subject sequences with at least 1 hit
    table(nhit_per_seq)
    which(nhit_per_seq == 37)  # 603
    sum(countPDict(pdict0, subject[[603]]))  # 37
  }

  ## ---------------------------------------------------------------------
  ## E. RELATIONSHIP BETWEEN vcountPDict(), countPDict() AND
  ## vcountPattern()
  ## ---------------------------------------------------------------------
  dict3 <- DNAStringSet(mkAllStrings(DNA_BASES, 3))  # all trinucleotides
  dict3
  pdict3 <- PDict(dict3)
  subject <- Dmelanogaster$upstream1000
  subject

  ## The 3 following calls are equivalent (from faster to slower):
  mat3a <- vcountPDict(pdict3, subject)
  mat3b <- sapply(dict3, function(pattern) vcountPattern(pattern, subject))
  mat3c <- sapply(unname(subject), function(x) countPDict(pdict3, x))
  stopifnot(identical(mat3a, t(mat3b)))
  stopifnot(identical(mat3a, mat3c))

  ## The 2 following calls are equivalent (from faster to slower):
  nhitpp3a <- vcountPDict(pdict3, subject, collapse=1)  # rowSums(mat3a)
  nhitpp3b <- sapply(dict3, function(pattern) sum(vcountPattern(pattern, subject)))
  stopifnot(identical(nhitpp3a, nhitpp3b))

  ## The 2 following calls are equivalent (from faster to slower):
  nhitps3a <- vcountPDict(pdict3, subject, collapse=2)  # colSums(mat3a)
  nhitps3b <- sapply(unname(subject), function(x) sum(countPDict(pdict3, x)))
  stopifnot(identical(nhitps3a, nhitps3b))

  ## ---------------------------------------------------------------------
  ## F. vwhichPDict()
  ## ---------------------------------------------------------------------
  ## The 3 following calls are equivalent (from faster to slower):
  vwp3a <- vwhichPDict(pdict3, subject)
  vwp3b <- lapply(seq_len(ncol(mat3a)), function(j) which(mat3a[ , j] != 0L))
  vwp3c <- lapply(unname(subject), function(x) whichPDict(pdict3, x))
  stopifnot(identical(vwp3a, vwp3b))
  stopifnot(identical(vwp3a, vwp3c))

  table(sapply(vwp3a, length))
  which.min(sapply(vwp3a, length))
  ## Get the trinucleotides not represented in reference sequence 9181:
  dict3[-vwp3a[[9181]]]  # 21 trinucleotides

  ## ---------------------------------------------------------------------
  ## G. MAPPING PROBE SET IDS BETWEEN CHIPS WITH vwhichPDict()
  ## ---------------------------------------------------------------------
  ## Here we show a simple (and very naive) algorithm for mapping probe
  ## set IDs between the hgu95av2 and hgu133a chips (Affymetrix).
  ## 2 probe set IDs are considered mapped iff they share at least one
  ## probe.
  ## WARNING: This example takes about 25 minutes to run.
  if (interactive()) {

    library(hgu95av2probe)
    library(hgu133aprobe)
    probes1 <- DNAStringSet(hgu95av2probe)
    probes2 <- DNAStringSet(hgu133aprobe)
    pdict2 <- PDict(probes2)

    ## Get the mapping from probes1 to probes2 (based on exact matching):
    map1to2 <- vwhichPDict(pdict2, probes1)  # takes about 10 minutes

    ## The following helper function uses the probe level mapping to induce
    ## the mapping at the probe set IDs level (from hgu95av2 to hgu133a).
    ## To keep things simple, 2 probe set IDs are considered mapped iff
    ## each of them contains at least one probe mapped to one probe of
    ## the other:
    mapProbeSetIDs1to2 <- function(psID)
      unique(hgu133aprobe$Probe.Set.Name[unlist(
        map1to2[hgu95av2probe$Probe.Set.Name == psID]
      )])

    ## Use the helper function to build the complete mapping:
    psIDs1 <- unique(hgu95av2probe$Probe.Set.Name)
    mapPSIDs1to2 <- lapply(psIDs1, mapProbeSetIDs1to2)  # about 3 min.
    names(mapPSIDs1to2) <- psIDs1

    ## Do some basic stats:
    table(sapply(mapPSIDs1to2, length))

    ## [ADVANCED USERS ONLY]
    ## An alternative that is slightly faster is to put all the probes
    ## (hgu95av2 + hgu133a) in a single PDict object and then query its
    ## 'dups0' slot directly. This slot is a Dups object containing the
    ## mapping between duplicated patterns.
    ## Note that we can do this only because all the probes have the
    ## same length (25) and because we are doing exact matching:

    probes12 <- DNAStringSet(c(hgu95av2probe$sequence, hgu133aprobe$sequence))
    pdict12 <- PDict(probes12)
    dups0 <- pdict12@dups0

    mapProbeSetIDs1to2alt <- function(psID)
    {
      ii1 <- unique(togroup(dups0, which(hgu95av2probe$Probe.Set.Name == psID)))
      ii2 <- members(dups0, ii1) - length(probes1)
      ii2 <- ii2[ii2 >= 1L]
      unique(hgu133aprobe$Probe.Set.Name[ii2])
    }

    mapPSIDs1to2alt <- lapply(psIDs1, mapProbeSetIDs1to2alt)  # about 10 min.
    names(mapPSIDs1to2alt) <- psIDs1

    ## 'mapPSIDs1to2alt' and 'mapPSIDs1to2' contain the same mapping:
    stopifnot(identical(lapply(mapPSIDs1to2alt, sort), lapply(mapPSIDs1to2, sort)))
  }
}

\keyword{methods}
