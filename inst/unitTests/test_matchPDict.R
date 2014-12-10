

# --- to be included in an upcoming biocDatasets package
randomDNASequences <- function(n, w)
{
  alphabet <- DNA_BASES
  w <- rep(w, length=n)
  sequences <- sapply(seq(1, n, length=n),
                      function(x) {
                        s <- sample(alphabet, w[x], replace=TRUE)
                        s <- paste(s, collapse="")
                        return(s)
                      })
  return(Biostrings::DNAStringSet(sequences))
}

msubseq <- function(x, ir)
{
  ## differs from subseq in the sense that several subsequences
  ## from the same sequence are extracted
  ## x:  XString
  ## ir: IRanges
  res <- vector("character", length = length(ir))
  for (i in seq(along=res)) {
    res[i] <- as.character(subseq(x, start=ir@start[i], width=width(ir)[i]))
    ## forced cast: chain of tools for DNAString seems interupted for
    ##              some use cases (or I missed something)
  }
  res <- DNAStringSet(res)
  return(res)
}


# ---

test_pdictConstantWidth <- function()
{
  set.seed(1)
  l <- 150
  target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6 
  ir <- successiveIRanges(rep(W, L), gapwidth = 1)
  short_sequences <- msubseq(target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))
  
  dna_short <- DNAStringSet(short_sequences[o])
  pdict <- PDict(dna_short)
  checkEquals(L, length(pdict))
  checkEquals(rep(W, L), width(pdict))
  checkEquals(NULL, head(pdict))
  checkEquals(W, tb.width(pdict))
  checkEquals(NULL, tail(pdict))
}

test_pdictVariableWidth <- function()
{
  set.seed(1)
  l <- 150
  target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  n_cut <- sample(0:5, L, replace=TRUE)
  ir <- successiveIRanges(rep(W, L) - n_cut, gapwidth = 1 + n_cut[-length(n_cut)])
  short_sequences <- msubseq(target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))
  
  dna_var_short <- DNAStringSet(short_sequences[o])
  
  pdict <- PDict(dna_var_short,
                 tb.start=1,                        # can't this be
                 tb.width=min(width(short_sequences)) # the default for
                                                    # variable width ?
                 )
  checkEquals(L, length(pdict))
  checkEquals( (rep(W, L) - n_cut)[o], width(pdict))
  checkEquals(NULL, head(pdict))
  shortest_seq_width <- min(width(dna_var_short))
  checkEquals(shortest_seq_width,
              tb.width(pdict))           # mostly a sanity check
  checkEquals(substring(short_sequences, shortest_seq_width+1)[o],
              as.character(tail(pdict)))
}


test_matchConstantWidth <- function()
{
  set.seed(1)
  l <- 150
  dna_target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6 
  ir <- successiveIRanges(rep(W, L), gapwidth = 1)
  short_sequences <- msubseq(dna_target, ir)
  # shuffle the sequences (so they are not in consecutive order)
  o <- sample(seq(along=short_sequences))
  
  dna_short <- DNAStringSet(short_sequences[o])
  pdict <- PDict(dna_short)
  
  res <- matchPDict(pdict, dna_target)

  # mostly a sanity check
  checkEquals(L, length(res))
  
  for (i in seq(along=res)) {
    m_start <- ir[o][i]@start
    checkEquals(m_start, start(res[[i]]))
    checkEquals(W, width(res[[i]]))
    checkEquals(m_start + W - 1, end(res[[i]]))  # mostly a sanity check
  }
  
  
}

test_matchVariableWidth <- function()
{
  set.seed(1)
  l <- 150
  dna_target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  n_cut <- sample(0:5, L, replace=TRUE)
  ir <- successiveIRanges(rep(W, L) - n_cut, gapwidth = 1 + n_cut[-length(n_cut)])
  short_sequences <- msubseq(dna_target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))
  
  dna_var_short <- DNAStringSet(short_sequences[o])
  
  pdict <- PDict(dna_var_short,
                 tb.start=1,                        # can't this be
                 tb.width=min(width(dna_var_short)) # the default for
                                                    # variable width ?
                 )

  res <- matchPDict(pdict, dna_target)

  # mostly a sanity check
  checkEquals(L, length(res))
  
  iro <- ir[o]
  for (i in seq(along=res)) {
    checkEquals(start(iro)[i], start(res[[i]]))
    checkEquals(width(iro)[i], width(res[[i]]))
    checkEquals(end(iro)[i], end(res[[i]]))  # mostly a sanity check
  }
    
}

