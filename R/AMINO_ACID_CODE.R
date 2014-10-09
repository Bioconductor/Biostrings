### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Single-Letter Amino Acid Code.
###

AMINO_ACID_CODE <- c(
    A="Ala", # Alanine
    R="Arg", # Arginine
    N="Asn", # Asparagine
    D="Asp", # Aspartic Acid
    C="Cys", # Cysteine
    Q="Gln", # Glutamine
    E="Glu", # Glutamic Acid
    G="Gly", # Glycine
    H="His", # Histidine
    I="Ile", # Isoleucine
    L="Leu", # Leucine
    K="Lys", # Lysine
    M="Met", # Methionine
    F="Phe", # Phenylalanine
    P="Pro", # Proline
    S="Ser", # Serine
    T="Thr", # Threonine
    W="Trp", # Tryptophan
    Y="Tyr", # Tyrosine
    V="Val", # Valine

  ## 21st and 22nd proteinogenic amino acids. Not coded for directly in
  ## the Standard Genetic Code:
    U="Sec", # Selenocysteine
    O="Pyl", # Pyrrolysine

  ## Ambiguities:
    B="Asx", # Asparagine or Aspartic Acid
    J="Xle", # Leucine or Isoleucine
    Z="Glx", # Glutamine or Glutamic Acid
    X="Xaa"  # Any amino acid
)

.AMINO_ACID_CODE_names <- names(AMINO_ACID_CODE)

### Amino Acid alphabet ("*" is a translation stop)
AA_ALPHABET <- c(.AMINO_ACID_CODE_names, "*", "-", "+", ".")

AA_STANDARD <- .AMINO_ACID_CODE_names[1:20]

AA_PROTEINOGENIC <- .AMINO_ACID_CODE_names[1:22]

