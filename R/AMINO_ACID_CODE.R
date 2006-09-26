# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The Single-Letter Amino Acid Code

AMINO_ACID_CODE <- c(
    "A"="Ala", # Alanine
    "R"="Arg", # Arginine
    "N"="Asn", # Asparagine
    "D"="Asp", # Aspartic Acid
    "C"="Cys", # Cysteine
    "Q"="Gln", # Glutamine
    "E"="Glu", # Glutamic Acid
    "G"="Gly", # Glycine
    "H"="His", # Histidine
    "I"="Ile", # Isoleucine
    "L"="Leu", # Leucine
    "K"="Lys", # Lysine
    "M"="Met", # Methionine
    "F"="Phe", # Phenylalanine
    "P"="Pro", # Proline
    "S"="Ser", # Serine
    "T"="Thr", # Threonine
    "W"="Trp", # Tryptophan
    "Y"="Tyr", # Tyrosine
    "V"="Val"  # Valine
)

# Amino Acid alphabet ("*" is a translation stop)
AA_ALPHABET <- c(names(AMINO_ACID_CODE), "*")