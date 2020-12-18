test_writeXStringSet <- function()
{
    tmpdir <- tempdir()
    tmpfilefa <- file.path(tmpdir,"whee.fa")
    tmpfilefq <- file.path(tmpdir,"whee.fq")
    x1 <- DNAStringSet(c("TTGA", "CTCN"))
    q1 <- PhredQuality(c("*+,-", "6789"))
    qdna1 <- QualityScaledDNAStringSet(x1, q1)
    names_dna <- c("A", "B")
    names(qdna1) <- names_dna

    writeXStringSet(x1, format="fasta", file=tmpfilefa)
    dna1 <- readDNAStringSet(tmpfilefa, format="fasta")
    checkTrue(is.null(mcols(dna1)))
    checkTrue(all(names(dna1) == c("","")))

    names(x1) <- names_dna
    writeXStringSet(x1, format="fastq", file=tmpfilefq)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    checkTrue(!all(quality(qdna2) == quality(qdna1)))
    checkTrue(all(all(as(quality(qdna2),"IntegerList") == 26L)))

    writeXStringSet(x1, format="fastq", file=tmpfilefq, qualities = q1)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    checkTrue(all(quality(qdna2) == quality(qdna1)))

    writeXStringSet(qdna1, format="fastq", file=tmpfilefq)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    checkTrue(all(names(qdna2) == names_dna))
    checkTrue(all(quality(qdna2) == quality(qdna1)))

    qdna2 <- readDNAStringSet(tmpfilefq, format="fastq")
    checkTrue(is.null(mcols(qdna2)))

    qdna2 <- readDNAStringSet(tmpfilefq, format="fastq", with.qualities = TRUE)
    checkTrue(all(mcols(qdna2)$qualities == quality(qdna1)))

    unlink(tmpfilefa)
    unlink(tmpfilefq)
}
