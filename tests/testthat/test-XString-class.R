dnastr <- paste(DNA_ALPHABET, collapse="")
rnastr <- paste(RNA_ALPHABET, collapse="")
aastr <- paste(AA_ALPHABET, collapse="")
bstr <- rawToChar(as.raw(32:126))

test_allstrings_for_error <- function(input, exp_error, ignore="") {
    if (!("DNA" %in% ignore)) expect_error2(DNAString(input), exp_error)
    if (!("RNA" %in% ignore)) expect_error2(RNAString(input), exp_error)
    if (!("AA" %in% ignore))  expect_error2(AAString(input), exp_error)
    if (!("B" %in% ignore))   expect_error2(BString(input), exp_error)
}

test_that("seqtype() correctly infers types", {
    expect_equal(seqtype(BString("ABC")), "B")
    expect_equal(seqtype(RNAString("AUG")), "RNA")
    expect_equal(seqtype(DNAString("ATG")), "DNA")
    expect_equal(seqtype(AAString("ARN")), "AA")
})


test_that("character, vector conversion works properly", {
    expect_equal(as.character(DNAString(dnastr)), dnastr)
    expect_equal(as.character(RNAString(rnastr)), rnastr)
    expect_equal(as.character(AAString(aastr)), aastr)
    expect_equal(as.character(BString(bstr)), bstr)
})

test_that("encode/decode tables work correctly", {
    expect_equal(DNAString(tolower(dnastr)), DNAString(dnastr))
    expect_equal(RNAString(tolower(rnastr)), RNAString(rnastr))
    expect_equal(AAString(tolower(aastr)), AAString(aastr))

    expect_equal(DNAString(as.factor(dnastr)), DNAString(dnastr))
    expect_equal(RNAString(as.factor(rnastr)), RNAString(rnastr))
    expect_equal(AAString(as.factor(aastr)), AAString(aastr))
    expect_equal(BString(as.factor(bstr)), BString(bstr))

    test_allstrings_for_error(bstr, "not in lookup table", ignore="B")
})

test_that("letter() works as expected", {
    expect_equal(letter(DNAString("ATGCATGC"), c(1,3,6)), "AGT")
    expect_equal(letter(RNAString("AUGCAUGC"), c(1,3,6)), "AGU")
    expect_equal( letter(AAString("ARNDARND"), c(1,3,6)), "ANR")
    expect_equal(  letter(BString("ABCDEFGH"), c(1,3,6)), "ACF")
    expect_equal(letter(bstr, seq(nchar(bstr), 1)),
                 rawToChar(rev(charToRaw(bstr))))

    ## TODO: Better error messages
    expect_error2(letter(DNAString(""), 10), "out of bounds")
    expect_error2(letter(RNAString(""), 10), "out of bounds")
    expect_error2(letter(AAString(""), 10), "out of bounds")
    expect_error2(letter(BString(""), 10), "out of bounds")
    expect_error2(letter("", 10), "out of bounds")
})

test_that("constructors handle invalid and non-char input correctly", {
    # note that this only displays for character input
    # meaning DNAString(1) returns a non-informative error
    test_allstrings_for_error(NA_character_,
                              "input must be a single non-NA string")
    test_allstrings_for_error(as.factor(1:4),
                              "input must be a single non-NA string")
    test_allstrings_for_error(c("A", "A"),
                              "input must be a single non-NA string")
})

test_that("conversion between XString seqtypes works properly", {
    # DNA <-> RNA
    expect_equal(RNAString(DNAString(dnastr)), RNAString(rnastr))
    expect_equal(DNAString(RNAString(rnastr)), DNAString(dnastr))

    # X -> B
    expect_equal(BString(DNAString(dnastr)), BString(dnastr))
    expect_equal(BString(RNAString(rnastr)), BString(rnastr))
    expect_equal(BString(AAString(aastr)), BString(aastr))

    # valid B -> X
    expect_equal(DNAString(BString(dnastr)), DNAString(dnastr))
    expect_equal(RNAString(BString(rnastr)), RNAString(rnastr))
    expect_equal(AAString(BString(aastr)), AAString(aastr))

    # invalid DNA,RNA <-> AA
    expect_error2(AAString(DNAString("ATGC")), "incompatible sequence types")
    expect_error2(AAString(RNAString("AUGC")), "incompatible sequence types")
    expect_error2(DNAString(AAString("ARND")), "incompatible sequence types")
    expect_error2(RNAString(AAString("ARND")), "incompatible sequence types")

    # invalid B -> X
    bbad <- BString(";")
    test_allstrings_for_error(bbad, "not in lookup table", ignore="B")
})

test_that("as.vector() methods work as expected", {
    ## TODO: sometimes returns factor, sometimes character
    ##                 note the BString case as an example
    dnafac <- factor(seq_len(nchar(dnastr)))
    attr(dnafac, "levels") <- strsplit(dnastr, "")[[1]]

    rnafac <- factor(seq_len(nchar(rnastr)))
    attr(rnafac, "levels") <- strsplit(rnastr, "")[[1]]

    aafac <- factor(seq_len(nchar(aastr)))
    attr(aafac, "levels") <- strsplit(aastr, "")[[1]]

    bfac <- factor(seq_len(nchar(bstr)))
    attr(bfac, "levels") <- strsplit(bstr, "")[[1]]

    expect_equal(as.vector(DNAString(dnastr)), dnafac)
    expect_equal(as.vector(RNAString(rnastr)), rnafac)
    expect_equal(as.vector(AAString(aastr)), aafac)
    #expect_equal(as.vector(BString(bstr)), bfac)
})

test_that("methods from the Compare group generic work as expected", {
    check_Compare_methods <- function(x, y, compare_code) {
        expect_identical(x == y, compare_code == 0L)
        expect_identical(x != y, compare_code != 0L)
        expect_identical(x <= y, compare_code <= 0L)
        expect_identical(x >= y, compare_code >= 0L)
        expect_identical(x <  y, compare_code <  0L)
        expect_identical(x >  y, compare_code >  0L)
    }

    pcompare_character <- function(x, y) {
        stopifnot(is.character(x), is.character(y))
        prev_locale <- Sys.getlocale("LC_COLLATE")
        Sys.setlocale("LC_COLLATE", "C")
        on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
        ifelse(x < y, -1L, ifelse(x > y, 1L, 0L))
    }

    BString_objects <- list(BString(), BString("AAA"), BString("aaa"),
                            BString("aaz"), BString("aaaa"),
                            BString("TGM"), BString("tgg"))
    ## Note that the objects in 'DNAString_objects' and 'RNAString_objects'
    ## are sorted in strictly increasing order. Yes, because of how the
    ## letters in an DNAString or RNAString object are encoded,
    ## DNAString("TGM") is considered < DNAString("TGG").
    DNAString_objects <- list(DNAString(), DNAString("TGM"), DNAString("TGG"))
    RNAString_objects <- list(RNAString(), RNAString("UGM"), RNAString("UGG"))
    AAString_objects <- list(AAString(), AAString("AAA"), AAString("TGM"))
    char_vecs <- c("", "aaa", "aaz", "aaaa", "tgg", "TGM")

    ## Between two BString objects.
    for (x in BString_objects) {
        for (y in BString_objects) {
            compare_code <- pcompare_character(as.character(x), as.character(y))
            check_Compare_methods(x, y, compare_code)
        }
    }

    ## Between a BString object and a DNAString object.
    for (x in BString_objects) {
        for (y in DNAString_objects) {
            compare_code <- pcompare_character(as.character(x), as.character(y))
            check_Compare_methods(x, y, compare_code)
            check_Compare_methods(y, x, -compare_code)
        }
    }

    ## Between a BString object and an RNAString object.
    for (x in BString_objects) {
        for (y in RNAString_objects) {
            compare_code <- pcompare_character(as.character(x), as.character(y))
            check_Compare_methods(x, y, compare_code)
            check_Compare_methods(y, x, -compare_code)
        }
    }

    ## Between a BString object and a AAString object.
    for (x in BString_objects) {
        for (y in AAString_objects) {
            compare_code <- pcompare_character(as.character(x), as.character(y))
            check_Compare_methods(x, y, compare_code)
            check_Compare_methods(y, x, -compare_code)
        }
    }

    ## Between a BString object and a character vector.
    for (x in BString_objects) {
        compare_codes <- pcompare_character(as.character(x), char_vecs)
        check_Compare_methods(x, char_vecs, compare_codes)
        check_Compare_methods(char_vecs, x, -compare_codes)
    }

    ## Between two DNAString objects.
    for (i in seq_along(DNAString_objects)) {
        x <- DNAString_objects[[i]]
        for (j in seq_along(DNAString_objects)) {
            y <- DNAString_objects[[j]]
            compare_code <- as.integer(sign(i - j))
            check_Compare_methods(x, y, compare_code)
        }
    }

    ## Between a DNAString object and an RNAString object.
    for (i in seq_along(DNAString_objects)) {
        x <- DNAString_objects[[i]]
        for (j in seq_along(RNAString_objects)) {
            y <- RNAString_objects[[j]]
            compare_code <- as.integer(sign(i - j))
            check_Compare_methods(x, y, compare_code)
            check_Compare_methods(y, x, -compare_code)
        }
    }

    ## Between a DNAString object and an AAString object.
    msg <- paste0("comparison between a DNAString object ",
                  "and a AAString object is not supported")
    expect_error2(DNAString() == AAString(), msg)
    expect_error2(DNAString() != AAString(), msg)
    expect_error2(DNAString() <= AAString(), msg)
    expect_error2(DNAString() >= AAString(), msg)
    expect_error2(DNAString() <  AAString(), msg)
    expect_error2(DNAString() >  AAString(), msg)
    msg <- paste0("comparison between a AAString object ",
                  "and a DNAString object is not supported")
    expect_error2(AAString() == DNAString(), msg)
    expect_error2(AAString() != DNAString(), msg)
    expect_error2(AAString() <= DNAString(), msg)
    expect_error2(AAString() >= DNAString(), msg)
    expect_error2(AAString() <  DNAString(), msg)
    expect_error2(AAString() >  DNAString(), msg)

    ## Between a DNAString object and a character vector.
    for (x in DNAString_objects) {
        compare_codes <- pcompare_character(as.character(x), char_vecs)
        check_Compare_methods(x, char_vecs, compare_codes)
        check_Compare_methods(char_vecs, x, -compare_codes)
    }

    ## Between two RNAString objects.
    for (i in seq_along(RNAString_objects)) {
        x <- RNAString_objects[[i]]
        for (j in seq_along(RNAString_objects)) {
            y <- RNAString_objects[[j]]
            compare_code <- as.integer(sign(i - j))
            check_Compare_methods(x, y, compare_code)
        }
    }

    ## Between an RNAString object and an AAString object.
    msg <- paste0("comparison between a RNAString object ",
                  "and a AAString object is not supported")
    expect_error2(RNAString() == AAString(), msg)
    expect_error2(RNAString() != AAString(), msg)
    expect_error2(RNAString() <= AAString(), msg)
    expect_error2(RNAString() >= AAString(), msg)
    expect_error2(RNAString() <  AAString(), msg)
    expect_error2(RNAString() >  AAString(), msg)
    msg <- paste0("comparison between a AAString object ",
                  "and a RNAString object is not supported")
    expect_error2(AAString() == RNAString(), msg)
    expect_error2(AAString() != RNAString(), msg)
    expect_error2(AAString() <= RNAString(), msg)
    expect_error2(AAString() >= RNAString(), msg)
    expect_error2(AAString() <  RNAString(), msg)
    expect_error2(AAString() >  RNAString(), msg)

    ## Between an RNAString object and a character vector.
    for (x in RNAString_objects) {
        compare_codes <- pcompare_character(as.character(x), char_vecs)
        check_Compare_methods(x, char_vecs, compare_codes)
        check_Compare_methods(char_vecs, x, -compare_codes)
    }

    ## Between two AAString objects.
    for (x in AAString_objects) {
        for (y in AAString_objects) {
            compare_code <- pcompare_character(as.character(x), as.character(y))
            check_Compare_methods(x, y, compare_code)
        }
    }

    ## Between an AAString object and a character vector.
    for (x in AAString_objects) {
        compare_codes <- pcompare_character(as.character(x), char_vecs)
        check_Compare_methods(x, char_vecs, compare_codes)
        check_Compare_methods(char_vecs, x, -compare_codes)
    }
})

test_that("output works correctly", {
    s <- paste(rep("A", 100L), collapse="")
    expect_output(show(DNAString(s)),
        "100-letter DNAString object\\nseq: .+\\.\\.\\..+$", width=80)
    expect_output(show(RNAString(s)),
        "100-letter RNAString object\\nseq: .+\\.\\.\\..+$", width=80)
    expect_output(show(AAString(s)),
        "100-letter AAString object\\nseq: .+\\.\\.\\..+$", width=80)
    expect_output(show(BString(s)),
        "100-letter BString object\\nseq: A{36}\\.\\.\\.A{36}$", width=80)

    # width of sequence is 5 less than that of 'width'
    expect_output(show(DNAString(s)),
        "100-letter DNAString object\\nseq: .+\\.\\.\\..+$", width=10)
    expect_output(show(RNAString(s)),
        "100-letter RNAString object\\nseq: .+\\.\\.\\..+$", width=10)
    expect_output(show(AAString(s)),
        "100-letter AAString object\\nseq: .+\\.\\.\\..+$", width=10)
    expect_output(show(BString(s)),
        "100-letter BString object\\nseq: AA\\.\\.\\.AA$", width=10)
})

test_that("substr() and substring() methods work as expected", {
    d <- DNAString(dnastr)
    r <- RNAString(rnastr)
    a <- AAString(aastr)
    b <- BString(bstr)

    expect_equal(as.character(substr(d, 1, 10)), substr(dnastr, 1, 10))
    expect_equal(as.character(substr(r, 1, 10)), substr(rnastr, 1, 10))
    expect_equal(as.character(substr(a, 1, 10)), substr(aastr, 1, 10))
    expect_equal(as.character(substr(b, 1, 10)), substr(bstr, 1, 10))

    expect_equal(as.character(substring(d, 5, 10)), substring(dnastr, 5, 10))
    expect_equal(as.character(substring(r, 5, 10)), substring(rnastr, 5, 10))
    expect_equal(as.character(substring(a, 5, 10)), substring(aastr, 5, 10))
    expect_equal(as.character(substring(b, 5, 10)), substring(bstr, 5, 10))

    expect_error2(substring(d, 10, 5), "Invalid sequence coordinates")
    expect_error2(substring(r, 10, 5), "Invalid sequence coordinates")
    expect_error2(substring(a, 10, 5), "Invalid sequence coordinates")
    expect_error2(substring(b, 10, 5), "Invalid sequence coordinates")

    # `[` dispatch
    expect_equal(as.character(d[1:10]), substr(dnastr, 1, 10))
    expect_equal(as.character(d[-1]), substr(dnastr, 2, nchar(dnastr)))
})

test_that("reverse(), complement(), reverseComplement() work as expected", {
    ## reverse tests
    .revString <- function(s) paste(rev(safeExplode(s)), collapse="")
    dna <- DNAString(dnastr)
    rna <- RNAString(rnastr)
    aaa <- AAString(aastr)
    bbb <- BString(bstr)
    d_comp <- "TGCAKYWSRMBDHVN-+."
    r_comp <- "UGCAKYWSRMBDHVN-+."

    ## example Views on strings
    d_v <- Views(dna, start=rep(1L,3L), end=rep(nchar(dnastr),3L))
    mdna <- dna
    mrna <- rna
    m1 <- Mask(nchar(dnastr), start=c(2,7), end=c(3,10))
    masks(mdna) <- masks(mrna) <- m1
    md_comp <- strsplit(d_comp, "")[[1]]
    mr_comp <- strsplit(r_comp, "")[[1]]
    md_comp[c(2:3,7:10)] <- mr_comp[c(2:3,7:10)] <- "#"
    md_comp <- paste(md_comp, collapse="")
    mr_comp <- paste(mr_comp, collapse="")


    ## reverse method
    expect_equal(reverse(dnastr), .revString(dnastr))
    expect_equal(as.character(reverse(dna)), .revString(dnastr))
    expect_equal(as.character(reverse(rna)), .revString(rnastr))
    expect_equal(as.character(reverse(aaa)), .revString(aastr))
    expect_equal(as.character(reverse(bbb)), .revString(bstr))
    expect_equal(as.character(reverse(mdna)), .revString(as.character(mdna)))
    expect_equal(as.character(reverse(mrna)), .revString(as.character(mrna)))

    ## complement method
    expect_equal(as.character(complement(dna)), d_comp)
    expect_equal(as.character(complement(rna)), r_comp)
    expect_error2(complement(AAString()), "unable to find an inherited method")
    expect_error2(complement(BString()), "unable to find an inherited method")
    expect_equal(as.character(complement(d_v)), rep(d_comp, 3L))
    expect_equal(as.character(complement(mdna)), md_comp)
    expect_equal(as.character(complement(mrna)), mr_comp)

    ## reverseComplement method
    expect_equal(as.character(reverseComplement(dna)), .revString(d_comp))
    expect_equal(as.character(reverseComplement(rna)), .revString(r_comp))
    expect_error2(reverseComplement(AAString()),
                  "unable to find an inherited method")
    expect_error2(reverseComplement(BString()),
                  "unable to find an inherited method")
    expect_equal(as.character(reverseComplement(d_v)),
                 rep(.revString(d_comp), 3L))
    expect_equal(as.character(reverseComplement(mdna)), .revString(md_comp))
    expect_equal(as.character(reverseComplement(mrna)), .revString(mr_comp))
})

## Porting RUnit tests
test_that("alphabet() finds the correct values", {
    expect_equal(alphabet(DNAString(dnastr)), strsplit(dnastr, "")[[1]])
    expect_equal(alphabet(RNAString(rnastr)), strsplit(rnastr, "")[[1]])
    expect_equal(alphabet(AAString(aastr)), strsplit(aastr, "")[[1]])
    expect_equal(alphabet(BString(bstr)), NULL)

    expect_equal(alphabet(DNAString(), baseOnly=TRUE), DNA_BASES)
    expect_equal(alphabet(RNAString(), baseOnly=TRUE), RNA_BASES)
    expect_error2(alphabet(DNAString(), baseOnly=1),
                  "'baseOnly' must be TRUE or FALSE")
})

