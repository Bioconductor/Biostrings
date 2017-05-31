### =========================================================================
### The Standard Genetic Code and its known variants
### -------------------------------------------------------------------------
###
### File: R/GENETIC_CODE.R
### Last update: January 28, 2014
###
### Reference document (last updated April 30, 2013):
###
###   http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
###
### The "official names" of the various codes ("Standard", "SGC0",
### "Vertebrate Mitochondrial", "SGC1", etc..) and their ids (1, 2, etc...)
### were taken from the print-form ASN.1 version of the above document
### (version 3.9 at the time of this writting):
###
###   ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
###
### By historical convention the codons (which are 3-nucleotide mRNA
### subsequences) are represented using the DNA alphabet (i.e. using T
### instead of U). Another way to look at this is to consider that the codons
### are coming from the "coding DNA strand" (a.k.a. "sense DNA strand" or
### "non-template DNA strand").
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Standard Code
###
### Official names:
###   - Standard
###   - SGC0
### Official id: 1
###

GENETIC_CODE <- c(
    TTT="F",
    TTC="F",
    TTA="L",
    TTG="L",

    TCT="S",
    TCC="S",
    TCA="S",
    TCG="S",

    TAT="Y",
    TAC="Y",
    TAA="*",
    TAG="*",

    TGT="C",
    TGC="C",
    TGA="*",
    TGG="W",

    CTT="L",
    CTC="L",
    CTA="L",
    CTG="L",

    CCT="P",
    CCC="P",
    CCA="P",
    CCG="P",

    CAT="H",
    CAC="H",
    CAA="Q",
    CAG="Q",

    CGT="R",
    CGC="R",
    CGA="R",
    CGG="R",

    ATT="I",
    ATC="I",
    ATA="I",
    ATG="M",

    ACT="T",
    ACC="T",
    ACA="T",
    ACG="T",

    AAT="N",
    AAC="N",
    AAA="K",
    AAG="K",

    AGT="S",
    AGC="S",
    AGA="R",
    AGG="R",

    GTT="V",
    GTC="V",
    GTA="V",
    GTG="V",

    GCT="A",
    GCC="A",
    GCA="A",
    GCG="A",

    GAT="D",
    GAC="D",
    GAA="E",
    GAG="E",

    GGT="G",
    GGC="G",
    GGA="G",
    GGG="G"
)

attr(GENETIC_CODE, "alt_init_codons") <- c("TTG", "CTG")

### From the "RNA transcript".
RNA_GENETIC_CODE <- GENETIC_CODE
names(RNA_GENETIC_CODE) <-
    chartr("T", "U", names(GENETIC_CODE))
attr(RNA_GENETIC_CODE, "alt_init_codons") <-
    chartr("T", "U", attr(GENETIC_CODE, "alt_init_codons"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All the known genetic codes
###

### All the known codes taken from
###   ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
### on May 25, 2017 (Version 4.0), and stored in one big list.
.genetic_code_table <- list(

    list(name   = "Standard",
         name2  = "SGC0",
         id     = 1,
         AAs    = c("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("---M---------------M------------",
                    "---M----------------------------")),

    list(name   = "Vertebrate Mitochondrial",
         name2  = "SGC1",
         id     = 2,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "MMMM---------------M------------")),

    list(name   = "Yeast Mitochondrial",
         name2  = "SGC2",
         id     = 3,
         AAs    = c("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRR",
                    "IIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "--MM----------------------------")),

    list(name   = paste0(c("Mold Mitochondrial",
                           "Protozoan Mitochondrial",
                           "Coelenterate Mitochondrial",
                           "Mycoplasma",
                           "Spiroplasma"), collapse="; "),
         name2  = "SGC3",
         id     = 4,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--MM---------------M------------",
                    "MMMM---------------M------------")),

    list(name   = "Invertebrate Mitochondrial",
         name2  = "SGC4",
         id     = 5,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"),
         Starts = c("---M----------------------------",
                    "MMMM---------------M------------")),

    list(name   = "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
         name2  = "SGC5",
         id     = 6,
         AAs    = c("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Echinoderm Mitochondrial; Flatworm Mitochondrial",
         name2  = "SGC8",
         id     = 9,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M---------------M------------")),

    list(name   = "Euplotid Nuclear",
         name2  = "SGC9",
         id     = 10,
         AAs    = c("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Bacterial, Archaeal and Plant Plastid",
         name2  = NA,
         id     = 11,
         AAs    = c("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("---M---------------M------------",
                    "MMMM---------------M------------")),

    list(name   = "Alternative Yeast Nuclear",
         name2  = NA,
         id     = 12,
         AAs    = c("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("-------------------M------------",
                    "---M----------------------------")),

    list(name   = "Ascidian Mitochondrial",
         name2  = NA,
         id     = 13,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"),
         Starts = c("---M----------------------------",
                    "--MM---------------M------------")),

    list(name   = "Alternative Flatworm Mitochondrial",
         name2  = NA,
         id     = 14,
         AAs    = c("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Blepharisma Macronuclear",
         name2  = NA,
         id     = 15,
         AAs    = c("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Chlorophycean Mitochondrial",
         name2  = NA,
         id     = 16,
         AAs    = c("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Trematode Mitochondrial",
         name2  = NA,
         id     = 21,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M---------------M------------")),

    list(name   = "Scenedesmus obliquus Mitochondrial",
         name2  = NA,
         id     = 22,
         AAs    = c("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "---M----------------------------")),

    list(name   = "Thraustochytrium Mitochondrial",
         name2  = NA,
         id     = 23,
         AAs    = c("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("--------------------------------",
                    "M--M---------------M------------")),

    list(name   = "Pterobranchia Mitochondrial",
         name2  = NA,
         id     = 24,
         AAs    = c("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
         Starts = c("---M---------------M------------",
                    "---M---------------M------------")),

    list(name   = "Candidate Division SR1 and Gracilibacteria",
         name2  = NA,
         id     = 25,
         AAs    = c("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("---M----------------------------",
                    "---M---------------M------------")),

    list(name   = "Pachysolen tannophilus Nuclear",
         name2  = NA,
         id     = 26,
         AAs    = c("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRR",
                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
         Starts = c("-------------------M------------",
                    "---M----------------------------"))

### The following genetic codes are listed on
###   https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
### but not on
###   ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
### They pose the problem that they use stop codons that are not always
### translated as such. For example, in code 27 below, TGA can mean STOP
### or be translated as Trp but I don't see any information on the
### https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi page about
### when one or the other occurs. So I've no idea how could one use these
### codes to predict translation. More precisely, translate() needs to be
### deterministic and I don't see how it could with these genetic codes.
 
#    list(name   = "Karyorelict Nuclear",
#         name2  = NA,
#         id     = 27,
#         AAs    = c("FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRR",
#                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
#         Starts = c("--------------------------------",
#                    "---M----------------------------")),

#    list(name   = "Condylostoma Nuclear",
#         name2  = NA,
#         id     = 28,
#         AAs    = c("FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRR",
#                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
#         Starts = c("--------------------------------",
#                    "---M----------------------------")),

#    list(name   = "Mesodinium Nuclear",
#         name2  = NA,
#         id     = 29,
#         AAs    = c("FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRR",
#                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
#         Starts = c("--------------------------------",
#                    "---M----------------------------")),

#    list(name   = "Peritrich Nuclear",
#         name2  = NA,
#         id     = 30,
#         AAs    = c("FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRR",
#                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
#         Starts = c("--------------------------------",
#                    "---M----------------------------")),

#    list(name   = "Blastocrithidia Nuclear",
#         name2  = NA,
#         id     = 31,
#         AAs    = c("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRR",
#                    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
#         Starts = c("--------------------------------",
#                    "---M----------------------------"))
)

### All the known codes summarized in one big data.frame:
GENETIC_CODE_TABLE <- do.call(rbind,
    lapply(.genetic_code_table,
        function(genetic_code) {
            data.frame(name=genetic_code$name,
                       name2=as.character(genetic_code$name2),
                       id=as.character(genetic_code$id),
                       AAs=paste0(genetic_code$AAs, collapse=""),
                       Starts=paste0(genetic_code$Starts, collapse=""),
                       stringsAsFactors=FALSE)
        }))

.solve_id_or_name2 <- function(id_or_name2, full.search=FALSE)
{
    if (!isSingleString(id_or_name2))
        stop("'id_or_name2' must be a single string")
    if (id_or_name2 == "")
        stop("'id_or_name2' must be a non empty string")
    if (!isTRUEorFALSE(full.search))
        stop("'full.search' must be TRUE or FALSE")
    ## Search the 'id' col.
    idx <- match(id_or_name2, GENETIC_CODE_TABLE$id)
    if (!is.na(idx))
        return(idx)
    ## Search the 'name2' col.
    idx <- match(id_or_name2, GENETIC_CODE_TABLE$name2)
    if (!is.na(idx))
        return(idx)
    if (!full.search)
        stop("Genetic code not found for 'id_or_name2=\"",
             id_or_name2, "\"'.\n  ",
             "Maybe try again with 'full.search=TRUE'.")
    ## Search the 'name' col.
    idx <- grep(id_or_name2, GENETIC_CODE_TABLE$name, ignore.case=TRUE)
    if (length(idx) == 0L)
        stop("Genetic code not found for 'id_or_name2=\"",
             id_or_name2, "\"'.")
    if (length(idx) >= 2L)
        stop("More than one genetic code found for 'id_or_name2=\"",
             id_or_name2, "\"'.\n  ",
             "Try to be more specific.")
    idx
}

getGeneticCode <- function(id_or_name2="1", full.search=FALSE,
                           as.data.frame=FALSE)
{
    if (!isTRUEorFALSE(as.data.frame))
        stop("'as.data.frame' must be TRUE or FALSE")
    idx <- .solve_id_or_name2(id_or_name2, full.search=full.search)
    ans <- safeExplode(GENETIC_CODE_TABLE[idx, "AAs"])
    Start <- safeExplode(GENETIC_CODE_TABLE[idx, "Starts"])
    if (as.data.frame) {
        ans <- data.frame(AA=ans, Start=Start, stringsAsFactors=FALSE)
        rownames(ans) <- names(GENETIC_CODE)
        return(ans)
    }
    names(ans) <- names(GENETIC_CODE)
    idx <- which(Start != "-" & Start != ans)
    stopifnot(all(Start[idx] == "M"))  # should never happen
    attr(ans, "alt_init_codons") <- names(ans)[idx]
    ans
}

