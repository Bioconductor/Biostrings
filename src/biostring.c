#include <limits.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef MAX
#undef MAX
#endif

#define MAX(i, j) ((i)>(j)?(i):(j))

#define R_assert(e) ((e) ? (void) 0 : error("assertion `%s' failed: file `%s', line %d\n", #e, __FILE__, __LINE__))

SEXP
IntegerBitOr(SEXP x)
{
    unsigned int ans = 0U;
    unsigned int* vec;
    int i, n;
    if (TYPEOF(x) != INTSXP)
        error("bitwise or can be done only for integers");
    vec = (unsigned int*) INTEGER(x);
    n = LENGTH(x);
    for (i = 0; i < n; i++)
        ans |= vec[i];
    return ScalarInteger((int)ans);
}

static int
isFromClass(SEXP x, const char* klass)
{
    SEXP call = PROTECT(lang3(install("is"),
                              x, mkString((char *) klass)));
    SEXP ans = eval(call, R_GlobalEnv);
    UNPROTECT(1);
    return asLogical(ans);
}

SEXP
BioStringValues(SEXP alphabet_length, SEXP string_length)
{
    int alphlen = asInteger(alphabet_length);
    int n = asInteger(string_length)+1; /* one extra value at begining */
    SEXP ans, storage;

    if (alphlen <= CHAR_BIT) {
        storage = allocString(n);
        memset(CHAR(storage), 0, n+1);
    } else if (alphlen <= sizeof(int)) {
        storage = allocVector(INTSXP, n);
        memset(INTEGER(storage), 0, n*sizeof(int));
    } else {
        error("unable to create string for alphabet with %d letters", alphlen);
        storage = R_NilValue; /* -Wall */
    }
    PROTECT(storage);
    ans = R_MakeExternalPtr(NULL, storage, R_NilValue);
    UNPROTECT(1);
    return ans;
}

static void
reserveBioString(SEXP values, unsigned long newlength)
{
    SEXP vec = R_ExternalPtrTag(values);
    int oldlength = LENGTH(vec);
    newlength++; /* one extra value at begining */
    if (oldlength < newlength) {
        PROTECT(values);
        if (TYPEOF(vec) == CHARSXP) {
            /* lengthgets does not work with CHARSXP */
            SEXP newvec = allocString(newlength);
            memcpy(CHAR(newvec), CHAR(vec), oldlength);
            vec = newvec;
        } else {
            vec = lengthgets(vec, newlength);
        }
        R_SetExternalPtrTag(values, vec);
        UNPROTECT(1);
    }
}

/*
 *  Add characters from element srcindex in src after the
 *  current_length element in values after encoding the characters.
 *  The encoding is done for the alphabet in letters.
 *
 *  This may increase the length of the character or integer vector
 *  inside values. The (possibly modified) values object is returned.
 */
static int
appendCharacterToBioString(SEXP alphMapping,
                           SEXP values,
                           unsigned long current_length,
                           SEXP src,
                           int srcindex)
{
    int i;
    unsigned long count, nletters;
    unsigned int mapping[256];
    unsigned int maxcode = 0;
    SEXP letters = GET_NAMES(alphMapping);
    if (TYPEOF(src) != STRSXP)
        error("source is not a character vector");
    if (TYPEOF(alphMapping) != INTSXP || TYPEOF(letters) != STRSXP)
        error("invalid mapping");
    nletters = LENGTH(letters);
    if (LENGTH(alphMapping) != nletters)
        error("invalid names for mapping");
    if (srcindex < 0 || srcindex >=  LENGTH(src))
        error("source index out of bounds");
    if (nletters > sizeof(int)*CHAR_BIT)
        error("alphabet is too large");
    count = LENGTH(STRING_ELT(src, srcindex));
    PROTECT(src);
    PROTECT(alphMapping);
    reserveBioString(values, current_length+count);
    UNPROTECT(2);
    memset(mapping, 0, 256*sizeof(unsigned int));
    for (i = 0; i < nletters; i++) {
        SEXP c = STRING_ELT(letters, i);
        unsigned int tmpcode = ((unsigned int*) INTEGER(alphMapping))[i];
        if (LENGTH(c) != 1)
            error("invalid names for mapping");
        if (tmpcode > maxcode)
            maxcode = tmpcode;
        mapping[*((unsigned char*) CHAR(c))] = tmpcode;
    }
    if (maxcode < (1U << CHAR_BIT)) {
        /* skip one element in front */
        unsigned char* destptr =
            (unsigned char*) CHAR(R_ExternalPtrTag(values))+1+current_length;
        SEXP src_i = STRING_ELT(src, srcindex);
        unsigned char* srcptr = (unsigned char*) CHAR(src_i);
        int j;

        for (j = 0; j < count; j++) {
            unsigned int val = mapping[srcptr[j]];
            if (!val)
                error("invalid character `%c` in source at string %d, position %d",
                      srcptr[j], srcindex, j);
            destptr[j] = (unsigned char) val;
        }
    } else {
        unsigned int* destptr =
            (unsigned int*) INTEGER(R_ExternalPtrTag(values))+1+current_length;
        SEXP src_i = STRING_ELT(src, srcindex);
        unsigned char* srcptr = (unsigned char*) CHAR(src_i);
        int j;
        for (j = 0; j < count; j++) {
            unsigned int val = mapping[srcptr[j]];
            if (!val)
                error("invalid character `%c` in source at string %d, position %d",
                      srcptr[j], srcindex, j);
            destptr[j] = val;
        }
    }
    return current_length+count;
}

SEXP
setBioString(SEXP biostring, SEXP src)
{
    int start;
    int end = 0;
    if (!isFromClass(biostring, "BioString"))
        error("first argument must be from BioString class");
    if (asLogical(GET_SLOT(biostring, install("initialized"))))
        error("can not modify initialized strings");
    if (length(src) != 1)
        error("can only set BioString value from single character string");
    start = INTEGER(GET_SLOT(biostring, install("end")))[0];
    end = appendCharacterToBioString(GET_SLOT(GET_SLOT(biostring,
                                                       install("alphabet")),
                                              install("mapping")),
                                     GET_SLOT(biostring, install("values")),
                                     asInteger(GET_SLOT(biostring,
                                                        install("end"))),
                                     src,
                                     0);
    if (start == 0)
        start = 1;
    INTEGER(GET_SLOT(biostring, install("start")))[0] = start;
    INTEGER(GET_SLOT(biostring, install("end")))[0] = end;
    return biostring;
}

/* Convert a range in a BioString object to a character string */
static SEXP
BioStringToCharacter(int nletters_base,
                     SEXP alphMapping,
                     SEXP values,
                     unsigned long start,
                     unsigned long end)
{
    int nvalues = LENGTH(values)-1;
    int nletters;
    SEXP ans;
    char* dest;
    unsigned long ndest, i;
    SEXP letters = GET_NAMES(alphMapping);
    unsigned int bad_element = 0;
    unsigned long bad_index = 0;

    if (TYPEOF(alphMapping) != INTSXP || TYPEOF(letters) != STRSXP)
        error("invalid mapping");
    nletters = LENGTH(letters);
    if (LENGTH(alphMapping) != nletters)
        error("invalid names for mapping");
    if (end > nvalues)
        end = nvalues;
    if (start > end)
        return ScalarString(R_BlankString);
    if (start <= 0)
        start = 1;
    if (nletters_base > sizeof(int)*CHAR_BIT)
        error("alphabet is too large");
    values = R_ExternalPtrTag(values);

    ndest = end-start+1;
    ans = allocVector(STRSXP, 1);
    PROTECT(ans);
    SET_STRING_ELT(ans, 0, allocString(ndest));
    dest = CHAR(STRING_ELT(ans, 0));

    if (nletters_base <= CHAR_BIT) {
        unsigned char* src;
        char mapping[256];
        if (TYPEOF(values) != CHARSXP)
            error("invalid type of storage in values");
        memset(mapping, 0, 256);
        for (i = 0; i < nletters; i++) {
            SEXP c = STRING_ELT(letters, i);
            if (LENGTH(c) != 1)
                error("invalid names for mapping");
            mapping[INTEGER(alphMapping)[i]] = CHAR(c)[0];
        }
        src = (unsigned char*) CHAR(values) + start;
        for (i = 0; i < ndest; i++) {
            char v = mapping[src[i]];
            if (!v) {
                bad_element = src[i];
                bad_index = i;
                goto invalid_element;
            }
            dest[i] = v;
        }
    } else {
#ifndef NOT_YET
        error("large alphabets are not supported");
#else
        char* alph = CHAR(letters);
        unsigned int* src;
        if (sizeof(unsigned int) != 4) {
            error("integer size is assumed to be 32 bits");
        }
        if (TYPEOF(values) != INTSXP)
            error("invalid type of storage in values");
        src = (unsigned int*) INTEGER(values) + start;
        for (i = 0; i < ndest; i++) {
            char indx = 0;
            switch(src[i]) {

#define SRC_I_CASE(n)                                 \
            case (1U << (n)):                     \
                indx = (n);                       \
                break

                SRC_I_CASE(0);
                SRC_I_CASE(1);
                SRC_I_CASE(2);
                SRC_I_CASE(3);
                SRC_I_CASE(4);
                SRC_I_CASE(5);
                SRC_I_CASE(6);
                SRC_I_CASE(7);
                SRC_I_CASE(8);
                SRC_I_CASE(9);
                SRC_I_CASE(10);
                SRC_I_CASE(11);
                SRC_I_CASE(12);
                SRC_I_CASE(13);
                SRC_I_CASE(14);
                SRC_I_CASE(15);
                SRC_I_CASE(16);
                SRC_I_CASE(17);
                SRC_I_CASE(18);
                SRC_I_CASE(19);
                SRC_I_CASE(20);
                SRC_I_CASE(21);
                SRC_I_CASE(22);
                SRC_I_CASE(23);
                SRC_I_CASE(24);
                SRC_I_CASE(25);
                SRC_I_CASE(26);
                SRC_I_CASE(27);
                SRC_I_CASE(28);
                SRC_I_CASE(29);
                SRC_I_CASE(30);
                SRC_I_CASE(31);

#undef SRC_I_CASE

            default:
                bad_element = src[i];
                bad_index = i;
                goto invalid_element;
            }
            if (indx >= nletters) {
                bad_element = src[i];
                bad_index = i;
                goto invalid_element;
            }
            dest[i] = alph[indx];
        }
#endif
    }
    UNPROTECT(1);
    return ans;

invalid_element:
    UNPROTECT(1);
    error("invalid %d-th element in values: %d",
          bad_index+start, bad_element);
    return R_NilValue; /* -Wall */
}

SEXP
BioStringToRString(SEXP x)
{
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP tmp = alph;
    while (isFromClass(tmp, "BioPatternAlphabet"))
        tmp = GET_SLOT(tmp, install("baseAlphabet"));
    return BioStringToCharacter(LENGTH(GET_SLOT(tmp, install("mapping"))),
                                GET_SLOT(alph, install("mapping")),
                                GET_SLOT(x, install("values")),
                                asInteger(GET_SLOT(x,
                                                   install("start"))),
                                asInteger(GET_SLOT(x,
                                                   install("end"))));
}

/*
 * Warning: This uses 1-based indexing rather than 0-based indexing.
 * Warning: pattern is assumed to have an extra character in front (which
 * will be overwritten by this function).
 * */
static SEXP
reverseFundamentalPreprocessing(char* pattern, int n, int anymatch)
{
    SEXP ans;
    int* N;

    ans = allocVector(INTSXP, n);
    /* Using one based indexing! */
    N = INTEGER(ans)-1;

    N[n] = n;
    if (n > 1) {
        int i, l, r, k;
        unsigned char prev, next;
        unsigned char mismatch;
        /* do an explicit calculation for n-1 */
        pattern[0] = (char)0xff;
        prev = pattern[n];
        i = k = n-1;
        next = pattern[i];
        if (anymatch) {
            while ((prev & next)) {
                prev = next;
                next = pattern[--i];
            }
        } else {
            while ((prev & next) == next) {
                prev = next;
                next = pattern[--i];
            }
        }
        if (i == k) {
            N[k] = 0;
            l = r = n;
        } else {
            N[k] = k-i;
            l = i+1;
            r = k;
        }
        for (k = n-2; k > 0; k--) {
            if (k < l) {
                int j;
                if (anymatch) {
                    for (j = n, i = k;
                         (pattern[i] & pattern[j]);
                         i--, --j) {
                    }
                } else {
                    for (j = n, i = k;
                         (pattern[i] & pattern[j]) == pattern[i];
                         i--, --j) {
                    }
                }
                if (i == k)
                    N[k] = 0;
                else {
                    N[k] = k-i;
                    l = i+1;
                    r = k;
                }
            } else {
                int k1 = n-r+k;
                int k_l = k-l;
                if (N[k1] <= k_l)
                    N[k] = N[k1];
                else if (N[k1] > k_l+1) {
                    N[k] = k_l+1;
                    r = k;
                } else {
                    int j;
                    if (anymatch) {
                        for (j = k_l, i = l-1;
                             (pattern[i] & pattern[j]);
                             i--, j--) {
                        }
                    } else {
                        for (j = k_l, i = l-1;
                             (pattern[i] & pattern[j]) == pattern[i];
                             i--, j--) {
                        }
                    }
                    N[k] = k-i;
                    l = i+1;
                    r = k;
                }
            }
        }
#ifdef DEBUG_BIOSTRINGS
        for (i = 1; i <= n; i++)
            Rprintf("%d, ", N[i]);
        Rprintf("\n");
#endif
    }
    return ans;
}

/*
 * Return an integer vector N such that (with 1 based subsetting)
 * N[j] is the length of the longest suffix of the substring
 * substr(x, 1, j) that is also a suffix of x.
 */
SEXP
do_reverseFundamentalPreprocessing(SEXP x, SEXP anymatch)
{
    SEXP ans;
    int n;
    char* pattern;
    if (TYPEOF(x) != CHARSXP) {
        if (TYPEOF(x) != STRSXP || LENGTH(x) != 1)
            error("argument not a character string");
        x = STRING_ELT(x, 0);
    }
    n = LENGTH(x);
    pattern = CHAR(PROTECT(allocString(n)));
    pattern[0] = 0;
    memcpy(pattern+1, CHAR(x), n);
    ans = reverseFundamentalPreprocessing(pattern, n, asLogical(anymatch));
    UNPROTECT(1);
    return ans;
}

typedef struct {
    int length;
    SEXP pattern;
    int start;
    int* good_suffix_shift;
    union {
        int R[256];
        SEXP letterIndex;
    } bad_char;
    int nletters;
    int usesChar;
} BoyerMoore_compiledPattern_t;

static void
BoyerMoore_preprocess(SEXP x,
                      BoyerMoore_compiledPattern_t* pattern)
{
    int xstart = asInteger(GET_SLOT(x, install("start")));
    int xend = asInteger(GET_SLOT(x, install("end")));
    int i, n;
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP alphMapping, letters;
    alphMapping = GET_SLOT(alph, install("mapping"));
    letters = GET_NAMES(alphMapping);
    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (xstart <= 0 || xstart > xend || xend > length(x)-1)
        error("BoyerMoore_preprocess: inconsistent BioString object");
    pattern->usesChar = TYPEOF(x) == CHARSXP;
    pattern->length = n = xend-xstart+1;
    pattern->nletters = LENGTH(alphMapping);
    if (pattern->usesChar) {
        SEXP tmp = PROTECT(allocVector(INTSXP,
                                       n+1));
        int* tmpptr = INTEGER(tmp);
        unsigned int* alphMappingptr = (unsigned int*) INTEGER(alphMapping);
        int* N;
        SEXP Nvec;
        char* patternptr;
        int k;

        pattern->pattern = x;
        pattern->start = xstart;
        patternptr = CHAR(x)+xstart-1;

        /* calculation of all occurance of each bit-pattern (used by
         * the bad character rule) */
        pattern->bad_char.letterIndex = allocVector(VECSXP, 256);
        PROTECT(pattern->bad_char.letterIndex);
        for (i = 0; i < pattern->nletters; i++) {
            unsigned int pat = alphMappingptr[i];
            int j;
            int lastj = n;
            int* indx;
            if (pat >= (1U << CHAR_BIT))
                error("invalid mapping with character storage");

            SET_VECTOR_ELT(pattern->bad_char.letterIndex, pat,
                           allocVector(INTSXP, n+1));
            indx = INTEGER(VECTOR_ELT(pattern->bad_char.letterIndex, pat));
            for (j = n-1; j > 0; j--) {
                /* does the j-th pattern contain pat? */
                if ((patternptr[j] & pat) == pat) {
                    while (lastj > j)
                        indx[lastj--] = j;
                }
            }
            while (lastj >= 0)
                indx[lastj--] = 0;
#ifdef DEBUG_BIOSTRINGS
            Rprintf("bad char index for pattern %d\n", pat);
            for (j = 1; j <= n; j++)
                Rprintf("%d,", indx[j]);
            Rprintf("\n");
#endif
        }

        /* calculation of the good suffix rule shift. */
        /* first calculate L' */
        Nvec = reverseFundamentalPreprocessing(patternptr, n, 0);
        PROTECT(Nvec);
        N = INTEGER(Nvec)-1;
        pattern->good_suffix_shift = (int*) R_alloc(n+1,
                                                    sizeof(int));
        memset(pattern->good_suffix_shift, 0, (1+n)*sizeof(int));
        for (i = 1; i < n; i++) {
            int j = n-N[i]+1;
            pattern->good_suffix_shift[j] = i;
        }
        /* then calculate l' */
        UNPROTECT(1);
        Nvec = reverseFundamentalPreprocessing(patternptr, n, 1);
        PROTECT(Nvec);
        N = INTEGER(Nvec)-1;
        memset(tmpptr, 0, (1+n)*sizeof(int));
        for (i = n, k=1; i > 0; i--) {
            if (N[i] == i) {
                int last = n-i+1;
                for (; k < last; k++) {
                    tmpptr[k] = i;
                }
                tmpptr[k++] = i;
            }
        }
        /* finally set shift[i] to n-L'[i] if L'[i] non-zero and to
         * n-l'[i] otherwise */
        for (i = 1; i <= n; i++) {
            if (pattern->good_suffix_shift[i])
                pattern->good_suffix_shift[i] =
                    n-pattern->good_suffix_shift[i];
            else pattern->good_suffix_shift[i] = n-tmpptr[i];
        }
        /* used when a match occurs */
        pattern->good_suffix_shift[0] = n-tmpptr[2];
#ifdef DEBUG_BIOSTRINGS
        for (i = 0; i <= n; i++)
            Rprintf("%d, ", pattern->good_suffix_shift[i]);
        Rprintf("\n");
#endif
        UNPROTECT(3);
    } else {
        error("non-character patterns and strings unimplemented");
    }
}

SEXP
BoyerMoore_exactMatch(SEXP origPattern, SEXP x)
{
    int xstart = asInteger(GET_SLOT(x, install("start")));
    int xend = asInteger(GET_SLOT(x, install("end")));
    BoyerMoore_compiledPattern_t pattern;
    SEXP matchIndex = R_NilValue;
    int* index = 0;
    int k, patlen;

    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (xstart <= 0 || xstart > xend || xend > length(x)-1)
        error("BoyerMoore_exactMatch: inconsistent BioString object");
    BoyerMoore_preprocess(origPattern, &pattern);
    
    patlen = pattern.length;
    if (pattern.usesChar) {
        unsigned char* str = (unsigned char*) CHAR(x)+xstart-1;
        char* patternptr;
        int m = xend-xstart+1;
        int nmatch, nmatchIndex;
        if (TYPEOF(x) != CHARSXP)
            error("type mismatch between pattern and string");

        PROTECT(pattern.bad_char.letterIndex);
        nmatchIndex = (m-patlen+1)*exp(patlen*
                                       log(-((double) pattern.nletters)));
        if (nmatchIndex > 64)
            nmatchIndex /= 4;
        else if (nmatchIndex > 16)
            nmatchIndex = 16;
        else nmatchIndex = 4;
        if (nmatchIndex > m-patlen+1)
            nmatchIndex = m-patlen+1;
#ifdef DEBUG_BIOSTRINGS
        Rprintf("nmatchIndex: %d\n", nmatchIndex);
#endif
        matchIndex = allocVector(INTSXP, nmatchIndex);
        PROTECT(matchIndex);
        index = INTEGER(matchIndex);
        nmatch = 0;
        patternptr = CHAR(pattern.pattern)+pattern.start-1;
        patternptr[0] = 0; /* make the first (dummy)
                            * element 0 */
        for (k = patlen; k <= m; ) {
            int i, h;
#ifdef DEBUG_BIOSTRINGS
            Rprintf("k: %d\n", k);
#endif
            for (i = patlen, h = k; 
                 patternptr[i] & str[h]; /* this is 0 when
                                          * i == 0 */
                 i--, h--) {
                /* empty body */
            }
            if (i == 0) {
                if (nmatchIndex == nmatch) {
                    double proportion = (nmatch+1)/(double)(k-patlen+1);
                    int estimate = proportion*(m-k);
                    if (estimate == 0 && m > k)
                        nmatchIndex += 2;
                    else nmatchIndex += estimate+1;
                    matchIndex = lengthgets(matchIndex, nmatchIndex);
                    UNPROTECT(1);
                    PROTECT(matchIndex);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = k-patlen+1;
                k += pattern.good_suffix_shift[0];
            } else {
                SEXP letterIndex =
                    VECTOR_ELT(pattern.bad_char.letterIndex, str[h]);
                int bad_character_shift;
                R_assert(letterIndex != R_NilValue);
                bad_character_shift = i-INTEGER(letterIndex)[i];
#ifdef DEBUG_BIOSTRINGS
                Rprintf("i:%d, h:%d, ch:%d\n", i, h, str[h]);
#endif
                if (i == patlen) {
                    /* good suffix rule gives 1 here - so we ignore it */
                    k += bad_character_shift;
#ifdef DEBUG_BIOSTRINGS
                    Rprintf("char shift:%d\n", bad_character_shift);
#endif
                } else {
                    int good_suffix_shift = pattern.good_suffix_shift[i+1];
                    k += MAX(good_suffix_shift, bad_character_shift);
#ifdef DEBUG_BIOSTRINGS
                    Rprintf("char shift:%d\n", bad_character_shift);
                    Rprintf("suff shift:%d\n", good_suffix_shift);
#endif
                }
            }
        }
        matchIndex = lengthgets(matchIndex, nmatch);
        UNPROTECT(2);
    } else {
        error("non-character patterns and strings unimplemented");
        if (TYPEOF(x) != INTSXP)
            error("type mismatch between pattern and string");
    }
    /* BoyerMoore_releasePattern(&pattern); */
    return matchIndex;
}

SEXP
LengthOne_exactMatch(SEXP pattern, SEXP x)
{
    int start = asInteger(GET_SLOT(pattern, install("start")));
    int end = asInteger(GET_SLOT(pattern, install("end")));
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP matchIndex = R_NilValue;
    int* index = NULL;
    int nmatch, nmatchIndex, nletters;
    int m;

    pattern = R_ExternalPtrTag(GET_SLOT(pattern, install("values")));
    if (start <= 0 || end-start != 1 || end > length(pattern)-1)
        error("LengthOne_exactMatch: invalid pattern");
    start = asInteger(GET_SLOT(x, install("start")));
    end = asInteger(GET_SLOT(x, install("end")));
    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (start <= 0 || start > end || end > length(x)-1)
        error("LengthOne_exactMatch: inconsistent BioString object");
    if (TYPEOF(x) != TYPEOF(pattern))
        error("pattern and text must be of same type");
    while (isFromClass(alph, "BioPatternAlphabet"))
        alph = GET_SLOT(alph, install("baseAlphabet"));
    nletters = LENGTH(GET_SLOT(alph, install("mapping")));

    m = end-start+1;
    nmatchIndex = m/((double) nletters);
    if (nmatchIndex > 64)
        nmatchIndex /= 4;
    else if (nmatchIndex > 16)
        nmatchIndex = 16;
    else if (4 <= m)
        nmatchIndex = 4;
    matchIndex = allocVector(INTSXP, nmatchIndex);

    nmatch = 0;
    if (TYPEOF(pattern) == CHARSXP) {
        unsigned char pat = ((unsigned char*) CHAR(pattern))[start];
        unsigned char* xptr = (unsigned char*) CHAR(x);
        int i;
        for (i = start; i <= end; i++) {
            if (pat & xptr[i]) {
                if (nmatchIndex == nmatch) {
                    double proportion = (nmatch+1)/(double)i;
                    int estimate = proportion*(m-i);
                    if (estimate == 0 && m > i)
                        nmatchIndex += 2;
                    else nmatchIndex += estimate+1;
                    PROTECT(matchIndex);
                    matchIndex = lengthgets(matchIndex, nmatchIndex);
                    UNPROTECT(1);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = i;
            }
        }
    } else {
        unsigned int pat = ((unsigned int*) INTEGER(pattern))[start];
        unsigned int* xptr = (unsigned int*) INTEGER(x);
        int i;
        for (i = start; i <= end; i++) {
            if (pat & xptr[i]) {
                if (nmatchIndex == nmatch) {
                    double proportion = (nmatch+1)/(double)i;
                    int estimate = proportion*(m-i);
                    if (estimate == 0 && m > i)
                        nmatchIndex += 2;
                    else nmatchIndex += estimate+1;
                    PROTECT(matchIndex);
                    matchIndex = lengthgets(matchIndex, nmatchIndex);
                    UNPROTECT(1);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = i;
            }
        }
    }
    matchIndex = lengthgets(matchIndex, nmatch);
    return matchIndex;
}

SEXP
reverseComplementBioString(SEXP x)
{
    int start = asInteger(GET_SLOT(x, install("start")));
    int end = asInteger(GET_SLOT(x, install("end")));
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP mapping;
    SEXP letters;
    SEXP ans;

    while (isFromClass(alph, "BioPatternAlphabet"))
        alph = GET_SLOT(alph, install("baseAlphabet"));
    mapping = GET_SLOT(alph, install("mapping"));
    letters = GET_NAMES(mapping);
    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(x) == CHARSXP) {
        int n = end-start+1;
        if (start <= 0 || n < 0 || end > length(x)-1)
            error("LengthOne_exactMatch: invalid pattern");
        if (n == 0) {
            ans = allocString(LENGTH(x));
        } else {
            unsigned char* src = (unsigned char*) CHAR(x)+start;
            unsigned char* dest;
            unsigned char A = 0;
            unsigned char C = 0;
            unsigned char G = 0;
            unsigned char T = 0;
            unsigned char gap = 0;
            unsigned char revmap[256];
            int i;

            if (TYPEOF(mapping) != INTSXP ||
                TYPEOF(letters) != STRSXP ||
                LENGTH(mapping) != 5 ||
                LENGTH(letters) != 5)
                error("incorrect mapping");
            for (i = 0; i < 5; i++) {
                SEXP tmp = STRING_ELT(letters, i);
                if (TYPEOF(tmp) != CHARSXP || LENGTH(tmp) != 1)
                    error("incorrect mapping");
                switch (CHAR(tmp)[0]) {
                case 'a': case 'A':
                    A = 1 << i;
                    break;
                case 'c': case 'C':
                    C = 1 << i;
                    break;
                case 'g': case 'G':
                    G = 1 << i;
                    break;
                case 't': case 'T': case 'u': case 'U':
                    T = 1 << i;
                    break;
                default:
                    gap = 1 << i;
                    break;
                }
            }
            if (!A || !G || !C || !T || !gap)
                error("Could not find some of the nucleotide letters");
            memset(revmap, 0, 256);
            revmap[A] = T;
            revmap[C] = G;
            revmap[G] = C;
            revmap[T] = A;
            revmap[gap] = gap;
            ans = allocString(length(x));
            dest = (unsigned char*) CHAR(ans)+start-1;
            for (i = 0; i < n; i++) {
                unsigned char pat = revmap[src[i]];
                if (!pat)
                    error("Unrecognized code in nucleotide sequence");
                dest[n-i] = pat;
            }
        }
    } else {
        error("Only character storage is supported now");
        ans = R_NilValue; /* -Wall */
    }
    PROTECT(ans);
    ans = R_MakeExternalPtr(NULL, ans, R_NilValue);
    UNPROTECT(1);
    return ans;
}
