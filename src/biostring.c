/* Copyright (C) 2003 by Saikat DebRoy */
#include <limits.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#define INTERRUPT_BIOSTRINGS
#ifdef INTERRUPT_BIOSTRINGS
#define INTERRUPTCHECK_AFTER (1U << 21)
#endif

#define MAX(i, j) ((i)>(j)?(i):(j))
#define MIN(i, j) ((i)<(j)?(i):(j))

#define R_assert(e) ((e) ? (void) 0 : error("assertion `%s' failed: file `%s', line %d\n", #e, __FILE__, __LINE__))

static SEXP
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
    SEXP Class = PROTECT(mkString((char *) klass));
    SEXP call = PROTECT(lang3(install("is"),
                              x, Class));
    SEXP ans = eval(call, R_GlobalEnv);
    UNPROTECT(2);
    return asLogical(ans);
}

static SEXP
expandIndex(SEXP index, int ndone, int nleft)
{
    int n = LENGTH(index);
    double proportion = (n+1)/(double)(ndone);
    int estimate = proportion*nleft+1;
    n = 2*(n+estimate);
    return lengthgets(index, n);
}

static int
getBioStringLength(SEXP x, int** startvecptr, int** endvecptr)
{
    SEXP offsets, dim;
    if (!isFromClass(x, "BioString"))
        error("argument must be of class BioString");
    offsets = GET_SLOT(x, install("offsets"));
    dim = GET_DIM(offsets);
    if (TYPEOF(offsets) != INTSXP || TYPEOF(dim) != INTSXP ||
        LENGTH(dim) != 2 || INTEGER(dim)[1] != 2)
        error("offsets slot of BioString must be integer matrix with two columns");
    if (startvecptr)
        *startvecptr = INTEGER(offsets);
    if (endvecptr)
        *endvecptr = INTEGER(offsets)+INTEGER(dim)[0];
    return INTEGER(dim)[0];
}

static SEXP
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
        if (TYPEOF(vec) == CHARSXP) {
            /* lengthgets does not work with CHARSXP */
            SEXP newvec = allocString(newlength);
            memcpy(CHAR(newvec), CHAR(vec), oldlength);
            vec = newvec;
        } else {
            vec = lengthgets(vec, newlength);
        }
        R_SetExternalPtrTag(values, vec);
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
    reserveBioString(values, current_length+count);
    memset(mapping, 0, 256*sizeof(unsigned int));
    for (i = 0; i < nletters; i++) {
        SEXP c = STRING_ELT(letters, i);
        unsigned int tmpcode = ((unsigned int*) INTEGER(alphMapping))[i];
        int cval;
        if (LENGTH(c) != 1)
            error("invalid names for mapping");
        if (tmpcode > maxcode)
            maxcode = tmpcode;
        cval = CHAR(c)[0];
        mapping[(unsigned char)toupper(cval)] = tmpcode;
        mapping[(unsigned char)tolower(cval)] = tmpcode;
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
                      srcptr[j], srcindex+1, j+1);
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

static SEXP
setBioString(SEXP biostring, SEXP src)
{
    SEXP offsets;
    int i, n, last;
    int* start;
    int* end;

    if (!isFromClass(biostring, "BioString"))
        error("first argument must be from BioString class");
    if (asLogical(GET_SLOT(biostring, install("initialized"))))
        error("can not modify initialized strings");
    n = length(src);
    biostring = duplicate(biostring);
    PROTECT(biostring);
    offsets = allocMatrix(INTSXP, n, 2);
    PROTECT(offsets);
    SET_SLOT(biostring, install("offsets"), offsets);
    UNPROTECT(1);
    offsets = GET_SLOT(biostring, install("offsets"));
    start = INTEGER(offsets);
    end = start+n;
    for (last = 0, i = 0; i < n; i++) {
        start[i] = last+1;
        last = appendCharacterToBioString(
            GET_SLOT(GET_SLOT(biostring,
                              install("alphabet")),
                     install("mapping")),
            GET_SLOT(biostring, install("values")),
            last, src, i);
        end[i] = last;
    }
    UNPROTECT(1);
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
    if (start <= 0)
        start = 1;
    if (end > nvalues)
        end = nvalues;
    if (start > end)
        return R_BlankString;
    if (nletters_base > sizeof(int)*CHAR_BIT)
        error("alphabet is too large");
    values = R_ExternalPtrTag(values);

    ndest = end-start+1;
    ans = allocString(ndest);
    PROTECT(ans);
    dest = CHAR(ans);

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

static SEXP
BioStringToRString(SEXP x)
{
    SEXP alph, mapping, values, tmp, dim, ans, offsets;
    int nletters_base;
    int* startvec;
    int* endvec;
    int i, n;

    if (!isFromClass(x, "BioString"))
        error("argument must be a BioString");
    offsets = GET_SLOT(x, install("offsets"));
    if (TYPEOF(offsets) != INTSXP)
        error("offsets must be integer");
    dim = GET_DIM(offsets);
    if (TYPEOF(dim) != INTSXP || LENGTH(dim) != 2)
        error("offsets must be a matrix");
    if (INTEGER(dim)[1] != 2)
        error("offsets must have two columns");
    alph = GET_SLOT(x, install("alphabet"));
    mapping = GET_SLOT(alph, install("mapping"));
    values = GET_SLOT(x, install("values"));
    tmp = alph;
    n = INTEGER(dim)[0];
    startvec = INTEGER(offsets);
    endvec = INTEGER(offsets)+n;
    
    while (isFromClass(tmp, "BioPatternAlphabet"))
        tmp = GET_SLOT(tmp, install("baseAlphabet"));
    nletters_base = LENGTH(GET_SLOT(tmp, install("mapping")));
    ans = allocVector(STRSXP, n);
    PROTECT(ans);
    for (i = 0; i < n; i++) {
        int start_i = startvec[i];
        int end_i = endvec[i];
        SET_STRING_ELT(ans, i,
                       BioStringToCharacter(nletters_base,
                                            mapping, values,
                                            start_i, end_i));
    }
    UNPROTECT(1);
    return ans;
}

#define mod_iterate(n,n1,n2,n3, i,i1,i2,i3) for (i=i1=i2=i3=0; i<n; \
	i1 = (++i1 == n1) ? 0 : i1,\
	i2 = (++i2 == n2) ? 0 : i2,\
	i3 = (++i3 == n3) ? 0 : i3,\
	++i)

static SEXP
BioString_substring(SEXP x, SEXP start, SEXP stop, SEXP doSubstring)
{
    int* startvec;
    int* stopvec;
    int* current_startvec;
    int* current_stopvec;
    int* ans_startvec;
    int* ans_stopvec;
    int i, icurrent, istart, istop;
    int n, ncurrent, nstart, nstop;
    int substring = asLogical(doSubstring);
    SEXP offsets, dim, ans;

    if (!isFromClass(x, "BioString"))
        error("invalid argument to substr for BioString");
    offsets = GET_SLOT(x, install("offsets"));
    dim = GET_DIM(offsets);
    if (TYPEOF(offsets) != INTSXP || TYPEOF(dim) != INTSXP ||
        LENGTH(dim) != 2 || INTEGER(dim)[1] != 2)
        error("offsets slot of BioString must be integer matrix with two columns");

    start = coerceVector(start, INTSXP);
    PROTECT(start);
    nstart = LENGTH(start);
    startvec = INTEGER(start);
    stop = coerceVector(stop, INTSXP);
    PROTECT(stop);
    nstop = LENGTH(stop);
    stopvec = INTEGER(stop);

    ncurrent = INTEGER(dim)[0];
    current_startvec = INTEGER(offsets);
    current_stopvec = INTEGER(offsets)+ncurrent;

    ans = duplicate(x);
    PROTECT(ans);
    n = ncurrent;
    if (substring) {
        n = (n>nstart)?n:nstart;
        n = (n>nstop)?n:nstop;
    }
    if (n != ncurrent) {
        SEXP dimnames = GET_DIMNAMES(offsets);
        SEXP ans_offsets = allocVector(INTSXP, 2*n);
        PROTECT(ans_offsets);
        memcpy(INTEGER(ans_offsets), current_startvec,
               sizeof(int)*ncurrent);
        memset(INTEGER(ans_offsets)+ncurrent, 0,
               sizeof(int)*(n-ncurrent));
        memcpy(INTEGER(ans_offsets)+n, current_stopvec,
               sizeof(int)*ncurrent);
        memset(INTEGER(ans_offsets)+n+ncurrent, 0,
               sizeof(int)*(n-ncurrent));
        dim = allocVector(INTSXP, 2);
        INTEGER(dim)[0] = n;
        INTEGER(dim)[1] = 2;
        PROTECT(dim);
        SET_DIM(ans_offsets, dim);
        if (TYPEOF(dimnames) == VECSXP && LENGTH(dimnames) == 2) {
            SEXP tmp = allocVector(VECSXP, 2);
            SET_VECTOR_ELT(tmp, 1, VECTOR_ELT(dimnames, 1));
            SET_DIMNAMES(ans_offsets, tmp);
        }
        SET_SLOT(ans, install("offsets"), ans_offsets);
        UNPROTECT(2);
    }

    offsets = GET_SLOT(ans, install("offsets"));
    ans_startvec = INTEGER(offsets);
    ans_stopvec = INTEGER(offsets)+n;

    mod_iterate(n, ncurrent, nstart, nstop,
                i, icurrent, istart, istop) {
        int current_first = current_startvec[icurrent];
        int current_last = current_stopvec[icurrent];
        int slen = current_last-current_first+1;
        if (slen > 0) {
            int first = startvec[istart];
            int last = stopvec[istop];
            if (first <= 0)
                first = 1;
            if (first > last || first > slen) {
                ans_startvec[i] = 1;
                ans_stopvec[i] = 0;
            } else {
                if (last < slen)
                    ans_stopvec[i] = current_first+last-1;
                else ans_stopvec[i] = current_last;
                ans_startvec[i] = current_first+first-1;
            }
        }
    }
    UNPROTECT(3);
    return ans;
}

#undef mod_iterate

/*
 * Warning: This uses 1-based indexing rather than 0-based indexing.
 * Warning: pattern is assumed to have an extra character in front (which
 * will be overwritten by this function).
 *
 * We are looking at a string of length n and finding N[j], the length
 * of the longest suffix of pattern[1:j] that is also a suffix of
 * pattern. However, we use two rules for matching pattern. In one, we
 * say pattern[i] matches pattern[j] (where i < j) if all the set bits
 * in pattern[i] are also set in pattern[j]. In the other rule we say
 * they match if there is any bit which is set in both. The N[j]
 * from the first rule are stored in Nsome[j] and the N[j] from the
 * second is stored in Nany[j].
 * 
 * */
static SEXP
reverseFundamentalPreprocessing(char* pattern, int n)
{
    int i;
    SEXP ans;
    int* Nsome;
    int* Nany;

    ans = allocVector(VECSXP, 2);
    PROTECT(ans);
    SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, n));
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));
    /* Using one based indexing! */
    Nsome = INTEGER(VECTOR_ELT(ans, 0))-1;
    Nany = INTEGER(VECTOR_ELT(ans, 1))-1;

    Nsome[n] = n;
    Nany[n] = n;
    if (n > 1) {
        int k;
        pattern[0] = (char) 0xff;
        for (k = n-1; k > 0; k--) {
            int j;
            for (j = n, i = k;
                 !(pattern[i] & ~pattern[j]);
                 i--, --j) {
            }
            Nsome[k] = k-i;
        }

        pattern[0] = (char) 0;
        for (k = n-1; k > 0; k--) {
            int j;
            for (j = n, i = k;
                 (pattern[i] & pattern[j]);
                 i--, --j) {
            }
            Nany[k] = k-i;
        }
    }
#ifdef DEBUG_BIOSTRINGS
    for (i = 1; i <= n; i++)
        Rprintf("%d,", Nsome[i]);
    Rprintf("\n");
    for (i = 1; i <= n; i++)
        Rprintf("%d,", Nany[i]);
    Rprintf("\n");
#endif
    UNPROTECT(1);
    return ans;
}
#if 0
static SEXP
reverseFundamentalPreprocessing(char* pattern, int n)
{
    SEXP ans;
    int* Nsome;
    int* Nany;

    ans = allocVector(VECSXP, 2);
    PROTECT(ans);
    SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, n));
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));
    /* Using one based indexing! */
    Nsome = INTEGER(VECTOR_ELT(ans, 0))-1;
    Nany = INTEGER(VECTOR_ELT(ans, 1))-1;

    Nsome[n] = n;
    Nany[n] = n;
    if (n > 1) {
        int i, lsome, rsome, lany, rany, k;
        unsigned char prev, next;
        char save_first = pattern[0];
        /* set the first character */

        /* do an explicit calculation for n-1 */
        pattern[0] = (char) 0xff;
        prev = pattern[n];
        i = k = n-1;
        next = pattern[i];
        while (!(~prev & next)) {
            prev = next;
            next = pattern[--i];
        }
        if (i == k) {
            Nsome[k] = 0;
            lsome = rsome = n;
        } else {
            Nsome[k] = k-i;
            lsome = i+1;
            rsome = k;
        }

        pattern[0] = (char) 0;
        while ((prev & next)) {
            prev = next;
            next = pattern[--i];
        }
        if (i == k) {
            Nany[k] = 0;
            lany = rany = n;
        } else {
            Nany[k] = k-i;
            lany = i+1;
            rany = k;
        }
        for (k = n-2; k > 0; k--) {
#if 0
            if (k < l) {
                int j;
#ifdef DEBUG_BIOSTRINGS
                Rprintf("%d<%d:", k, l);
#endif
                if (anymatch) {
                    for (j = n, i = k;
                         (pattern[i] & pattern[j]);
                         i--, --j) {
                    }
                } else {
                    for (j = n, i = k;
                         !(pattern[i] & ~pattern[j]);
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
#ifdef DEBUG_BIOSTRINGS
                Rprintf("%d\n", N[k]);
#endif
            } else {
                int k1 = n-r+k;
                int k_l = k-l;
                if (N[k1] <= k_l) {
                    N[k] = N[k1];
#ifdef DEBUG_BIOSTRINGS
                    Rprintf("N[%d]<%d-%d+1:%d\n", k1, k, l, N[k]);
#endif
                } else if (N[k1] > k_l+1) {
                    N[k] = k_l+1;
                    r = k;
#ifdef DEBUG_BIOSTRINGS
                    Rprintf("N[%d]>%d-%d+1:%d\n", k1, k, l, N[k]);
#endif
                } else {
                    int j;
#ifdef DEBUG_BIOSTRINGS
                    Rprintf("N[%d]==%d-%d+1:", k1, k, l);
#endif
                    if (anymatch) {
                        for (j = k_l, i = l-1;
                             (pattern[i] & pattern[j]);
                             i--, j--) {
                        }
                    } else {
                        for (j = k_l, i = l-1;
                             (pattern[i] & ~pattern[j]);
                             i--, j--) {
                        }
                    }
                    N[k] = k-i;
                    l = i+1;
                    r = k;
#ifdef DEBUG_BIOSTRINGS
                Rprintf("%d\n", N[k]);
#endif
                }
            }
#else
            int k_l = k-lsome;
            int j;
            pattern[0] = (char) 0xff;
            if (k_l < 0) {
                for (j = n, i = k;
                     !(pattern[i] & ~pattern[j]);
                     i--, --j) {
                }
                if (i == k)
                    Nsome[k] = 0;
                else {
                    Nsome[k] = k-i;
                    lsome = i+1;
                    rsome = k;
                }
            } else {
                int k1 = n-rsome+k;
                if (Nsome[k1] <= k_l) {
                    Nsome[k] = Nsome[k1];
                } else if (Nsome[k1] > k_l+1) {
                    Nsome[k] = k_l+1;
                    rsome = k;
                } else {
                    for (j = k_l, i = lsome-1;
                         (pattern[i] & ~pattern[j]);
                         i--, j--) {
                    }
                    Nsome[k] = k-i;
                    lsome = i+1;
                    rsome = k;
                }
            }

            i = k-Nsome[k];
            j = n-Nsome[k];
            pattern[0] = (char) 0;
            for (;
                 (pattern[i] & pattern[j]);
                 i--, j--) {
            }
            if (i == k)
                Nany[k] = 0;
            else {
                Nany[k] = k-i;
                lany = i+1;
                rany = k;
            }
#endif
        }
        /* restore the first character */
        pattern[0] = save_first;
#ifdef DEBUG_BIOSTRINGS
        for (i = 1; i <= n; i++)
            Rprintf("%d,", Nsome[i]);
        Rprintf("\n");
        for (i = 1; i <= n; i++)
            Rprintf("%d,", Nany[i]);
        Rprintf("\n");
#endif
    }
    UNPROTECT(1);
    return ans;
}
#endif

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
getLengthOneBioStringRange(SEXP x, int*start, int* end)
{
    int xstart, xend;
    SEXP offsets, dim;

    if (!isFromClass(x, "BioString"))
        error("x must be a BioString");
    offsets = GET_SLOT(x, install("offsets"));
    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    dim = GET_DIM(offsets);
    if (TYPEOF(offsets) != INTSXP || TYPEOF(dim) != INTSXP ||
        LENGTH(dim) != 2 || INTEGER(dim)[1] != 2)
        error("offsets slot of BioString must be integer matrix with two columns");
    if (INTEGER(dim)[0] != 1)
        error("not a single BioString");
    xstart = INTEGER(offsets)[0];
    xend = INTEGER(offsets)[1];
    if (xstart <= 0)
        xstart = 1;
    if (xstart > xend || xstart > length(x)-1) {
        *start = 1;
        *end = 0;
        return;
    }
    if (xend > length(x)-1)
        xend = length(x)-1;
    *start = xstart;
    *end = xend;
}

/*
 *  This may trash the matchIndex argument.
 */
static SEXP
matchIndexToBioString(SEXP x, SEXP matchIndex, int nmatch, int patlen)
{
    int nmatchIndex = LENGTH(matchIndex);
    int* index = INTEGER(matchIndex);
    x = duplicate(x);
    PROTECT(x);
#ifdef DEBUG_BIOSTRINGS
    Rprintf("In matchIndexToBioString\nnmatch: %d\n", nmatch);
#endif
    if (nmatch == 0) {
        SEXP offsets = allocMatrix(INTSXP, 0, 2);
        PROTECT(offsets);
        SET_SLOT(x, install("offsets"), offsets);
        UNPROTECT(1);
    } else if (nmatch == 1) {
        int* offsets = INTEGER(GET_SLOT(x, install("offsets")));
        offsets[1] = index[0];
        offsets[0] = index[0]-patlen+1;
    } else if (nmatchIndex == 2*nmatch) {
        SEXP dim = GET_DIM(GET_SLOT(x, install("offsets")));
        int i;

        INTEGER(dim)[0] = nmatch;
        INTEGER(dim)[1] = 2;
        SET_DIM(matchIndex, dim);
        memcpy(index+nmatch, index, sizeof(int)*nmatch);
        for (i = 0; i < nmatch; i++)
            index[i] -= patlen-1;
        SET_SLOT(x, install("offsets"), matchIndex);
    } else {
        SEXP offsets = allocMatrix(INTSXP, nmatch, 2);
        int* tmp = INTEGER(offsets);
        int i;

        PROTECT(offsets);
        for (i = 0; i < nmatch; i++)
            tmp[i] = index[i]-patlen+1;
        tmp += nmatch;
        memcpy(tmp, index, sizeof(int)*nmatch);
        SET_SLOT(x, install("offsets"), offsets);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return x;
}

static SEXP
ForwardSearch_exactMatch(SEXP pattern, SEXP x)
{
    int pstart, pend, xstart, xend, patlen = 0;
    SEXP alph, vec;
    SEXP matchIndex = R_NilValue;
    PROTECT_INDEX matchIndex_pi;
    int* index = NULL;
    int nmatchIndex, nletters;
    int nmatch = 0;
#ifdef INTERRUPT_BIOSTRINGS
    unsigned long interruptcheck;
#endif
    int m;

    getLengthOneBioStringRange(pattern, &pstart, &pend);
    PROTECT_WITH_INDEX(matchIndex, &matchIndex_pi);
    if (pstart > pend)
        goto finished_match;
    getLengthOneBioStringRange(x, &xstart, &xend);
    patlen = pend-pstart+1;
    if (xstart > xend)
        goto finished_match;
    alph = GET_SLOT(x, install("alphabet"));
    pattern = R_ExternalPtrTag(GET_SLOT(pattern, install("values")));
    vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(vec) != TYPEOF(pattern))
        error("pattern and text must be of same type");
    while (isFromClass(alph, "BioPatternAlphabet"))
        alph = GET_SLOT(alph, install("baseAlphabet"));
    nletters = LENGTH(GET_SLOT(alph, install("mapping")));

    m = xend-xstart+1;
    nmatchIndex = (m-patlen+1)*exp(patlen*
                                   log(-((double) nletters)));
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
    REPROTECT(matchIndex, matchIndex_pi);
    index = INTEGER(matchIndex);

    nmatch = 0;
#ifdef INTERRUPT_BIOSTRINGS
    interruptcheck = 0UL;
#endif
    if (TYPEOF(pattern) == CHARSXP) {
        unsigned char* pptr = ((unsigned char*) CHAR(pattern))+pstart-1;
        unsigned char* xptr = (unsigned char*) CHAR(vec)+xstart-1;
        int k;
        unsigned char save_first = pptr[0];
        pptr[0] = 0;
        for (k = patlen; k <= m; k++) {
            int i, j;
            for (i = patlen, j = k; pptr[i] & xptr[j]; i--, j--) {
            }
            if (i == 0) {
                if (nmatchIndex == nmatch) {
                    pptr[0] = save_first;
                    matchIndex = expandIndex(matchIndex, k-patlen+1, m-k);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                    pptr[0] = 0;
                }
                index[nmatch++] = k;
            }
#ifdef INTERRUPT_BIOSTRINGS
            interruptcheck += patlen-i;
            if (interruptcheck > INTERRUPTCHECK_AFTER) {
                pptr[0] = save_first;
                R_CheckUserInterrupt();
                interruptcheck = 0UL;
                pptr[0] = 0;
            }
#endif
        }
        pptr[0] = save_first;
    } else {
        unsigned int* pptr = ((unsigned int*) INTEGER(pattern))+pstart-1;
        unsigned int* xptr = (unsigned int*) INTEGER(vec)+xstart-1;
        int k;
        unsigned int save_first = pptr[0];
        pptr[0] = 0;
        for (k = patlen; k <= m; k++) {
            int i, j;
            for (i = patlen, j = k; pptr[i] & xptr[j]; i--, j--) {
            }
            if (i == 0) {
                if (nmatchIndex == nmatch) {
                    pptr[0] = save_first;
                    matchIndex = expandIndex(matchIndex, k-patlen+1, m-k);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                    pptr[0] = 0;
                }
                index[nmatch++] = k;
            }
#ifdef INTERRUPT_BIOSTRINGS
            interruptcheck += patlen-i;
            if (interruptcheck > INTERRUPTCHECK_AFTER) {
                pptr[0] = save_first;
                R_CheckUserInterrupt();
                interruptcheck = 0UL;
                pptr[0] = 0;
            }
#endif
        }
        pptr[0] = save_first;
    }
finished_match:
    x = matchIndexToBioString(x, matchIndex, nmatch, patlen);
    UNPROTECT(1);
    return x;
}

static SEXP
LengthOne_exactMatch(SEXP pattern, SEXP x)
{
    int start, end;
    SEXP alph;
    SEXP matchIndex = R_NilValue;
    PROTECT_INDEX matchIndex_pi;
    int* index = NULL;
    int nmatchIndex, nletters;
    int nmatch = 0;
    int m;

    getLengthOneBioStringRange(pattern, &start, &end);
    if (start <= 0 || end != start)
        error("not a length one pattern");
    getLengthOneBioStringRange(x, &start, &end);
    if (start > end)
        goto finished_match;
    alph = GET_SLOT(x, install("alphabet"));
    pattern = R_ExternalPtrTag(GET_SLOT(pattern, install("values")));
    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
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
    PROTECT_WITH_INDEX(matchIndex, &matchIndex_pi);

    nmatch = 0;
    if (TYPEOF(pattern) == CHARSXP) {
        unsigned char pat = ((unsigned char*) CHAR(pattern))[start];
        unsigned char* xptr = (unsigned char*) CHAR(x);
        int i;
        for (i = start; i <= end; i++) {
            if (pat & xptr[i]) {
                if (nmatchIndex == nmatch) {
                    matchIndex = expandIndex(matchIndex, i, m-i);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = i;
            }
#ifdef INTERRUPT_BIOSTRINGS
            if (i % INTERRUPTCHECK_AFTER) {
                R_CheckUserInterrupt();
            }
#endif
        }
    } else {
        unsigned int pat = ((unsigned int*) INTEGER(pattern))[start];
        unsigned int* xptr = (unsigned int*) INTEGER(x);
        int i;
        for (i = start; i <= end; i++) {
            if (pat & xptr[i]) {
                if (nmatchIndex == nmatch) {
                    matchIndex = expandIndex(matchIndex, i, m-i);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = i;
            }
#ifdef INTERRUPT_BIOSTRINGS
            if (i % INTERRUPTCHECK_AFTER) {
                R_CheckUserInterrupt();
            }
#endif
        }
    }
    UNPROTECT(1);
finished_match:
    return matchIndexToBioString(x, matchIndex, nmatch, 1);
}

static void
BoyerMoore_preprocess(SEXP x, BoyerMoore_compiledPattern_t* pattern)
{
    int i, n, xstart, xend;
    SEXP alph;
    SEXP alphMapping, letters;

    getLengthOneBioStringRange(x, &xstart, &xend);
    if (xstart > xend) {
        pattern->length = 0;
        return;
    }

    alph = GET_SLOT(x, install("alphabet"));
    alphMapping = GET_SLOT(alph, install("mapping"));
    letters = GET_NAMES(alphMapping);

    x = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    pattern->usesChar = TYPEOF(x) == CHARSXP;
    pattern->length = n = xend-xstart+1;
    pattern->nletters = LENGTH(alphMapping);

    if (pattern->usesChar) {
        SEXP tmp = PROTECT(allocVector(INTSXP,
                                       n+1));
        int* tmpptr = INTEGER(tmp);
        unsigned int* alphMappingptr = (unsigned int*) INTEGER(alphMapping);
        int* N;
        SEXP Nvecs;
        char* patternptr;
        int k;

        pattern->pattern = x;
        pattern->start = xstart;
        patternptr = CHAR(x)+xstart-1;

        /* calculation of the good suffix rule shift. */
        /* first calculate L' */
        Nvecs = reverseFundamentalPreprocessing(patternptr, n);
        PROTECT(Nvecs);

        /* calculation of all occurance of each bit-pattern (used by
         * the bad character rule) */
        pattern->bad_char.letterIndex = allocVector(VECSXP, 256);
        PROTECT(pattern->bad_char.letterIndex);
#ifdef NOT_YET
        N = INTEGER(VECTOR_ELT(Nvecs, 1))-1;
#endif
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
#ifdef NOT_YET
            indx[n] = 0;
            lastj = n+1;
            for (j = n; j > 0; j--) {
                /* does the j-th pattern contain pat? */
                if ((patternptr[j] & pat) == pat) {
                    for (; lastj > j; lastj--)
                        indx[lastj] = j-lastj;
                }
            }
#else
            for (j = n-1; j > 0; j--) {
                /* does the j-th pattern contain pat? */
                if ((patternptr[j] & pat) == pat) {
                    for (; lastj > j; lastj--)
                        indx[lastj] = lastj-j;
                }
            }
#endif
            for (; lastj >= 0; lastj--)
                indx[lastj] = lastj;
#ifdef DEBUG_BIOSTRINGS
            Rprintf("bad char index for pattern %d\n", pat);
            for (j = 1; j <= n; j++)
                Rprintf("%d,", indx[j]);
            Rprintf("\n");
#endif
        }

        N = INTEGER(VECTOR_ELT(Nvecs, 0))-1;
        memset(tmpptr, 0, (1+n)*sizeof(int));
        for (i = 1; i < n; i++) {
            int j = n-N[i];
            tmpptr[j] = i;
        }
        N = INTEGER(VECTOR_ELT(Nvecs, 1))-1;
        pattern->good_suffix_shift = (int*) R_alloc(n+1,
                                                    sizeof(int));
        memset(pattern->good_suffix_shift, 0, (1+n)*sizeof(int));
        for (i = 1; i < n; i++) {
            int j = n-N[i]+1;
            pattern->good_suffix_shift[j] = i;
        }
        for (i = 1; i <= n; i++) {
            int L0 = tmpptr[i];
            int L1 = pattern->good_suffix_shift[i];
            if (L0 > L1)
                pattern->good_suffix_shift[i] = L0;
        }
        
        /* then calculate l' */
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

static SEXP
BoyerMoore_exactMatch(SEXP origPattern, SEXP x)
{
#ifdef INTERRUPT_BIOSTRINGS
    unsigned long interruptcheck = 0UL;
#endif
    int xstart, xend;
    BoyerMoore_compiledPattern_t pattern;
    SEXP matchIndex = R_NilValue;
    PROTECT_INDEX matchIndex_pi;
    SEXP vec;
    int* index = 0;
    int k, patlen = 0, nmatch = 0;

    getLengthOneBioStringRange(x, &xstart, &xend);
    if (xstart > xend)
        goto finished_match;
    BoyerMoore_preprocess(origPattern, &pattern);
    patlen = pattern.length;
    if (patlen == 0)
        goto finished_match;
    if (patlen == 1)
        return LengthOne_exactMatch(origPattern, x);
    if (pattern.usesChar) {
        int* shiftTable[256];
        unsigned char* str;
        char* patternptr;
        char save_first;
        int m = xend-xstart+1;
        int nmatchIndex;
        PROTECT(pattern.bad_char.letterIndex);
        vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
        if (TYPEOF(vec) != CHARSXP)
            error("type mismatch between pattern and string");

        str = (unsigned char*) CHAR(vec)+xstart-1;
        nmatchIndex = 20;
        if (nmatchIndex > m-patlen+1)
            nmatchIndex = 2;
#ifdef DEBUG_BIOSTRINGS
        Rprintf("nmatchIndex: %d\n", nmatchIndex);
#endif
        matchIndex = allocVector(INTSXP, nmatchIndex);
        PROTECT_WITH_INDEX(matchIndex, &matchIndex_pi);
        index = INTEGER(matchIndex);
        nmatch = 0;
        patternptr = CHAR(pattern.pattern)+pattern.start-1;
        save_first = patternptr[0];
        patternptr[0] = 0; /* make the first (dummy)
                            * element 0 */
        for (k = 0; k < 256; k++) {
            SEXP tmp = VECTOR_ELT(pattern.bad_char.letterIndex, k);
            if (tmp == R_NilValue)
                shiftTable[k] = NULL;
            else if (TYPEOF(tmp) != INTSXP || LENGTH(tmp) != patlen+1)
                error("invalid letterIndex, %d, %d, %d", TYPEOF(tmp),
                      LENGTH(tmp), patlen);
            else shiftTable[k] = INTEGER(tmp);
        }
        for (k = patlen; k <= m; ) {
            int i, h;
            for (i = patlen, h = k; 
                 patternptr[i] & str[h]; /* this is 0 when i == 0 */
                 i--, h--) {
                /* empty body */
            }
            if (i == 0) {
                if (nmatchIndex == nmatch) {
                    patternptr[0] = save_first;
                    matchIndex = expandIndex(matchIndex, k-patlen+1, m-k);
                    REPROTECT(matchIndex, matchIndex_pi);
                    index = INTEGER(matchIndex);
                    nmatchIndex = LENGTH(matchIndex);
                    patternptr[0] = 0;
                }
                index[nmatch++] = k;
                k += pattern.good_suffix_shift[0];
            } else {
                int* letterIndex = shiftTable[str[h]];
                int bad_character_shift;
                R_assert(letterIndex != NULL);
#if 0
                bad_character_shift = letterIndex[i];
                if (i == patlen) {
                    /* good suffix rule gives 1 here - so we ignore it */
                    k += bad_character_shift;
                } else {
                    int good_suffix_shift = pattern.good_suffix_shift[i+1];
                    k += MAX(good_suffix_shift, bad_character_shift);
                }
#else
                k += letterIndex[i];
#endif
            }
#ifdef INTERRUPT_BIOSTRINGS
            interruptcheck += patlen-i;
            if (interruptcheck > INTERRUPTCHECK_AFTER) {
                patternptr[0] = save_first;
                R_CheckUserInterrupt();
                interruptcheck = 0UL;
                patternptr[0] = 0;
            }
#endif
        }
        patternptr[0] = save_first;
        UNPROTECT(2);
    } else {
        error("non-character patterns and strings unimplemented");
    }
    /* BoyerMoore_releasePattern(&pattern); */
finished_match:
    return matchIndexToBioString(x, matchIndex, nmatch, patlen);
}

static SEXP
reverseComplementBioString(SEXP x)
{
    SEXP alph, mapping, letters, xvec;
    int n;

    if (!isFromClass(x, "BioString"))
        error("argument must be of class BioString");
    alph = GET_SLOT(x, install("alphabet"));
    while (isFromClass(alph, "BioPatternAlphabet"))
        alph = GET_SLOT(alph, install("baseAlphabet"));
    mapping = GET_SLOT(alph, install("mapping"));
    letters = GET_NAMES(mapping);
    x = duplicate(x);
    PROTECT(x);
    xvec = R_ExternalPtrProtected(GET_SLOT(x, install("values")));
    if (xvec != R_NilValue) {
        SET_SLOT(x, install("values"),
                 R_MakeExternalPtr(NULL, xvec,
                                   R_ExternalPtrTag(GET_SLOT(x,
                                                             install("values")))));
        UNPROTECT(1);
        return x;
    }
    xvec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(xvec) != CHARSXP)
        error("Only character storage is supported now");

    n = LENGTH(xvec)-1;
    if (n > 0) {
        unsigned char revmap[256];
        unsigned char* src = (unsigned char*) CHAR(xvec)+1;
        unsigned char* dest;
        unsigned char A = 0;
        unsigned char C = 0;
        unsigned char G = 0;
        unsigned char T = 0;
        unsigned char gap = 0;
        SEXP ansvec, offsets, dim, xptr;
        int* start;
        int* end;
        int i, noffsets;

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

        offsets = GET_SLOT(x, install("offsets"));
        dim = GET_DIM(offsets);
        if (TYPEOF(offsets) != INTSXP || TYPEOF(dim) != INTSXP ||
            LENGTH(dim) != 2 || INTEGER(dim)[1] != 2)
            error("offsets slot of BioString must be integer matrix with two columns");
        noffsets = INTEGER(dim)[0];

        memset(revmap, 0, 256);
        for (i = 1; i < 32; i++) {
            if (i & A)
                revmap[i] |= T;
            if (i & C)
                revmap[i] |= G;
            if (i & T)
                revmap[i] |= A;
            if (i & G)
                revmap[i] |= C;
            if (i & gap)
                revmap[i] |= gap;
        }
        ansvec = allocString(LENGTH(xvec));
        PROTECT(ansvec);
        dest = (unsigned char*) CHAR(ansvec);
        for (i = 0; i < n; i++) {
            unsigned char v = revmap[src[i]];
            if (!v)
                error("unrecognized code: %d", src[i]);
            /* not dest[n-i+1] - skip one character in front */
            dest[n-i] = v;
        }
        start = INTEGER(offsets);
        end = start+noffsets;
        for (i = 0; i < noffsets; i++) {
            int tmp = end[i];
            if (tmp) {
                end[i] = n-start[i]+1;
                start[i] = n-tmp+1;
            }
        }
        R_SetExternalPtrProtected(GET_SLOT(x, install("values")),
                                  ansvec);
        xptr = R_MakeExternalPtr(NULL, ansvec, xvec);
        UNPROTECT(1);
        PROTECT(xptr);
        SET_SLOT(x, install("values"), xptr);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return x;
}

typedef Rboolean (BioStringCall_t)(unsigned char* str, int slen,
                                   int i, void* user_data);

static void
foreach_BioStringC(SEXP x, BioStringCall_t* f, void* user_data)
{
    SEXP xvec;
    int i;
    int* startvec;
    int* endvec;
    int len = getBioStringLength(x, &startvec, &endvec);
    unsigned char* seq;

    xvec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(xvec) != CHARSXP)
        error("Only character storage is supported now");
    seq = (unsigned char*) CHAR(xvec);

    for (i = 0; i < len ; i++) {
        int start = startvec[i];
        int end = endvec[i];
        if (!f(seq+start, end-start+1, i, user_data))
            break;
    }
}

typedef struct {
    int* ans;
    unsigned char c;
} allSameLetter_t;

static Rboolean
allSameLetter_func(unsigned char* str, int slen,
                   int j, allSameLetter_t* data)
{
    unsigned char c = data->c;
    int i;
    char savefirst = str[-1];
    str[-1] = ~c;
    for (i = slen-1; str[i] == c; i--) {
    }
    str[-1] = savefirst;
    data->ans[j] = (i == -1);
    return 1;
}

static SEXP
allSameLetter(SEXP x, SEXP pattern)
{
    int len = getBioStringLength(x, NULL, NULL);
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP mapping = GET_SLOT(alph, install("mapping"));
    SEXP letters = GET_NAMES(mapping);
    int nletters = LENGTH(letters);
    int start, end;
    SEXP ans;
    allSameLetter_t data;

    if (TYPEOF(mapping) != INTSXP || TYPEOF(letters) != STRSXP ||
        nletters == 0)
        error("invalid mapping");
    getLengthOneBioStringRange(pattern, &start, &end);
    if (start != end)
        error("pattern is not a single letter");
    pattern = R_ExternalPtrTag(GET_SLOT(pattern, install("values")));
    if (TYPEOF(pattern) != CHARSXP)
        error("can only handle character storage");
    data.c = ((unsigned char*)CHAR(pattern))[start];
    ans = allocVector(LGLSXP, len);
    PROTECT(ans);
    data.ans = INTEGER(ans);
    foreach_BioStringC(x, (BioStringCall_t*) allSameLetter_func,
                       &data);
    UNPROTECT(1);
    return ans;
}

static Rboolean
baseFrequency_func(unsigned char* str, int slen,
                   int j, int* counts)
{
    int i;
    for (i = 0; i < slen; i++) {
        counts[str[i]]++;
    }
    return 1;
}

static SEXP
baseFrequency(SEXP x)
{
    SEXP alph = GET_SLOT(x, install("alphabet"));
    SEXP mapping = GET_SLOT(alph, install("mapping"));
    SEXP letters = GET_NAMES(mapping);
    int nletters = LENGTH(letters);
    SEXP ans;
    int counts[256];
    int i;

    if (TYPEOF(mapping) != INTSXP || TYPEOF(letters) != STRSXP ||
        nletters == 0)
        error("invalid mapping");

    memset(counts, 0, 256*sizeof(int));
    foreach_BioStringC(x, (BioStringCall_t*) baseFrequency_func,
                       counts);

    ans = allocVector(INTSXP, nletters);
    PROTECT(ans);
    SET_NAMES(ans, letters);
    memset(INTEGER(ans), 0, nletters*sizeof(int));
    for (i = 0; i < nletters; i++) {
        unsigned int pat = ((unsigned int*)INTEGER(mapping))[i];
        if (pat >= (1U << CHAR_BIT))
            error("invalid mapping with character storage");
        INTEGER(ans)[i] = counts[pat];
        counts[pat] = 0;
    }
    for (i = 0; i < 256; i++) {
        if (counts[i] > 0)
            error("found pattern (%d) not in mapping", i);
    }
    UNPROTECT(1);
    return ans;
}

static Rboolean
suffixLess(unsigned char* start1, unsigned char* end1,
           unsigned char* start2, unsigned char* end2,
           int prefixLength)
{
    size_t len1 = end1-start1;
    size_t len2 = end2-start2;
    int cmplen = MIN(len1, len2);
    int j;
    if (cmplen >= prefixLength)
        cmplen = prefixLength-1;
    for (j = 0; j <= cmplen; j++) {
        if (start1[j] != start2[j]) {
            return start1[j] < start2[j];
        }
    }
    return (len1 < prefixLength && len1 < len2);
}

/*
 *  A simple insertion sort of suffix indexes. In each suffix, we only
 *  use the cmpStart to cmpEnd bytes.
 */
static void
insertionSortDNASuffixByPrefix(unsigned char* str, int len,
                               int* suffixIndex, int* startvec,
                               int* endvec, int prefixLength)
{
    int i;
    unsigned long interruptcheck = 0UL;
    for (i = len-1; i > 0; i--) {
        int s_i = suffixIndex[i];
        int s_im1 = suffixIndex[i-1];
        if (suffixLess(str+startvec[s_i], str+endvec[s_i],
                       str+startvec[s_im1], str+endvec[s_im1],
                       prefixLength)) {
            suffixIndex[i-1] = s_i;
            suffixIndex[i] = s_im1;
        }
    }
    for (i = 2; i < len; i++) {
        int j;
        int s_i = suffixIndex[i];
        int s_jm1 = suffixIndex[i-1];
        for (j = i;
             suffixLess(str+startvec[s_i], str+endvec[s_i],
                        str+startvec[s_jm1], str+endvec[s_jm1],
                        prefixLength);
             s_jm1 = suffixIndex[--j-1]) {
            suffixIndex[j] = s_jm1;
        }
        if (j != i) {
            suffixIndex[j] = s_i;
        }
#ifdef INTERRUPT_BIOSTRINGS
        interruptcheck += i;
        if (interruptcheck > (1UL << 10)) {
            R_CheckUserInterrupt();
            interruptcheck = 0UL;
        }
#endif
    }
}

/*
 *  Given a string of encoded patterns (str) of length len, return an
 *  integer vector of length len with starting index (with 1 based
 *  indexing) of suffixes of str which are sorted according to their
 *  first prefixLength elements in increasing order. For the purpose
 *  of sorting, any suffix with less than prefixLength is assumed to
 *  be padded with patterns which have no bits set (ie. they are zero).
 *
 *  We assume that the encoding uses 5 least significant bits to do
 *  the encoding and each byte has at least one of the five bits
 *  set. We also assume that the rest of the bits are not set.
 *  
 */
static void
sortDNASuffixByPrefix(unsigned char* str, int len, int* suffixIndex,
                      int* startvec, int* endvec, int prefixLength)
{
    int M, p, ncounts, i, K;
    int* counts;
    int* oldIndex;
    int* newIndex;

    if (len <= 0 || prefixLength <= 0)
        return;

    /*
     *  We take a two step approach. In the first step, we do an LSD
     *  sort on the most significant bytes. If the prefixLength is
     *  too large, we then do an insertion sort to finish the
     *  sorting.
     *  
     *  Knuth describes this type of sorting algorithm as LSD-first
     *  sort on most significant digits. See the end of secion 5.2.5
     *  of The Art of Compter Prgramming (volume 3, second edition)
     *  for details.
     */

    /*
     *  We use a simple formula that gives a close approximation to
     *  the optimum number of words to be combined. This is based on
     *  the table in Knuth's book in the previously mentioned section.
     * 
     *  We use at most two passes. The actual value of M=K/p and p are
     *  chosen so that M*p is approximately prefixLength. If
     *  prefixLength is too large, this is not possible and we settle
     *  on an insertion sort to finish the sorting.
     */

    M = (.6*log((double) len)+.7)/log(32.0)+.5;
    if (M <= 0)
        M = 1;
    if (M > 3)
        M = 3;
    if (M >= prefixLength) {
        M = prefixLength;
        p = 1;
    } else {
        p = 2.0*log((double)len)/(1.4*M)-12.0/M+0.5;
        if (p <= 1)
            p = 2;
        if (p*M > prefixLength)
            M = (prefixLength+p-1)/p;
    }
    K = p*M;

    if (suffixIndex) {
        PROTECT(R_NilValue);
        oldIndex = suffixIndex;
    } else {
        oldIndex = INTEGER(PROTECT(allocVector(INTSXP, len)));
        for (i = 0; i < len; i++)
            oldIndex[i] = i;
    }
    newIndex = INTEGER(PROTECT(allocVector(INTSXP, len)));

    ncounts = 1 << 5*M;
    counts = INTEGER(PROTECT(allocVector(INTSXP, ncounts)));
    for (p--; p >= 0; p--) {
        int j;
        int* tmp;

        memset(counts, 0, ncounts*sizeof(int));
        switch (M) {
        case 1:
            for (j = 0; j < len; j++) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                if (endvec[oi_j] < s_j) {
                    counts[0]++;
                } else {
                    counts[str[s_j]]++;
                }
            }
            break;
        case 2:
            for (j = 0; j < len; j++) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                int end = endvec[oi_j]-s_j;
                if (end > 0) {
                    counts[str[s_j+1]+(str[s_j] << 5)]++;
                } else if (end == 0) {
                    counts[str[s_j] << 5]++;
                } else {
                    counts[0]++;
                }
            }
            break;
        case 3:
            for (j = 0; j < len; j++) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                int end = endvec[oi_j]-s_j;
                if (end > 1) {
                    counts[str[s_j+2]+(str[s_j+1] << 5)+(str[s_j] << 10)]++;
                } else {
                    switch (end) {
                    case 1:
                        counts[(str[s_j+1] << 5)+(str[s_j] << 10)]++;
                        break;
                    case 0:
                        counts[str[s_j] << 10]++;
                        break;
                    default:
                        counts[0]++;
                        break;
                    }
                }
            }
            break;
        default:
            error("invalid value for M %d", M);
        }
        for (j = 1; j < ncounts; j++)
            counts[j] += counts[j-1];
        switch (M) {
        case 1:
            for (j = len-1; j >= 0; j--) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                if (endvec[oi_j] < s_j) {
                    newIndex[--counts[0]] = oldIndex[j];
                } else {
                    newIndex[--counts[str[s_j]]] = oldIndex[j];
                }
            }
            break;
        case 2:
            for (j = len-1; j >= 0; j--) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                int end = endvec[oi_j]-s_j;
                if (end > 0) {
                    newIndex[--counts[str[s_j+1]+(str[s_j] << 5)]] =
                        oldIndex[j];
                } else if (end == 0) {
                    newIndex[--counts[str[s_j] << 5]] = oldIndex[j];
                } else {
                    newIndex[--counts[0]] = oldIndex[j];
                }
            }
            break;
        case 3:
            for (j = len-1; j >= 0; j--) {
                int oi_j = oldIndex[j];
                int s_j = startvec[oi_j]+p;
                int end = endvec[oi_j]-s_j;
                if (end > 1) {
                    newIndex[--counts[str[s_j+2]+(str[s_j+1] << 5)+
                                      (str[s_j] << 10)]]
                        = oldIndex[j];
                } else {
                    switch (end) {
                    case 1:
                        newIndex[--counts[(str[s_j+1] << 5)+(str[s_j] << 10)]]
                            = oldIndex[j];
                        break;
                    case 0:
                        newIndex[--counts[str[s_j] << 10]] = oldIndex[j];
                        break;
                    default:
                        newIndex[--counts[0]] = oldIndex[j];
                        break;
                    }
                }
            }
            break;
        default:
            error("invalid value for M %d", M);
        }
        tmp = oldIndex;
        oldIndex = newIndex;
        newIndex = tmp;
    }
    if (prefixLength > K)
        insertionSortDNASuffixByPrefix(str, len, oldIndex,
                                       startvec, endvec, prefixLength);
    if (suffixIndex == NULL) {
        memcpy(newIndex, startvec, len*sizeof(int));
        for (i = 0; i < len; i++)
            startvec[i] = newIndex[oldIndex[i]];
        memcpy(newIndex, endvec, len*sizeof(int));
        for (i = 0; i < len; i++)
            endvec[i] = newIndex[oldIndex[i]];
    } else if (suffixIndex != oldIndex)
        memcpy(suffixIndex, oldIndex, len*sizeof(int));
    UNPROTECT(3);
}

static SEXP
DNASuffixArray(SEXP x, SEXP prefixLength)
{
    int* startvec;
    int* endvec;
    int len;
    int plen = asInteger(prefixLength);
    SEXP vec;

    vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(vec) != CHARSXP)
        error("values must be a CHARSXP");
    x = duplicate(x);
    PROTECT(x);
    len = getBioStringLength(x, &startvec, &endvec);
    if (len > 0) {
        int* offsets;
        int totlen = len;
        int i, k;

        for (i = 0; i < len; i++) {
            totlen += endvec[i]-startvec[i];
        }
        SET_SLOT(x, install("offsets"), allocMatrix(INTSXP, totlen, 2));
        offsets = INTEGER(GET_SLOT(x, install("offsets")));
        for (i = 0, k = 0; i < len; i++) {
            int j;
            int start = startvec[i];
            int end = endvec[i];
            int len_i = end-start+1;
            for (j = 0; j < len_i; j++, k++) {
                offsets[k] = end-j;
                offsets[k+totlen] = end;
            }
        }
        sortDNASuffixByPrefix(CHAR(vec),
                              totlen, NULL, offsets,
                              offsets+totlen, plen);
    }
    UNPROTECT(1);
    return x;
}

static SEXP
SortDNAString(SEXP x, SEXP prefixLength)
{
    int* startvec;
    int* endvec;
    int len;
    SEXP vec;
    int plen = asInteger(prefixLength);

    vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(vec) != CHARSXP)
        error("values must be a CHARSXP");
    x = duplicate(x);
    PROTECT(x);
    len = getBioStringLength(x, &startvec, &endvec);
    if (len > 0 && plen > 0) {
        sortDNASuffixByPrefix(CHAR(vec),
                              len, NULL, startvec,
                              endvec, plen);
    }
    UNPROTECT(1);
    return x;
}

static SEXP
longestCommonPrefixSuffixArray(SEXP x)
{
    int* startvec;
    int* endvec;
    int len;
    SEXP vec;
    SEXP ans;

    len = getBioStringLength(x, &startvec, &endvec);
    vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(vec) != CHARSXP)
        error("values must be a CHARSXP");
    ans = allocVector(INTSXP, len);
    if (len > 0) {
        int* lcp = INTEGER(ans);
        unsigned char* str;
        unsigned char* str_ip1;
        int len_ip1;
        int i;

        PROTECT(ans);
        lcp[len-1] = 0;
        str = (unsigned char*) CHAR(vec);
        str_ip1 = str+startvec[len-1];
        len_ip1 = endvec[len-1]-startvec[len-1]+1;
        for (i = len-2; i >= 0; i--) {
            int j = 0;
            unsigned char* str_i = str+startvec[i];
            int len_i = endvec[i]-startvec[i]+1;
            int minlen = MIN(len_ip1, len_i);

            while (j < minlen && str_ip1[j] == str_i[j])
                j++;
            lcp[i] = j;
            len_ip1 = len_i;
            str_ip1 = str_i;
        }
        UNPROTECT(1);
    }
    return ans;
}

#include "common.h"
#include <R.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef R_CallDef  [] = {
    CALL_DEF(IntegerBitOr, 1),
    CALL_DEF(BioStringValues, 2),
    CALL_DEF(setBioString, 2),
    CALL_DEF(BioString_substring, 4),
    CALL_DEF(BioStringToRString, 1),
    CALL_DEF(BoyerMoore_exactMatch, 2),
    CALL_DEF(ForwardSearch_exactMatch, 2),
    CALL_DEF(reverseComplementBioString, 1),
    CALL_DEF(allSameLetter, 2),
    CALL_DEF(baseFrequency, 1),
    CALL_DEF(DNASuffixArray, 2),
    CALL_DEF(SortDNAString, 2),
    CALL_DEF(longestCommonPrefixSuffixArray, 1),
    {NULL, NULL, 0},
};

void R_init_Biostrings(DllInfo *info)
{
    R_registerRoutines(info, NULL, R_CallDef, NULL, NULL);
}
