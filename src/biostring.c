/* Copyright (C) 2003 by Saikat DebRoy */
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

static SEXP
expandIndex(SEXP index, int ndone, int nleft)
{
    int n = LENGTH(index);
    double proportion = (n+1)/(double)(ndone);
    int estimate = proportion*nleft+1;
    n = 2*(n+estimate);
    return lengthgets(index, n);
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

SEXP
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
    if (NAMED(biostring))
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

SEXP
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

SEXP
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

    if (NAMED(x))
        ans = duplicate(x);
    else ans = x;
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
    PROTECT(matchIndex);
    if (NAMED(x))
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
    UNPROTECT(2);
    return x;
}

SEXP
ForwardSearch_exactMatch(SEXP pattern, SEXP x)
{
    int pstart, pend, xstart, xend, patlen = 0;
    SEXP alph, vec;
    SEXP matchIndex = R_NilValue;
    PROTECT_INDEX matchIndex_pi;
    int* index = NULL;
    int nmatchIndex, nletters;
    int nmatch = 0;
#ifdef PROGRESS_BIOSTRINGS
    int twopercent, progresscheck;
#endif
    int m;

    getLengthOneBioStringRange(pattern, &pstart, &pend);
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
    PROTECT_WITH_INDEX(matchIndex, &matchIndex_pi);
    index = INTEGER(matchIndex);

    nmatch = 0;
#ifdef PROGRESS_BIOSTRINGS
    progresscheck = twopercent = (m-patlen+1)/50;
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
                    matchIndex = expandIndex(matchIndex, k-patlen+1, m-k);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = k;
            }
#ifdef PROGRESS_BIOSTRINGS
            if (k > progresscheck) {
                R_CheckUserInterrupt();
                Rprintf(".");
                progresscheck += twopercent;
            }
#endif
        }
#ifdef PROGRESS_BIOSTRINGS
        Rprintf("\n");
#endif
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
                    matchIndex = expandIndex(matchIndex, k-patlen+1, m-k);
                    REPROTECT(matchIndex, matchIndex_pi);
                    nmatchIndex = LENGTH(matchIndex);
                    index = INTEGER(matchIndex);
                }
                index[nmatch++] = k;
            }
#ifdef PROGRESS_BIOSTRINGS
            if (k > progresscheck) {
                R_CheckUserInterrupt();
                Rprintf(".");
                progresscheck += twopercent;
            }
#endif
        }
#ifdef PROGRESS_BIOSTRINGS
        Rprintf("\n");
#endif
        pptr[0] = save_first;
    }
    UNPROTECT(1);
finished_match:
    return matchIndexToBioString(x, matchIndex, nmatch, patlen);
}

SEXP
LengthOne_exactMatch(SEXP pattern, SEXP x)
{
    int start, end;
    SEXP alph;
    SEXP matchIndex = R_NilValue;
    PROTECT_INDEX matchIndex_pi;
    int* index = NULL;
    int nmatchIndex, nletters;
    int nmatch = 0;
#ifdef PROGRESS_BIOSTRINGS
    int twopercent, progresscheck;
#endif
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
#ifdef PROGRESS_BIOSTRINGS
    progresscheck = twopercent = m/50;
#endif
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
#ifdef PROGRESS_BIOSTRINGS
            if (i > progresscheck) {
                R_CheckUserInterrupt();
                Rprintf(".");
                progresscheck += twopercent;
            }
#endif
        }
#ifdef PROGRESS_BIOSTRINGS
        Rprintf("\n");
#endif
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
#ifdef PROGRESS_BIOSTRINGS
            if (i > progresscheck) {
                R_CheckUserInterrupt();
                Rprintf(".");
                progresscheck += twopercent;
            }
#endif
        }
#ifdef PROGRESS_BIOSTRINGS
        Rprintf("\n");
#endif
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
        Nvecs = reverseFundamentalPreprocessing(patternptr, n);
        PROTECT(Nvecs);
        N = INTEGER(VECTOR_ELT(Nvecs, 0))-1;
        memset(tmpptr, 0, (1+n)*sizeof(int));
        for (i = 1; i < n; i++) {
            int j = n-N[i]+1;
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

SEXP
BoyerMoore_exactMatch(SEXP origPattern, SEXP x)
{
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
    vec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (pattern.usesChar) {
        unsigned char* str;
        char* patternptr;
        char save_first;
        int m = xend-xstart+1;
        int nmatchIndex;
#ifdef PROGRESS_BIOSTRINGS
        int twopercent = (m-patlen+1)/50;
        int progresscheck = twopercent;
#endif
        if (TYPEOF(vec) != CHARSXP)
            error("type mismatch between pattern and string");

        str = (unsigned char*) CHAR(vec)+xstart-1;
        PROTECT(pattern.bad_char.letterIndex);
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
                SEXP letterIndex =
                    VECTOR_ELT(pattern.bad_char.letterIndex, str[h]);
                int bad_character_shift;
                R_assert(letterIndex != R_NilValue);
                bad_character_shift = i-INTEGER(letterIndex)[i];
#if 0
                if (i == patlen) {
                    /* good suffix rule gives 1 here - so we ignore it */
                    k += bad_character_shift;
                } else {
                    int good_suffix_shift = pattern.good_suffix_shift[i+1];
                    k += MAX(good_suffix_shift, bad_character_shift);
                }
#else
                k += bad_character_shift;
#endif
            }
#ifdef PROGRESS_BIOSTRINGS
            if (k > progresscheck) {
                R_CheckUserInterrupt();
                Rprintf(".");
                progresscheck += twopercent;
            }
#endif
        }
#ifdef PROGRESS_BIOSTRINGS
        Rprintf("\n");
#endif
        patternptr[0] = save_first;
        UNPROTECT(2);
    } else {
        error("non-character patterns and strings unimplemented");
    }
    /* BoyerMoore_releasePattern(&pattern); */
finished_match:
    return matchIndexToBioString(x, matchIndex, nmatch, patlen);
}

SEXP
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
    xvec = R_ExternalPtrTag(GET_SLOT(x, install("values")));
    if (TYPEOF(xvec) != CHARSXP)
        error("Only character storage is supported now");

    n = LENGTH(xvec)-1;
    if (NAMED(x))
        x = duplicate(x);
    PROTECT(x);
    if (n > 0) {
        unsigned char revmap[256];
        unsigned char* src = (unsigned char*) CHAR(xvec)+1;
        unsigned char* dest;
        unsigned char A = 0;
        unsigned char C = 0;
        unsigned char G = 0;
        unsigned char T = 0;
        unsigned char gap = 0;
        SEXP ansvec, offsets, dim;
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
        R_SetExternalPtrTag(GET_SLOT(x, install("values")), ansvec);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return x;
}
