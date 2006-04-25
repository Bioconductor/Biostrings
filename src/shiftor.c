/****************************************************************************
                           THE SHIFT-OR ALGORITHM
                                   aka
                            THE BITAP ALGORITHM
 On Wikipedia:
   http://en.wikipedia.org/wiki/Shift-or_algorithm
 Other resources:
   http://www-igm.univ-mlv.fr/~lecroq/string/node6.html
 For all kind of string algorithms with animation in Java, see also:
   http://www-igm.univ-mlv.fr/~lecroq/string/index.html

 The agrep homepage (DOS, Windows and OS/2 only):
   http://www.tgries.de/agrep

 ****************************************************************************/
#include "Biostrings.h"
#include <limits.h>
#include <Rinternals.h>

#define CHAR_SIZE               (sizeof(char))
#define LONG_SIZE               (sizeof(long))
#define BITS_PER_CHAR           CHAR_BIT
#define BITS_PER_SIZEOF_UNIT    (BITS_PER_CHAR / (int) CHAR_SIZE)
#define BITS_PER_LONG           ((int) LONG_SIZE * BITS_PER_SIZEOF_UNIT)

/*
 * Expected to be 32-bit on 32-bit machines and 64-bit on 64-bit machines
 */
typedef unsigned long ShiftOrWord_t;
int shiftor_maxbits = sizeof(ShiftOrWord_t) * CHAR_BIT;

/****************************************************************************/
static int debug = 0;

SEXP shiftor_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'shiftor.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'shiftor.c'\n");
#endif
	return R_NilValue;
}

#ifdef DEBUG_BIOSTRINGS
static void debug_printULBits(unsigned long bits)
{
	unsigned long current_bit = 1UL << (BITS_PER_LONG-1);
	int i;

	for (i = 0; i < BITS_PER_LONG; i++) {
		printf("%d", (bits & current_bit) != 0UL);
		if ((i % 8) == 7) {
			printf(" ");
		}
		current_bit >>= 1;
	}
	printf("-> %lu\n", bits);
	return;
}
#endif


/****************************************************************************/
/* A naive estimate */
static int estimateMatchNumber(int pat_length, int subj_length)
{
	int matchpos_length, tmp;

	matchpos_length = 1;
	if (pat_length >= shiftor_maxbits)
		pat_length--;
	tmp = 1 << pat_length;
	if (tmp < subj_length)
		matchpos_length = subj_length / tmp;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] estimateMatchNumber(): ");
		Rprintf("matchpos_length=%d\n", matchpos_length);
	}
#endif
	return matchpos_length;
}

static SEXP expandIndex(SEXP index, int ndone, int nleft)
{
	int n = LENGTH(index);
	int n1;
	double proportion = (n+1)/(double)(ndone);
	int estimate = proportion*nleft+1;
	SEXP temp;

	n1 = 2*(n+estimate);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] expandIndex(): ");
		Rprintf("ndone=%d nleft=%d n=%d n1=%d\n", ndone, nleft, n, n1);
	}
#endif
	temp = allocVector(INTSXP, n1);
	memcpy(INTEGER(temp), INTEGER(index), n*sizeof(int));
	return temp;
}

static void _set_pmaskmap(
		int is_fixed,
		int pmaskmap_length,
		ShiftOrWord_t *pmaskmap,
		int pat_length,
		char *pat)
{
	ShiftOrWord_t pmask;
	int nncode, i;

	/* Why go to 255? Only pmaskmap[nncode] will be used,
	where nncode is a numerical nucleotide code.
	nncode can only have 16 possible values: 1, 2, 4, 6, ..., 30.
	Not even all values <= 30 are used!
	*/
	for (nncode = 0; nncode < pmaskmap_length; nncode++) {
		pmask = 0LU;
		for (i = 0; i < pat_length; i++) {
			pmask <<= 1;
			if (is_fixed) {
				if (((unsigned char) pat[i]) != nncode)
					pmask |= 1UL;
			} else {
				if ((((unsigned char) pat[i]) & nncode) == 0)
					pmask |= 1UL;
			}
		}
		pmaskmap[nncode] = pmask;
	}
	return;
}

static void _update_PMmasks(
		int PMmask_length,
		ShiftOrWord_t *PMmask,
		ShiftOrWord_t pmask)
{
	static ShiftOrWord_t PMmaskA, PMmaskB;
	static int e;

	PMmaskA = PMmask[0] >> 1;
	PMmask[0] = PMmaskA | pmask;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _update_PMmasks(): PMmask[%d]=", 0);
		debug_printULBits(PMmask[0]);
	}
#endif
	for (e = 1; e < PMmask_length; e++) {
		PMmaskB = PMmaskA;
		PMmaskA = PMmask[e] >> 1;
		PMmask[e] = (PMmaskA | pmask) & PMmaskB & PMmask[e-1];
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] _update_PMmasks(): PMmask[%d]=", e);
			debug_printULBits(PMmask[e]);
		}
#endif
	}
	return;
}

/*
 * Returns -1 if no match is found.
 * Returns nb of mismatches (>= 0) if a fuzzy match is found.
 */
static int _next_match(
		int *Lpos,
		int *Rpos,
		int subj_length,
		char *subj,
		ShiftOrWord_t *pmaskmap,
		int PMmask_length, /* PMmask_length = kerr+1 */
		ShiftOrWord_t *PMmask)
{
	static ShiftOrWord_t pmask;
	static int nncode;
	static int e;

	while (*Lpos < subj_length) {
		if (*Rpos < subj_length) {
			nncode = (unsigned char) subj[*Rpos];
			pmask = pmaskmap[nncode];
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] _next_match(): ");
				Rprintf("pmaskmap[%d]=", nncode);
				debug_printULBits(pmask);
			}
#endif
		} else {
			pmask = ~0UL;
		}
		_update_PMmasks(PMmask_length, PMmask, pmask);
		(*Lpos)++;
		(*Rpos)++;
		for (e = 0; e < PMmask_length; e++) {
			if ((PMmask[e] & 1UL) == 0UL) {
				return e;
			}
		}
	}
	return -1;
}

static void _shiftor(
		int is_fixed,
		int PMmask_length,
		int pat_length,
		char *pat,
		int subj_length,
		char *subj,
		int *p_nmatch,
		SEXP *p_matchpos,
		PROTECT_INDEX matchpos_pi)
{
	ShiftOrWord_t *PMmask, pmaskmap[256];
	int i, e, Lpos, Rpos, ret;
	int *index;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _shiftor(): BEGIN\n");
	}
#endif
	_set_pmaskmap(is_fixed, 256, pmaskmap, pat_length, pat);
	PMmask = (ShiftOrWord_t *)
			R_alloc(PMmask_length, sizeof(ShiftOrWord_t));
	PMmask[0] = 1UL;
	for (i = 1; i < pat_length; i++) {
		PMmask[0] <<= 1;
		PMmask[0] |= 1UL;
	}
	for (e = 1; e < PMmask_length; e++) {
		PMmask[e] = PMmask[e-1] >> 1;
	}
	if (*p_matchpos != R_NilValue) {
		index = INTEGER(*p_matchpos);
	}
	*p_nmatch = 0;
	Lpos = 1 - pat_length;
	Rpos = 0;
	while (1) {
		ret = _next_match(
			&Lpos,
			&Rpos,
			subj_length,
			subj,
			pmaskmap,
			PMmask_length,
			PMmask);
		if (ret == -1) {
			break;
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] _shiftor(): ");
			Rprintf("match found for Lpos=%d Rpos=%d\n",
				Lpos-1, Rpos-1);
		}
#endif
		if (*p_matchpos != R_NilValue) {
			if (*p_nmatch == LENGTH(*p_matchpos)) {
				*p_matchpos = expandIndex(*p_matchpos, Rpos,
							  subj_length - Lpos);
				REPROTECT(*p_matchpos, matchpos_pi);
				index = INTEGER(*p_matchpos);
			}
			index[*p_nmatch] = Lpos - 1;
		}
		(*p_nmatch)++;
	}
	/* No need to free PMmask, R does that for us */
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _shiftor(): END\n");
	}
#endif
	return;
}

/*
 * Arguments:
 *   'p_xp': pattern@data@xp
 *   'p_offset': pattern@offset
 *   'p_length': pattern@length
 *   's_xp': subject@data@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *   'mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': single logical
 *   'count_only': single logical
 * Return an integer vector containing the relative pos of the matches.
 * All matches have the length of the pattern.
 */
SEXP shiftor(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP mismatch, SEXP fixed, SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    kerr, is_fixed, is_count_only, ans_length;
	char *pat, *subj;
	SEXP ans;
	int matchpos_length;
	SEXP matchpos = R_NilValue;
	PROTECT_INDEX matchpos_pi;

	pat_length = INTEGER(p_length)[0];
	if (pat_length > shiftor_maxbits)
		error("pattern is too long");
	pat_offset = INTEGER(p_offset)[0];
	pat = CHAR(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_length = INTEGER(s_length)[0];
	subj_offset = INTEGER(s_offset)[0];
	subj = CHAR(R_ExternalPtrTag(s_xp)) + subj_offset;
	kerr = INTEGER(mismatch)[0];
	is_fixed = LOGICAL(fixed)[0];
	is_count_only = LOGICAL(count_only)[0];
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] shiftor(): BEGIN\n");
		Rprintf("[DEBUG] shiftor(): ");
		Rprintf("pat_offset=%d pat_length=%d pat[0]=%d\n",
			pat_offset, pat_length, (unsigned char) pat[0]);
		Rprintf("subj_offset=%d subj_length=%d subj[0]=%d\n",
			subj_offset, subj_length, (unsigned char) subj[0]);
		Rprintf("kerr=%d is_fixed=%d is_count_only=%d\n",
			kerr, is_fixed, is_count_only);
	}
#endif
	if (!is_count_only) {
		matchpos_length = estimateMatchNumber(pat_length, subj_length);
		PROTECT_WITH_INDEX(matchpos, &matchpos_pi);
		matchpos = allocVector(INTSXP, matchpos_length);
		REPROTECT(matchpos, matchpos_pi);
	}
	_shiftor(
		is_fixed,
		kerr+1,
		pat_length,
		pat,
		subj_length,
		subj,
		&ans_length,
		&matchpos,
		matchpos_pi
	);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] shiftor(): END\n");
	}
#endif
	if (is_count_only) {
		PROTECT(ans = allocVector(INTSXP, 1));
		INTEGER(ans)[0] = ans_length;
		UNPROTECT(1);
	} else {
		PROTECT(ans = allocVector(INTSXP, ans_length));
		memcpy(INTEGER(ans), INTEGER(matchpos),
			sizeof(int) * ans_length);
		UNPROTECT(2);
	}
	return ans;
}
