/****************************************************************************
                           THE SHIFT-OR ALGORITHM
                                   aka
                            THE BITAP ALGORITHM

                             Author: H. Pages

 This is a complete rewrite from Biostrings 1.

 References:
   - Wikipedia:
       http://en.wikipedia.org/wiki/Shift-or_algorithm
   - other:
       http://www-igm.univ-mlv.fr/~lecroq/string/node6.html
   - For all kind of string algorithms with animation in Java, see also:
       http://www-igm.univ-mlv.fr/~lecroq/string/index.html
   - The agrep homepage (DOS, Windows and OS/2 only):
       http://www.tgries.de/agrep

 Note that in Biostrings 1 (1.4.0) the "shift-or" algo was almost always
 returning wrong results for 'max_mismatch' != 0 ("inexact searching").
 The bug was in the C function ShiftOr_matchInternal() which was using
 the following formula for computing the "pattern bitmasks":
 
         ShiftOrWord_t tmp = M_k[0];
	 M_k[0] = (tmp << 1) | U[xptr[k]];
	 M_k[1] = ((M_k[1] << 1) | U[xptr[k]]) & tmp & (tmp << 1) & (M_k[0] << 1);
 
 where the correct formula would have been to replace the last line by:
 
         M_k[1] = ((M_k[1] << 1) | U[xptr[k]]) & (tmp << 1) & (M_k[0]);

 Also note that using the formula given in D. Gusfield's book (Algorithms on
 strings, trees, and sequences), p. 73, would result in doing:
 
         M_k[1] = ((M_k[1] << 1) | U[xptr[k]]) & (tmp) & (M_k[0]);

 which is wrong too. 
 
 In this complete rewrite of the "shift-or" algo, we use a more general
 version of the correct formula. The code doing this is in the
 _update_PMmasks() function:
 
	 for (e = 1; e < PMmask_length; e++) {
	     PMmaskB = PMmaskA;
	     PMmaskA = PMmask[e] >> 1;
	     PMmask[e] = (PMmaskA | pmask) & PMmaskB & PMmask[e-1];
	 }

 Note that the bit order in those bitmasks is reverted relatively to the
 order used in Biostrings 1. In Biostrings 1 the first position in the
 pattern (the most left) was mapped to the most right bit (the weak bit)
 of the bitmask. Now it is the last position in the pattern (the most right)
 that is mapped to the most right bit of the bitmask which is more intuitive.
 
 In addition to the "inexact matching" bug previously described, the following
 problems have been addressed:
   - The 'kerr' arg of ShiftOr_matchInternal() (nb of mismatches) is no
     longer required to be <= 3. Now it can be whatever positive int. 
   - The "border problem" occuring with "inexact matching" is fixed.
     This new version of the "shift-or" algo should find ALL matches, even
     those that are out of bounds (i.e. that start before the first or end
     after the last character of the subject).
   - The bus-error that occured on Solaris 2.9 and Mac OS X for:
       > pattern <- DNAString("AAAA")
       > subject <- DNAString("AAAAC")
       > matchDNAPattern(pattern, subject, mis=2)
     is gone.

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

SEXP match_shiftor_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_shiftor.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_shiftor.c'\n");
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

SEXP bits_per_long()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = BITS_PER_LONG;
	UNPROTECT(1);
	return ans;
}


/****************************************************************************/

static void _set_pmaskmap(
		int is_fixed,
		int pmaskmap_length,
		ShiftOrWord_t *pmaskmap,
		int nP,
		const char *pat)
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
		for (i = 0; i < nP; i++) {
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
 * Returns nb of mismatches (>= 0) if an inexact match is found.
 */
static int _next_match(
		int *Lpos,
		int *Rpos,
		int nS,
		const char *subj,
		ShiftOrWord_t *pmaskmap,
		int PMmask_length, /* PMmask_length = kerr+1 */
		ShiftOrWord_t *PMmask)
{
	static ShiftOrWord_t pmask;
	static int nncode;
	static int e;

	while (*Lpos < nS) {
		if (*Rpos < nS) {
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

static void _match_shiftor(const char *P, int nP, const char *S, int nS, 
		int PMmask_length, int is_fixed)
{
	ShiftOrWord_t *PMmask, pmaskmap[256];
	int i, e, Lpos, Rpos, ret;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _match_shiftor(): BEGIN\n");
	}
#endif
	if (nP <= 0)
		error("empty pattern");
	_set_pmaskmap(is_fixed, 256, pmaskmap, nP, P);
	PMmask = (ShiftOrWord_t *)
			R_alloc(PMmask_length, sizeof(ShiftOrWord_t));
	PMmask[0] = 1UL;
	for (i = 1; i < nP; i++) {
		PMmask[0] <<= 1;
		PMmask[0] |= 1UL;
	}
	for (e = 1; e < PMmask_length; e++) {
		PMmask[e] = PMmask[e-1] >> 1;
	}
	Lpos = 1 - nP;
	Rpos = 0;
	while (1) {
		ret = _next_match(
			&Lpos,
			&Rpos,
			nS,
			S,
			pmaskmap,
			PMmask_length,
			PMmask);
		if (ret == -1) {
			break;
		}
		_Biostrings_report_match(Lpos - 1, 0);
	}
	/* No need to free PMmask, R does that for us */
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _match_shiftor(): END\n");
	}
#endif
	return;
}

/*
 * Arguments:
 *   'p_xp': pattern@xdata@xp
 *   'p_offset': pattern@offset
 *   'p_length': pattern@length
 *   's_xp': subject@xdata@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *   'max_mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 * Return an integer vector containing the relative pos of the matches.
 * All matches have the length of the pattern.
 */
SEXP match_shiftor(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP max_mismatch, SEXP fixed, SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    kerr, fixedP, fixedS, is_count_only;
	const Rbyte *pat, *subj;
	SEXP ans;

	pat_length = INTEGER(p_length)[0];
	if (pat_length > shiftor_maxbits)
		error("pattern is too long");
	pat_offset = INTEGER(p_offset)[0];
	pat = RAW(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_length = INTEGER(s_length)[0];
	subj_offset = INTEGER(s_offset)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	kerr = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (fixedP != fixedS)
		error("fixedP != fixedS not yet supported");
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	_match_shiftor((char *) pat, pat_length, (char *) subj, subj_length,
		       kerr+1, fixedP);
	if (is_count_only)
		PROTECT(ans = _Biostrings_viewsbuf_count_asINTEGER());
	else
		PROTECT(ans = _Biostrings_viewsbuf_start_asINTEGER());
	UNPROTECT(1);
	return ans;
}
