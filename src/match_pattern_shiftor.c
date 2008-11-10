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
 update_PMmasks() function:
 
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

SEXP debug_match_pattern_shiftor()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pattern_shiftor.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pattern_shiftor.c'\n");
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

static void set_pmaskmap(
		int is_fixed,
		int pmaskmap_length,
		ShiftOrWord_t *pmaskmap,
		const RoSeq *P)
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
		for (i = 0; i < P->nelt; i++) {
			pmask <<= 1;
			if (is_fixed) {
				if (((unsigned char) P->elts[i]) != nncode)
					pmask |= 1UL;
			} else {
				if ((((unsigned char) P->elts[i]) & nncode) == 0)
					pmask |= 1UL;
			}
		}
		pmaskmap[nncode] = pmask;
	}
	return;
}

static void update_PMmasks(
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
		Rprintf("[DEBUG] update_PMmasks(): PMmask[%d]=", 0);
		debug_printULBits(PMmask[0]);
	}
#endif
	for (e = 1; e < PMmask_length; e++) {
		PMmaskB = PMmaskA;
		PMmaskA = PMmask[e] >> 1;
		PMmask[e] = (PMmaskA | pmask) & PMmaskB & PMmask[e-1];
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] update_PMmasks(): PMmask[%d]=", e);
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
static int next_match(
		int *Lpos,
		int *Rpos,
		const RoSeq *S,
		ShiftOrWord_t *pmaskmap,
		int PMmask_length, /* PMmask_length = kerr+1 */
		ShiftOrWord_t *PMmask)
{
	static ShiftOrWord_t pmask;
	static int nncode;
	static int e;

	while (*Lpos < S->nelt) {
		if (*Rpos < S->nelt) {
			nncode = (unsigned char) S->elts[*Rpos];
			pmask = pmaskmap[nncode];
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] next_match(): ");
				Rprintf("pmaskmap[%d]=", nncode);
				debug_printULBits(pmask);
			}
#endif
		} else {
			pmask = ~0UL;
		}
		update_PMmasks(PMmask_length, PMmask, pmask);
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

static void shiftor(const RoSeq *P, const RoSeq *S, int PMmask_length, int is_fixed)
{
	ShiftOrWord_t *PMmask, pmaskmap[256];
	int i, e, Lpos, Rpos, ret;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] shiftor(): BEGIN\n");
	}
#endif
	if (P->nelt <= 0)
		error("empty pattern");
	set_pmaskmap(is_fixed, 256, pmaskmap, P);
	PMmask = (ShiftOrWord_t *)
			R_alloc(PMmask_length, sizeof(ShiftOrWord_t));
	PMmask[0] = 1UL;
	for (i = 1; i < P->nelt; i++) {
		PMmask[0] <<= 1;
		PMmask[0] |= 1UL;
	}
	for (e = 1; e < PMmask_length; e++) {
		PMmask[e] = PMmask[e-1] >> 1;
	}
	Lpos = 1 - P->nelt;
	Rpos = 0;
	while (1) {
		ret = next_match(
			&Lpos,
			&Rpos,
			S,
			pmaskmap,
			PMmask_length,
			PMmask);
		if (ret == -1) {
			break;
		}
		_report_match(Lpos, P->nelt);
	}
	/* No need to free PMmask, R does that for us */
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] shiftor(): END\n");
	}
#endif
	return;
}

void _match_pattern_shiftor(const RoSeq *P, const RoSeq *S,
		int max_mm, int fixedP, int fixedS)
{
	if (P->nelt > shiftor_maxbits)
		error("pattern is too long");
	if (fixedP != fixedS)
		error("fixedP != fixedS not supported by shift-or algo");
	shiftor(P, S, max_mm + 1, fixedP);
}

