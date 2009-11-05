/****************************************************************************
                      A BOYER-MOORE-LIKE MATCHING ALGO
		            Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>

/* I've turned off the "MWshift feature". Doesn't seem to give very good
 * results :-( To turn it on, set MWSHIFT_NPMAX to a positive integer
 * (typically 128, 256 or 512, doesn't have to be a power of 2).
 * The original idea of this feature was to try to boost the boyermoore()
 * function when the pattern is small.
 */
#define MWSHIFT_NPMAX	0


/****************************************************************************/
static int debug = 0;

SEXP debug_match_pattern_boyermoore()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pattern_boyermoore.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pattern_boyermoore.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * The 'ppP' (Preprocessed Pattern) struct holds a copy of the current
 * pattern (eventually reverted if init_ppP_seq() was called with
 * walk_backward = 1) + the results of some preprocessing operations on it.
 * IMPORTANT: All members in 'ppP' that point to dynamically allocated
 * memory must be *persistent* buffers so they must point to user-controlled
 * memory (i.e. memory that is not reclaimed by R at the end of the .Call()
 * call). Hence the use of malloc()/free() instead of Salloc() for memory
 * allocation.
 * Members of 'ppP' are:
 *   buflength: the size of the buffer pointed by the 'seq' member, which, in
 *              the current implemenation, is also the length of the longest
 *              pattern seen so far (i.e. since the beginning of the current
 *              R session);
 *   seq: the letters of the current pattern (eventually in reverse order
 *              if init_ppP_seq() was called with walk_backward = 1);
 *   seqlength: the length of the current pattern (must be <= 'buflength');
 *   LCP: see init_ppP_seq() below;
 *   j0, shift0: see "j0/shift0" section below;
 *   VSGSshift_table: see "The Very Strong Good Suffix shifts" section below;
 *   MWshift_table: see "The Matching Window shifts" section below.
 */
static struct {
	int buflength;
	char *seq;
	int seqlength;
	int LCP;
	int j0, shift0;
	int *VSGSshift_table;
	int *MWshift_table;
} ppP = {0, NULL, 0, -1, 0, 0, NULL, NULL};

/* The 'LCP' member:
 *     -1: init_ppP_seq() changed the value of ppP.buflength.
 *   >= 0: init_ppP_seq() didn't change the value of ppP.buflength.
 *         The non-negative integer is the length of the Longest Common
 *         Prefix between old and new current pattern (LCP will always be <=
 *         min(P->length, ppP.seqlength)).
 */
static void init_ppP_seq(const cachedCharSeq *P, int walk_backward)
{
	int LCP, j1, j2;
	char c;

	if (P->length == 0) { /* should never happen but safer anyway... */
		ppP.LCP = 0;
		return;
	}
	if (P->length > 20000)
		error("pattern is too long");
	if (P->length > ppP.buflength) {
		/* We need to extend the size of 'ppP'. In that case, we
		   don't need to compute the LCP and we set it to -1. */
		if (ppP.seq != NULL)
			free(ppP.seq);
		ppP.buflength = 0;
		ppP.seq = (char *) malloc(P->length * sizeof(char));
		if (ppP.seq == NULL)
			error("can't allocate memory for ppP.seq");
		ppP.buflength = P->length;
		LCP = -1;
	} else {
		/* We don't need to extend the size of 'ppP'. In that case,
		   we compute the LCP. */
		LCP = 0;
	}
	for (j1 = 0, j2 = P->length - 1; j1 < P->length; j1++, j2--) {
		c = P->seq[walk_backward ? j2 : j1];
		if (LCP != -1 && j1 < ppP.seqlength && c == ppP.seq[j1])
			LCP++;
		else
			ppP.seq[j1] = c;
	}
	ppP.seqlength = P->length;
	ppP.LCP = LCP;
	return;
}

/****************************************************************************
 * j0/shift0
 * =========
 * 
 * For any j1 and j2 such that 0 <= j1 < j2 <= nP (nP being the length of
 * the pattern P), we note P(j1, j2) the subpattern of P starting at
 * position j1 and ending at position j2-1.
 * Note that with this definition: P(0, nP) = P.
 *
 * Definition of j0 and shift0
 * ---------------------------
 *   j0 = start position in P of the smallest suffix of P that is not a
 *        substring of P(1, nP-1)
 * Note that j0 can also be defined by using the MWshift function (see
 * section "The Matching Window shifts" below):
 *   j0 = max{j1 | MWshift(j1, nP) >= j}
 *
 *  shift0 =  MWshift(j0, nP)
 *
 * Properties
 * ----------
 *
 *   (a) If j1 <= j0, then MWshift(j1, nP) = shift0 >= j0
 *   (b) If j1 > j0, then MWshift(j1, nP) < j1
 *   (c) MWshift(j1, j2) <= shift0
 *   (d) If j < j0, VSGSshift(c, j) = shift0
 *   (e) VSGSshift(P[0], 0) = shift0
 */

static void init_ppP_j0shift0()
{
	int j0, shift0, length, j;

	length = 1;
	j0 = ppP.seqlength - 1;
	for (j = j0 - 1; j >= 1; j--) {
		if (memcmp(ppP.seq + j, ppP.seq + j0, length) == 0) {
			length++;
			j0--;
		}
	}
	for (shift0 = j0 - j; shift0 < ppP.seqlength; shift0++, length--) {
		if (memcmp(ppP.seq, ppP.seq + shift0, length) == 0)
			break;
	}
	ppP.j0 = j0;
	ppP.shift0 = shift0;
	/*Rprintf("j0=%d shift0=%d\n", j0, shift0);*/
}


/****************************************************************************
 * The Very Strong Good Suffix shifts
 * ==================================
 *
 * ppP.VSGSshift_table is a 256 x ppP.buflength matrix.
 * Its layout is (only the values marked with an "x" will be potentially
 * used):
 *
 *           0 1 2 3 4 5 j
 *         0 x x x x - -
 *         1 x x x x - - 
 *         2 x x x x - -    ppP.seqlength = 4 <= ppP.buflength = 6
 *         .............
 *       256 x x x x - -
 *         c
 *
 * The "x" region is defined by 0 <= j < ppP.seqlength
 */

#define VSGS_SHIFT(c, j) (ppP.VSGSshift_table[ppP.buflength * ((unsigned char) (c)) + (j)])

static int get_VSGSshift(char c, int j)
{
	int shift, k, k1, k2, length;
	const char *tmp;

	if (j < ppP.j0)
		return ppP.shift0;
	shift = VSGS_SHIFT(c, j);
	if (shift != 0)
		return shift;
	for (shift = 1; shift < ppP.seqlength; shift++) {
		if (shift <= j) {
			k = j - shift;
			if (ppP.seq[k] != c)
				continue;
			k1 = k + 1;
		} else {
			k1 = 0;
		}
		k2 = ppP.seqlength - shift;
		if (k1 == k2)
			break;
		length = k2 - k1;
		tmp = ppP.seq + k1;
		if (memcmp(tmp, tmp + shift, length) == 0)
			break;
	}
	/* shift is ppP.seqlength when the "for" loop is not interrupted by "break" */
	/*Rprintf("VSGSshift(c=%c, j=%d) = %d\n", c, j, shift);*/
	return VSGS_SHIFT(c, j) = shift;
}

static void init_ppP_VSGSshift_table()
{
	int u, j;
	char c;

	if (ppP.LCP == -1 && ppP.VSGSshift_table != NULL) {
		free(ppP.VSGSshift_table);
		ppP.VSGSshift_table = NULL;
	}
	if (ppP.buflength != 0 && ppP.VSGSshift_table == NULL) {
		ppP.VSGSshift_table = (int *)
			malloc(256 * ppP.buflength * sizeof(int));
		if (ppP.VSGSshift_table == NULL)
			error("can't allocate memory for ppP.VSGSshift_table");
	}
	for (u = 0; u < 256; u++) {
		for (j = 0; j < ppP.seqlength; j++) {
			c = (char) u;
			VSGS_SHIFT(c, j) = 0;
		}
	}
}


/****************************************************************************
 * The Matching Window shifts
 * ==========================
 *
 * Definition
 * ----------
 *
 * MWshift(j1, j2), the "Matching Window shift" for P(j1, j2), is defined by:
 *   Imagine that, for a given alignement of P with the subject S, all
 *   letters in subpattern P(j1, j2) are matching S (we say that the current
 *   Matching Window is (j1, j2)). Then MWshift(j1, j2) is the smallest
 *   (non-zero) amount of letters that we can shift P to the right without
 *   introducing an immediate mismatch.
 *
 * How to compute the MWshift function
 * -----------------------------------
 * 
 * The MWshift function depends only on the pattern and therefore can be
 * preprocessed.
 * To compute MWshift(j1, j2) we need to solve a matching problem again but
 * within P itself i.e. we must find the smallest (non-zero) amount of letters
 * that we need to shift P(j1, j2) to the _left_ until it matches P again.
 * IMPORTANT: This match must be a "full" match (when we are still within the
 * limits of P) and a "partial" match (when we are out of limits i.e. we've
 * moved to far to the left). For a "partial" match, all letters that are
 * still within P limits must match. When we have moved so far that all
 * letters are out of limits, then we have a "partial" match too.
 * This ensure that MWshift(j1, j2) is always <= j2.
 * 
 * Example: P = acbaba
 *              012345
 * 
 *   MWshift(0,1) = 1
 *   MWshift(1,2) = 2
 *   MWshift(2,3) = 3
 *   MWshift(3.4) = MWshift(2,4) = MWshift(1,4) = MWshift(0,4) = 3
 *   ...
 *   MWshift(4,5) = MWshift(5,6) = MWshift(4,6) = 2
 *   MWshift(0,6) = 5 (MWshift max)
 * 
 * A nice property of the MWshift function is that:
 *     if j1 <= j1' < j2' <= j2, then MWshift(j1', j2') <= MWshift(j1, j2)
 * Therefore the biggest value for MWshift is reached for j1 = 0 and j2 = nP.
 * To summarize:
 *     1 <= MWshift(j1, j2) <= min(j2, MWshift(0, nP))
 * Another property of the MWshift function is that if P1 and P2 are 2
 * patterns such that P1 is a prefix of P2 then MWshift1 and MWshift2 (their
 * respective MWshift functions) are equal on the set of points (j1, j2) that
 * are valid for MWshift1.
 * This last property is used by the init_MWshift_table() C function below.
 * 
 * No need to preprocess the MWshift function
 * ------------------------------------------
 *
 * However this preprocessing would be too expensive since it requires the
 * evaluation of nP * (nP + 1) / 2 values, and each evaluation itself is a
 * string matching sub-problem with a cost of its own. So the trick is to
 * delay those evaluations until they are needed, and the fact is that, in
 * practise, very few of them are actually needed compared to the total
 * number of possible MWshift(j1, j2) values.
 *
 * ppP.MWshift_table is a 2-dim array with nrow = ncol = ppP.buflength.
 * The layout of ppP.MWshift_table is (only the values marked with an "x"
 * will be potentially used):
 *
 *           1 2 3 4 5 6 j2
 *         0 x x x x - -
 *         1 - x x x - - 
 *         2 - - x x - -    ppP.seqlength = 4 <= ppP.buflength = 6
 *         3 - - - x - -
 *         4 - - - - - -
 *         5 - - - - - -
 *        j1
 *
 * The "x" region is defined by 0 <= j1 < j2 <= ppP.seqlength
 */

#define MWSHIFT(j1, j2) (ppP.MWshift_table[ppP.buflength * (j1) + (j2) - 1])

static int get_MWshift(int j1, int j2)
{
	int shift, k1, k2, length;
	const char *tmp;

	shift = MWSHIFT(j1, j2);
	if (shift != 0)
		return shift;
	for (shift = 1; shift < j2; shift++) {
		if (shift < j1) k1 = j1 - shift; else k1 = 0;
		k2 = j2 - shift;
		length = k2 - k1;
		tmp = ppP.seq + k1;
		if (memcmp(tmp, tmp + shift, length) == 0)
			break;
	}
	/* shift is j2 when the "for" loop is not interrupted by "break" */
	return MWSHIFT(j1, j2) = shift;
}

static void init_ppP_MWshift_table()
{
	int j1, j2 = 1;

	if (ppP.LCP == -1 && ppP.MWshift_table != NULL) {
		free(ppP.MWshift_table);
		ppP.MWshift_table = NULL;
	}
	if (ppP.buflength != 0 && ppP.MWshift_table == NULL) {
		ppP.MWshift_table = (int *)
			malloc(ppP.buflength * ppP.buflength * sizeof(int));
		if (ppP.MWshift_table == NULL)
			error("can't allocate memory for ppP.MWshift_table");
	}
	if (ppP.LCP != -1)
		j2 = ppP.LCP + 1;
	for ( ; j2 <= ppP.seqlength; j2++) {
		for (j1 = 0; j1 < j2; j1++) {
			MWSHIFT(j1, j2) = 0;
		}
	}
}


/****************************************************************************
 * The boyermoore() function
 * =========================
 * 
 * My poor understanding of the original Boyer-Moore algo as explained in Dan
 * Gusfield book "Algorithms on strings, trees, and sequences", is that it
 * will be inefficent in this situation:
 *   S = "Xbabababababab"
 *   P = "abababababab"
 * because after it has applied the "strong good suffix rule" and shifted P 2
 * letters to the right, it will start comparing P and S suffixes again which
 * is a waste of time.
 * The original idea behind our "MWshift feature" is to address this problem.
 * By using this "feature", our boyermoore() function below never compares
 * twice the letters that are in the current Matching Window (defined by
 * (j1, j2) in P and (i1, i2) in S).
 */

#define GET_S_LETTER(S, n, walk_backward) \
	((walk_backward) ? (S)->seq[(n)] : (S)->seq[(S)->length - 1 - (n)])

#define ADJUST_MW(i, j, shift) \
{ \
	if ((shift) <= (j)) \
		(j) -= (shift); \
	else { \
		(i) += (shift) - (j); \
		(j) = 0; \
	} \
}

/* Return 1-based end of last match or -1 if no match */
int _match_pattern_boyermoore(const cachedCharSeq *P, const cachedCharSeq *S,
		int nfirstmatches, int walk_backward)
{
	int nmatches, last_match_end, n, i1, i2, j1, j2, shift, shift1, i, j;
	char ppP_rmc, c; /* ppP_rmc is 'ppP.seq' right-most char */

	if (P->length <= 0)
		error("empty pattern");
	nmatches = 0;
	last_match_end = -1;
	init_ppP_seq(P, walk_backward);
	init_ppP_j0shift0();
	init_ppP_VSGSshift_table();
	if (ppP.seqlength <= MWSHIFT_NPMAX)
		init_ppP_MWshift_table();
	n = ppP.seqlength - 1;
	ppP_rmc = ppP.seq[n];
	j2 = 0;
	while (n < S->length) {
		if (j2 == 0) {
			/* No Matching Window yet, we need to find one */
			c = GET_S_LETTER(S, n, walk_backward);
			if (c != ppP_rmc) {
				shift = get_VSGSshift(c, ppP.seqlength - 1);
				n += shift;
				continue;
			}
			i1 = n;
			i2 = i1 + 1;
			j2 = ppP.seqlength;
			j1 = j2 - 1;
			/* Now we have a Matching Window (1-letter suffix) */
		}
		/* We try to extend the current Matching Window... */
		if (j1 > 0) {
			/* ... to the left */
			for (i = i1-1, j = j1-1; j >= 0; i--, j--)
				if ((c = GET_S_LETTER(S, i, walk_backward)) != ppP.seq[j])
					break;
			i1 = i + 1;
			j1 = j + 1;
		}
		if (j2 < ppP.seqlength) {
			/* ... to the right */
			for ( ; j2 < ppP.seqlength; i2++, j2++)
				if (GET_S_LETTER(S, i2, walk_backward) != ppP.seq[j2])
					break;
		}
		if (j2 == ppP.seqlength) { /* the Matching Window is a suffix */
			if (j1 == 0) {
				/* we have a full match! */
				_report_match(i1 + 1, ppP.seqlength);
				nmatches++;
				last_match_end = i1 + ppP.seqlength;
				if (nfirstmatches >= 0 && nmatches >= nfirstmatches)
					break;
				shift = ppP.shift0;
			} else {
				shift = get_VSGSshift(c, j1 - 1);
			}
		} else {
			shift = get_MWshift(j1, j2);
			c = GET_S_LETTER(S, n, walk_backward);
			if (c != ppP_rmc) {
				shift1 = get_VSGSshift(c, ppP.seqlength - 1);
				if (shift1 > shift)
					shift = shift1;
			}
		}
		n += shift;
		if (ppP.seqlength <= MWSHIFT_NPMAX) {
			ADJUST_MW(i1, j1, shift)
			ADJUST_MW(i2, j2, shift)
		} else {
			j2 = 0; /* forget the current Matching Window */
		}
	}
	return last_match_end;
}

