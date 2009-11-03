#include "Biostrings.h"

#include <stdio.h>


/****************************************************************************/
static int debug = 0;

SEXP debug_find_palindromes()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'find_palindromes.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'find_palindromes.c'\n");
#endif
	return R_NilValue;
}

/****************************************************************************/

static void naive_palindrome_search(const char *S, int nS, int armlen_min, int ngaps_max)
{
	int n1, n2, ngaps, armlen, Lpos, Rpos, all_letter0;
	char letter0;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] naive_palindrome_search(): nS=%d armlen_min=%d ngaps_max=%d\n",
			nS, armlen_min, ngaps_max);
	}
#endif
	for (n1 = armlen_min, n2 = 2 * armlen_min; n2 <= nS; n1++, n2++) {
		for (ngaps = 0; ngaps <= ngaps_max; ngaps++) {
			armlen = 0;
			Lpos = n1 - 1;
			Rpos = n1 + ngaps;
			while (0 <= Lpos && Rpos < nS && S[Lpos] == S[Rpos]) {
				if (ngaps == 0) {
					if (armlen == 0) {
						letter0 = S[Rpos];
						all_letter0 = 1;
					} else {
						if (S[Rpos] != letter0)
							all_letter0 = 0;
					} 
				}
				armlen++;
				Lpos--;
				Rpos++;
			}
			Lpos++;
			if (ngaps == 0 && armlen != 0 && all_letter0) {
				// The current palindrome is in fact the left part of a region where the
				// same letter (letter0) is repeated. We move to the right end of this region.
				while (Rpos < nS && S[Rpos] == letter0)
					Rpos++;
				if (Rpos - Lpos < 2 * armlen_min)
					continue;
				Rpos--;
				n1 = Rpos;
				n2 = Rpos + armlen_min;
			} else {
				if (armlen < armlen_min)
					continue;
				Rpos--;
			}
			_report_match(Lpos + 1, Rpos - Lpos + 1);
			break;
		}
	}
	return;
}

static void naive_antipalindrome_search(const char *S, int nS, int armlen_min, int ngaps_max,
		const int *lkup, int lkup_length)
{
	int n1, n2, ngaps, armlen, Lpos, Rpos, all_letter0, lkup_key, lkup_val;
	char letter0;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] naive_antipalindrome_search(): nS=%d armlen_min=%d ngaps_max=%d\n",
			nS, armlen_min, ngaps_max);
	}
#endif
	for (n1 = armlen_min, n2 = 2 * armlen_min; n2 <= nS; n1++, n2++) {
		for (ngaps = 0; ngaps <= ngaps_max; ngaps++) {
			armlen = 0;
			Lpos = n1 - 1;
			Rpos = n1 + ngaps;
			while (0 <= Lpos && Rpos < nS) {
				lkup_key = (unsigned char) S[Lpos];
				if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
					error("key %d not in lookup table", lkup_key);
				}
				if (((char) lkup_val) != S[Rpos])
					break;
				if (ngaps == 0) {
					if (armlen == 0) {
						letter0 = S[Rpos];
						// Will be 1 iff S[Rpos] is its own complementary (only
						// IUPAC letter N and the gap letter - have this property)
						all_letter0 = S[Lpos] == S[Rpos];
					} else {
						if (S[Rpos] != letter0)
							all_letter0 = 0;
					} 
				}
				armlen++;
				Lpos--;
				Rpos++;
			}
			Lpos++;
			if (ngaps == 0 && armlen != 0 && all_letter0) {
				// The current palindrome is in fact the left part of a region where the
				// same letter (letter0) is repeated. We move to the right end of this region.
				while (Rpos < nS && S[Rpos] == letter0)
					Rpos++;
				if (Rpos - Lpos < 2 * armlen_min)
					continue;
				Rpos--;
				n1 = Rpos;
				n2 = Rpos + armlen_min;
			} else {
				if (armlen < armlen_min)
					continue;
				Rpos--;
			}
			_report_match(Lpos + 1, Rpos - Lpos + 1);
			break;
		}
	}
	return;
}

/*
 *
 */
SEXP find_palindromes(SEXP s_xp, SEXP s_offset, SEXP s_length,
                SEXP min_armlength, SEXP max_ngaps, SEXP L2R_lkup)
{
	int subj_offset, subj_length, armlen_min, ngaps_max;
	const Rbyte *subj;

	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	armlen_min = INTEGER(min_armlength)[0];
	ngaps_max = INTEGER(max_ngaps)[0];
	_init_match_reporting("ASIRANGES");
	if (L2R_lkup == R_NilValue)
		naive_palindrome_search((char *) subj, subj_length,
			armlen_min, ngaps_max);
	else
		naive_antipalindrome_search((char *) subj, subj_length,
			armlen_min, ngaps_max,
			INTEGER(L2R_lkup), LENGTH(L2R_lkup));
	return _reported_matches_asSEXP();
}

