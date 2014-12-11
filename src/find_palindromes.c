#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"

#include <stdio.h>


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

static void find_palindromes_at(const char *x, int x_len,
	int i1, int i2, int max_loop_len1, int min_arm_len,
	const int *lkup, int lkup_len)
{
	int arm_len, valid_indices, key, val, is_match;
	char x1, x2;

	arm_len = 0;
	while (((valid_indices = i1 >= 0 && i2 < x_len) &&
                i2 - i1 <= max_loop_len1) || arm_len != 0)
	{
		if (valid_indices) {
			x1 = x[i1];
			x2 = x[i2];
			if (lkup == NULL) {
				is_match = x1 == x2;
			} else {
				key = (unsigned char) x1;
				if (key >= lkup_len
				 || (val = lkup[key]) == NA_INTEGER)
					is_match = FALSE;
				else
					is_match = val == x2;
			}
			if (is_match) {
				arm_len++;
				goto next;
			}
		}
		if (arm_len >= min_arm_len)
			_report_match(i1 + 2, i2 - i1 - 1);
		arm_len = 0;
	next:
		i1--;
		i2++;
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP find_palindromes(SEXP x, SEXP min_armlength, SEXP max_looplength,
		      SEXP L2R_lkup)
{
	Chars_holder x_holder;
	int x_len, min_arm_len, max_loop_len1, lkup_len, n;
	const int *lkup;

	x_holder = hold_XRaw(x);
	x_len = x_holder.length;
	min_arm_len = INTEGER(min_armlength)[0];
	max_loop_len1 = INTEGER(max_looplength)[0] + 1;
	if (L2R_lkup == R_NilValue) {
		lkup = NULL;
		lkup_len = 0;
	} else {
		lkup = INTEGER(L2R_lkup);
		lkup_len = LENGTH(L2R_lkup);
	}
	_init_match_reporting("MATCHES_AS_RANGES", 1);
	for (n = 0; n < x_len; n++) {
		/* Find palindromes centered on n. */
		find_palindromes_at(x_holder.ptr, x_len, n - 1, n + 1,
				    max_loop_len1, min_arm_len,
				    lkup, lkup_len);
		/* Find palindromes centered on n + 0.5. */
		find_palindromes_at(x_holder.ptr, x_len, n, n + 1,
				    max_loop_len1, min_arm_len,
				    lkup, lkup_len);
	}
	return _reported_matches_asSEXP();
}

