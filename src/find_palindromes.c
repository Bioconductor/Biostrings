#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"

#include <stdio.h>


static int is_match(char c1, char c2, const int *lkup, int lkup_len)
{
	int key, val;

	if (lkup != NULL) {
		key = (unsigned char) c1;
		if (key >= lkup_len || (val = lkup[key]) == NA_INTEGER)
			return 0;
		c1 = (char) val;
	}
	return c1 == c2;
}

static void get_find_palindromes_at(const char *x, int x_len,
	int i1, int i2, int max_loop_len1, int min_arm_len, int max_nmis,
	const int *lkup, int lkup_len)
{
	int arm_len, valid_indices;
	char c1, c2;

	arm_len = 0;
	while (((valid_indices = i1 >= 0 && i2 < x_len) &&
                i2 - i1 <= max_loop_len1) || arm_len != 0)
	{
		if (valid_indices) {
			c1 = x[i1];
			c2 = x[i2];
			if (is_match(c1, c2, lkup, lkup_len) ||
			    max_nmis-- > 0) {
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

static int get_palindrome_arm_length(const char *x, int x_len, int max_nmis,
	const int *lkup, int lkup_len)
{
	int i1, i2;
	char c1, c2;

	for (i1 = 0, i2 = x_len - 1; i1 < i2; i1++, i2--) {
		c1 = x[i1];
		c2 = x[i2];
		if (!(is_match(c1, c2, lkup, lkup_len) ||
		      max_nmis-- > 0))
			break;
	}
	return i1;
}

/* --- .Call ENTRY POINT --- */
SEXP find_palindromes(SEXP x, SEXP min_armlength, SEXP max_looplength,
		      SEXP max_mismatch, SEXP L2R_lkup)
{
	Chars_holder x_holder;
	int x_len, min_arm_len, max_loop_len1, max_nmis, lkup_len, n;
	const int *lkup;

	x_holder = hold_XRaw(x);
	x_len = x_holder.length;
	min_arm_len = INTEGER(min_armlength)[0];
	max_loop_len1 = INTEGER(max_looplength)[0] + 1;
	max_nmis = INTEGER(max_mismatch)[0];
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
		get_find_palindromes_at(x_holder.ptr, x_len, n - 1, n + 1,
					max_loop_len1, min_arm_len, max_nmis,
					lkup, lkup_len);
		/* Find palindromes centered on n + 0.5. */
		get_find_palindromes_at(x_holder.ptr, x_len, n, n + 1,
					max_loop_len1, min_arm_len, max_nmis,
					lkup, lkup_len);
	}
	return _reported_matches_asSEXP();
}

/* --- .Call ENTRY POINT --- */
SEXP palindrome_arm_length(SEXP x, SEXP max_mismatch, SEXP L2R_lkup)
{
	Chars_holder x_holder;
	int x_len, max_nmis, lkup_len, arm_len;
	const int *lkup;

	x_holder = hold_XRaw(x);
	x_len = x_holder.length;
	max_nmis = INTEGER(max_mismatch)[0];
	if (L2R_lkup == R_NilValue) {
		lkup = NULL;
		lkup_len = 0;
	} else {
		lkup = INTEGER(L2R_lkup);
		lkup_len = LENGTH(L2R_lkup);
	}
	arm_len = get_palindrome_arm_length(x_holder.ptr, x_len, max_nmis,
					    lkup, lkup_len);
	return ScalarInteger(arm_len);
}

