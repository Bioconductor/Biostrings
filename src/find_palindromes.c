#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"

#include <stdio.h>


static int is_match(char c1, char c2, int allow_wobble1,
	const int *lkup, int lkup_len)
{
	int key, val;

	if (lkup != NULL) {
		key = (unsigned char) c1;
		if (key >= lkup_len || (val = lkup[key]) == NA_INTEGER)
			return 0;
		c1 = (char) val;
	}
	if (allow_wobble1) {
		return c1 == c2 ||
			(c1 == 2 && c2 == 8) || // G and T/U
			(c1 == 1 && c2 == 4); // T/U and G
	} else {
		return c1 == c2;
	}
}

static void get_find_palindromes_at(const char *x, int x_len,
	int i1, int i2, int max_loop_len1, int min_arm_len,
	int max_nmis, int allow_wobble1,
	const int *lkup, int lkup_len)
{
	int stack = max_nmis + 1; // size of stack (assumes max_nmis is small)
	int arms[stack]; // lengths of matched arms
	int arm_len, valid_indices, nmis;
	char c1, c2;
	static int ismatch;

	arm_len = 0;
	nmis = 0;
	int count = 0; // position in stack
	while (((valid_indices = i1 >= 0 && i2 < x_len) &&
		i2 - i1 <= max_loop_len1) || arm_len != 0)
	{
		if (valid_indices) {
			c1 = x[i1];
			c2 = x[i2];
			ismatch = is_match(c1, c2, allow_wobble1,
					   lkup, lkup_len);
			if (!ismatch) {
				arms[count] = arm_len;
				count++;
			}
			if (ismatch || nmis++ < max_nmis) {
				arm_len++;
				goto next;
			}
		}
		if (arm_len >= min_arm_len)
			_report_match(i1 + 2, i2 - i1 - 1);
		if ((i1 > 0 && i2 + 1 < x_len) || nmis > max_nmis) {
			nmis--;
			arm_len -= arms[0];
			/* Pop off bottom of stack */
			for (int i = 1; i < count; i++)
				arms[i] -= arms[0];
			for (int i = 1; i < count; i++)
				arms[i - 1] = arms[i];
			count--;
		} else {
			arm_len = 0;
			nmis = 0;
		}
	    next:
		i1--;
		i2++;
	}
	return;
}

static int get_palindrome_arm_length(const char *x, int x_len, int max_nmis,
	int allow_wobble1, const int *lkup, int lkup_len)
{
	int i1, i2;
	char c1, c2;

	for (i1 = 0, i2 = x_len - 1; i1 < i2; i1++, i2--) {
		c1 = x[i1];
		c2 = x[i2];
		if (!(is_match(c1, c2, allow_wobble1, lkup, lkup_len) ||
		      max_nmis-- > 0))
			break;
	}
	return i1;
}

/* --- .Call ENTRY POINT --- */
SEXP find_palindromes(SEXP x, SEXP min_armlength, SEXP max_looplength,
		      SEXP max_mismatch, SEXP min_looplength, SEXP allow_wobble,
		      SEXP L2R_lkup)
{
	Chars_holder x_holder;
	int x_len, min_arm_len, max_loop_len1, max_nmis, lkup_len, n,
		min_loop_len1, min_loop_len2, allow_wobble1;
	const int *lkup;

	x_holder = hold_XRaw(x);
	x_len = x_holder.length;
	min_arm_len = INTEGER(min_armlength)[0];
	max_loop_len1 = INTEGER(max_looplength)[0] + 1;
	max_nmis = INTEGER(max_mismatch)[0];
	min_loop_len1 = INTEGER(min_looplength)[0];
	min_loop_len2 = (min_loop_len1 + 1) / 2; // 0 1 1 2 2 ...
	min_loop_len1 /= 2;                      // 0 0 1 1 2 2 ...
	allow_wobble1 = INTEGER(allow_wobble)[0];
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
		get_find_palindromes_at(x_holder.ptr, x_len,
					n - 1 - min_loop_len1,
					n + 1 + min_loop_len1,
					max_loop_len1, min_arm_len,
					max_nmis, allow_wobble1,
					lkup, lkup_len);
		/* Find palindromes centered on n + 0.5. */
		get_find_palindromes_at(x_holder.ptr, x_len,
					n - min_loop_len2,
					n + 1 + min_loop_len2,
					max_loop_len1, min_arm_len,
					max_nmis, allow_wobble1,
					lkup, lkup_len);
	}
	return _reported_matches_asSEXP();
}

/* --- .Call ENTRY POINT --- */
SEXP palindrome_arm_length(SEXP x, SEXP max_mismatch, SEXP allow_wobble,
	SEXP L2R_lkup)
{
	XStringSet_holder x_holder;
	int x_len, max_nmis, allow_wobble1, lkup_len;
	const int *lkup;
	SEXP ans;
	int *ans_p;
	Chars_holder x_elt_holder;

	x_holder = _hold_XStringSet(x);
	x_len = _get_XStringSet_length(x);
	max_nmis = INTEGER(max_mismatch)[0];
	allow_wobble1 = INTEGER(allow_wobble)[0];
	if (L2R_lkup == R_NilValue) {
		lkup = NULL;
		lkup_len = 0;
	} else {
		lkup = INTEGER(L2R_lkup);
		lkup_len = LENGTH(L2R_lkup);
	}
	PROTECT(ans = NEW_INTEGER(x_len));
	ans_p = INTEGER(ans);
	for (int i = 0; i < x_len; i++) {
		x_elt_holder = _get_elt_from_XStringSet_holder(&x_holder, i);
		ans_p[i] = get_palindrome_arm_length(x_elt_holder.ptr,
						     x_elt_holder.length,
						     max_nmis, allow_wobble1,
						    lkup, lkup_len);
	}
	UNPROTECT(1);
	return ans;
}

