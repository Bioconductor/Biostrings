#include "Biostrings.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"

#define REPLACE_IFNOTEXTEND	1
#define SKIP_IFNOTEXTEND	2
#define MERGE_IFNOTEXTEND	3
#define ERROR_IFNOTEXTEND	4

static int notextend_action;
static int skip_or_merge_count;
static char errmsg_buf[200];
static ByteTrTable byte2code;

static void set_notextend_action(const char *if_not_extending)
{
	if (strcmp(if_not_extending, "replace") == 0)
		notextend_action = REPLACE_IFNOTEXTEND;
	else if (strcmp(if_not_extending, "skip") == 0)
		notextend_action = SKIP_IFNOTEXTEND;
	else if (strcmp(if_not_extending, "merge") == 0)
		notextend_action = MERGE_IFNOTEXTEND;
	else if (strcmp(if_not_extending, "error") == 0)
		notextend_action = ERROR_IFNOTEXTEND;
	else
		error("invalid 'if_not_extending' value %s", if_not_extending);
	return;
}

static int replace_letter_at(char *dest, int dest_length,
		const int *at, int at_length, const char *src, int use_byte2code)
{
	int i, j, byte, are_IUPAC;
	char old_letter, new_letter;

	for (j = 0; j < at_length; j++) {
		i = at[j] - 1;
		if (i == NA_INTEGER || i < 0 || i >= dest_length) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'at' contains NAs or \"out of limits\" locations");
			return -1;
		}
		new_letter = src[j];
		if (use_byte2code) {
			byte = (unsigned char) new_letter;
			if ((new_letter = byte2code.byte2code[byte]) == NA_INTEGER) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "'letter' contains invalid letters "
					 "(first found has code %d)",
					 byte);
				return -1;
			}
		}
		old_letter = dest[i];
		if (old_letter == new_letter)
			continue;
		if (notextend_action == REPLACE_IFNOTEXTEND) {
			dest[i] = new_letter;
			continue;
		}
		are_IUPAC = ((unsigned char) old_letter) < 16
			 && ((unsigned char) new_letter) < 16;
		if (are_IUPAC && (old_letter & ~new_letter) == 0) {
			// new_letter extends old_letter
			dest[i] = new_letter;
			continue;
		}
		if (notextend_action == ERROR_IFNOTEXTEND) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "new letter (code %d) does not extend old letter (code %d) "
				 "at location %d", new_letter, old_letter, i + 1);
			return -1;
		}
		skip_or_merge_count++;
		if (notextend_action == SKIP_IFNOTEXTEND)
			continue;
		// only IUPAC letters are mergeable
		if (!are_IUPAC) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot merge non IUPAC letters at location %d", i + 1);
			return -1;
		}
		dest[i] |= new_letter;
	}
	return 0;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_replace_letter_at(SEXP x, SEXP at, SEXP letter, SEXP lkup,
		SEXP if_not_extending, SEXP verbose)
{
	const char *x_classname;
	cachedCharSeq X;
	SEXP tag, letter_elt, ans;
	int at_length, letter_length, letter_elt_length, letter_ncharsum, i;
	const int *at_p;

	x_classname = get_classname(x);
	X = cache_XRaw(x);
	at_length = LENGTH(at);
	letter_length = LENGTH(letter);
	if (lkup != R_NilValue)
		_init_ByteTrTable_with_lkup(&byte2code, lkup);
	set_notextend_action(CHAR(STRING_ELT(if_not_extending, 0)));

	PROTECT(tag = NEW_RAW(X.length));
	memcpy((char *) RAW(tag), X.seq, X.length);
	skip_or_merge_count = letter_ncharsum = 0;
	at_p = INTEGER(at);
	for (i = 0; i < letter_length; i++) {
		letter_elt = STRING_ELT(letter, i);
		if (letter_elt == NA_STRING) {
			UNPROTECT(1);
			error("'letter' contains NAs");
		}
		letter_ncharsum += letter_elt_length = LENGTH(letter_elt);
		if (letter_ncharsum > at_length)
			break;
		if (replace_letter_at((char *) RAW(tag), LENGTH(tag),
				      at_p, letter_elt_length, CHAR(letter_elt),
				      lkup != R_NilValue) != 0) {
			UNPROTECT(1);
			error("%s", errmsg_buf);
		}
		at_p += letter_elt_length;
	}
	if (letter_ncharsum != at_length) {
		UNPROTECT(1);
		error("total nb of letters in 'letter' must be the same as nb of locations");
	}
	if (skip_or_merge_count != 0
         && notextend_action != REPLACE_IFNOTEXTEND
         && LOGICAL(verbose)[0])
		warning("%s %d letter(s)",
			notextend_action == SKIP_IFNOTEXTEND ? "skipped" : "merged",
			skip_or_merge_count);
	PROTECT(ans = new_XRaw_from_tag(x_classname, tag));
	UNPROTECT(2);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_inplace_replace_letter_at(SEXP x, SEXP at, SEXP letter, SEXP lkup)
{
	SEXP tag, letter_elt;
	int at_length, letter_length, letter_elt_length, letter_ncharsum, i;
	const int *at_p;

	at_length = LENGTH(at);
	letter_length = LENGTH(letter);
	if (lkup != R_NilValue)
		_init_ByteTrTable_with_lkup(&byte2code, lkup);
	notextend_action = MERGE_IFNOTEXTEND;

	tag = get_XVector_tag(x);
	skip_or_merge_count = letter_ncharsum = 0;
	at_p = INTEGER(at);
	for (i = 0; i < letter_length; i++) {
		letter_elt = STRING_ELT(letter, i);
		if (letter_elt == NA_STRING)
			error("'letter' contains NAs");
		letter_ncharsum += letter_elt_length = LENGTH(letter_elt);
		if (letter_ncharsum > at_length)
			break;
		if (replace_letter_at((char *) RAW(tag), LENGTH(tag),
				      at_p, letter_elt_length, CHAR(letter_elt),
				      lkup != R_NilValue) != 0) {
			error("%s", errmsg_buf);
		}
		at_p += letter_elt_length;
	}
	if (letter_ncharsum != at_length)
		error("total nb of letters in 'letter' must be the same as nb of locations");
	return x;
}

