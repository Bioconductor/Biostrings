#include "Biostrings.h"

#define REPLACE_IFNOTEXTEND	1
#define SKIP_IFNOTEXTEND	2
#define MERGE_IFNOTEXTEND	3
#define ERROR_IFNOTEXTEND	4

static int notextend_action;
static int skip_or_merge_count;
static char errmsg_buf[200];
static int chrtrtable[CHRTRTABLE_LENGTH];

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

static void init_chrtrtable_with_lkup(SEXP lkup)
{
	int trkey, trval;

	if (LENGTH(lkup) > CHRTRTABLE_LENGTH)
		error("Biostrings internal error in init_chrtrtable_with_lkup(): "
		      "LENGTH(lkup) > CHRTRTABLE_LENGTH");
	for (trkey = 0; trkey < LENGTH(lkup); trkey++) {
		trval = INTEGER(lkup)[trkey];
		chrtrtable[trkey] = trval == NA_INTEGER ? -1 : trval;
	}
	for ( ; trkey < CHRTRTABLE_LENGTH; trkey++)
		chrtrtable[trkey] = -1;
	return;
}

static int replace_locs(char *dest, int dest_length,
		const int *loc, int loc_length, const char *src, int use_chrtrtable)
{
	int i, j, trkey, are_IUPAC;
	char old_letter, new_letter;

	for (j = 0; j < loc_length; j++) {
		i = loc[j] - 1;
		if (i == NA_INTEGER || i < 0 || i >= dest_length) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'loc' contains NAs or out of limits locations");
			return -1;
		}
		new_letter = src[j];
		if (use_chrtrtable) {
			trkey = (unsigned char) new_letter;
			if ((new_letter = chrtrtable[trkey]) == -1) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "'letter' contains invalid letters "
					 "(first found has code %d)",
					 trkey);
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
SEXP XString_replace_locs_bySTRSXP(SEXP x, SEXP loc, SEXP letter, SEXP lkup,
		SEXP if_not_extending, SEXP verbose)
{
	const char *x_class;
	RoSeq x_seq;
	SEXP tag, letter_elt, data, ans;
	int loc_length, letter_length, letter_elt_length, letter_ncharsum, i;
	const int *loc_p;

	x_class = _get_class(x);
	x_seq = _get_XString_asRoSeq(x);
	loc_length = LENGTH(loc);
	letter_length = LENGTH(letter);
	if (lkup != R_NilValue)
		init_chrtrtable_with_lkup(lkup);
	set_notextend_action(CHAR(STRING_ELT(if_not_extending, 0)));

	PROTECT(tag = NEW_RAW(x_seq.nelt));
	memcpy((char *) RAW(tag), x_seq.elts, x_seq.nelt);
	skip_or_merge_count = letter_ncharsum = 0;
	loc_p = INTEGER(loc);
	for (i = 0; i < LENGTH(letter); i++) {
		letter_elt = STRING_ELT(letter, i);
		if (letter_elt == NA_STRING) {
			UNPROTECT(1);
			error("'letter' contains NAs");
		}
		letter_ncharsum += letter_elt_length = LENGTH(letter_elt);
		if (letter_ncharsum > loc_length)
			break;
		if (replace_locs((char *) RAW(tag), x_seq.nelt,
                		 loc_p, letter_elt_length, CHAR(letter_elt),
				 lkup != R_NilValue) != 0) {
			UNPROTECT(1);
			error("%s", errmsg_buf);
		}
		loc_p += letter_elt_length;
	}
	if (letter_ncharsum != loc_length) {
		UNPROTECT(1);
		error("total nb of letters in 'letter' must be the same as nb of locations");
	}
	if (skip_or_merge_count != 0
         && notextend_action != REPLACE_IFNOTEXTEND
         && LOGICAL(verbose)[0])
		warning("%s %d letter(s)",
			notextend_action == SKIP_IFNOTEXTEND ? "skipped" : "merged",
			skip_or_merge_count);
	PROTECT(data = _new_XRaw(tag));
	PROTECT(ans = _new_XString(x_class, data, 0, LENGTH(tag)));
	UNPROTECT(3);
	return ans;
}

