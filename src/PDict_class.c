/****************************************************************************
 *                   Basic manipulation of PDict objects                    *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc()/Srealloc() */

static int debug = 0;

SEXP debug_Pdict_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'Pdict_class.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'Pdict_class.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * The cropped dictionary buffer.
 *
 * The 'cropped_dict' object is used for storing the pointers to the patterns
 * contained in the input dictionary provided by the user.
 * These patterns can only be accessed for reading!
 */

typedef struct croppeddict {
	// 'start' and 'end' define the cropping operation
	int start;
	int end;
	int width;
	char cropping_type; // 'a', 'b', 'c', 'd'
	int head_min_width, head_max_width;
	int tail_min_width, tail_max_width;
	// subpatterns resulting from the cropping operation are stored in
	// the (const char *) array below (of length 'length')
	const char **patterns;
	int length;
} CroppedDict;

static CroppedDict cropped_dict;

/*
 * 'start' and 'end' are user-specified values that describe the cropping
 * operation that needs to be applied to an arbitrary set of input sequences
 * to turn it into a constant width dictionary. In the most general case, this
 * cropping operation consists of removing a prefix, a suffix or both from each
 * input sequence. The set of prefixes to remove is called "the head", and the
 * set of suffixes "the tail".
 * 'start' and 'end' must be integers (eventually NA) giving the relative
 * starting and ending position of the substring to extract from each input
 * sequence. Negative values indicate positions relative to the end of the
 * input sequences.
 * 4 types of cropping operations are supported:
 *   a) 1 <= start <= end:
 *        o the croppping operation will work on any variable width input
 *          where the shortest sequence has at least 'end' letters,
 *        o the head is rectangular (constant width),
 *        o the tail has a variable width (tail_min_width, tail_max_width);
 *   b) start <= end <= -1:
 *        o the croppping operation will work on any variable width input
 *          where the shortest sequence has at least '-start' letters,
 *        o the tail is rectangular (constant width),
 *        o the head has a variable width (head_min_width, head_max_width).
 *   c) 1 <= start and end == NA:
 *        o the set of input sequences must be already of constant width or
 *          the cropping operation will raise an error,
 *        o the head is rectangular (constant width),
 *        o no tail;
 *   d) start == NA and end <= -1:
 *        o the set of input sequences must be already of constant width or
 *          the cropping operation will raise an error,
 *        o the tail is rectangular (constant width),
 *        o no head;
 */
static char cropping_type(int start, int end)
{
	if (start == NA_INTEGER && end == NA_INTEGER)
		error("'start' and 'end' cannot both be NA");
	if (start == 0)
		error("'start' must be a single >= 1, <= -1 or NA integer");
	if (end == 0)
		error("'end' must be a single >= 1, <= -1 or NA integer");
	if (end == NA_INTEGER) {
		if (start < 0)
			error("'start' must be positive when 'end' is NA");
		return 'c';
	}
	if (start == NA_INTEGER) {
		if (end > 0)
			error("'end' must be negative when 'start' is NA");
		return 'd';
	}
	if ((start > 0) != (end > 0))
		error("'start' and 'end' must have the same sign");
	if (end < start)
		error("'end' must be >= 'start'");
	return start > 0 ? 'a' : 'b';
}

static void alloc_cropped_dict(int length, int start, int end)
{
	cropped_dict.cropping_type = cropping_type(start, end);
	cropped_dict.start = start;
	cropped_dict.end = end;
	if (cropped_dict.cropping_type == 'a') {
		cropped_dict.width = cropped_dict.end - cropped_dict.start + 1;
		cropped_dict.tail_min_width = cropped_dict.tail_max_width = -1;
	}
	if (cropped_dict.cropping_type == 'b') {
		cropped_dict.width = cropped_dict.end - cropped_dict.start + 1;
		cropped_dict.head_min_width = cropped_dict.head_max_width = -1;
	}
	cropped_dict.patterns = Salloc((long) length, const char *);
	cropped_dict.length = length;
	return;
}

static void add_subpattern_to_cropped_dict(int poffset,
		const char *pattern, int pattern_length)
{
	int head_width, tail_width;

	if (pattern_length == 0)
		error("'dict' contains empty patterns");

	switch (cropped_dict.cropping_type) {

	    case 'a':
		tail_width = pattern_length - cropped_dict.end;
		if (tail_width < 0)
			error("'dict' contains patterns with less than %d characters",
			      cropped_dict.end);
		cropped_dict.patterns[poffset] = pattern + cropped_dict.start - 1;
		if (cropped_dict.tail_min_width == -1) {
			cropped_dict.tail_min_width = cropped_dict.tail_max_width = tail_width;
			break;
		}
		if (tail_width < cropped_dict.tail_min_width)
			cropped_dict.tail_min_width = tail_width;
		if (tail_width > cropped_dict.tail_max_width)
			cropped_dict.tail_max_width = tail_width;
		break;

	    case 'b':
		head_width = pattern_length + cropped_dict.start;
		if (head_width < 0)
			error("'dict' contains patterns with less than %d characters",
			      -cropped_dict.start);
		cropped_dict.patterns[poffset] = pattern + head_width;
		if (cropped_dict.head_min_width == -1) {
			cropped_dict.head_min_width = cropped_dict.head_max_width = head_width;
			break;
		}
		if (head_width < cropped_dict.head_min_width)
			cropped_dict.head_min_width = head_width;
		if (head_width > cropped_dict.head_max_width)
			cropped_dict.head_max_width = head_width;
		break;

	    case 'c':
		if (cropped_dict.end == NA_INTEGER) {
			cropped_dict.end = pattern_length;
			cropped_dict.width = cropped_dict.end - cropped_dict.start + 1;
			if (cropped_dict.width < 1)
				error("'dict' contains patterns with less than %d characters",
				      cropped_dict.start);
		} else {
			if (pattern_length != cropped_dict.end)
				error("all patterns in 'dict' must have the same length");
		}
		cropped_dict.patterns[poffset] = pattern + cropped_dict.start - 1;
		break;

	    case 'd':
		if (cropped_dict.start == NA_INTEGER) {
			cropped_dict.start = -pattern_length;
			cropped_dict.width = cropped_dict.end - cropped_dict.start + 1;
			if (cropped_dict.width < 1)
				error("'dict' contains patterns with less than %d characters",
				      -cropped_dict.end);
		} else {
			if (pattern_length != -cropped_dict.start)
				error("all patterns in 'dict' must have the same length");
		}
		cropped_dict.patterns[poffset] = pattern;
		break;
	}
	return;
}

int _init_CroppedDict_with_CHARACTER(SEXP dict, int tb_start, int tb_end)
{
	int length, poffset;
	SEXP dict_elt;

	length = LENGTH(dict);
	alloc_cropped_dict(length, tb_start, tb_end);
	for (poffset = 0; poffset < length; poffset++) {
		dict_elt = STRING_ELT(dict, poffset);
		if (dict_elt == NA_STRING)
			error("'dict' contains NAs");
		add_subpattern_to_cropped_dict(poffset, CHAR(dict_elt), LENGTH(dict_elt));
	}
	return length;
}

int _init_CroppedDict_with_XStringSet(SEXP dict, int tb_start, int tb_end)
{
	int length, poffset;
	CachedXStringSet cached_dict;
	RoSeq pattern;

	length = _get_XStringSet_length(dict);
	alloc_cropped_dict(length, tb_start, tb_end);
	cached_dict = _new_CachedXStringSet(dict);
	for (poffset = 0; poffset < length; poffset++) {
		pattern = _get_CachedXStringSet_elt_asRoSeq(&cached_dict, poffset);
		add_subpattern_to_cropped_dict(poffset, pattern.elts, pattern.nelt);
	}
	return length;
}

int _CroppedDict_length()
{
	return cropped_dict.length;
}

int _CroppedDict_width()
{
	return cropped_dict.width;
}

const char *_CroppedDict_pattern(int poffset)
{
	return cropped_dict.patterns[poffset];
}


/****************************************************************************
 * Buffer of duplicates.
 */

static IntBuf dup2unq_buf;

void init_dup2unq_buf(int length)
{
	dup2unq_buf = _new_IntBuf(length, length, NA_INTEGER);
	return;
}

void report_dup(int poffset, int P_id)
{
	dup2unq_buf.elts[poffset] = P_id;
	return;
}


/****************************************************************************
 * Turning our local data structures into an R list (SEXP).
 *
 * Note that none of the functions below is a .Call() entry point!
 */

static SEXP oneint_asINTEGER(int oneint)
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = oneint;
	UNPROTECT(1);
	return ans;
}

/* NOT a .Call() entry point! */
SEXP _CroppedDict_geom_asLIST()
{
	SEXP ans, ans_names, ans_elt;
	const char *min_width_name, *max_width_name;
	int min_width_val, max_width_val;

	if (cropped_dict.cropping_type != 'a' && cropped_dict.cropping_type != 'b') {
		PROTECT(ans_names = NEW_CHARACTER(1));
		PROTECT(ans = NEW_LIST(1));
	} else {
		if (cropped_dict.cropping_type == 'a') {
			min_width_name = "tail.min.width";
			min_width_val = cropped_dict.tail_min_width;
			max_width_name = "tail.max.width";
			max_width_val = cropped_dict.tail_max_width;
		} else {
			min_width_name = "head.min.width";
			min_width_val = cropped_dict.head_min_width;
			max_width_name = "head.max.width";
			max_width_val = cropped_dict.head_max_width;
		}

		PROTECT(ans_names = NEW_CHARACTER(3));
		PROTECT(ans = NEW_LIST(3));

		/* set the "min.width" element */
		SET_STRING_ELT(ans_names, 1, mkChar(min_width_name));
		PROTECT(ans_elt = oneint_asINTEGER(min_width_val));
		SET_ELEMENT(ans, 1, ans_elt);
		UNPROTECT(1);

		/* set the "max.width" element */
		SET_STRING_ELT(ans_names, 2, mkChar(max_width_name));
		PROTECT(ans_elt = oneint_asINTEGER(max_width_val));
		SET_ELEMENT(ans, 2, ans_elt);
		UNPROTECT(1);
	}

	/* set the "width" element */
	SET_STRING_ELT(ans_names, 0, mkChar("width"));
	PROTECT(ans_elt = oneint_asINTEGER(cropped_dict.width));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}

/* NOT a .Call() entry point! */
SEXP _dup2unq_asINTEGER()
{
	return _IntBuf_asINTEGER(&dup2unq_buf);
}

