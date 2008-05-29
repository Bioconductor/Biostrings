/****************************************************************************
 *               Preprocessing a constant width DNA dictionary              *
 *                                   i.e.                                   *
 *                 storing it in an Aho-Corasick 4-ary tree                 *
 *                                                                          *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc()/Srealloc() */

#include <stdlib.h>
#include <string.h>
#include <limits.h>

static int debug = 0;


/****************************************************************************
 * The input_uldict object is used for storing the pointers to the patterns
 * contained in the input dictionary provided by the user.
 * These patterns can only be accessed for reading!
 */

typedef struct uldict {
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
} ULdict;

static ULdict input_uldict;

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

static void alloc_input_uldict(int length, int start, int end)
{
	input_uldict.cropping_type = cropping_type(start, end);
	input_uldict.start = start;
	input_uldict.end = end;
	if (input_uldict.cropping_type == 'a') {
		input_uldict.width = input_uldict.end - input_uldict.start + 1;
		input_uldict.tail_min_width = input_uldict.tail_max_width = -1;
	}
	if (input_uldict.cropping_type == 'b') {
		input_uldict.width = input_uldict.end - input_uldict.start + 1;
		input_uldict.head_min_width = input_uldict.head_max_width = -1;
	}
	input_uldict.patterns = Salloc((long) length, const char *);
	input_uldict.length = length;
	return;
}

static void add_subpattern_to_input_uldict(int poffset,
		const char *pattern, int pattern_length)
{
	int head_width, tail_width;

	if (pattern_length == 0)
		error("'dict' contains empty patterns");

	switch (input_uldict.cropping_type) {

	    case 'a':
		tail_width = pattern_length - input_uldict.end;
		if (tail_width < 0)
			error("'dict' contains patterns with less than %d characters",
			      input_uldict.end);
		input_uldict.patterns[poffset] = pattern + input_uldict.start - 1;
		if (input_uldict.tail_min_width == -1) {
			input_uldict.tail_min_width = input_uldict.tail_max_width = tail_width;
			break;
		}
		if (tail_width < input_uldict.tail_min_width)
			input_uldict.tail_min_width = tail_width;
		if (tail_width > input_uldict.tail_max_width)
			input_uldict.tail_max_width = tail_width;
		break;

	    case 'b':
		head_width = pattern_length + input_uldict.start;
		if (head_width < 0)
			error("'dict' contains patterns with less than %d characters",
			      -input_uldict.start);
		input_uldict.patterns[poffset] = pattern + head_width;
		if (input_uldict.head_min_width == -1) {
			input_uldict.head_min_width = input_uldict.head_max_width = head_width;
			break;
		}
		if (head_width < input_uldict.head_min_width)
			input_uldict.head_min_width = head_width;
		if (head_width > input_uldict.head_max_width)
			input_uldict.head_max_width = head_width;
		break;

	    case 'c':
		if (input_uldict.end == NA_INTEGER) {
			input_uldict.end = pattern_length;
			input_uldict.width = input_uldict.end - input_uldict.start + 1;
			if (input_uldict.width < 1)
				error("'dict' contains patterns with less than %d characters",
				      input_uldict.start);
		} else {
			if (pattern_length != input_uldict.end)
				error("all patterns in 'dict' must have the same length");
		}
		input_uldict.patterns[poffset] = pattern + input_uldict.start - 1;
		break;

	    case 'd':
		if (input_uldict.start == NA_INTEGER) {
			input_uldict.start = -pattern_length;
			input_uldict.width = input_uldict.end - input_uldict.start + 1;
			if (input_uldict.width < 1)
				error("'dict' contains patterns with less than %d characters",
				      -input_uldict.end);
		} else {
			if (pattern_length != -input_uldict.start)
				error("all patterns in 'dict' must have the same length");
		}
		input_uldict.patterns[poffset] = pattern;
		break;
	}
	return;
}


/****************************************************************************
 * Buffer of duplicates.
 */

static SEXP dup2unq_env0;
static SEXP unq2dup_env0;

static void init_dups_buf(int length, SEXP dup2unq_env, SEXP unq2dup_env)
{
	dup2unq_env0 = dup2unq_env;
	unq2dup_env0 = unq2dup_env;
	return;
}

static void report_dup(int poffset, int P_id)
{
	SEXP symbol, old_value, new_value;
	int old_length, new_length;

	poffset++;
	/* add a new entry to the 'dup2unq_env0' map */
	PROTECT(symbol = _SparseList_int2symb(poffset));
	PROTECT(new_value = ScalarInteger(P_id));
	defineVar(install(translateChar(symbol)), new_value, dup2unq_env0);
	UNPROTECT(2);
	/* update the 'unq2dup_env0' map */
	PROTECT(symbol = _SparseList_int2symb(P_id));
	old_value = _get_val_from_env(symbol, unq2dup_env0, 0);
	if (old_value == R_UnboundValue) {
		PROTECT(new_value = ScalarInteger(poffset));
	} else {
		old_length = LENGTH(old_value);
		new_length = old_length + 1;
		PROTECT(new_value = NEW_INTEGER(new_length));
		memcpy(INTEGER(new_value), INTEGER(old_value), old_length * sizeof(int));
		INTEGER(new_value)[old_length] = poffset;
	}
	defineVar(install(translateChar(symbol)), new_value, unq2dup_env0);
	UNPROTECT(2);
	return;
}


/****************************************************************************
 * Building the Aho-Corasick 4-ary tree
 * ====================================
 *
 * For this Aho-Corasick implementation, we take advantage of 2
 * important specifities of the dictionary (aka pattern set):
 *   1. it's a constant width dictionary (all words have the same length)
 *   2. it's based on a 4-letter alphabet
 * Because of this, the Aho-Corasick tree (which is in fact a graph if we
 * consider the failure links) can be stored in an array of ACNode elements.
 * This has the following advantages:
 *   - Speed: no need to call alloc for each new node (in the other hand
 *     memory usage is no optimal, there will be a lot of unused space in the
 *     array, this would need to be quantified though).
 *   - Can be stored in an integer vector in R: this is because the size of
 *     a node (ACNode struct) is a multiple of the size of an int (1 node = 6
 *     ints). If this was not the case, we could still use a raw vector.
 *   - Easy to serialize.
 *   - Easy to reallocate.
 * Note that the id of an ACNode element is just its offset in the array.
 */

static int actree_base_codes_buf[MAX_CHILDREN_PER_ACNODE];

#define INTS_PER_ACNODE (sizeof(ACNode) / sizeof(int))
#define MAX_ACNODEBUF_LENGTH (INT_MAX / INTS_PER_ACNODE)

static ACNode *actree_nodes_buf = NULL;
static int actree_nodes_buf_count;

SEXP debug_ACtree_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'ACtree_utils.c'\n",
		debug ? "on" : "off");
	if (debug) {
		Rprintf("[DEBUG] debug_ACtree_utils(): INTS_PER_ACNODE=%d\n",
			INTS_PER_ACNODE);
		Rprintf("[DEBUG] debug_ACtree_utils(): MAX_ACNODEBUF_LENGTH=%d\n",
			MAX_ACNODEBUF_LENGTH);
	}
#else
	Rprintf("Debug mode not available in 'ACtree_utils.c'\n");
#endif
	return R_NilValue;
}

static void init_actree_base_codes_buf()
{
	int childslot;

	for (childslot = 0; childslot < MAX_CHILDREN_PER_ACNODE; childslot++)
		actree_base_codes_buf[childslot] = -1;
	return;
}

SEXP CWdna_free_actree_nodes_buf()
{
	if (actree_nodes_buf != NULL) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] CWdna_free_actree_nodes_buf(): freeing actree_nodes_buf ... ");
		}
#endif
		free(actree_nodes_buf);
		actree_nodes_buf = NULL;
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("OK\n");
		}
#endif
	}
	return R_NilValue;
}

/*
 * We want to avoid using reallocation for actree_nodes_buf (the buffer where
 * we are going to build the Aho-Corasick tree) so we need to know what's the
 * maximum number of nodes that is needed for storing a dictionary of a given
 * width and length. First some notations:
 *   L: dictionary length i.e. number of patterns
 *   W: dictionary width i.e. number of chars per pattern
 *   A: length of the alphabet i.e. max number of children per node
 *   maxnodes: maximum number of nodes needed
 * Now here is how to get maxnodes:
 *   maxnodes = sum(from: i=0, to: i=W, of: min(A^i, L))
 */
static void alloc_actree_nodes_buf(int length, int width)
{
	int maxnodes, pow, depth;
	size_t bufsize;

	if (actree_nodes_buf != NULL) {
		// We use the on.exit() mechanism to call CWdna_free_actree_nodes_buf()
		// to free the buffer so if this mechanism is reliable we should
		// never come here. Anyway just in case...
		warning("actree_nodes_buf was not previously freed, this is anormal, please report");
		CWdna_free_actree_nodes_buf();
	}
	maxnodes = 0;
	for (depth = 0, pow = 1; depth <= width; depth++) {
		if (pow >= length)
			break;
		maxnodes += pow;
		pow *= MAX_CHILDREN_PER_ACNODE;
	}
	maxnodes += (width - depth + 1) * length;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] alloc_actree_nodes_buf(): length=%d width=%d maxnodes=%d\n",
			length, width, maxnodes);
	}
#endif
	if (maxnodes >= MAX_ACNODEBUF_LENGTH)
		error("'dict' is too big");
	bufsize = sizeof(ACNode) * maxnodes;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] alloc_actree_nodes_buf(): allocating actree_nodes_buf (bufsize=%lu) ... ",
			bufsize);
	}
#endif
	actree_nodes_buf = (ACNode *) malloc(bufsize);
	if (actree_nodes_buf == NULL)
		error("alloc_actree_nodes_buf(): failed to alloc actree_nodes_buf");
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("OK\n");
	}
#endif
	actree_nodes_buf_count = 0;
	return;
}

static void init_acnode(ACNode *node, int parent_id)
{
	ACNode *parent_node;
	int childslot;

	node->parent_id = parent_id;
	parent_node = actree_nodes_buf + parent_id;
	if (parent_node == node)	
		node->depth = 0;
	else
		node->depth = parent_node->depth + 1;
	for (childslot = 0; childslot < MAX_CHILDREN_PER_ACNODE; childslot++)
		node->child_id[childslot] = -1;
	node->flink = -1;
	node->P_id = -1;
	return;
}

static int append_acnode(int parent_id)
{
	ACNode *node;

	node = actree_nodes_buf + actree_nodes_buf_count;
	init_acnode(node, parent_id);
	return actree_nodes_buf_count++;
}

static int try_moving_to_acnode_child(int node_id, int childslot, char c)
{
	int child_id;

	if (actree_base_codes_buf[childslot] == -1) {
		actree_base_codes_buf[childslot] = c;
	} else {
		if (c != actree_base_codes_buf[childslot])
			return -1;
	}
	child_id = actree_nodes_buf[node_id].child_id[childslot];
	if (child_id != -1)
		return child_id;
	child_id = append_acnode(node_id);
	return actree_nodes_buf[node_id].child_id[childslot] = child_id;
}

static void pp_pattern(int poffset)
{
	const char *pattern;
	int n, node_id, child_id, childslot;
	char c;
	ACNode *node;

	pattern = input_uldict.patterns[poffset];
	for (n = 0, node_id = 0; n < input_uldict.width; n++, node_id = child_id) {
		c = pattern[n];
		for (childslot = 0; childslot < MAX_CHILDREN_PER_ACNODE; childslot++) {
			child_id = try_moving_to_acnode_child(node_id, childslot, c);
			if (child_id != -1)
				break;
		}
		if (child_id == -1)
			error("'dict' contains more than %d distinct letters", MAX_CHILDREN_PER_ACNODE);
	}
	node = actree_nodes_buf + node_id;
	if (node->P_id == -1)
		node->P_id = poffset + 1;
	else
		report_dup(poffset, node->P_id);
	return;
}

static void build_actree(SEXP dup2unq_env, SEXP unq2dup_env)
{
	int poffset;

	init_dups_buf(input_uldict.length, dup2unq_env, unq2dup_env);
	init_actree_base_codes_buf();
	alloc_actree_nodes_buf(input_uldict.length, input_uldict.width);
	append_acnode(0);
	for (poffset = 0; poffset < input_uldict.length; poffset++)
		pp_pattern(poffset);
	return;
}


/****************************************************************************
 * Turning our local data structures into R objects (SEXP)
 * =======================================================
 */

/*
 * Turn the dictionary width into an R integer vector of length 1.
 */
static SEXP oneint_asINTEGER(int oneint)
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = oneint;
	UNPROTECT(1);
	return ans;
}

/*
 * Turn the Aho-Corasick 4-ary tree stored in actree_nodes_buf into an R eXternal
 * Pointer.
 */
static SEXP actree_nodes_buf_asXP()
{
	SEXP ans, tag;
	int tag_length;

	PROTECT(ans = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	tag_length = actree_nodes_buf_count * INTS_PER_ACNODE;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), actree_nodes_buf, tag_length * sizeof(int));
	R_SetExternalPtrTag(ans, tag);
	UNPROTECT(2);
	return ans;
}

/*
 * Turn the actree_base_codes_buf array into an R integer vector of length
 * MAX_CHILDREN_PER_ACNODE.
 */
static SEXP actree_base_codes_buf_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(MAX_CHILDREN_PER_ACNODE));
	memcpy(INTEGER(ans), actree_base_codes_buf, MAX_CHILDREN_PER_ACNODE * sizeof(int));
	UNPROTECT(1);
	return ans;
}

static SEXP stats_asLIST()
{
	SEXP ans, ans_names, ans_elt;
	const char *min_width_name, *max_width_name;
	int min_width_val, max_width_val;

	if (input_uldict.cropping_type != 'a' && input_uldict.cropping_type != 'b')
		return NEW_LIST(0);

	if (input_uldict.cropping_type == 'a') {
		min_width_name = "tail.min.width";
		max_width_name = "tail.max.width";
		min_width_val = input_uldict.tail_min_width;
		max_width_val = input_uldict.tail_max_width;
	} else {
		min_width_name = "head.min.width";
		max_width_name = "head.max.width";
		min_width_val = input_uldict.head_min_width;
		max_width_val = input_uldict.head_max_width;
	}
	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar(min_width_name));
	SET_STRING_ELT(ans_names, 1, mkChar(max_width_name));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "min.width" element */
	PROTECT(ans_elt = oneint_asINTEGER(min_width_val));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "max.width" element */
	PROTECT(ans_elt = oneint_asINTEGER(max_width_val));
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

/*
 * Return an R list with the following elements:
 *   - width: dictionary width (single integer);
 *   - actree_nodes_xp: "externalptr" object pointing to the Aho-Corasick 4-ary
 *         tree built from 'dict';
 *   - actree_base_codes: integer vector containing the MAX_CHILDREN_PER_ACNODE character
 *         codes (ASCII) attached to the MAX_CHILDREN_PER_ACNODE child slots of any
 *         node in the tree pointed by actree_nodes_xp;
 *   - stats: a list containing some stats about the input data.
 */
static SEXP uldna_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(4));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(4));
	SET_STRING_ELT(ans_names, 0, mkChar("width"));
	SET_STRING_ELT(ans_names, 1, mkChar("actree_nodes_xp"));
	SET_STRING_ELT(ans_names, 2, mkChar("actree_base_codes"));
	SET_STRING_ELT(ans_names, 3, mkChar("stats"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "width" element */
	PROTECT(ans_elt = oneint_asINTEGER(input_uldict.width));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "actree_nodes_xp" element */
	PROTECT(ans_elt = actree_nodes_buf_asXP());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	/* set the "actree_base_codes" element */
	PROTECT(ans_elt = actree_base_codes_buf_asINTEGER());
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	/* set the "stats" element */
	PROTECT(ans_elt = stats_asLIST());
	SET_ELEMENT(ans, 3, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry points for building the Aho-Corasick 4-ary tree from the input
 * dictionary (preprocessing). The input dictionary must be:
 *   - of length >= 1
 *   - of constant width (i.e. all words have the same length)
 ****************************************************************************/

/*
 * .Call entry point: "CWdna_pp_STRSXP"
 *
 * Argument:
 *   'dict': a string vector (aka character vector) containing the input
 *           sequences
 *   'start': single >= 1, <= -1 or NA integer
 *   'end': single >= 1, <= -1 or NA integer
 *
 * See uldna_asLIST() for a description of the returned SEXP.
 */
SEXP CWdna_pp_STRSXP(SEXP dict, SEXP start, SEXP end, SEXP dup2unq_env, SEXP unq2dup_env)
{
	int dict_len, poffset;
	SEXP dict_elt;

	dict_len = LENGTH(dict);
	alloc_input_uldict(dict_len, INTEGER(start)[0], INTEGER(end)[0]);
	for (poffset = 0; poffset < dict_len; poffset++) {
		dict_elt = STRING_ELT(dict, poffset);
		if (dict_elt == NA_STRING)
			error("'dict' contains NAs");
		add_subpattern_to_input_uldict(poffset, CHAR(dict_elt), LENGTH(dict_elt));
	}
	build_actree(dup2unq_env, unq2dup_env);
	return uldna_asLIST();
}

/*
 * .Call entry point: "CWdna_pp_XStringSet"
 *
 * Argument:
 *   'dict': a DNAStringSet object containing the input sequences
 *   'start': single >= 1, <= -1 or NA integer
 *   'end': single >= 1, <= -1 or NA integer
 *
 * See uldna_asLIST() for a description of the returned SEXP.
 */
SEXP CWdna_pp_XStringSet(SEXP dict, SEXP start, SEXP end, SEXP dup2unq_env, SEXP unq2dup_env)
{
	int dict_len, poffset;
	CachedXStringSet cached_dict;
	RoSeq pattern;

	dict_len = _get_XStringSet_length(dict);
	alloc_input_uldict(dict_len, INTEGER(start)[0], INTEGER(end)[0]);
	cached_dict = _new_CachedXStringSet(dict);
	for (poffset = 0; poffset < dict_len; poffset++) {
		pattern = _get_CachedXStringSet_elt_asRoSeq(&cached_dict, poffset);
		add_subpattern_to_input_uldict(poffset, pattern.elts, pattern.nelt);
	}
	build_actree(dup2unq_env, unq2dup_env);
	return uldna_asLIST();
}

