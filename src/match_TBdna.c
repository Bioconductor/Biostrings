/****************************************************************************
 *              Aho-Corasick for constant width DNA dictionary              *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc()/Srealloc() */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


/****************************************************************************/
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

static IntBuf dups_buf;

static void init_dups_buf(int length)
{
	dups_buf = _new_IntBuf(length, length);
	return;
}

static void report_dup(int poffset, int P_id)
{
	dups_buf.elts[poffset] = P_id;
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

#define ALPHABET_LENGTH 4

static int actree_base_codes_buf[ALPHABET_LENGTH];

typedef struct acnode {
	int parent_id;
	int depth;
	int child_id[ALPHABET_LENGTH];
	int flink;
	int P_id;
} ACNode;

#define INTS_PER_ACNODE (sizeof(ACNode) / sizeof(int))
#define MAX_ACNODEBUF_LENGTH (INT_MAX / INTS_PER_ACNODE)

static ACNode *actree_nodes_buf = NULL;
static int actree_nodes_buf_count;

SEXP debug_match_TBdna()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_TBdna.c'\n", debug ? "on" : "off");
	if (debug) {
		Rprintf("[DEBUG] debug_match_TBdna(): INTS_PER_ACNODE=%d\n",
			INTS_PER_ACNODE);
		Rprintf("[DEBUG] debug_match_TBdna(): MAX_ACNODEBUF_LENGTH=%d\n",
			MAX_ACNODEBUF_LENGTH);
	}
#else
	Rprintf("Debug mode not available in 'match_TBdna.c'\n");
#endif
	return R_NilValue;
}

static void init_actree_base_codes_buf()
{
	int childslot;

	for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++)
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
		pow *= ALPHABET_LENGTH;
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
	for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++)
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
		for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++) {
			child_id = try_moving_to_acnode_child(node_id, childslot, c);
			if (child_id != -1)
				break;
		}
		if (child_id == -1)
			error("'dict' contains more than %d distinct letters", ALPHABET_LENGTH);
	}
	node = actree_nodes_buf + node_id;
	if (node->P_id == -1)
		node->P_id = poffset + 1;
	else
		report_dup(poffset, node->P_id);
	return;
}

static void build_actree()
{
	int poffset;

	init_dups_buf(input_uldict.length);
	init_actree_base_codes_buf();
	alloc_actree_nodes_buf(input_uldict.length, input_uldict.width);
	append_acnode(0);
	for (poffset = 0; poffset < input_uldict.length; poffset++)
		pp_pattern(poffset);
	return;
}


/****************************************************************************
 * Match reporting
 * ===============
 */

static int match_reporting_mode; // 0, 1 or 2
static IntBuf match_count; // used when mode == 0 and initialized when mode == 2
static IntBBuf ends_bbuf;  // used when mode >= 1

static void init_match_reporting(int not_tail, int is_count_only, int length)
{
	match_reporting_mode = is_count_only ? (not_tail ? 0 : 2) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		match_count = _new_IntBuf(length, length);
	if (match_reporting_mode >= 1)
		ends_bbuf = _new_IntBBuf(length, length);
	return;
}

static void report_match(int poffset, int end)
{
	IntBuf *ends_buf;

	if (match_reporting_mode == 0) {
		match_count.elts[poffset]++;
		return;
	}
	ends_buf = ends_bbuf.elts + poffset;
	_IntBuf_insert_at(ends_buf, ends_buf->nelt, end);
	return;
}

static void report_matches_for_dups(const int *dups, int length)
{
	int poffset, *val;
	IntBuf *ends_buf;

	if (match_reporting_mode == 0) {
		for (poffset = 0, val = match_count.elts;
		     poffset < length;
		     poffset++, val++, dups++) {
			if (*dups == 0)
				continue;
			*val = *(match_count.elts + *dups - 1);
		}
		return;
	}
	for (poffset = 0, ends_buf = ends_bbuf.elts;
	     poffset < length;
	     poffset++, ends_buf++, dups++) {
		if (*dups == 0)
			continue;
		*ends_buf = *(ends_bbuf.elts + *dups - 1);
	}
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

static int code2childoffset_chrtrtable[CHRTRTABLE_LENGTH];

static int get_child_id(const ACNode *node, char c)
{
	int offset;

	offset = code2childoffset_chrtrtable[(unsigned char) c];
	if (offset == -1)
		return -1;
	return node->child_id[offset];
}

/*
 * We use the child slots for storing the shortcuts so it's important to
 * remember that a child slot doesn't necessary contain the ID of a child
 * node anymore: it can be any other node in the tree. In fact, it can't be
 * whatever other node either: its depth can't be greater than the depth of
 * the referring node. This property provides an efficient way to know whether
 * N1 -> N2 is a parent-to-child link or a shortcut:
 *   parent-to-child: depth(N2) == depth(N1) + 1
 *   shortcut: depth(N2) <= depth(N1)
 * Note that this trick is not needed by the current implementation of the
 * follow_string() function.
 */
static void set_shortcut(ACNode *basenode, char c, int node_id)
{
	int offset, *slot;

	offset = code2childoffset_chrtrtable[(unsigned char) c];
	if (offset == -1)
		return;
	slot = basenode->child_id + offset;
	if (*slot == -1)
		*slot = node_id;
	return;
}

static int follow_string(ACNode *node0, const int *base_codes, const char *S, int nS)
{
	int basenode_id, node_id, child_id, n, subcall_nS;
	ACNode *basenode, *node;
	static int rec_level = 0;
#ifdef DEBUG_BIOSTRINGS
	char format[20], pathbuf[2000];
#endif

	basenode_id = 0;
	basenode = node0;
	for (n = 0; n < nS; n++, S++) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] follow_string():");
			sprintf(format, "%%%ds", 1 + 2*rec_level);
			Rprintf(format, " ");
			snprintf(pathbuf, basenode->depth + 1, "%s", S - basenode->depth);
			Rprintf("On basenode_id=%d (basepath=%s), reading S[%d]=%c\n", basenode_id, pathbuf, n, *S);
		}
#endif
		node_id = basenode_id;
		node = basenode;
		while (1) {
			child_id = get_child_id(node, *S);
			if (child_id != -1) {
				node_id = child_id;
				node = node0 + node_id;
				break;
			}
			if (node_id == 0) {
				node = node0; /* node == node0 */
				break;
			}
			if (node->flink == -1) {
				rec_level++;
				subcall_nS = node->depth - 1;
				node->flink = follow_string(node0, base_codes, S - subcall_nS, subcall_nS);
				rec_level--;
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					Rprintf("[DEBUG] follow_string():");
					Rprintf(format, " ");
					Rprintf("setting failure link %d -> %d\n", node_id, node->flink);
				}
#endif
			}
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] follow_string():");
				Rprintf(format, " ");
				Rprintf("following failure link %d -> %d\n", node_id, node->flink);
			}
#endif
			node_id = node->flink;
			node = node0 + node_id;
		}
		set_shortcut(basenode, *S, node_id);
		basenode_id = node_id;
		basenode = node0 + basenode_id;
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] follow_string():");
			Rprintf(format, " ");
			Rprintf("moving to basenode %d\n", basenode_id);
		}
#endif
		// Finding a match cannot happen during a nested call to follow_string()
		// so there is no need to check that rec_level is 0
		if (basenode->P_id != -1)
			report_match(basenode->P_id - 1, n + 1);
	}
	return basenode_id;
}

static void CWdna_exact_search(ACNode *node0, const int *base_codes, const char *S, int nS)
{
	_init_chrtrtable(base_codes, ALPHABET_LENGTH, code2childoffset_chrtrtable);
	follow_string(node0, base_codes, S, nS);
	return;
}


/****************************************************************************
 * Inexact matching on the tails of the TBdna_PDict object
 * =======================================================
 */

static void TBdna_match_pattern_tail(RoSeq pattern_tail, RoSeq S,
		int poffset, const int *dups, int dups_len,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int dup0, i, end;
	IntBuf *ends_buf, *ends_buf0;

	ends_buf = ends_buf0 = ends_bbuf.elts + poffset;
	dup0 = dups[poffset];
	if (dup0 != 0)
		ends_buf0 = ends_bbuf.elts + dup0 - 1;
	for (i = 0; i < ends_buf0->nelt; i++) {
		end = ends_buf0->elts[i];
		if (_is_matching(pattern_tail, S, end, max_mm, fixedP, fixedS)) {
			/* Match */
			if (is_count_only) {
				match_count.elts[poffset]++;
				continue;
			}
			if (dup0 == 0) {
				ends_buf0->elts[i] += pattern_tail.nelt;
				continue;
			}
			end += pattern_tail.nelt;
			_IntBuf_insert_at(ends_buf, ends_buf->nelt, end);
			continue;
		}
		/* Mismatch */
		if (is_count_only)
			continue;
		if (dup0 != 0)
			continue;
		/* We need to shrink the buffer we are walking on! This is safe
		 * because shrinking a IntBuf object should never trigger reallocation.
		 */
		_IntBuf_delete_at(ends_buf0, i--);
	}
	return;
}

static void TBdna_match_tail(SEXP tail, RoSeq S,
		const int *dups, int dups_len,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int poffset;
	CachedXStringSet cached_tail;
	RoSeq pattern_tail;

	cached_tail = _new_CachedXStringSet(tail);
	// The duplicated must be treated BEFORE the first pattern they
	// duplicate, hence we must walk from last to first.
	for (poffset = dups_len - 1; poffset >= 0; poffset--) {
		pattern_tail = _get_CachedXStringSet_elt_asRoSeq(&cached_tail, poffset);
		TBdna_match_pattern_tail(pattern_tail, S,
			poffset, dups, dups_len,
			max_mm, fixedP, fixedS, is_count_only);
	}
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
 * ALPHABET_LENGTH.
 */
static SEXP actree_base_codes_buf_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ALPHABET_LENGTH));
	memcpy(INTEGER(ans), actree_base_codes_buf, ALPHABET_LENGTH * sizeof(int));
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
 *   - actree_base_codes: integer vector containing the ALPHABET_LENGTH character
 *         codes (ASCII) attached to the ALPHABET_LENGTH child slots of any
 *         node in the tree pointed by actree_nodes_xp;
 *   - dups: an integer vector of the same length as 'dict' containing the
 *         duplicate info;
 *   - stats: a list containing some stats about the input data.
 */
static SEXP uldna_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(5));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(5));
	SET_STRING_ELT(ans_names, 0, mkChar("width"));
	SET_STRING_ELT(ans_names, 1, mkChar("actree_nodes_xp"));
	SET_STRING_ELT(ans_names, 2, mkChar("actree_base_codes"));
	SET_STRING_ELT(ans_names, 3, mkChar("dups"));
	SET_STRING_ELT(ans_names, 4, mkChar("stats"));
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

	/* set the "dups" element */
	PROTECT(ans_elt = _IntBuf_asINTEGER(&dups_buf));
	SET_ELEMENT(ans, 3, ans_elt);
	UNPROTECT(1);

	/* set the "stats" element */
	PROTECT(ans_elt = stats_asLIST());
	SET_ELEMENT(ans, 4, ans_elt);
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
SEXP CWdna_pp_STRSXP(SEXP dict, SEXP start, SEXP end)
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
	build_actree();
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
SEXP CWdna_pp_XStringSet(SEXP dict, SEXP start, SEXP end)
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
	build_actree();
	return uldna_asLIST();
}


/****************************************************************************
 * .Call entry point: "match_TBdna"
 *
 * Arguments:
 *   'actree_nodes_xp': uldna_pdict@actree@nodes@xp
 *   'actree_base_codes': uldna_pdict@actree@base_codes
 *   'pdict_dups': pdict@dups
 *   'pdict_head_XStringSet': pdict@head (can be NULL)
 *   'pdict_tail_XStringSet': pdict@tail (can be NULL)
 *   'subject_XString': subject
 *   'max_mismatch': max.mismatch (max nb of mismatches in the tail)
 *   'fixed': logical vector of length 2
 *   'count_only': TRUE or FALSE
 *   'envir': environment to be populated with the matches
 *
 * Return an R object that will be assigned to the 'ends' slot of the
 * MIndex object returned by matchPDict(). Refer to the description
 * of this slot in the matchPDict.R file for the details.
 *
 ****************************************************************************/

SEXP match_TBdna(SEXP actree_nodes_xp, SEXP actree_base_codes,
		SEXP pdict_dups, SEXP pdict_head_XStringSet, SEXP pdict_tail_XStringSet,
		SEXP subject_XString,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	ACNode *actree_nodes;
	RoSeq S;
	int is_count_only, no_head, no_tail;

	actree_nodes = (ACNode *) INTEGER(R_ExternalPtrTag(actree_nodes_xp));
	S = _get_XString_asRoSeq(subject_XString);
	is_count_only = LOGICAL(count_only)[0];
	no_head = pdict_head_XStringSet == R_NilValue;
	no_tail = pdict_tail_XStringSet == R_NilValue;

	if (!no_head)
		error("matchPDict() doesn't support PDict objects with a head yet, sorry!");
	init_match_reporting(no_tail, is_count_only, LENGTH(pdict_dups));
	CWdna_exact_search(actree_nodes, INTEGER(actree_base_codes), S.elts, S.nelt);
	if (no_tail) {
		report_matches_for_dups(INTEGER(pdict_dups), LENGTH(pdict_dups));
	} else {
		TBdna_match_tail(pdict_tail_XStringSet, S,
			INTEGER(pdict_dups), LENGTH(pdict_dups),
			INTEGER(max_mismatch)[0], LOGICAL(fixed)[0], LOGICAL(fixed)[1],
			is_count_only);
	}
	if (is_count_only)
		return _IntBuf_asINTEGER(&match_count);
	if (envir == R_NilValue)
		return _IntBBuf_asLIST(&ends_bbuf, 1);
	return _IntBBuf_toEnvir(&ends_bbuf, envir, 1);
}


/****************************************************************************
 * Some additional utility functions used fast data extraction from the
 * MIndex object returned by matchPDict().
 */

/* 'symbol' must be a CHARSXP */
static SEXP getSymbolVal(SEXP symbol, SEXP envir)
{
	SEXP ans;

	/* The following code was inspired by R's do_get() code.
	 * Note that do_get() doesn't use PROTECT at all and so do we...
	 */
	ans = findVar(install(translateChar(symbol)), envir);
	if (ans == R_UnboundValue)
		error("Biostrings internal error in getSymbolVal(): unbound value");
	if (TYPEOF(ans) == PROMSXP)
		ans = eval(ans, envir);
	if (ans != R_NilValue && NAMED(ans) == 0)
		SET_NAMED(ans, 1);
	return ans;
}

/*
 * 'e1' must be an INTSXP.
 * addInt() must ALWAYS duplicate 'e1', even when e2 = 0!
 */
static SEXP addInt(SEXP e1, int e2)
{
	SEXP ans;
	int i, *val;

	PROTECT(ans = duplicate(e1));
	for (i = 0, val = INTEGER(ans); i < LENGTH(ans); i++, val++)
		*val += e2;
	UNPROTECT(1);
	return ans;
}

/*
 * Only elements in 'x' that are integer vectors are shifted.
 */
SEXP shiftListOfInts(SEXP x, SEXP shift)
{
	SEXP ans, ans_elt;
	int shiftval, i, j, *val;

	PROTECT(ans = duplicate(x));
	shiftval = INTEGER(shift)[0];
	for (i = 0; i < LENGTH(ans); i++) {
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt))
			continue;
		for (j = 0, val = INTEGER(ans_elt); j < LENGTH(ans_elt); j++, val++)
			*val += shiftval;
	}
	UNPROTECT(1);
	return ans;
}

/* All the keys in ends_envir must be representing integers left-padded with 0s
 * so they have the same length. This works properly:
     library(Biostrings)
     ends_envir <- new.env(parent = emptyenv())
     ends_envir[['0000000010']] <- -2:1
     ends_envir[['0000000004']] <- 9:6
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], TRUE, PACKAGE="Biostrings")
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 */
SEXP extract_endIndex(SEXP ends_envir, SEXP shift, SEXP names, SEXP all_names)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IntBuf poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(ends_envir, 1));
	poffsets = _CHARACTER_asIntBuf(symbols, -1);
	if (LOGICAL(all_names)[0]) {
		PROTECT(ans = NEW_LIST(LENGTH(names)));
		for (i = 0; i < poffsets.nelt; i++) {
			end = getSymbolVal(STRING_ELT(symbols, i), ends_envir);
			PROTECT(ans_elt = addInt(end, INTEGER(shift)[0]));
			SET_ELEMENT(ans, poffsets.elts[i], ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(names));
		UNPROTECT(1);
	} else {
		//poffsets_order = _new_IntBuf(poffsets.nelt, 0);
		//get_intorder(poffsets.nelt, poffsets.elts, poffsets_order.elts);
		//poffsets_order.nelt = poffsets.nelt; /* = poffsets_order.buflength */
		PROTECT(ans = NEW_LIST(poffsets.nelt));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.nelt));
		for (i = 0; i < poffsets.nelt; i++) {
			//j = poffsets_order.elts[i];
			j = i;
			end = getSymbolVal(STRING_ELT(symbols, j), ends_envir);
			SET_ELEMENT(ans, i, addInt(end, INTEGER(shift)[0]));
			SET_STRING_ELT(ans_names, i, duplicate(STRING_ELT(names, poffsets.elts[j])));
		}
		SET_NAMES(ans, ans_names);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

