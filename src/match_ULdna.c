/****************************************************************************
 *              Aho-Corasick for uniform-length DNA dictionary              *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a uniform-length dictionary is a non-empty set of non-empty        *
 * strings of the same length.                                              *
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
 * contained in the uniform length dictionary provided by the user.
 * These patterns can only be accessed for reading!
 */

typedef struct uldict {
	const char **patterns;	// array of pointers to the read-only patterns
	int length;		// number of patterns
	int width;		// number of chars per pattern
	int truncate_on_add;	// should patterns be truncated when added?
	int min_width;
	int max_width;
} ULdict;

static ULdict input_uldict;

static void alloc_input_uldict(int length, SEXP width)
{
	input_uldict.patterns = Salloc((long) length, const char *);
	input_uldict.length = length;
	input_uldict.truncate_on_add = width != R_NilValue;
	if (input_uldict.truncate_on_add) {
		input_uldict.width = INTEGER(width)[0];
		input_uldict.min_width = input_uldict.max_width = -1;
	} else {
		input_uldict.width = -1;
	}
	return;
}

static void add_pattern_to_input_uldict(int poffset, const char *pattern, int pattern_length)
{
	if (input_uldict.truncate_on_add) {
		if (pattern_length < input_uldict.width)
			error("'dict' contains patterns with less than %d characters",
			      input_uldict.width);
		if (input_uldict.min_width == -1) {
			input_uldict.min_width = input_uldict.max_width = pattern_length;
		} else {
			if (pattern_length < input_uldict.min_width)
				input_uldict.min_width = pattern_length;
			if (pattern_length > input_uldict.max_width)
				input_uldict.max_width = pattern_length;
		}	
	} else if (input_uldict.width == -1) {
		if (pattern_length == 0)
			error("first pattern in 'dict' is empty");
		input_uldict.width = pattern_length;
	} else {
		if (pattern_length != input_uldict.width)
			error("all patterns in 'dict' must have the same length");
	}
	input_uldict.patterns[poffset] = pattern;
	return;
}


/****************************************************************************
 * Buffer of duplicates.
 */

static IBuf dups_buf;

static void init_dups_buf(int length)
{
	_IBuf_init(&dups_buf, length, length);
	return;
}

static void report_dup(int poffset, int P_id)
{
	dups_buf.vals[poffset] = P_id;
	return;
}


/****************************************************************************
 * Building the Aho-Corasick 4-ary tree
 * ====================================
 *
 * For this Aho-Corasick implementation, we take advantage of 2
 * important specifities of the dictionary (aka pattern set):
 *   1. it's a uniform length dictionary (all words have the same length)
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

SEXP match_ULdna_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_ULdna.c'\n", debug ? "on" : "off");
	if (debug) {
		Rprintf("[DEBUG] match_ULdna_debug(): INTS_PER_ACNODE=%d\n",
			INTS_PER_ACNODE);
		Rprintf("[DEBUG] match_ULdna_debug(): MAX_ACNODEBUF_LENGTH=%d\n",
			MAX_ACNODEBUF_LENGTH);
	}
#else
	Rprintf("Debug mode not available in 'match_ULdna.c'\n");
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

SEXP ULdna_free_actree_nodes_buf()
{
	if (actree_nodes_buf != NULL) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] ULdna_free_actree_nodes_buf(): freeing actree_nodes_buf ... ");
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
		// We use the on.exit() mechanism to call ULdna_free_actree_nodes_buf()
		// to free the buffer so if this mechanism is reliable we should
		// never come here. Anyway just in case...
		warning("actree_nodes_buf was not previously freed, this is anormal, please report");
		ULdna_free_actree_nodes_buf();
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
static IBuf match_count; // used when mode == 0 and initialized when mode == 2
static IBBuf ends_bbuf;  // used when mode >= 1

static void init_match_reporting(int not_tail, int is_count_only, int length)
{
	match_reporting_mode = is_count_only ? (not_tail ? 0 : 2) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		_IBuf_init(&match_count, length, length);
	if (match_reporting_mode >= 1)
		_IBBuf_init(&ends_bbuf, length, length);
	return;
}

static void report_match(int poffset, int end)
{
	IBuf *ends_buf;

	if (match_reporting_mode == 0) {
		match_count.vals[poffset]++;
		return;
	}
	ends_buf = ends_bbuf.ibufs + poffset;
	_IBuf_insert_at(ends_buf, ends_buf->count, end);
	return;
}

static void report_matches_for_dups(const int *dups, int length)
{
	int poffset, *val;
	IBuf *ends_buf;

	if (match_reporting_mode == 0) {
		for (poffset = 0, val = match_count.vals;
		     poffset < length;
		     poffset++, val++, dups++) {
			if (*dups == 0)
				continue;
			*val = *(match_count.vals + *dups - 1);
		}
		return;
	}
	for (poffset = 0, ends_buf = ends_bbuf.ibufs;
	     poffset < length;
	     poffset++, ends_buf++, dups++) {
		if (*dups == 0)
			continue;
		*ends_buf = *(ends_bbuf.ibufs + *dups - 1);
	}
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

static int code2childoffset_lkup[256];

static int get_child_id(const ACNode *node, char c)
{
	int offset;

	offset = code2childoffset_lkup[(unsigned char) c];
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

	offset = code2childoffset_lkup[(unsigned char) c];
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

static void ULdna_exact_search(ACNode *node0, const int *base_codes, const char *S, int nS)
{
	_init_code2offset_lkup(base_codes, ALPHABET_LENGTH, code2childoffset_lkup);
	follow_string(node0, base_codes, S, nS);
	return;
}


/****************************************************************************
 * Inexact matching on the tails of the TailedULdna_PDict object
 * =============================================================
 */

static int inexact_match(const char *P, int nP, const char *S, int nS, int max_mm)
{
	int min_pm, mm, pm, i;

	min_pm = nP - max_mm;
	if (min_pm <= 0)
		return 1;
	if (nP > nS) {
		max_mm -= nP - nS;
		if (max_mm < 0)
			return 0;
		nP = nS;
	}
	mm = pm = 0;
	/* From here we have:
	 *   0 = mm <= max_mm < nP <= nS
	 *   0 = pm < min_pm <= nP <= nS
         *   min_pm + max_mm = nP
	 */
	for (i = 0; i < nP; i++, P++, S++) {
		if (*P == *S) {
			if (++pm >= min_pm)
				return 1;
		} else {
			if (++mm > max_mm)
				return 0;
		}
	}
	error("Biostrings internal error in inexact_match(): we should never be here");
	return -1;
}

static void ULdna_match_tail(const char *tail, int tail_len, const char *S, int nS,
			int poffset, const int *dups, int dups_len,
			int max_mm, int is_count_only)
{
	int dup0, i, end;
	IBuf *ends_buf, *ends_buf0;

	ends_buf = ends_buf0 = ends_bbuf.ibufs + poffset;
	dup0 = dups[poffset];
	if (dup0 != 0)
		ends_buf0 = ends_bbuf.ibufs + dup0 - 1;
	for (i = 0; i < ends_buf0->count; i++) {
		end = ends_buf0->vals[i];
		if (inexact_match(tail, tail_len, S + end, nS - end, max_mm)) {
			/* Match */
			if (is_count_only) {
				match_count.vals[poffset]++;
				continue;
			}
			if (dup0 == 0) {
				ends_buf0->vals[i] += tail_len;
				continue;
			}
			end += tail_len;
			_IBuf_insert_at(ends_buf, ends_buf->count, end);
			continue;
		}
		/* Mismatch */
		if (is_count_only)
			continue;
		if (dup0 != 0)
			continue;
		/* We need to shrink the buffer we are walking on! This is safe
		 * because shrinking a IBuf object will never trigger reallocation.
		 */
		_IBuf_delete_at(ends_buf0, i--);
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

	if (!input_uldict.truncate_on_add)
		return NEW_LIST(0);

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("min.width"));
	SET_STRING_ELT(ans_names, 1, mkChar("max.width"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "min.width" element */
	PROTECT(ans_elt = oneint_asINTEGER(input_uldict.min_width));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "max.width" element */
	PROTECT(ans_elt = oneint_asINTEGER(input_uldict.max_width));
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
	PROTECT(ans_elt = _IBuf_asINTEGER(&dups_buf));
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
 *   - uniform-length (i.e. all words have the same length)
 ****************************************************************************/

/*
 * .Call entry point: "ULdna_pp_StrVect"
 *
 * Argument:
 *   'dict': a string vector (aka character vector) containing the
 *       uniform-length dictionary
 *   'width': if not NULL then truncate the input patterns to keep the first
 *       'width' nucleotides ('width' must be a single integer)
 *
 * See uldna_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_pp_StrVect(SEXP dict, SEXP width)
{
	int poffset;
	SEXP dict_elt;

	alloc_input_uldict(LENGTH(dict), width);
	for (poffset = 0; poffset < LENGTH(dict); poffset++) {
		dict_elt = STRING_ELT(dict, poffset);
		add_pattern_to_input_uldict(poffset, CHAR(dict_elt), LENGTH(dict_elt));
	}
	build_actree();
	return uldna_asLIST();
}

/*
 * .Call entry point: "ULdna_pp_BStringList"
 *
 * Argument:
 *   'dict': a list of (pattern@data@xp, pattern@offset, pattern@length)
 *           triplets containing the uniform-length dictionary
 *
 * See uldna_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_pp_BStringList(SEXP dict, SEXP width)
{
	error("Not ready yet!\n");
	return uldna_asLIST();
}

/*
 * .Call entry point: "ULdna_pp_views"
 *
 * Arguments:
 *   'dict_subj_xp': subject(dict)@data@xp
 *   'dict_subj_offset': subject(dict)@offset
 *   'dict_start': start(dict)
 *   'dict_end': end(dict)
 * Note that 'dict_start' and 'dict_end' describe a set of views on subject(dict)
 * with no "out of limits" views.
 *
 * See uldna_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_pp_views(SEXP dict_subj_xp, SEXP dict_subj_offset, SEXP dict_subj_length,
		SEXP dict_start, SEXP dict_end, SEXP width)
{
	int subj_offset, subj_length, dict_length;
	const Rbyte *subj;
	int poffset, *view_start, *view_end, view_offset;

	subj_offset = INTEGER(dict_subj_offset)[0];
	subj = RAW(R_ExternalPtrTag(dict_subj_xp)) + subj_offset;
	subj_length = INTEGER(dict_subj_length)[0];
	dict_length = LENGTH(dict_start); // must be the same as LENGTH(dict_end)
	alloc_input_uldict(dict_length, width);
	for (poffset = 0, view_start = INTEGER(dict_start), view_end = INTEGER(dict_end);
	     poffset < dict_length;
	     poffset++, view_start++, view_end++) {
		view_offset = *view_start - 1;
		if (view_offset < 0 || *view_end > subj_length)
			error("'dict' has out of limits views");
		add_pattern_to_input_uldict(poffset,
			(char *) (subj + view_offset), *view_end - view_offset);
	}
	build_actree();
	return uldna_asLIST();
}


/****************************************************************************
 * .Call entry point: "match_TailedULdna"
 *
 * Arguments:
 *   'actree_nodes_xp': uldna_pdict@actree@nodes@xp
 *   'actree_base_codes': uldna_pdict@actree@base_codes
 *   'pdict_dups': pdict@dups
 *   'pdict_tail_bstrings': pdict@tail@bstrings or NULL
 *   'subject_BString': subject
 *   'max_mismatch': max.mismatch (max nb of mismatches in the tail)
 *   'count_only': TRUE or FALSE
 *   'envir': environment to be populated with the matches
 *
 * Return an R object that will be assigned to the 'ends' slot of the
 * P2Views object returned by matchPDict(). Refer to the description
 * of this slot in the matchPDict.R file for the details.
 *
 ****************************************************************************/

SEXP match_TailedULdna(SEXP actree_nodes_xp, SEXP actree_base_codes,
                SEXP pdict_dups, SEXP pdict_tail_bstrings,
		SEXP subject_BString,
		SEXP max_mismatch, SEXP count_only, SEXP envir)
{
	ACNode *actree_nodes;
	const char *S, *tail;
	int nS, no_tail, is_count_only, tail_len, poffset;
	SEXP tail_BString;

	actree_nodes = (ACNode *) INTEGER(R_ExternalPtrTag(actree_nodes_xp));
	S = get_BString_seq(subject_BString, &nS);
	no_tail = pdict_tail_bstrings == R_NilValue;
	is_count_only = LOGICAL(count_only)[0];

	init_match_reporting(no_tail, is_count_only, LENGTH(pdict_dups));
	ULdna_exact_search(actree_nodes, INTEGER(actree_base_codes), S, nS);
	if (no_tail) {
		report_matches_for_dups(INTEGER(pdict_dups), LENGTH(pdict_dups));
	} else {
		/* The duplicated must be treated BEFORE the first pattern they
		 * duplicate, hence we must walk from last to first.
		 */
		for (poffset = LENGTH(pdict_tail_bstrings) - 1; poffset >= 0; poffset--) {
			tail_BString = VECTOR_ELT(pdict_tail_bstrings, poffset);
			tail = get_BString_seq(tail_BString, &tail_len);
			ULdna_match_tail(tail, tail_len, S, nS,
				poffset, INTEGER(pdict_dups), LENGTH(pdict_dups),
				INTEGER(max_mismatch)[0], is_count_only);
		}
	}
	if (is_count_only)
		return _IBuf_asINTEGER(&match_count);
	if (envir == R_NilValue)
		return _IBBuf_asLIST(&ends_bbuf, 1);
	return _IBBuf_toEnvir(&ends_bbuf, envir, 1);
}


/****************************************************************************
 * Some additional utility functions used fast data extraction from the
 * P2Views object returned by matchPDict().
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
     .Call("extract_p2end", ends_envir, 0L, letters[1:10], TRUE, PACKAGE="Biostrings")
     .Call("extract_p2end", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("extract_p2end", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 */
SEXP extract_p2end(SEXP ends_envir, SEXP shift, SEXP pids, SEXP all_pids)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IBuf poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(ends_envir, 1));
	poffsets = _CHARACTER_asIBuf(symbols, -1);
	if (LOGICAL(all_pids)[0]) {
		PROTECT(ans = NEW_LIST(LENGTH(pids)));
		for (i = 0; i < poffsets.count; i++) {
			end = getSymbolVal(STRING_ELT(symbols, i), ends_envir);
			PROTECT(ans_elt = addInt(end, INTEGER(shift)[0]));
			SET_ELEMENT(ans, poffsets.vals[i], ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(pids));
		UNPROTECT(1);
	} else {
		//_IBuf_init(&poffsets_order, poffsets.count, 0);
		//get_intorder(poffsets.vals, poffsets_order.vals, poffsets.count);
		//poffsets_order.count = poffsets.count; /* = poffsets_order.maxcount */
		PROTECT(ans = NEW_LIST(poffsets.count));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.count));
		for (i = 0; i < poffsets.count; i++) {
			//j = poffsets_order.vals[i];
			j = i;
			end = getSymbolVal(STRING_ELT(symbols, j), ends_envir);
			SET_ELEMENT(ans, i, addInt(end, INTEGER(shift)[0]));
			SET_STRING_ELT(ans_names, i, duplicate(STRING_ELT(pids, poffsets.vals[j])));
		}
		SET_NAMES(ans, ans_names);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

