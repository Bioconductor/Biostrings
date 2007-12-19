/****************************************************************************
 *              Aho-Corasick for uniform-length DNA dictionary              *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a uniform-length dictionary is a non-empty set of non-empty        *
 * strings of the same length.                                              *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Srealloc() */

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
	int width;		// number of chars per pattern
	int length;		// number of patterns
	const char **patterns;	// array of pointers to the read-only patterns
} ULdict;

static ULdict input_uldict;

static void alloc_input_uldict(int length)
{
	input_uldict.patterns = Salloc((long) length, const char *);
	input_uldict.length = length;
	input_uldict.width = -1;
	return;
}

static void add_pattern_to_input_uldict(int poffset, const char *pattern, int pattern_length)
{
	if (input_uldict.width == -1) {
		if (pattern_length == 0)
			error("dictionary contains empty patterns");
		input_uldict.width = pattern_length;
	} else {
		if (pattern_length != input_uldict.width)
			error("all patterns in dictionary must have the same length");
	}
	input_uldict.patterns[poffset] = pattern;
	return;
}


/****************************************************************************
 * Buffer of duplicates.
 */

IBuf dups_buf;

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

typedef struct ac_node {
	int parent_id;
	int child_id[ALPHABET_LENGTH];
	int flink;
	int P_id;
} ACNode;

#define INTS_PER_ACNODE (sizeof(ACNode) / sizeof(int))
#define MAX_ACNODEBUF_LENGTH (INT_MAX / INTS_PER_ACNODE)

static ACNode *ACnodebuf;
static int ACnodebuf_count;
static int AC_base_codes[ALPHABET_LENGTH];

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

/*
 * We want to avoid using reallocation for ACnodebuf (the buffer where we are
 * going to build the Aho-Corasick tree) so we need to know what's the maximum
 * number of nodes that is needed for storing a dictionary of a given width and
 * length. First some notations:
 *   W: dictionary width i.e. number of chars per pattern
 *   L: dictionary length i.e. number of patterns
 *   A: length of the alphabet i.e. max number of children per node
 *   maxnodes: maximum number of nodes needed
 * Now here is how to get maxnodes:
 *   maxnodes = sum(from: i=0, to: i=W, of: min(A^i, L))
 */
static void alloc_ACnodebuf(int width, int length)
{
	int maxnodes, pow, depth;

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
		Rprintf("[DEBUG] alloc_ACnodebuf(): width=%d length=%d maxnodes=%d\n",
			width, length, maxnodes);
	}
#endif
	if (maxnodes >= MAX_ACNODEBUF_LENGTH)
		error("dictionary is too big");
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] alloc_ACnodebuf(): allocating ACnodebuf ... ");
	}
#endif
	ACnodebuf = (ACNode *) malloc(sizeof(ACNode) * maxnodes);
	if (ACnodebuf == NULL)
		error("alloc_ACnodebuf(): failed to alloc ACnodebuf");
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("OK\n");
	}
#endif
	ACnodebuf_count = 0;
	return;
}

static void init_AC_base_codes()
{
	int childslot;

	for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++)
		AC_base_codes[childslot] = -1;
	return;
}

SEXP ULdna_free_ACnodebuf()
{
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] ULdna_free_ACnodebuf(): freeing ACnodebuf ... ");
	}
#endif
	free(ACnodebuf);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("OK\n");
	}
#endif
	return R_NilValue;
}

static void ACnodebuf_init_node(ACNode *node, int parent_id)
{
	int childslot;

	node->parent_id = parent_id;
	for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++)
		node->child_id[childslot] = -1;
	node->flink = -1;
	node->P_id = -1;
	return;
}

static int ACnodebuf_append_node(int parent_id)
{
	ACNode *node;

	node = ACnodebuf + ACnodebuf_count;
	ACnodebuf_init_node(node, parent_id);
	return ACnodebuf_count++;
}

static int ACnodebuf_tryMovingToChild(int node_id, int childslot, char c)
{
	int child_id;

	if (AC_base_codes[childslot] == -1) {
		AC_base_codes[childslot] = c;
	} else {
		if (c != AC_base_codes[childslot])
			return -1;
	}
	child_id = ACnodebuf[node_id].child_id[childslot];
	if (child_id != -1)
		return child_id;
	child_id = ACnodebuf_append_node(node_id);
	return ACnodebuf[node_id].child_id[childslot] = child_id;
}

static void ACnodebuf_addPattern(int poffset)
{
	const char *pattern;
	int n, node_id, child_id, childslot;
	char c;
	ACNode *node;

	pattern = input_uldict.patterns[poffset];
	for (n = 0, node_id = 0; n < input_uldict.width; n++, node_id = child_id) {
		c = pattern[n];
		for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++) {
			child_id = ACnodebuf_tryMovingToChild(node_id, childslot, c);
			if (child_id != -1)
				break;
		}
		if (child_id == -1)
			error("dictionary contains more than %d distinct letters", ALPHABET_LENGTH);
	}
	node = ACnodebuf + node_id;
	if (node->P_id == -1)
		node->P_id = poffset + 1;
	else
		report_dup(poffset, node->P_id);
	return;
}

static void build_ACtree()
{
	int poffset;

	init_dups_buf(input_uldict.length);
	init_AC_base_codes();
	alloc_ACnodebuf(input_uldict.width, input_uldict.length);
	ACnodebuf_append_node(0);
	for (poffset = 0; poffset < input_uldict.length; poffset++)
		ACnodebuf_addPattern(poffset);
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

IBBuf ends_bbuf;

static void init_ends_bbuf(int length)
{
	_IBBuf_init(&ends_bbuf, length, length);
	return;
}

static void report_match(int poffset, int end)
{
	IBuf *ends_buf;

	ends_buf = ends_bbuf.ibufs + poffset;
	_IBuf_insert_at(ends_buf, ends_buf->count, end);
	return;
}

static int get_child_id(ACNode *node, char c)
{
	int childslot;

	for (childslot = 0; childslot < ALPHABET_LENGTH; childslot++)
		if (c == AC_base_codes[childslot])
			return node->child_id[childslot];
	return -1;
}

static int get_node_depth(ACNode *node0, ACNode *node)
{
	int depth;

	for (depth = 0; node != node0; depth++)
		node = node0 + node->parent_id;
	return depth;
}

static int follow_string(ACNode *node0, const char *S, int nS)
{
	int node_id, child_id, depth, new_depth, n;
	const char *path;
	ACNode *node;
	static int rec_level = 0;
#ifdef DEBUG_BIOSTRINGS
	char format[20], pathbuf[2000];
#endif

	node = node0;
	node_id = 0;
	path = S;
	depth = 0;
	for (n = 0; n < nS; n++, S++) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] follow_string():");
			sprintf(format, "%%%ds", 1 + 2*rec_level);
			Rprintf(format, " ");
			snprintf(pathbuf, depth + 1, "%s", path);
			Rprintf("On node_id=%d (path=%s), reading S[%d]=%c\n", node_id, pathbuf, n, *S);
		}
#endif
		while (1) {
			child_id = get_child_id(node, *S);
			if (child_id != -1) {
				node_id = child_id;
				depth++;
				break;
			}
			if (node_id == 0) {
				path++;
				break;
			}
			if (node->flink == -1) {
				rec_level++;
				node->flink = follow_string(node0, path + 1,  depth - 1);
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
			new_depth = get_node_depth(node0, node);
			path += depth - new_depth;
			depth = new_depth;
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] follow_string():");
			Rprintf(format, " ");
			Rprintf("moving to node %d\n", node_id);
		}
#endif
		node = node0 + node_id;
		// Finding a match cannot happen during a nested call to follow_string()
		if (node->P_id != -1)
			report_match(node->P_id - 1, n + 1);
	}
	return node_id;
}

static void ULdna_exact_search(int uldna_len, ACNode *ACtree, const char *S, int nS)
{
	init_ends_bbuf(uldna_len);
	follow_string(ACtree, S, nS);
	return;
}


/****************************************************************************
 * Turning our local data structures into R objects
 * ================================================
 */

/*
 * Turn the Aho-Corasick 4-ary tree stored in ACnodebuf into an eXternal Pointer
 */
static SEXP ACnodebuf_asXP()
{
	SEXP ans, tag;
	int tag_length;

	PROTECT(ans = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	tag_length = ACnodebuf_count * INTS_PER_ACNODE;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), ACnodebuf, tag_length * sizeof(int));
	R_SetExternalPtrTag(ans, tag);
	UNPROTECT(2);
	return ans;
}

/*
 * Turn the AC_base_codes array into an R vector of ALPHABET_LENGTH integers
 */
static SEXP base_codes_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ALPHABET_LENGTH));
	memcpy(INTEGER(ans), AC_base_codes, ALPHABET_LENGTH * sizeof(int));
	UNPROTECT(1);
	return ans;
}

/*
 * Return an R list with the following elements:
 *   - ACtree_xp: "externalptr" object pointing to the Aho-Corasick 4-ary
 *         tree built from 'dict'.
 *   - AC_base_codes: integer vector containing the ALPHABET_LENGTH character
 *         codes (ASCII) attached to the ALPHABET_LENGTH child slots of any
 *         node in the tree pointed by ACtree_xp.
 *   - dups: an integer vector of the same length as 'dict' containing the
 *         duplicate info.
 */
static SEXP ACtree_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(3));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("ACtree_xp"));
	SET_STRING_ELT(ans_names, 1, mkChar("AC_base_codes"));
	SET_STRING_ELT(ans_names, 2, mkChar("dups"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "ACtree_xp" element */
	PROTECT(ans_elt = ACnodebuf_asXP());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "AC_base_codes" element */
	PROTECT(ans_elt = base_codes_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	/* set the "dups" element */
	PROTECT(ans_elt = _IBuf_asINTEGER(&dups_buf));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry points for building the Aho-Corasick 4-ary tree from the input
 * dictionary. The input dictionary must be:
 *   - of length >= 1
 *   - uniform-length (i.e. all words have the same length)
 ****************************************************************************/

/*
 * .Call entry point: "ULdna_init_with_StrVect"
 *
 * Argument:
 *   'dict': a string vector (aka character vector) containing the
 *           uniform-length dictionary
 *
 * See ACtree_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_init_with_StrVect(SEXP dict)
{
	int poffset;
	SEXP dict_elt;

	alloc_input_uldict(LENGTH(dict));
	for (poffset = 0; poffset < LENGTH(dict); poffset++) {
		dict_elt = STRING_ELT(dict, poffset);
		add_pattern_to_input_uldict(poffset, CHAR(dict_elt), LENGTH(dict_elt));
	}
	build_ACtree();
	return ACtree_asLIST();
}

/*
 * .Call entry point: "ULdna_init_with_BStringList"
 *
 * Argument:
 *   'dict': a list of (pattern@data@xp, pattern@offset, pattern@length)
 *           triplets containing the uniform-length dictionary
 *
 * See ACtree_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_init_with_BStringList(SEXP dict)
{
	error("Not ready yet!\n");
	return ACtree_asLIST();
}

/*
 * .Call entry point: "ULdna_init_with_views"
 *
 * Arguments:
 *   'dict_subj_xp': subject(dict)@data@xp
 *   'dict_subj_offset': subject(dict)@offset
 *   'dict_start': start(dict)
 *   'dict_end': end(dict)
 * Note that 'dict_start' and 'dict_end' describe a set of views on subject(dict)
 * with no "out of limits" views.
 *
 * See ACtree_asLIST() for a description of the returned SEXP.
 */
SEXP ULdna_init_with_views(SEXP dict_subj_xp, SEXP dict_subj_offset,
		SEXP dict_start, SEXP dict_end)
{
	SEXP tag;
	int subj_offset, dict_length, poffset, view_offset, view_end;
	const Rbyte *subj;

	PROTECT(tag = R_ExternalPtrTag(dict_subj_xp));
	subj_offset = INTEGER(dict_subj_offset)[0];
	subj = RAW(R_ExternalPtrTag(dict_subj_xp)) + subj_offset;
	dict_length = LENGTH(dict_start); // must be the same as LENGTH(dict_end)
	alloc_input_uldict(dict_length);
	for (poffset = 0; poffset < dict_length; poffset++) {
		view_offset = INTEGER(dict_start)[poffset] - 1;
		view_end = INTEGER(dict_end)[poffset];
		add_pattern_to_input_uldict(poffset, (char *) (subj + view_offset), view_end - view_offset);
	}
	build_ACtree();
	UNPROTECT(1);
	return ACtree_asLIST();
}


/****************************************************************************
 * .Call entry point: "match_ULdna_exact"
 *
 * Arguments:
 *   'uldna_length': uldna_pdict@length (this is the number of words in the
 *                   uniform-length dictionary)
 *   'uldna_dups': uldna_pdict@dups
 *   'ACtree_xp': uldna_pdict@ACtree@xp
 *   'AC_base_codes': uldna_pdict@AC_base_codes
 *   's_xp': subject@data@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *
 * For now return a list of length 'uldna_length' where each element is a
 * vector of integers containing matches (start position). More precisely,
 * for the i-th pattern in the dictionary, all the matches found in the
 * subject are stored in the i-th element of this list.
 *
 ****************************************************************************/

SEXP match_ULdna_exact(SEXP uldna_length, SEXP uldna_dups,
		SEXP ACtree_xp, SEXP AC_base_codes,
		SEXP s_xp, SEXP s_offset, SEXP s_length)
{
	int uldna_len, ACtree_length, subj_offset, subj_length;
	ACNode *ACtree;
	const Rbyte *subj;
	SEXP tag;

	uldna_len = INTEGER(uldna_length)[0];
	dups_buf = _INTEGER_asIBuf(uldna_dups);
	tag = R_ExternalPtrTag(ACtree_xp);
	ACtree = (ACNode *) INTEGER(tag);
	ACtree_length = LENGTH(tag) / INTS_PER_ACNODE;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;

	ULdna_exact_search(uldna_len, ACtree, (char *) subj, subj_length);

	return _IBBuf_asLIST(&ends_bbuf);
}

