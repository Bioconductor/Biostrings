/****************************************************************************
 *           A fast implementation of the Aho-Corasick algorithm            *
 *                   for constant width DNA dictionaries                    *
 *                                                                          *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#include <stdlib.h>
#include <limits.h>

static int debug = 0;


#define MAX_CHILDREN_PER_ACNODE 4
typedef struct acnode {
	int parent_id;
	int depth;
	int child_id[MAX_CHILDREN_PER_ACNODE];
	int flink;
	int P_id;
} ACNode;



/****************************************************************************
 *                                                                          *
 *                             A. PREPROCESSING                             *
 *                                                                          *
 ****************************************************************************/

/****************************************************************************
 * Building the Aho-Corasick 4-ary tree
 * ------------------------------------
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

#define INTS_PER_ACNODE (sizeof(ACNode) / sizeof(int))
#define MAX_ACNODEBUF_LENGTH (INT_MAX / INTS_PER_ACNODE)

static ACNode *actree_nodes_buf = NULL;
static int actree_nodes_buf_count;

SEXP debug_match_pdict_ACtree()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pdict_ACtree.c'\n",
		debug ? "on" : "off");
	if (debug) {
		Rprintf("[DEBUG] debug_match_pdict_ACtree(): "
			"INTS_PER_ACNODE=%d\n", INTS_PER_ACNODE);
		Rprintf("[DEBUG] debug_match_pdict_ACtree(): "
			"MAX_ACNODEBUF_LENGTH=%d\n", MAX_ACNODEBUF_LENGTH);
	}
#else
	Rprintf("Debug mode not available in 'match_pdict_ACtree.c'\n");
#endif
	return R_NilValue;
}

SEXP free_actree_nodes_buf()
{
	if (actree_nodes_buf != NULL) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] free_actree_nodes_buf(): "
				"freeing actree_nodes_buf ... ");
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
 *   maxnodes = sum(from i=0; to i=W; of min(A^i, L))
 */
static void alloc_actree_nodes_buf(int length, int width)
{
	int maxnodes, pow, depth;
	size_t bufsize;

	if (actree_nodes_buf != NULL) {
		// We use the on.exit() mechanism to call
		// free_actree_nodes_buf() to free the buffer so if this
		// mechanism is reliable we should never come here.
		// Anyway just in case...
		warning("actree_nodes_buf was not previously freed, "
			"this is anormal, please report");
		free_actree_nodes_buf();
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
		Rprintf("[DEBUG] alloc_actree_nodes_buf(): "
			"length=%d width=%d maxnodes=%d\n",
			length, width, maxnodes);
	}
#endif
	if (maxnodes >= MAX_ACNODEBUF_LENGTH)
		error("the length or the width of the Trusted Band is too big");
	bufsize = sizeof(ACNode) * maxnodes;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] alloc_actree_nodes_buf(): "
			"allocating actree_nodes_buf (bufsize=%lu) ... ",
			bufsize);
	}
#endif
	actree_nodes_buf = (ACNode *) malloc(bufsize);
	if (actree_nodes_buf == NULL)
		error("alloc_actree_nodes_buf(): "
		      "failed to alloc actree_nodes_buf");
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

static int try_moving_to_acnode_child(int node_id, const int *base_codes,
		int childslot, char c)
{
	int child_id;

	if (base_codes[childslot] != c)
		return -1;
	child_id = actree_nodes_buf[node_id].child_id[childslot];
	if (child_id != -1)
		return child_id;
	child_id = append_acnode(node_id);
	return actree_nodes_buf[node_id].child_id[childslot] = child_id;
}

static void pp_pattern(const RoSeq *pattern, const int *base_codes, int poffset)
{
	int n, node_id, child_id, childslot;
	char c;
	ACNode *node;

	for (n = 0, node_id = 0; n < pattern->nelt; n++, node_id = child_id) {
		c = pattern->elts[n];
		for (childslot = 0; childslot < MAX_CHILDREN_PER_ACNODE; childslot++) {
			child_id = try_moving_to_acnode_child(node_id,
					base_codes, childslot, c);
			if (child_id != -1)
				break;
		}
		if (child_id == -1)
			error("non-base DNA letter found in Trusted Band "
			      "for pattern %d", poffset + 1);
	}
	node = actree_nodes_buf + node_id;
	if (node->P_id == -1)
		node->P_id = poffset + 1;
	else
		_report_dup(poffset, node->P_id);
	return;
}


/****************************************************************************
 * Turning our local data structures into an R list (SEXP)
 * -------------------------------------------------------
 */

/*
 * Turn the Aho-Corasick 4-ary tree stored in actree_nodes_buf into an R
 * external pointer (EXTPTRSXP).
 */
static SEXP new_ExternalPtr_from_actree_nodes_buf()
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
 * ACtree_asLIST() returns an R list with the following elements:
 *   - actree_nodes_xp: "externalptr" object pointing to the Aho-Corasick 4-ary
 *         tree built from the Trusted Band;
 *   - dup2unq: an integer vector containing the mapping between duplicated and
 *         primary reads.
 */
static SEXP ACtree_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("actree_nodes_xp"));
	SET_STRING_ELT(ans_names, 1, mkChar("dup2unq"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "actree_nodes_xp" element */
	PROTECT(ans_elt = new_ExternalPtr_from_actree_nodes_buf());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "dup2unq" element */
	PROTECT(ans_elt = _dup2unq_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry point for preprocessing
 * -----------------------------------
 *
 * Arguments:
 *   tb:         the Trusted Band extracted from the original dictionary as a
 *               DNAStringSet object;
 *   dup2unq0:   NULL or an integer vector of the same length as 'tb'
 *               containing the mapping from duplicated to unique patterns in
 *               the original dictionary;
 *   base_codes: the internal codes for A, C, G and T.
 *
 * See ACtree_asLIST() for a description of the returned SEXP.
 */

SEXP build_ACtree(SEXP tb, SEXP dup2unq0, SEXP base_codes)
{
	int tb_length, tb_width, poffset;
	CachedXStringSet cached_tb;
	RoSeq pattern;

	if (LENGTH(base_codes) != MAX_CHILDREN_PER_ACNODE)
		error("Biostrings internal error in build_ACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_ACNODE");
	tb_length = _get_XStringSet_length(tb);
	_init_dup2unq_buf(tb_length);
	tb_width = -1;
	cached_tb = _new_CachedXStringSet(tb);
	for (poffset = 0; poffset < tb_length; poffset++) {
		/* Skip duplicated patterns */
		if (dup2unq0 != R_NilValue
		 && INTEGER(dup2unq0)[poffset] != NA_INTEGER)
			continue;
		pattern = _get_CachedXStringSet_elt_asRoSeq(&cached_tb,
				poffset);
		if (pattern.nelt == 0)
			error("empty trusted region for pattern %d",
			      poffset + 1);
		if (tb_width == -1) {
			tb_width = pattern.nelt;
			alloc_actree_nodes_buf(tb_length, tb_width);
			append_acnode(0);
		} else if (pattern.nelt != tb_width) {
			error("all the trusted regions must have "
			      "the same length");
		}
		pp_pattern(&pattern, INTEGER(base_codes), poffset);
	}
	return ACtree_asLIST();
}



/****************************************************************************
 *                                                                          *
 *                             B. MATCH FINDING                             *
 *                                                                          *
 ****************************************************************************/

static int slotno_chrtrtable[CHRTRTABLE_LENGTH];

static int get_child_node_id(const ACNode *node, char c)
{
	int slotno;

	slotno = slotno_chrtrtable[(unsigned char) c];
	if (slotno == -1)
		return -1;
	return node->child_id[slotno];
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
 * walk_subject() function.
 */
static void set_shortcut(ACNode *node, char c, int next_node_id)
{
	int slotno, *slot;

	slotno = slotno_chrtrtable[(unsigned char) c];
	if (slotno == -1)
		return;
	slot = node->child_id + slotno;
	if (*slot == -1)
		*slot = next_node_id;
	return;
}

/*
 * We use indirect recursion for walking the Aho-Corasick tree.
 * This indirect recursion involves the 2 core functions path_to_node_id() and
 * get_next_node_id(): the latter calls the former which in turn calls the
 * latter.
 */
static int get_next_node_id(ACNode *node0, const int *base_codes,
		int node_id, const char *Stail, char c);

static int path_to_node_id(ACNode *node0, const int *base_codes,
		const char *path, int path_len)
{
	int node_id, n;
	ACNode *node;

	node_id = 0;
	for (n = 0; n < path_len; n++, path++) {
		node = node0 + node_id;
		node_id = get_next_node_id(node0, base_codes,
				node_id, path, *path);
	}
	return node_id;
}

/*
 * An important trick here is that the chars located _before_ 'Stail' will
 * always describe the path that goes from the root node to 'node_id'.
 * More precisely, if d is the depth of 'node_id', then this path is made
 * of the path[i] chars where path is 'Stail' - d and 0 <= i < d.
 */
static int get_next_node_id(ACNode *node0, const int *base_codes,
		int node_id, const char *Stail, char c)
{
	ACNode *node, *next_node;
	int next_node_id, child_node_id, subpath_len;
	const char *subpath;
#ifdef DEBUG_BIOSTRINGS
	static int rec_level = 0;
	char format[20], pathbuf[2000];
#endif

	node = node0 + node_id;
	next_node_id = node_id;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] ENTERING get_next_node_id():");
		sprintf(format, "%%%ds", 1 + 2*rec_level);
		Rprintf(format, " ");
		snprintf(pathbuf, node->depth + 1, "%s", Stail - node->depth);
		Rprintf("node_id=%d path=%s c=%c\n",
			node_id, pathbuf, c);
	}
#endif
	while (1) {
		next_node = node0 + next_node_id;
		child_node_id = get_child_node_id(next_node, c);
		if (child_node_id != -1) {
			next_node_id = child_node_id;
			break;
		}
		if (next_node_id == 0) {
			//next_node->flink = 0;
			break;
		}
		if (next_node->flink == -1) {
			rec_level++;
			subpath_len = next_node->depth - 1;
			subpath = Stail - subpath_len;
			next_node->flink = path_to_node_id(node0, base_codes,
						subpath, subpath_len);
			rec_level--;
		}
		next_node_id = next_node->flink;
	}
	set_shortcut(node, c, next_node_id);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] LEAVING get_next_node_id(): ");
		Rprintf(format, " ");
		Rprintf("next_node_id=%d\n", next_node_id);
	}
#endif
	return next_node_id;
}

static int walk_subject(ACNode *node0, const int *base_codes,
		const char *S, int nS)
{
	int basenode_id, node_id, child_id, n, subwalk_nS;
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
			Rprintf("[DEBUG] walk_subject():");
			sprintf(format, "%%%ds", 1 + 2*rec_level);
			Rprintf(format, " ");
			snprintf(pathbuf, basenode->depth + 1, "%s", S - basenode->depth);
			Rprintf("On basenode_id=%d (basepath=%s), reading S[%d]=%c\n", basenode_id, pathbuf, n, *S);
		}
#endif
		node_id = basenode_id;
		node = basenode;
		while (1) {
			child_id = get_child_node_id(node, *S);
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
				subwalk_nS = node->depth - 1;
				node->flink = walk_subject(node0, base_codes, S - subwalk_nS, subwalk_nS);
				rec_level--;
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					Rprintf("[DEBUG] walk_subject():");
					Rprintf(format, " ");
					Rprintf("setting failure link %d -> %d\n", node_id, node->flink);
				}
#endif
			}
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] walk_subject():");
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
			Rprintf("[DEBUG] walk_subject():");
			Rprintf(format, " ");
			Rprintf("moving to basenode %d\n", basenode_id);
		}
#endif
		// Finding a match cannot happen during a nested call to
		// walk_subject() so there is no need to check that rec_level
		// is 0
		if (basenode->P_id != -1)
			_MIndex_report_match(basenode->P_id - 1, n + 1);
	}
	return basenode_id;
}

static void walk_nonfixed_subject(ACNode *node0, const int *base_codes,
		const RoSeq *S)
{
	IntAE cnode_ids; // buffer of current node ids
	int n, npointers, i, node_id, next_node_id, is_first, j, base, P_id;
	const char *S_tail;
	char c;

	cnode_ids = new_IntAE(256, 0, 0);
	IntAE_insert_at(&cnode_ids, 0, 0);
	for (n = 1, S_tail = S->elts; n <= S->nelt; n++, S_tail++) {
		c = *S_tail;
		npointers = cnode_ids.nelt;
		// move and split pointers
		for (i = 0; i < npointers; i++) {
			node_id = cnode_ids.elts[i];
			is_first = 1;
			for (j = 0, base = 1; j < 4; j++, base *= 2) {
				if ((((unsigned char) c) & base) != 0) {
					next_node_id = get_next_node_id(node0,
							base_codes,
							node_id, S_tail, base);
					if (is_first) {
						cnode_ids.elts[i] = next_node_id;
						is_first = 0;
					} else {
						IntAE_insert_at(&cnode_ids,
							cnode_ids.nelt, next_node_id);
					}
				}
			}
		}
		// merge pointers and report matches
		for (i = 0; i < cnode_ids.nelt; i++) {
			node_id = cnode_ids.elts[i];
			// FIXME: This merging algo is dumb and inefficient!
			// There must be a way to do something better.
			for (j = i + 1; j < cnode_ids.nelt; j++) {
				if (cnode_ids.elts[j] == node_id)
					IntAE_delete_at(&cnode_ids, j--);
			}
			P_id = node0[node_id].P_id;
			if (P_id != -1)
				_MIndex_report_match(P_id - 1, n);
		}
		// error if too many remaining pointers
		if (cnode_ids.nelt > 4096)
			error("too many IUPAC ambiguity letters in 'subject'");
	}
	return;
}

void _match_ACtree(SEXP pdict_pptb, const RoSeq *S, int fixedS)
{
	ACNode *node0;
	SEXP base_codes;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_ACtree()\n");
#endif
	node0 = (ACNode *) INTEGER(R_ExternalPtrTag(VECTOR_ELT(pdict_pptb, 0)));
	base_codes = VECTOR_ELT(pdict_pptb, 1);
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_ACNODE)
		error("Biostrings internal error in _match_ACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_ACNODE");
	_init_chrtrtable(INTEGER(base_codes), MAX_CHILDREN_PER_ACNODE,
			 slotno_chrtrtable);
	if (fixedS)
		walk_subject(node0, INTEGER(base_codes), S->elts, S->nelt);
	else
		walk_nonfixed_subject(node0, INTEGER(base_codes), S);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING _match_ACtree()\n");
#endif
	return;
}

