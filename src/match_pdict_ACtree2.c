/****************************************************************************
 *     A fast and compact implementation of the Aho-Corasick algorithm      *
 *                   for constant width DNA dictionaries                    *
 *                                                                          *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#include <S.h> /* for Salloc() */
#include <stdlib.h>
#include <limits.h>

static int debug = 0;

#define MAX_CHILDREN_PER_NODE 4



/****************************************************************************
 *    A. GENERAL UTILITIES FOR COMPUTING THE NB OF NODES AND EXTENSIONS     *
 ****************************************************************************/

/*
 * nleaves = number of leaves in the tree (it's also the number of unique
 *           patterns so 'nleaves' <= 'tb_length');
 * depth = depth of the tree (it's also the length of each pattern so 'depth'
 *           = 'tb_width').
 */
static int get_max_needed_nnodes(int nleaves, int depth)
{
	int max_nn, nnodes_inc, d;

	max_nn = 0;
	nnodes_inc = 1;
	for (d = 0; d <= depth; d++) {
		if (nnodes_inc >= nleaves)
			return max_nn + (depth + 1 - d) * nleaves;
		max_nn += nnodes_inc;
		nnodes_inc *= MAX_CHILDREN_PER_NODE;
	}
	return max_nn;
}

static int get_min_needed_nnodes(int nleaves, int depth)
{
	int min_nn, nnodes_inc, d;
	div_t q;

	min_nn = 0;
	nnodes_inc = nleaves;
	for (d = depth; d >= 0; d--) {
		if (nnodes_inc == 1)
			return min_nn + d + 1;
		min_nn += nnodes_inc;
		q = div(nnodes_inc, MAX_CHILDREN_PER_NODE);
		nnodes_inc = q.quot;
		if (q.rem != 0)
			nnodes_inc++;
	}
	return min_nn;
}

/*
 * OptMaxNN: Optimistic Max Needed nb of Nodes.
 * get_OptMaxNN() is a decreasing function of 'max_needed_nnodes' bounded by
 * 'max_nn' and 'min_nn2' (max_nn >= get_OptMaxNN() > min_nn2).
 */
static int get_OptMaxNN(int max_needed_nnodes, int min_needed_nnodes)
{
	double max_nn, min_nn2, x;

	max_nn = max_needed_nnodes;
	min_nn2 = min_needed_nnodes + 0.70 * (max_nn - min_needed_nnodes);
	x = max_needed_nnodes / 300000000.00;
	return (int) min_nn2 + (max_nn - min_nn2) / (1.00 + x);
}

/*
 * nleaves = number of leaves in the tree (it's also the number of unique
 *           patterns so 'nleaves' <= 'tb_length');
 * depth = depth of the tree (it's also the length of each pattern so 'depth'
 *           = 'tb_width').
 */
static int get_max_needed_nextensions_at_pp_time(int nleaves, int depth)
{
	int nextensions, d;

	nextensions = 1;
	for (d = 0; d < depth; d++) {
		if (nextensions > nleaves)
			break;
		nextensions *= 2;
	}
	return nextensions - 1;
}

/*
 * MEER (Maximum Expected Extension Ratio) is the maximum nextensions/nnodes
 * ratio that is expected during the lifetime of an ACtree2 object.
 * The nextensions/nnodes ratio can only increase during the lifetime of an
 * ACtree2 object (more precisely, it can only increase when the object is
 * used to walk on new subjects). Since the cost of reallocating is high, we
 * want to preallocate the buffer of extensions (stored in the 'extensions'
 * slot of the object) once for all with the hope that its length will be
 * enough for the entire lifetime of the object. The current approach is to
 * set its length to 'nnodes * MEER'. Some simulations with real data tend to
 * indicate that most of the times this will be enough.
 * Note that when 'extensions_buflength == 0.4 * nnodes', the 'extensions'
 * slot of the ACtree2 object has the same size as its 'nodes' slot. Therefore
 * the total size for these 2 slots ('nodes' + 'extensions') is half the size
 * of the 'nodes' slot of the corresponding old ACtree object.
 * We define MEER as an increasing function of the number of nodes ('nnodes').
 */
static double nnodes2MEER(int nnodes)
{
	double x;

	x = nnodes / 1000000.00;
	if (x <= 1)
		return 1.00;
	if (x <= 5)
		return 0.75;
	if (x <= 25)
		return 0.60;
	if (x <= 50)
		return 0.50;
	if (x <= 150)
		return 0.50 - 0.10 * (x / 50 - 1);
	if (x <= 300)
		return 0.30 - 0.05  * (x / 150 - 1);
	return 0.25;
}



/****************************************************************************
 *                                B. DEFINES                                *
 ****************************************************************************/

typedef struct acnode {
	int attribs;
	int nid_or_eid;
} ACnode;

typedef struct acnode_extension {
	int link_nid[MAX_CHILDREN_PER_NODE];
	int flink_nid;
} ACnodeExtension;

#define INTS_PER_NODE (sizeof(ACnode) / sizeof(int))
#define MAX_NNODES (INT_MAX / INTS_PER_NODE)

#define INTS_PER_EXTENSION (sizeof(ACnodeExtension) / sizeof(int))
#define MAX_NEXTENSIONS (INT_MAX / INTS_PER_EXTENSION)

#define LINKTAG_BITSHIFT 28
#define LINKTAG_BITMASK (3 << LINKTAG_BITSHIFT)
#define MAX_DEPTH ((1 << LINKTAG_BITSHIFT) - 1)
#define ISLEAF_BIT (1 << 30)
#define ISEXTENDED_BIT (ISLEAF_BIT << 1) /* strongest bit for 32-bit integers */
#define MAX_P_ID (ISLEAF_BIT - 1)  /* P_id values are encoded on 30 bits */

SEXP debug_match_pdict_ACtree2()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
	if (debug) {
		Rprintf("[DEBUG] debug_match_pdict_ACtree2():\n"
			"  INTS_PER_NODE=%d MAX_NNODES=%d\n"
			"  INTS_PER_EXTENSION=%d MAX_NEXTENSIONS=%d\n"
			"  LINKTAG_BITSHIFT=%d LINKTAG_BITMASK=%d\n"
			"  MAX_DEPTH=%d\n"
			"  ISLEAF_BIT=%d ISEXTENDED_BIT=%d\n"
			"  MAX_P_ID=%d\n",
			INTS_PER_NODE, MAX_NNODES,
			INTS_PER_EXTENSION, MAX_NEXTENSIONS,
			LINKTAG_BITSHIFT, LINKTAG_BITMASK,
			MAX_DEPTH,
			ISLEAF_BIT, ISEXTENDED_BIT,
			MAX_P_ID);
	}
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

#define ISEXTENDED(node) ((node)->attribs & ISEXTENDED_BIT)
#define ISLEAF(node) ((node)->attribs & ISLEAF_BIT)

typedef struct actree {
	int depth;  /* this is the depth of all leaf nodes */
	ACnode *nodes;
	int nodes_buflength, nnodes;
	ACnodeExtension *extensions;
	int extensions_buflength, nextensions, *nextensions_ptr;
	ByteTrTable char2linktag;
} ACtree;

#define TREE_DEPTH(tree) ((tree)->depth)
#define TREE_NNODES(tree) ((tree)->nnodes)
#define GET_NODE(tree, nid) ((tree)->nodes + (nid))
#define TREE_NEXTENSIONS(tree) (*((tree)->nextensions_ptr))
#define GET_EXTENSION(tree, eid) ((tree)->extensions + (eid))
#define CHAR2LINKTAG(tree, c) ((tree)->char2linktag[(unsigned char) (c)])

#define NODE_DEPTH(tree, node) \
	(ISLEAF(node) ? TREE_DEPTH(tree) : ((node)->attribs & MAX_DEPTH))
/* result of NODE_P_ID() is undefined on a non-leaf node */
#define NODE_P_ID(node) ((node)->attribs & MAX_P_ID)



/****************************************************************************
 *                          C. LOW-LEVEL UTILITIES                          *
 ****************************************************************************/

/* extends the 'nodes' slot */
static void extend_nodes_buffer(ACtree *tree)
{
	error("extend_nodes_buffer(): implement me");
	return;
}

/* extends the 'extensions' slot */
static void extend_extensions_buffer(ACtree *tree)
{
	error("extend_extensions_buffer(): implement me");
	return;
}

static void extend_ACnode(ACtree *tree, ACnode *node)
{
	ACnodeExtension *extension;
	int eid, i, linktag;

	if (TREE_NEXTENSIONS(tree) >= tree->extensions_buflength)
		extend_extensions_buffer(tree);
	eid = TREE_NEXTENSIONS(tree)++;
	extension = GET_EXTENSION(tree, eid);
	for (i = 0; i < MAX_CHILDREN_PER_NODE; i++)
		extension->link_nid[i] = -1;
	extension->flink_nid = -1;
	if (node->nid_or_eid != -1) {
		/* this is correct because 'node' cannot be a leaf node */
		linktag = node->attribs >> LINKTAG_BITSHIFT;
		extension->link_nid[linktag] = node->nid_or_eid;
	}
	node->nid_or_eid = eid;
	/* sets the "ISEXTENDED" bit to 1 */
	node->attribs |= ISEXTENDED_BIT;
	return;
}

static SEXP ACtree_nodes_asXInteger(ACtree *tree)
{
	int tag_length;
	SEXP ans, tag;

	tag_length = TREE_NNODES(tree) * INTS_PER_NODE;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), tree->nodes,
			TREE_NNODES(tree) * sizeof(ACnode));
	PROTECT(ans = new_XInteger_from_tag("XInteger", tag));
	UNPROTECT(2);
	return ans;
}

static SEXP ACtree_extensions_asXInteger(ACtree *tree)
{
	int extensions_buflength, tag_length;
	SEXP ans, tag;

	extensions_buflength = TREE_NNODES(tree)
			     * nnodes2MEER(TREE_NNODES(tree));
	/* very unlikely to happen */
	if (extensions_buflength < TREE_NEXTENSIONS(tree))
		extensions_buflength = TREE_NEXTENSIONS(tree);
	tag_length = extensions_buflength * INTS_PER_EXTENSION;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), tree->extensions,
			TREE_NEXTENSIONS(tree) * sizeof(ACnodeExtension));
	PROTECT(ans = new_XInteger_from_tag("XInteger", tag));
	UNPROTECT(2);
	return ans;
}

static SEXP ACtree_nextensions_asXInteger(ACtree *tree)
{
	SEXP ans, tag;

	PROTECT(tag = ScalarInteger(TREE_NEXTENSIONS(tree)));
	PROTECT(ans = new_XInteger_from_tag("XInteger", tag));
	UNPROTECT(2);
	return ans;
}



/****************************************************************************
 *                           D. Aho-Corasick tree API                       *
 ****************************************************************************/

static int new_ACnode(ACtree *tree, int depth)
{
	ACnode *node;
	int nid;

	if (depth >= TREE_DEPTH(tree))
		error("new_ACnode(): depth >= TREE_DEPTH(tree)");
	if (TREE_NNODES(tree) >= tree->nodes_buflength)
		extend_nodes_buffer(tree);
	nid = tree->nnodes++;
	node = GET_NODE(tree, nid);
	/* this sets the "ISEXTENDED" and "ISLEAF" bits to 0 */
	node->attribs = depth;
	node->nid_or_eid = -1;
	return nid;
}

static int new_leafACnode(ACtree *tree, int P_id)
{
	ACnode *node;
	int nid;

	if (TREE_NNODES(tree) >= tree->nodes_buflength)
		extend_nodes_buffer(tree);
	nid = tree->nnodes++;
	node = GET_NODE(tree, nid);
	/* this sets the "ISEXTENDED" bit to 0 and "ISLEAF" bit to 1 */
	node->attribs = ISLEAF_BIT | P_id;
	node->nid_or_eid = -1;
	return nid;
}

static int get_ACnode_link(ACtree *tree, ACnode *node, int linktag)
{
	ACnodeExtension *extension;

	if (node->nid_or_eid == -1)
		return -1;
	if (ISEXTENDED(node)) {
		extension = GET_EXTENSION(tree, node->nid_or_eid);
		return extension->link_nid[linktag];
	}
	/* the node has no extension and is not a leaf node */
	if (linktag == (node->attribs >> LINKTAG_BITSHIFT))
		return node->nid_or_eid;
	return -1;
}

/*
 * set_ACnode_flink() should always have been called before set_ACnode_link()
 * on a leaf node. So we can assume that, if set_ACnode_link() is called on a
 * leaf node, then this node is already extended.
 */
static void set_ACnode_link(ACtree *tree, ACnode *node, int linktag, int nid)
{
	ACnodeExtension *extension;

	if (node->nid_or_eid == -1) {
		/* cannot be a leaf node (see assumption above) and
		   no need to extend it */
		node->attribs |= linktag << LINKTAG_BITSHIFT;
		node->nid_or_eid = nid;
		return;
	}
	if (!ISEXTENDED(node)) {
		/* again, cannot be a leaf node (see assumption above) */
		extend_ACnode(tree, node);
	}
	extension = GET_EXTENSION(tree, node->nid_or_eid);
	extension->link_nid[linktag] = nid;
	return;
}

static int get_ACnode_flink(ACtree *tree, ACnode *node)
{
	ACnodeExtension *extension;

	if (!ISEXTENDED(node))
		return -1;
	extension = GET_EXTENSION(tree, node->nid_or_eid);
	return extension->flink_nid;
}

static void set_ACnode_flink(ACtree *tree, ACnode *node, int nid)
{
	ACnodeExtension *extension;

	if (!ISEXTENDED(node))
		extend_ACnode(tree, node);
	extension = GET_EXTENSION(tree, node->nid_or_eid);
	extension->flink_nid = nid;
	return;
}

static ACtree new_ACtree(int tb_length, int tb_width, SEXP base_codes)
{
	ACtree tree;
	int max_nn, min_nn, n1, n2;

	max_nn = get_max_needed_nnodes(tb_length, tb_width);
	min_nn = get_min_needed_nnodes(tb_length, tb_width);
	n1 = get_OptMaxNN(max_nn, min_nn);
	n2 = get_max_needed_nextensions_at_pp_time(tb_length, tb_width);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] new_ACtree():\n"
			"  tb_length=%d tb_width=%d\n"
			"  max_needed_nnodes=%d min_needed_nnodes=%d\n"
			"  nodes_buflength=%d (=OptMaxNN)\n"
			"  extensions_buflength=%d (=max_needed_nextensions)\n",
			tb_length, tb_width,
			max_nn, min_nn,
			n1,
			n2);
	}
#endif
	if (tb_length > MAX_P_ID)
		error("new_ACtree(): tb_length > MAX_P_ID");
	if (tb_width > MAX_DEPTH)
		error("new_ACtree(): tb_width > MAX_DEPTH");
	if (n1 >= MAX_NNODES || n2 >= MAX_NEXTENSIONS)
		error("Trusted Band is too big (please reduce its "
		      "width or its length)");
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_NODE)
		error("Biostrings internal error in new_ACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_NODE");

	tree.depth = tb_width;
	tree.nodes = Salloc((long) n1, ACnode);
	tree.nodes_buflength = n1;
	tree.nnodes = 0;
	tree.extensions = Salloc((long) n2, ACnodeExtension);
	tree.extensions_buflength = n2;
	tree.nextensions = 0;
	tree.nextensions_ptr = &(tree.nextensions);
	_init_byte2offset_with_INTEGER(tree.char2linktag, base_codes, 1);
	new_ACnode(&tree, 0);  /* create the root node */
	return tree;
}

static SEXP ACtree_asLIST(ACtree *tree)
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(3));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("nodes"));
	SET_STRING_ELT(ans_names, 1, mkChar("extensions"));
	SET_STRING_ELT(ans_names, 2, mkChar("nextensions"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "nodes" element */
	PROTECT(ans_elt = ACtree_nodes_asXInteger(tree));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "extensions" element */
	PROTECT(ans_elt = ACtree_extensions_asXInteger(tree));
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	/* set the "nextensions" element */
	PROTECT(ans_elt = ACtree_nextensions_asXInteger(tree));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

static ACtree pptb_asACtree(SEXP pptb)
{
	ACtree tree;
	SEXP tag, base_codes;

	tree.depth = _get_PreprocessedTB_width(pptb);
	tag = _get_ACtree2_nodes_tag(pptb);
	tree.nodes = (ACnode *) INTEGER(tag);
	tree.nodes_buflength = LENGTH(tag) / INTS_PER_NODE;
	tree.nnodes = tree.nodes_buflength;
	tag = _get_ACtree2_extensions_tag(pptb);
	tree.extensions = (ACnodeExtension *) INTEGER(tag);
	tree.extensions_buflength = LENGTH(tag) / INTS_PER_EXTENSION;
	// tree.nextensions is not set because it's not used
	tree.nextensions_ptr = _get_ACtree2_nextensions(pptb);
	base_codes = _get_ACtree2_base_codes(pptb);
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_NODE)
		error("Biostrings internal error in pptb_asACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_NODE");
	_init_byte2offset_with_INTEGER(tree.char2linktag, base_codes, 1);
	return tree;
}

static void print_ACnode(ACtree *tree, ACnode *node)
{
	error("print_ACnode(): implement me");
	return;
}

static int get_ACnode_nlinks(ACtree *tree, ACnode *node)
{
	int nlinks, linktag;

	nlinks = get_ACnode_flink(tree, node) != -1 ? 1 : 0;
	for (linktag = 0; linktag < MAX_CHILDREN_PER_NODE; linktag++)
		if (get_ACnode_link(tree, node, linktag) != -1)
			nlinks++;
	return nlinks;
}



/****************************************************************************
 *                       E. STATS AND DEBUG UTILITIES                       *
 ****************************************************************************/

SEXP ACtree2_print_nodes(SEXP pptb)
{
	ACtree tree;
	ACnode *node;
	int nnodes, nid;

	tree = pptb_asACtree(pptb);
	nnodes = TREE_NNODES(&tree);
	for (nid = 0; nid < nnodes; nid++) {
		node = GET_NODE(&tree, nid);
		print_ACnode(&tree, node);
	}
	return R_NilValue;
}

SEXP ACtree2_summary(SEXP pptb)
{
	ACtree tree;
	ACnode *node;
	int nnodes, nlinks_table[MAX_CHILDREN_PER_NODE+2], nleaves,
	    max_nn, min_nn, n1, i, nid, nlinks;

	tree = pptb_asACtree(pptb);
	nnodes = TREE_NNODES(&tree);
	Rprintf("  Total nb of nodes = %d\n", nnodes);
	for (i = 0; i < MAX_CHILDREN_PER_NODE+2; i++)
		nlinks_table[i] = 0;
	nleaves = 0;
	for (nid = 0; nid < nnodes; nid++) {
		node = GET_NODE(&tree, nid);
		nlinks = get_ACnode_nlinks(&tree, node);
		nlinks_table[nlinks]++;
		if (ISLEAF(node))
			nleaves++;
	}
	for (i = 0; i < MAX_CHILDREN_PER_NODE+2; i++)
		Rprintf("  - %d nodes (%.2f%) with %d links\n",
			nlinks_table[i],
			100.00 * nlinks_table[i] / nnodes,
			i);
	Rprintf("  Nb of leaf nodes (nleaves) = %d\n", nleaves);
	max_nn = get_max_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	min_nn = get_min_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	n1 = get_OptMaxNN(max_nn, min_nn);
	Rprintf("  - max_needed_nnodes(nleaves, TREE_DEPTH) = %d\n", max_nn);
	Rprintf("  - min_needed_nnodes(nleaves, TREE_DEPTH) = %d\n", min_nn);
	Rprintf("  - OptMaxNN(nleaves, TREE_DEPTH) = %d\n", n1);
	return R_NilValue;
}



/****************************************************************************
 *                             F. PREPROCESSING                             *
 ****************************************************************************/

static void add_pattern(ACtree *tree, const RoSeq *P, int P_offset)
{
	int P_id, depth, dmax, nid1, nid2, linktag;
	ACnode *node1, *node2;

	P_id = P_offset + 1;
	dmax = TREE_DEPTH(tree) - 1;
	for (depth = 0, nid1 = 0; depth <= dmax; depth++, nid1 = nid2) {
		node1 = GET_NODE(tree, nid1);
		linktag = CHAR2LINKTAG(tree, P->elts[depth]);
		if (linktag == NA_INTEGER)
			error("non base DNA letter found in Trusted Band "
			      "for pattern %d", P_id);
		nid2 = get_ACnode_link(tree, node1, linktag);
		if (depth < dmax) {
			if (nid2 != -1)
				continue;
			nid2 = new_ACnode(tree, depth + 1);
			set_ACnode_link(tree, node1, linktag, nid2);
			continue;
		}
		if (nid2 != -1) {
			node2 = GET_NODE(tree, nid2);
			_report_dup(P_offset, NODE_P_ID(node2));
		} else {
			nid2 = new_leafACnode(tree, P_id);
			set_ACnode_link(tree, node1, linktag, nid2);
		}
	}
	return;
}

/****************************************************************************
 * .Call entry point for preprocessing
 * -----------------------------------
 *
 * Arguments:
 *   tb:         the Trusted Band extracted from the original dictionary as a
 *               constant width DNAStringSet object;
 *   dup2unq0:   NULL or an integer vector of the same length as 'tb'
 *               containing the mapping from duplicated to unique patterns in
 *               the original dictionary;
 *   base_codes: the internal codes for A, C, G and T.
 *
 * See ACtree_asLIST() for a description of the returned SEXP.
 */

SEXP ACtree2_build(SEXP tb, SEXP dup2unq0, SEXP base_codes)
{
	ACtree tree;
	int tb_length, tb_width, P_offset;
	CachedXStringSet cached_tb;
	RoSeq P;
	SEXP ans, ans_names, ans_elt;

	tb_length = _get_XStringSet_length(tb);
	if (tb_length == 0)
		error("Trusted Band is empty");
	_init_dup2unq_buf(tb_length);
	tb_width = -1;
	cached_tb = _new_CachedXStringSet(tb);
	for (P_offset = 0; P_offset < tb_length; P_offset++) {
		/* skip duplicated patterns */
		if (dup2unq0 != R_NilValue
		 && INTEGER(dup2unq0)[P_offset] != NA_INTEGER)
			continue;
		P = _get_CachedXStringSet_elt_asRoSeq(&cached_tb, P_offset);
		if (tb_width == -1) {
			if (P.nelt == 0)
				error("first element in Trusted Band "
				      "is of length 0");
			tb_width = P.nelt;
			tree = new_ACtree(tb_length, tb_width, base_codes);
		} else if (P.nelt != tb_width) {
			error("element %d in Trusted Band has a different "
			      "length than first element", P_offset + 1);
		}
		add_pattern(&tree, &P, P_offset);
	}

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("ACtree"));
	SET_STRING_ELT(ans_names, 1, mkChar("dup2unq"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "ACtree" element */
	PROTECT(ans_elt = ACtree_asLIST(&tree));
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
 *                             G. MATCH FINDING                             *
 ****************************************************************************/

/*
 * We use indirect recursion for walking a sequence.
 * This indirect recursion involves the 2 core functions walk_shortseq() and
 * transition(): the former calls the latter which in turn calls the former.
 */
static int walk_shortseq(ACtree *tree, const char *seq, int seq_len);

/*
 * A trick is to have the path from the root node to the current node 'node'
 * stored at *negative* indices in 'seq_tail' i.e., if d is the depth of
 * 'node', then the path is seq_tail[-d], seq_tail[-d+1], seq_tail[-d+2], ...,
 * seq_tail[-1].
 */
static int transition(ACtree *tree, ACnode *node, int linktag, const char *seq_tail)
{
	int link, flink, newpath_len;
	const char *newpath;
/*
#ifdef DEBUG_BIOSTRINGS
	static int rec_level = -1;
	int node_depth;
	char format[20], pathbuf[2000];

	rec_level++;
	if (debug) {
		Rprintf("[DEBUG] ENTERING transition():");
		sprintf(format, "%%%ds", 1 + 2*rec_level);
		Rprintf(format, " ");
		node_depth = NODE_DEPTH(tree, node);
		snprintf(pathbuf, node_depth + 1, "%s", seq_tail - node_depth);
		Rprintf("nid=%d node_depth=%d linktag=%d path=%s\n",
			node - tree->nodes, node_depth, linktag, pathbuf);
	}
#endif
*/

	if (linktag == NA_INTEGER) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%d\n", 0);
		}
		rec_level--;
#endif
*/
		return 0;
	}
	link = get_ACnode_link(tree, node, linktag);
	if (link != -1) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%d\n", link);
		}
		rec_level--;
#endif
*/
		return link;
	}
	if (node == tree->nodes) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%d\n", 0);
		}
		rec_level--;
#endif
*/
		return 0;
	}
	flink = get_ACnode_flink(tree, node);
	if (flink == -1) {
		newpath_len = NODE_DEPTH(tree, node) - 1;
		newpath = seq_tail - newpath_len;
		flink = walk_shortseq(tree, newpath, newpath_len);
		set_ACnode_flink(tree, node, flink);
	}
	link = transition(tree, GET_NODE(tree, flink), linktag, seq_tail);
	set_ACnode_link(tree, node, linktag, link); /* sets a shortcut */
/*
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG]  LEAVING transition():");
		Rprintf(format, " ");
		Rprintf("link=%d\n", link);
	}
	rec_level--;
#endif
*/
	return link;
}

/*
 * Does NOT report matches. Would not find matches anyway because it's used
 * to walk on short sequences i.e. sequences with a length < TREE_DEPTH(tree).
 */
static int walk_shortseq(ACtree *tree, const char *seq, int seq_len)
{
	ACnode *node;
	int n, linktag, nid;
	const char *seq_tail;

	node = GET_NODE(tree, 0);
	for (n = 0, seq_tail = seq; n < seq_len; n++, seq_tail++) {
		linktag = CHAR2LINKTAG(tree, *seq_tail);
		nid = transition(tree, node, linktag, seq_tail);
		node = GET_NODE(tree, nid);
	}
	return nid;
}

/* Does report matches */
static void walk_subject(ACtree *tree, const RoSeq *S)
{
	ACnode *node;
	int n, linktag, nid;
	const char *S_tail;

	node = GET_NODE(tree, 0);
	for (n = 1, S_tail = S->elts; n <= S->nelt; n++, S_tail++) {
		linktag = CHAR2LINKTAG(tree, *S_tail);
		nid = transition(tree, node, linktag, S_tail);
		node = GET_NODE(tree, nid);
		if (ISLEAF(node))
			_MIndex_report_match(NODE_P_ID(node) - 1, n);
	}
	return;
}

/* Does report matches */
static void walk_nonfixed_subject(ACtree *tree, const RoSeq *S)
{
	error("walk_nonfixed_subject(): implement me");
	return;
}

/* Entry point for the MATCH FINDING section */
void _match_ACtree2(SEXP pptb, const RoSeq *S, int fixedS)
{
	ACtree tree;

	tree = pptb_asACtree(pptb);
	if (fixedS)
		walk_subject(&tree, S);
	else
		walk_nonfixed_subject(&tree, S);
	return;
}

