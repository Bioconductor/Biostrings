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

SEXP debug_match_pdict_ACtree2()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}



/****************************************************************************
 *                    A. DEFINES AND LOW-LEVEL UTILITIES                    *
 ****************************************************************************/

#define MAX_CHILDREN_PER_NODE 4
#define LINKTAG_BITSHIFT 28
#define LINKTAG_BITMASK (3 << LINKTAG_BITSHIFT)
#define MAX_DEPTH ((1 << LINKTAG_BITSHIFT) - 1)
#define ISLEAF_BIT (1 << 30)
#define ISEXTENDED_BIT (ISLEAF_BIT << 1) /* strongest bit for 32-bit integers */
#define MAX_P_ID (ISLEAF_BIT - 1)  /* P_id values are encoded on 30 bits */

typedef struct acnode {
	int attribs;
	int nid_or_eid;
} ACnode;

#define INTS_PER_NODE (sizeof(ACnode) / sizeof(int))
#define ISEXTENDED(node) ((node)->attribs & ISEXTENDED_BIT)
#define ISLEAF(node) ((node)->attribs & ISLEAF_BIT)
#define GET_DEPTH(node) ((node)->attribs & MAX_DEPTH)
#define GET_P_ID(node) ((node)->attribs & MAX_P_ID)

typedef struct acnode_extension {
	int link_nid[MAX_CHILDREN_PER_NODE];
	int flink_nid;
} ACnodeExtension;

#define INTS_PER_EXTENSION (sizeof(ACnodeExtension) / sizeof(int))

typedef struct actree {
	int depth;  /* this is the depth of all leaf nodes */
	ACnode *nodes;
	int nnodes, nodes_buflength;
	ACnodeExtension *extensions;
	int nextensions, extensions_buflength;
	ByteTrTable char2linktag;
} ACtree;

#define TREE_DEPTH(tree) ((tree)->depth)
#define TREE_NNODES(tree) ((tree)->nnodes)
#define GET_NODE(tree, nid) ((tree)->nodes + (nid))
#define CHAR2LINKTAG(tree, c) ((tree)->char2linktag[(unsigned char) (c)])

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

	if (tree->nextensions >= tree->extensions_buflength)
		extend_extensions_buffer(tree);
	extension = tree->extensions + (eid = tree->nextensions++);
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

	tag_length = tree->nnodes * INTS_PER_NODE;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), tree->nodes, tag_length * sizeof(int));
	PROTECT(ans = new_XInteger_from_tag("XInteger", tag));
	UNPROTECT(2);
	return ans;
}

static SEXP ACtree_extensions_asXInteger(ACtree *tree)
{
	int tag_length;
	SEXP ans, tag;

	tag_length = tree->nextensions * INTS_PER_EXTENSION;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy(INTEGER(tag), tree->extensions, tag_length * sizeof(int));
	PROTECT(ans = new_XInteger_from_tag("XInteger", tag));
	UNPROTECT(2);
	return ans;
}



/****************************************************************************
 *                      B. Aho-Corasick tree low-level API                  *
 ****************************************************************************/

static int new_ACnode(ACtree *tree, int depth)
{
	ACnode *node;
	int nid;

	if (depth >= TREE_DEPTH(tree))
		error("new_ACnode(): depth >= TREE_DEPTH(tree)");
	if (tree->nnodes >= tree->nodes_buflength)
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

	if (tree->nnodes >= tree->nodes_buflength)
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
		extension = tree->extensions + node->nid_or_eid;
		return extension->link_nid[linktag];
	}
	/* the node has no extension and is not a leaf node */
	if (linktag == (node->attribs >> LINKTAG_BITSHIFT))
		return node->nid_or_eid;
	return -1;
}

/* on a leaf node, set_ACnode_flink() should always have been called first
   so we can assume that the node is already extended */
static void set_ACnode_link(ACtree *tree, ACnode *node, int linktag, int nid)
{
	ACnodeExtension *extension;

	if (node->nid_or_eid == -1) {
		/* cannot a leaf node (see assumption above) and
                   no need to extend it */
		node->attribs |= linktag << LINKTAG_BITSHIFT;
		node->nid_or_eid = nid;
		return;
	}
	if (!ISEXTENDED(node)) {
		/* again, cannot be a leaf node (see assumption above) */
		extend_ACnode(tree, node);
	}
	extension = tree->extensions + node->nid_or_eid;
	extension->link_nid[linktag] = nid;
	return;
}

static int get_ACnode_flink(ACtree *tree, ACnode *node)
{
	ACnodeExtension *extension;

	if (!ISEXTENDED(node))
		return -1;
	extension = tree->extensions + node->nid_or_eid;
	return extension->flink_nid;
}

static void set_ACnode_flink(ACtree *tree, ACnode *node, int nid)
{
	error("set_ACnode_flink(): implement me");
	return;
}

static void init_ACtree(ACtree *tree, int tb_length, int tb_width,
		int max_needed_nodes, SEXP base_codes)
{
	if (tb_length > MAX_P_ID)
		error("init_ACtree(): tb_length > MAX_P_ID");
	if (tb_width > MAX_DEPTH)
		error("init_ACtree(): tb_width > MAX_DEPTH");
	tree->depth = tb_width;

	tree->nodes = NULL;
	tree->nnodes = tree->nodes_buflength = 0;
	/* for testing */
	tree->nodes = Salloc((long) max_needed_nodes, ACnode);
	tree->nnodes = 0;
	tree->nodes_buflength = max_needed_nodes;

	tree->extensions = NULL;
	tree->nextensions = tree->extensions_buflength = 0;
	/* for testing */
	max_needed_nodes = max_needed_nodes / 4;
	tree->extensions = Salloc((long) max_needed_nodes, ACnodeExtension);
	tree->nextensions = 0;
	tree->extensions_buflength = max_needed_nodes;

	_init_byte2offset_with_INTEGER(tree->char2linktag, base_codes, 1);
	new_ACnode(tree, 0);  /* create the root node */
	return;
}

static SEXP ACtree_asLIST(ACtree *tree)
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("nodes"));
	SET_STRING_ELT(ans_names, 1, mkChar("extensions"));
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
	tree.nnodes = LENGTH(tag) / INTS_PER_NODE;
	tree.nodes_buflength = -1;
	tag = _get_ACtree2_extensions_tag(pptb);
	tree.extensions = (ACnodeExtension *) INTEGER(tag);
	tree.nextensions = LENGTH(tag) / INTS_PER_EXTENSION;
	tree.extensions_buflength = -1;
	base_codes = _get_ACtree2_base_codes(pptb);
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
 *                       C. STATS AND DEBUG UTILITIES                       *
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
	    i, nid, nlinks;

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
	Rprintf("  Nb of leaf nodes = %d\n", nleaves);
	return R_NilValue;
}



/****************************************************************************
 *                             D. PREPROCESSING                             *
 ****************************************************************************/

static int max_needed_nodes(int tb_length, int tb_width)
{
	int nnodes, depth, pow;

	nnodes = 0;
	for (depth = 0, pow = 1; depth <= tb_width; depth++) {
		if (pow >= tb_length)
			return nnodes + (tb_width - depth + 1) * tb_length;
		nnodes += pow;
		pow *= MAX_CHILDREN_PER_NODE;
	}
	return nnodes;
}

static void pp_pattern(ACtree *tree, const RoSeq *P, int P_offset)
{
	int P_id, depth, dmax, nid1, nid2, linktag;;
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
			_report_dup(P_offset, GET_P_ID(node2));
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
	int tb_length, tb_width, P_offset, m;
	CachedXStringSet cached_tb;
	RoSeq P;
	SEXP ans, ans_names, ans_elt;

	tb_length = _get_XStringSet_length(tb);
	if (tb_length == 0)
		error("Trusted Band is empty");
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_NODE)
		error("Biostrings internal error in build_ACtree2(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_NODE");
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
			m = max_needed_nodes(tb_length, tb_width);
			init_ACtree(&tree, tb_length, tb_width, m, base_codes);
		} else if (P.nelt != tb_width) {
			error("element %d in Trusted Band has a different "
			      "length than first element", P_offset + 1);
		}
		pp_pattern(&tree, &P, P_offset);
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
 *                             E. MATCH FINDING                             *
 ****************************************************************************/

