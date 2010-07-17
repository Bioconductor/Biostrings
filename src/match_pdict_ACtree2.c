/****************************************************************************
 *     A fast and compact implementation of the Aho-Corasick algorithm      *
 *                     for rectangular DNA dictionaries                     *
 *                                                                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#include <stdlib.h> /* for div() */
#include <limits.h> /* for UINT_MAX */

static int debug = 0;


/*
 * Internal representation of the Aho-Corasick tree
 * ------------------------------------------------
 *
 * For this Aho-Corasick implementation, we take advantage of 2 important
 * properties of the input dictionary (aka pattern set):
 *   1. It's rectangular (i.e. all patterns have the same length).
 *   2. It's based on a 4-letter alphabet (4-ary tree). Note that this tree
 *      becomes an oriented graph when we start adding the failure links (or
 *      the shortcut links) to it.
 * Failure/shortcut links are not precomputed, but computed on-the-fly when
 * the tree is used to walk along a subject.
 * A node is represented with either 2 ints (8 bytes) before extension,
 * or 7 ints (28 bytes) once it has been extended. Extending a node is done
 * by linking its 2-int part to a 5-int extension (so the 2-int part doesn't
 * need to be reallocated).
 * Two separate buffers are used to store the nodes: one for the 2-int parts
 * of the nodes (every node has one, whether it's extended or not) and one for
 * the 5-int extensions. Since the number of nodes doesn't change during the
 * life of the tree, the first buffer will grow only while we are building the
 * tree (preprocessing) and then it will not change anymore (i.e. when the tree
 * is used to walk along a subject). However, since failure/shortcut links are
 * not precomputed, some nodes will need to be extended (and therefore the
 * second buffer will grow) in order to store the links that are computed
 * on-the-fly.
 * The two node buffers are made of R integer vectors so they can be
 * serialized.
 * Some testing with real data shows that, typically, less than 10% of the
 * nodes are already extended right after preprocessing (i.e. before any use
 * of the PDict object), and that this percentage grows up to 30% or 40%
 * during the typical life of the PDict object (e.g. after it has been used
 * to walk along a full genome).
 * Using node extensions makes the tree more compact in memory than with
 * fixed-size nodes. The latter tend to waste memory when a high percentage
 * of them have 1 child only. Also the two buffer approach allows the tree to
 * grow smoothly without any reallocations while links are added to it.
 * We use unsigned ints for the node ids so, on Intel i386/x86_64 platforms,
 * the maximum number of nodes in a tree is 2^32-1 nodes (UINT_MAX is used as
 * a special value).
 */

#define MAX_CHILDREN_PER_NODE 4  /* do NOT change this */



/****************************************************************************
 *  A. GENERAL UTILITIES FOR COMPUTING THE NB OF NODES AND NODE EXTENSIONS  *
 ****************************************************************************/

/*
 * nleaves = number of leaves in the tree (it's also the number of unique
 *           patterns so 'nleaves' <= 'tb_length');
 * depth = depth of the tree (it's also the length of each pattern so 'depth'
 *           = 'tb_width').
 * 'nleaves' is assumed to be '<= 4^depth'. The functions below do not check
 * this and will return invalid results if this is not true.
 */
static unsigned int count_max_needed_nnodes(int nleaves, int depth)
{
	unsigned int res;
	int inc, d;

	res = 0U;
	inc = 1;
	for (d = 0; d <= depth; d++) {
		if (inc >= nleaves) {
			res += (unsigned int) (depth + 1 - d) * nleaves;
			break;
		}
		res += (unsigned int) inc;
		inc *= MAX_CHILDREN_PER_NODE;
	}
	return res;
}

static unsigned int count_min_needed_nnodes(int nleaves, int depth)
{
	unsigned int res;
	int inc, d;
	div_t q;

	res = 0U;
	inc = nleaves;
	for (d = depth; d >= 0; d--) {
		if (inc == 1) {
			res += (unsigned int) (d + 1);
			break;
		}
		res += (unsigned int) inc;
		q = div(inc, MAX_CHILDREN_PER_NODE);
		inc = q.quot;
		if (q.rem != 0)
			inc++;
	}
	return res;
}

static unsigned int count_max_needed_nnodeexts_at_pp_time(int nleaves, int depth)
{
	unsigned int res;
	int inc, d, four_power_d;
	div_t q;

	res = 0U;
	for (d = depth - 1; d >= 0; d--) {
		q = div(nleaves, 2);
		inc = q.quot;
		nleaves = inc + q.rem;
		if (d >= 16) {
			res += (unsigned int) inc;
			continue;
		}
		four_power_d = 1 << (2 * d);  /* = MAX_CHILDREN_PER_NODE ^ d */
		if (nleaves <= four_power_d) {
			res += (unsigned int) inc;
			continue;
		}
		res += count_max_needed_nnodes(four_power_d, d);
		break;
	}
	return res;
}

#ifdef DEBUG_BIOSTRINGS
/* I've checked by hand the output of this function up to depth=3 nleaves=17
 * and it looked correct... pffff!!! :-b */
static void debug_node_counting_functions(int maxdepth)
{
	int depth, four_power_d, nleaves, delta;
	unsigned int max_nn, min_nn, n2;

	Rprintf("[DEBUG] debug_node_counting_functions():\n");
	for (depth = 1; depth <= maxdepth; depth++) {
		four_power_d = 1 << (2 * depth);  /* = MAX_CHILDREN_PER_NODE ^ depth */
		for (nleaves = 1; nleaves <= four_power_d; nleaves++) {
			max_nn = count_max_needed_nnodes(nleaves, depth);
			min_nn = count_min_needed_nnodes(nleaves, depth);
			n2 = count_max_needed_nnodeexts_at_pp_time(nleaves, depth);
			delta = max_nn - nleaves - n2;  /* should always be >= 0 */
			Rprintf("  depth=%d nleaves=%d --> ", depth, nleaves);
			Rprintf("max_nn=%u min_nn=%u n2=%u max_nn-nleaves-n2=%d\n",
				max_nn, min_nn, n2, delta);
			if (delta < 0)
				error("max_nn-nleaves-n2 < 0");
		}
	}
	return;
}
#endif



/****************************************************************************
 *                         B. ACnode AND ACnodeBuf                          *
 ****************************************************************************/

#define NOT_AN_ID UINT_MAX

/*
 * A node id (nid) is represented by an unsigned int. The ANSI C standard only
 * guarantees that an (unsigned) int will take at least 2 bytes in memory i.e.
 * that UINT_MAX will always be >= 2^16-1 even though, on most modern
 * platforms (including PC and Mac), an int will take at least 4 bytes in
 * memory i.e. UINT_MAX will be 2^32-1 or more.
 * The same apply to nodeext ids (eid).
 */
typedef struct acnode {
	int attribs;
	unsigned int nid_or_eid;
} ACnode;

#define INTS_PER_NODE (sizeof(ACnode) / sizeof(int))

/*
 * We must have:
 *   (a) ACNODEBUF_MAX_NBLOCK * ACNODEBUF_MAX_NELT_PER_BLOCK <= UINT_MAX + 1
 *   (b) ACNODEBUF_MAX_NELT_PER_BLOCK * INTS_PER_NODE <= INT_MAX
 * On 64-bit Unix/Linux (4 bytes per int):
 *   UINT_MAX = 2^32-1 = 4294967295U
 *   INT_MAX = 2^31-1 = 2147483647
 * The following settings are assuming UINT_MAX >= 2^32-1 and
 * INT_MAX >= 2^31-1. They result in blocks of size 32 MB.
 */
#define ACNODEBUF_MAX_NBLOCK 1024  /* = 2^10 */
#define ACNODEBUF_MAX_NELT_PER_BLOCK 4194304U  /* = 2^22 */

typedef struct acnodebuf {
	SEXP bab;  /* Big Atomic Buffer */
	int *nblock;
	int *lastblock_nelt;
	ACnode *block[ACNODEBUF_MAX_NBLOCK];
} ACnodeBuf;

#define _IS_ROOTNODE(nodebuf, node) ((nodebuf)->block[0] == (node))

static ACnode *get_node_from_buf(ACnodeBuf *buf, unsigned int nid)
{
	unsigned int b, i;

	b = nid >> 22U;
	i = nid & (ACNODEBUF_MAX_NELT_PER_BLOCK - 1U);
	return buf->block[b] + i;
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_nodebuf_max_nblock()
{
	return ScalarInteger(ACNODEBUF_MAX_NBLOCK);
}

static int ACnodeBuf_is_full(ACnodeBuf *buf)
{
	return *(buf->nblock) == 0
	       || *(buf->lastblock_nelt) >= ACNODEBUF_MAX_NELT_PER_BLOCK;
}

static unsigned int get_ACnodeBuf_nelt(ACnodeBuf *buf)
{
	int nblock;

	nblock = *(buf->nblock);
	if (nblock == 0)
		return 0U;
	return (unsigned int) (nblock - 1) * ACNODEBUF_MAX_NELT_PER_BLOCK
	       + *(buf->lastblock_nelt);
}

static ACnodeBuf new_ACnodeBuf(SEXP bab)
{
	ACnodeBuf buf;
	SEXP bab_blocks;
	int nblock, b;

	buf.bab = bab;
	nblock = *(buf.nblock = _get_BAB_nblock_ptr(bab));
	buf.lastblock_nelt = _get_BAB_lastblock_nelt_ptr(bab);
	bab_blocks = _get_BAB_blocks(bab);
	for (b = 0; b < nblock; b++)
		buf.block[b] = (ACnode *) INTEGER(VECTOR_ELT(bab_blocks, b));
	return buf;
}

static void extend_ACnodeBuf(ACnodeBuf *buf)
{
	int length;
	SEXP bab_block;

	length = ACNODEBUF_MAX_NELT_PER_BLOCK * INTS_PER_NODE;
	bab_block = _IntegerBAB_addblock(buf->bab, length);
	/* sync 'buf->block' with 'buf->bab' */
	buf->block[*(buf->nblock) - 1] = (ACnode *) INTEGER(bab_block);
	return;
}

static unsigned int new_nid(ACnodeBuf *buf)
{
	unsigned int nid;

	if (ACnodeBuf_is_full(buf))
		extend_ACnodeBuf(buf);
	nid = get_ACnodeBuf_nelt(buf);
	if (nid == NOT_AN_ID)
		error("reached max number of nodes (%u)", NOT_AN_ID);
	(*(buf->lastblock_nelt))++;
	return nid;
}



/****************************************************************************
 *                       C. ACnodeext AND ACnodeextBuf                      *
 ****************************************************************************/

typedef struct acnodeext {
	unsigned int link_nid[MAX_CHILDREN_PER_NODE];
	unsigned int flink_nid;
} ACnodeext;

#define INTS_PER_NODEEXT (sizeof(ACnodeext) / sizeof(int))

/*
 * We must have:
 *   (a) ACNODEEXTBUF_MAX_NBLOCK * ACNODEEXTBUF_MAX_NELT_PER_BLOCK <=
 *       UINT_MAX + 1
 *   (b) ACNODEEXTBUF_MAX_NELT_PER_BLOCK * INTS_PER_NODEEXT <= INT_MAX
 * The following settings are assuming UINT_MAX >= 2^32-1 and
 * INT_MAX >= 2^31-1. They result in blocks of size 80 MB.
 */
#define ACNODEEXTBUF_MAX_NBLOCK 1024  /* = 2^10 */
#define ACNODEEXTBUF_MAX_NELT_PER_BLOCK 4194304U  /* = 2^22 */

typedef struct acnodeextbuf {
	SEXP bab;  /* Big Atomic Buffer */
	int *nblock;
	int *lastblock_nelt;
	ACnodeext *block[ACNODEEXTBUF_MAX_NBLOCK];
} ACnodeextBuf;

static ACnodeext *get_nodeext_from_buf(ACnodeextBuf *buf, unsigned int eid)
{
	unsigned int b, i;

	b = eid >> 22U;
	i = eid & (ACNODEEXTBUF_MAX_NELT_PER_BLOCK - 1U);
	return buf->block[b] + i;
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_nodeextbuf_max_nblock()
{
	return ScalarInteger(ACNODEEXTBUF_MAX_NBLOCK);
}

static int ACnodeextBuf_isfull(ACnodeextBuf *buf)
{
	return *(buf->nblock) == 0
	       || *(buf->lastblock_nelt) >= ACNODEEXTBUF_MAX_NELT_PER_BLOCK;
}

static unsigned int get_ACnodeextBuf_nelt(ACnodeextBuf *buf)
{
	int nblock;

	nblock = *(buf->nblock);
	if (nblock == 0)
		return 0U;
	return (unsigned int) (nblock - 1) * ACNODEEXTBUF_MAX_NELT_PER_BLOCK
	       + *(buf->lastblock_nelt);
}

static ACnodeextBuf new_ACnodeextBuf(SEXP bab)
{
	ACnodeextBuf buf;
	SEXP bab_blocks;
	int nblock, b;

	buf.bab = bab;
	nblock = *(buf.nblock = _get_BAB_nblock_ptr(bab));
	buf.lastblock_nelt = _get_BAB_lastblock_nelt_ptr(bab);
	bab_blocks = _get_BAB_blocks(bab);
	for (b = 0; b < nblock; b++)
		buf.block[b] = (ACnodeext *) INTEGER(VECTOR_ELT(bab_blocks, b));
	return buf;
}

static void extend_ACnodeextBuf(ACnodeextBuf *buf)
{
	int length;
	SEXP bab_block;

	length = ACNODEEXTBUF_MAX_NELT_PER_BLOCK * INTS_PER_NODEEXT;
	bab_block = _IntegerBAB_addblock(buf->bab, length);
	/* sync 'buf->block' with 'buf->bab' */
	buf->block[*(buf->nblock) - 1] = (ACnodeext *) INTEGER(bab_block);
	return;
}

static unsigned int new_eid(ACnodeextBuf *buf)
{
	unsigned int eid;

	if (ACnodeextBuf_isfull(buf))
		extend_ACnodeextBuf(buf);
	eid = get_ACnodeextBuf_nelt(buf);
	(*(buf->lastblock_nelt))++;
	return eid;
}



/****************************************************************************
 *                                 D. ACtree                                *
 ****************************************************************************/

#define LINKTAG_BITSHIFT 28
#define MAX_DEPTH ((1 << LINKTAG_BITSHIFT) - 1)
#define ISLEAF_BIT (1 << 30)
#define ISEXTENDED_BIT (ISLEAF_BIT << 1) /* strongest bit for 32-bit integers */
#define MAX_P_ID (ISLEAF_BIT - 1)  /* P_id values are encoded on 30 bits */

#define IS_EXTENDEDNODE(node) ((node)->attribs & ISEXTENDED_BIT)
#define _NODE_DEPTH(node) ((node)->attribs & MAX_DEPTH)
/* result of NODE_P_ID() is undefined on a non-leaf node */
#define NODE_P_ID(node) ((node)->attribs & MAX_P_ID)

/* --- .Call ENTRY POINT --- */
SEXP debug_match_pdict_ACtree2()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
	if (debug) {
		Rprintf("[DEBUG] debug_match_pdict_ACtree2():\n");
		Rprintf("  INTS_PER_NODE=%d INTS_PER_NODEEXT=%d\n",
			INTS_PER_NODE, INTS_PER_NODEEXT);
		Rprintf("  LINKTAG_BITSHIFT=%d\n"
			"  MAX_DEPTH=%d\n"
			"  ISLEAF_BIT=%d ISEXTENDED_BIT=%d\n"
			"  MAX_P_ID=%d\n",
			LINKTAG_BITSHIFT,
			MAX_DEPTH,
			ISLEAF_BIT, ISEXTENDED_BIT,
			MAX_P_ID);
		debug_node_counting_functions(3);
	}
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

/*
 * Always set 'max_nodeextbuf_nelt' to 0U (no max) and 'dont_extend_nodes' to
 * 0 during preprocessing.
 */
typedef struct actree {
	int depth;  /* this is the depth of all leaf nodes */
	ACnodeBuf nodebuf;
	ACnodeextBuf nodeextbuf;
	ByteTrTable char2linktag;
	unsigned int max_nodeextbuf_nelt;  /* 0U means "no max" */
	int dont_extend_nodes;  /* always at 0 during preprocessing */
} ACtree;

#define GET_NODEEXT(tree, eid) get_nodeext_from_buf(&((tree)->nodeextbuf), eid)

static void extend_ACnode(ACtree *tree, ACnode *node)
{
	ACnodeext *nodeext;
	unsigned int eid;
	int i, linktag;

	eid = new_eid(&(tree->nodeextbuf));
	if (eid + 1U == tree->max_nodeextbuf_nelt) {
		tree->dont_extend_nodes = 1;
		warning("Reached max nb of node extensions (%u) so I will\n"
			"stop extending the nodes of this ACtree2 object.\n"
			"As a consequence not all new links and failure\n"
			"links will be set. This might (slightly) affect\n"
                        "speed but not the results.",
			tree->max_nodeextbuf_nelt);
	}
	nodeext = GET_NODEEXT(tree, eid);
	for (i = 0; i < MAX_CHILDREN_PER_NODE; i++)
		nodeext->link_nid[i] = NOT_AN_ID;
	nodeext->flink_nid = NOT_AN_ID;
	if (node->nid_or_eid != NOT_AN_ID) {
		/* this is correct because 'node' cannot be a leaf node */
		linktag = node->attribs >> LINKTAG_BITSHIFT;
		nodeext->link_nid[linktag] = node->nid_or_eid;
	}
	node->nid_or_eid = eid;
	/* sets the ISEXTENDED_BIT bit to 1 */
	node->attribs |= ISEXTENDED_BIT;
	return;
}



/****************************************************************************
 *                           E. Aho-Corasick tree API                       *
 ****************************************************************************/

/*
 * Formal API
 */
#define TREE_SIZE(tree) get_ACnodeBuf_nelt(&((tree)->nodebuf)) /* nb nodes */
#define TREE_DEPTH(tree) ((tree)->depth)
#define GET_NODE(tree, nid) get_node_from_buf(&((tree)->nodebuf), nid)
#define IS_ROOTNODE(tree, node) _IS_ROOTNODE(&((tree)->nodebuf), node)
#define IS_LEAFNODE(node) ((node)->attribs & ISLEAF_BIT)
#define NODE_DEPTH(tree, node) \
		(IS_LEAFNODE(node) ? TREE_DEPTH(tree) : _NODE_DEPTH(node))
#define CHAR2LINKTAG(tree, c) ((tree)->char2linktag[(unsigned char) (c)])
#define NEW_NODE(tree, depth) new_ACnode(tree, depth)
#define NEW_LEAFNODE(tree, P_id) new_leafACnode(tree, P_id)
#define GET_NODE_LINK(tree, node, linktag) \
		get_ACnode_link(tree, node, linktag)
#define SET_NODE_LINK(tree, node, linktag, nid) \
		set_ACnode_link(tree, node, linktag, nid)
#define GET_NODE_FLINK(tree, node) get_ACnode_flink(tree, node)
#define SET_NODE_FLINK(tree, node, nid) set_ACnode_flink(tree, node, nid)

/*
 * API implementation
 */

static unsigned int new_ACnode(ACtree *tree, int depth)
{
	ACnodeBuf *nodebuf;
	unsigned int nid;
	ACnode *node;

	if (depth >= TREE_DEPTH(tree))
		error("new_ACnode(): depth >= TREE_DEPTH(tree)");
	nodebuf = &(tree->nodebuf);
	nid = new_nid(nodebuf);
	node = get_node_from_buf(nodebuf, nid);
	/* this sets the ISEXTENDED_BIT and ISLEAF_BIT bits to 0 */
	node->attribs = depth;
	node->nid_or_eid = NOT_AN_ID;
	return nid;
}

static unsigned int new_leafACnode(ACtree *tree, int P_id)
{
	ACnodeBuf *nodebuf;
	unsigned int nid;
	ACnode *node;

	nodebuf = &(tree->nodebuf);
	nid = new_nid(nodebuf);
	node = get_node_from_buf(nodebuf, nid);
	/* this sets the ISEXTENDED_BIT bit to 0 and ISLEAF_BIT bit to 1 */
	node->attribs = ISLEAF_BIT | P_id;
	node->nid_or_eid = NOT_AN_ID;
	return nid;
}

static unsigned int get_ACnode_link(ACtree *tree, ACnode *node, int linktag)
{
	ACnodeext *nodeext;

	if (node->nid_or_eid == NOT_AN_ID)
		return NOT_AN_ID;
	if (IS_EXTENDEDNODE(node)) {
		nodeext = GET_NODEEXT(tree, node->nid_or_eid);
		return nodeext->link_nid[linktag];
	}
	/* the node has no extension and is not a leaf node */
	if (linktag == (node->attribs >> LINKTAG_BITSHIFT))
		return node->nid_or_eid;
	return NOT_AN_ID;
}

/*
 * set_ACnode_flink() should always have been called before set_ACnode_link()
 * on a leaf node. So we can assume that, if set_ACnode_link() is called on a
 * leaf node, then this node is already extended.
 */
static void set_ACnode_link(ACtree *tree, ACnode *node, int linktag, unsigned int nid)
{
	ACnodeext *nodeext;

	if (node->nid_or_eid == NOT_AN_ID) {
		/* cannot be a leaf node (see assumption above) and
		   no need to extend it */
		node->attribs |= linktag << LINKTAG_BITSHIFT;
		node->nid_or_eid = nid;
		return;
	}
	if (!IS_EXTENDEDNODE(node)) {
		if (tree->dont_extend_nodes) {
			/* NEVER during preprocessing */
			return;
		}
		/* again, cannot be a leaf node (see assumption above) */
		extend_ACnode(tree, node);
	}
	nodeext = GET_NODEEXT(tree, node->nid_or_eid);
	nodeext->link_nid[linktag] = nid;
	return;
}

static unsigned int get_ACnode_flink(ACtree *tree, const ACnode *node)
{
	ACnodeext *nodeext;

	if (!IS_EXTENDEDNODE(node))
		return NOT_AN_ID;
	nodeext = GET_NODEEXT(tree, node->nid_or_eid);
	return nodeext->flink_nid;
}

static void set_ACnode_flink(ACtree *tree, ACnode *node, unsigned int nid)
{
	ACnodeext *nodeext;

	if (!IS_EXTENDEDNODE(node)) {
		if (tree->dont_extend_nodes)
			return;
		extend_ACnode(tree, node);
	}
	nodeext = GET_NODEEXT(tree, node->nid_or_eid);
	//Rprintf("set flink: %d --> %d\n", node - GET_NODE(tree, 0U), nid);
	nodeext->flink_nid = nid;
	return;
}

/*
 * Not part of the API
 */

static ACtree new_ACtree(int tb_length, int tb_width, SEXP base_codes,
		SEXP nodebuf_ptr, SEXP nodeextbuf_ptr)
{
	ACtree tree;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] new_ACtree():\n"
			"  tb_length=%d tb_width=%d\n",
			tb_length, tb_width);
	}
#endif
	if (tb_length > MAX_P_ID)
		error("new_ACtree(): tb_length > MAX_P_ID");
	if (tb_width > MAX_DEPTH)
		error("new_ACtree(): tb_width > MAX_DEPTH");
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_NODE)
		error("Biostrings internal error in new_ACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_NODE");

	tree.depth = tb_width;
	tree.nodebuf = new_ACnodeBuf(nodebuf_ptr);
	tree.nodeextbuf = new_ACnodeextBuf(nodeextbuf_ptr);
	_init_byte2offset_with_INTEGER(tree.char2linktag, base_codes, 1);
	tree.max_nodeextbuf_nelt = 0U;
	tree.dont_extend_nodes = 0;
	NEW_NODE(&tree, 0);  /* create the root node */
	return tree;
}

static unsigned int a_nice_max_nodeextbuf_nelt(int nnodes)
{
	unsigned int n, rem;

	/* when nb of nodeexts has reached 0.40 * nb of nodes, then the size
	   in memory of the 2 buffers (nodebuf and nodeextbuf) is about the
	   same */
	n = (unsigned int) (0.40 * nnodes);
	/* then we round up to the closer multiple of
	   ACNODEBUF_MAX_NELT_PER_BLOCK so we don't waste space in the last
	   block */
	rem = n % ACNODEBUF_MAX_NELT_PER_BLOCK;
	if (rem != 0U)
		n += ACNODEBUF_MAX_NELT_PER_BLOCK - rem;
	return n;
}

static ACtree pptb_asACtree(SEXP pptb)
{
	ACtree tree;
	SEXP base_codes;
	unsigned int max_nelt, nelt;

	tree.depth = _get_PreprocessedTB_width(pptb);
	tree.nodebuf = new_ACnodeBuf(_get_ACtree2_nodebuf_ptr(pptb));
	tree.nodeextbuf = new_ACnodeextBuf(_get_ACtree2_nodeextbuf_ptr(pptb));
	base_codes = _get_PreprocessedTB_base_codes(pptb);
	if (LENGTH(base_codes) != MAX_CHILDREN_PER_NODE)
		error("Biostrings internal error in pptb_asACtree(): "
		      "LENGTH(base_codes) != MAX_CHILDREN_PER_NODE");
	_init_byte2offset_with_INTEGER(tree.char2linktag, base_codes, 1);
/*
  Using max_nelt = 0U will turn off the "dont_extend_nodes" feature
  for now. Seems like having dont_extend_nodes at 1 causes segfaults.
  To reproduce, put
    tree.dont_extend_nodes = 1;
  just before the return statement below and reinstall (i.e. recompile)
  Biostrings. Then run the following code:
    library(Biostrings)
    dict0 <- DNAStringSet(c("TACCNG", "TAGT", "CGGNT", "AGTAG", "TAGT"))
    pdict <- PDict(dict0, tb.end=3)
    subject <- DNAString("TAGTACCAGTTTCGGG")
    m0 <- matchPDict(pdict, subject)
    m1 <- matchPDict(pdict, subject, max.mismatch=1)
  Is the "dont_extend_nodes" feature reasonable i.e. is it safe to
  not set all new links and failure links? Still need to think about it.
*/
	//max_nelt = a_nice_max_nodeextbuf_nelt(TREE_SIZE(&tree));
	max_nelt = 0U;
	tree.max_nodeextbuf_nelt = max_nelt;
	nelt = get_ACnodeextBuf_nelt(&(tree.nodeextbuf));
	tree.dont_extend_nodes = max_nelt != 0U && nelt >= max_nelt;
	return tree;
}



/****************************************************************************
 *                       F. STATS AND DEBUG UTILITIES                       *
 ****************************************************************************/

static void print_ACnode(ACtree *tree, ACnode *node)
{
	error("print_ACnode(): implement me");
	return;
}

static int get_ACnode_nlink(ACtree *tree, ACnode *node)
{
	int nlink, linktag;

	nlink = get_ACnode_flink(tree, node) != NOT_AN_ID;
	for (linktag = 0; linktag < MAX_CHILDREN_PER_NODE; linktag++)
		if (get_ACnode_link(tree, node, linktag) != NOT_AN_ID)
			nlink++;
	return nlink;
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_nnodes(SEXP pptb)
{
	ACtree tree;

	tree = pptb_asACtree(pptb);
	return ScalarInteger(TREE_SIZE(&tree));
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_print_nodes(SEXP pptb)
{
	ACtree tree;
	unsigned int nnodes, nid;
	ACnodeBuf *nodebuf;
	ACnode *node;

	tree = pptb_asACtree(pptb);
	nnodes = TREE_SIZE(&tree);
	nodebuf = &(tree.nodebuf);
	for (nid = 0U; nid < nnodes; nid++) {
		node = get_node_from_buf(nodebuf, nid);
		print_ACnode(&tree, node);
	}
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_summary(SEXP pptb)
{
	ACtree tree;
	unsigned int nnodes, nlink_table[MAX_CHILDREN_PER_NODE+2],
		     nid, max_nn, min_nn;
	ACnodeBuf *nodebuf;
	ACnode *node;
	int nleaves, nlink;

	tree = pptb_asACtree(pptb);
	nnodes = TREE_SIZE(&tree);
	nodebuf = &(tree.nodebuf);
	Rprintf("| Total nb of nodes = %u\n", nnodes);
	for (nlink = 0; nlink < MAX_CHILDREN_PER_NODE+2; nlink++)
		nlink_table[nlink] = 0U;
	nleaves = 0;
	for (nid = 0U; nid < nnodes; nid++) {
		node = get_node_from_buf(nodebuf, nid);
		nlink = get_ACnode_nlink(&tree, node);
		nlink_table[nlink]++;
		if (IS_LEAFNODE(node))
			nleaves++;
	}
	for (nlink = 0; nlink < MAX_CHILDREN_PER_NODE+2; nlink++)
		Rprintf("| - %u nodes (%.2f%) with %d links\n",
			nlink_table[nlink],
			100.00 * nlink_table[nlink] / nnodes,
			nlink);
	Rprintf("| Nb of leaf nodes (nleaves) = %d\n", nleaves);
	max_nn = count_max_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	min_nn = count_min_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	Rprintf("| - max_needed_nnodes(nleaves, TREE_DEPTH) = %u\n", max_nn);
	Rprintf("| - min_needed_nnodes(nleaves, TREE_DEPTH) = %u\n", min_nn);
	return R_NilValue;
}



/****************************************************************************
 *                             G. PREPROCESSING                             *
 ****************************************************************************/

static void add_pattern(ACtree *tree, const cachedCharSeq *P, int P_offset)
{
	int P_id, depth, dmax, linktag;
	unsigned int nid1, nid2;
	ACnode *node1, *node2;

	P_id = P_offset + 1;
	dmax = TREE_DEPTH(tree) - 1;
	for (depth = 0, nid1 = 0U; depth <= dmax; depth++, nid1 = nid2) {
		node1 = GET_NODE(tree, nid1);
		linktag = CHAR2LINKTAG(tree, P->seq[depth]);
		if (linktag == NA_INTEGER)
			error("non base DNA letter found in Trusted Band "
			      "for pattern %d", P_id);
		nid2 = GET_NODE_LINK(tree, node1, linktag);
		if (depth < dmax) {
			if (nid2 != NOT_AN_ID)
				continue;
			nid2 = NEW_NODE(tree, depth + 1);
			SET_NODE_LINK(tree, node1, linktag, nid2);
			continue;
		}
		if (nid2 != NOT_AN_ID) {
			node2 = GET_NODE(tree, nid2);
			_report_ppdup(P_offset, NODE_P_ID(node2));
		} else {
			nid2 = NEW_LEAFNODE(tree, P_id);
			SET_NODE_LINK(tree, node1, linktag, nid2);
		}
	}
	return;
}

/* --- .Call ENTRY POINT ---
 * Arguments:
 *   tb:         the Trusted Band extracted from the input dictionary as a
 *               rectangular DNAStringSet object;
 *   pp_exclude: NULL or an integer vector of the same length as 'tb' where
 *               non-NA values indicate the elements to exclude from
 *               preprocessing;
 *   base_codes: the internal codes for A, C, G and T.
 */
SEXP ACtree2_build(SEXP tb, SEXP pp_exclude, SEXP base_codes,
		SEXP nodebuf_ptr, SEXP nodeextbuf_ptr)
{
	ACtree tree;
	int tb_length, tb_width, P_offset;
	cachedXStringSet cached_tb;
	cachedCharSeq P;
	SEXP ans, ans_names, ans_elt;

	tb_length = _get_XStringSet_length(tb);
	if (tb_length == 0)
		error("Trusted Band is empty");
	_init_ppdups_buf(tb_length);
	tb_width = -1;
	cached_tb = _cache_XStringSet(tb);
	for (P_offset = 0; P_offset < tb_length; P_offset++) {
		/* skip duplicated patterns */
		if (pp_exclude != R_NilValue
		 && INTEGER(pp_exclude)[P_offset] != NA_INTEGER)
			continue;
		P = _get_cachedXStringSet_elt(&cached_tb, P_offset);
		if (tb_width == -1) {
			if (P.length == 0)
				error("first element in Trusted Band "
				      "is of length 0");
			tb_width = P.length;
			tree = new_ACtree(tb_length, tb_width, base_codes,
					nodebuf_ptr, nodeextbuf_ptr);
		} else if (P.length != tb_width) {
			error("element %d in Trusted Band has a different "
			      "length than first element", P_offset + 1);
		}
		add_pattern(&tree, &P, P_offset);
	}

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("ACtree"));
	SET_STRING_ELT(ans_names, 1, mkChar("high2low"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "ACtree" element */
	SET_ELEMENT(ans, 0, R_NilValue);

	/* set the "high2low" element */
	PROTECT(ans_elt = _get_ppdups_buf_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}



/****************************************************************************
 *                   H. CORE AHO-CORASICK WALKING HELPERS                   *
 ****************************************************************************/

/*
 * The core Aho-Corasick walking helpers use recursion and indirect
 * recursion:
 *   o transition() calls compute_flink() and transition();
 *   o compute_flink() calls transition();
 * A trick is to have the path from the root node to the current node 'node'
 * stored at *negative* indices in 'node_path' i.e., if d is the depth of
 * 'node', then the path is node_path[-d], node_path[-d+1], node_path[-d+2],
 * ..., node_path[-1].
 */
static unsigned int compute_flink(ACtree *tree,
		const ACnode *node, const char *node_path);

/*
 * 'node_path' will only be used to compute failure links so it's safe to not
 * provide it (i.e. NULL) if all the nodes already have one.
 */
static unsigned int transition(ACtree *tree,
		ACnode *node, const char *node_path, int linktag)
{
	unsigned int link, flink;

	if (linktag == NA_INTEGER)
		return 0U;
	link = GET_NODE_LINK(tree, node, linktag);
	if (link != NOT_AN_ID)
		return link;
	if (IS_ROOTNODE(tree, node))
		return 0U;
	flink = GET_NODE_FLINK(tree, node);
	if (flink == NOT_AN_ID) {
		flink = compute_flink(tree, node, node_path);
		SET_NODE_FLINK(tree, node, flink);
	}
	link = transition(tree, GET_NODE(tree, flink), node_path, linktag);
	SET_NODE_LINK(tree, node, linktag, link); /* sets a shortcut */
	return link;
}

static unsigned int compute_flink(ACtree *tree,
		const ACnode *node, const char *node_path)
{
	int seq_len, n, linktag;
	const char *seq, *node1_path;
	unsigned int nid;
	ACnode *node1;

	/* 'seq' is obtained by trimming the first letter from the node path */
	seq_len = NODE_DEPTH(tree, node) - 1;
	seq = node_path - seq_len;
	/* walk on 'seq' */
	nid = 0U;
	for (n = 0, node1_path = seq; n < seq_len; n++, node1_path++) {
		node1 = GET_NODE(tree, nid);
		linktag = CHAR2LINKTAG(tree, *node1_path);
		nid = transition(tree, node1, node1_path, linktag);
	}
	return nid;
}

static int has_all_flinks(ACtree *tree)
{
	unsigned int nnodes, nid, flink;
	const ACnode *node;

	nnodes = TREE_SIZE(tree);
	for (nid = 1U; nid < nnodes; nid++) {
		node = GET_NODE(tree, nid);
		flink = GET_NODE_FLINK(tree, node);
		if (flink == NOT_AN_ID)
			return 0;
	}
	return 1;
}

static void compute_flinks_along_pattern(ACtree *tree, const cachedCharSeq *P)
{
	ACnode *node;
	int n, linktag;
	const char *node_path;
	unsigned int nid, flink;

	node = GET_NODE(tree, 0U);
	node_path = P->seq;
	for (n = 1; n <= P->length; n++) {
		linktag = CHAR2LINKTAG(tree, *node_path);
		nid = transition(tree, node, node_path, linktag);
		node = GET_NODE(tree, nid);
		node_path++;
		flink = GET_NODE_FLINK(tree, node);
		if (flink == NOT_AN_ID) {
			flink = compute_flink(tree, node, node_path);
			SET_NODE_FLINK(tree, node, flink);
		}
	}
	return;
}

static void compute_all_flinks(ACtree *tree, const cachedXStringSet *tb)
{
	unsigned int nnodes, nid;
	ACnode *node;
	int P_offset;
	cachedCharSeq P;

	nnodes = TREE_SIZE(tree);
	for (nid = 1U; nid < nnodes; nid++) {
		node = GET_NODE(tree, nid);
		if (!IS_LEAFNODE(node))
			continue;
		P_offset = NODE_P_ID(node) - 1;
		P = _get_cachedXStringSet_elt(tb, P_offset);
		compute_flinks_along_pattern(tree, &P);
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_has_all_flinks(SEXP pptb)
{
	ACtree tree;

	tree = pptb_asACtree(pptb);
	return ScalarLogical(has_all_flinks(&tree));
}

/* --- .Call ENTRY POINT --- */
SEXP ACtree2_compute_all_flinks(SEXP pptb)
{
	ACtree tree;
	SEXP tb;
	cachedXStringSet cached_tb;

	tree = pptb_asACtree(pptb);
	tb = _get_PreprocessedTB_tb(pptb);
	cached_tb = _cache_XStringSet(tb);
	compute_all_flinks(&tree, &cached_tb);
	return R_NilValue;
}



/****************************************************************************
 *                             I. MATCH FINDING                             *
 ****************************************************************************/

/* Does report matches */
static void walk_tb_subject(ACtree *tree, const cachedCharSeq *S,
		TBMatchBuf *tb_matches)
{
	ACnode *node;
	int n, linktag;
	const char *node_path;
	unsigned int nid;

	node = GET_NODE(tree, 0U);
	node_path = S->seq;
	for (n = 1; n <= S->length; n++) {
		linktag = CHAR2LINKTAG(tree, *node_path);
		nid = transition(tree, node, node_path, linktag);
		node = GET_NODE(tree, nid);
		node_path++;
		if (IS_LEAFNODE(node))
			_TBMatchBuf_report_match(tb_matches,
					NODE_P_ID(node) - 1, n);
	}
	return;
}

/* 1st helper function for walk_tb_nonfixed_subject() */
#define	NODE_SUBSET_MAXSIZE	5000000 /* 5 million node pointers */
static ACnode *node_subset[NODE_SUBSET_MAXSIZE];
static int node_subset_size = 0;

static void split_and_move_pointers(ACtree *tree, unsigned char c)
{
	int node_subset_size0, i, is_first, j, linktag;
	ACnode *node0, *node1, *node2;
	unsigned char base;
	unsigned int nid;

	node0 = GET_NODE(tree, 0U);
	node_subset_size0 = node_subset_size;
	for (i = 0; i < node_subset_size0; i++) {
		node1 = node_subset[i];
		is_first = 1;
		for (j = 0, base = 1; j < 4; j++, base *= 2) {
			if ((c & base) == 0)
				continue;
			linktag = CHAR2LINKTAG(tree, base);
			nid = transition(tree, node1, NULL, linktag);
			//Rprintf("%d --[%d]--> %d\n", node1 - node0, base, nid);
			node2 = GET_NODE(tree, nid);
			if (is_first) {
				node_subset[i] = node2;
				is_first = 0;
			} else {
				if (node_subset_size >= NODE_SUBSET_MAXSIZE) {
					node_subset_size = 0;
					error("too many IUPAC ambiguity "
					      "letters in 'subject'");
				}
				node_subset[node_subset_size++] = node2;
			}
		}
	}
	return;
}

/* 2nd helper function for walk_tb_nonfixed_subject() */
static int compar_node_pointers_for_sort(const void *p1, const void *p2)
{
	return *((const ACnode * const *) p1) - *((const ACnode * const *) p2);
}
static void sort_node_pointer_array(ACnode **x, int nelt)
{
	qsort(x, nelt, sizeof(ACnode *), compar_node_pointers_for_sort);
	return;
}
static void merge_pointers(ACtree *tree, int n)
{
	int i1, i2;
	ACnode *node0, *node1, *node2;

	node0 = GET_NODE(tree, 0U);
/*
	Rprintf("n=%d nodes before merging: ", n);
	for (i1 = 0; i1 < node_subset_size; i1++) {
		node1 = node_subset[i1];
		Rprintf(" %d", node1 - node0);
	}
	Rprintf("\n");
*/
	sort_node_pointer_array(node_subset, node_subset_size);
	i1 = 0;
	node1 = node_subset[i1];
	for (i2 = 1; i2 < node_subset_size; i2++) {
		node2 = node_subset[i2];
		if (node2 == node1)
			continue;
		i1++;
		node1 = node_subset[i1] = node2;
	}
	node_subset_size = i1 + 1;
/*
	Rprintf("n=%d nodes after merging: ", n);
	for (i1 = 0; i1 < node_subset_size; i1++) {
		node1 = node_subset[i1];
		Rprintf(" %d", node1 - node0);
	}
	Rprintf("\n");
*/
	return;
}

static void report_matches(TBMatchBuf *tb_matches, int n)
{
	int i;
	ACnode *node;

	for (i = 0; i < node_subset_size; i++) {
		node = node_subset[i];
		if (IS_LEAFNODE(node))
			_TBMatchBuf_report_match(tb_matches,
					NODE_P_ID(node) - 1, n);
	}
	return;
}

/* Does report matches */
static void walk_tb_nonfixed_subject(ACtree *tree, const cachedCharSeq *S,
		TBMatchBuf *tb_matches)
{
	int max_size, n;
	const unsigned char *c;

	if (node_subset_size != 0)
		error("Biostrings internal error in "
		      "walk_tb_nonfixed_subject(): node_subset_size != 0... "
		      "PLEASE REPORT THIS! THANKS.\n");
	node_subset_size = max_size = 1;
	node_subset[0] = GET_NODE(tree, 0U);
	for (n = 1, c = (unsigned char *) S->seq; n <= S->length; n++, c++) {
		if (*c >= 16) {
			/* '*c' is not an IUPAC (base or extended) code */
			node_subset[0] = GET_NODE(tree, 0U);
			node_subset_size = 1;
			continue;
		}
		split_and_move_pointers(tree, *c);
		merge_pointers(tree, n);
/*
		if (node_subset_size > max_size) {
			max_size = node_subset_size;
			Rprintf("walk_tb_nonfixed_subject(): "
				"n=%d max_size=%d\n",
				n, max_size);
		}
*/
		report_matches(tb_matches, n);
	}
	node_subset_size = 0;
	return;
}

/* Entry point for the MATCH FINDING section */
void _match_tbACtree2(SEXP pptb, const cachedCharSeq *S, int fixedS,
		TBMatchBuf *tb_matches)
{
	ACtree tree;
	SEXP tb;
	cachedXStringSet cached_tb;

	tree = pptb_asACtree(pptb);
	if (fixedS) {
		walk_tb_subject(&tree, S, tb_matches);
		return;
	}
	if (!has_all_flinks(&tree)) {
		tb = _get_PreprocessedTB_tb(pptb);
		cached_tb = _cache_XStringSet(tb);
		//Rprintf("computing all flinks... ");
		compute_all_flinks(&tree, &cached_tb);
		//Rprintf("OK\n");
	}
	walk_tb_nonfixed_subject(&tree, S, tb_matches);
	return;
}



/****************************************************************************
 *                          J. MORE MATCH FINDING                           *
 ****************************************************************************/

/* Does report matches */
static void walk_pdict_subject(ACtree *tree,
		SEXP low2high, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_nmis, int min_nmis, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	ACnode *node;
	int n, linktag;
	unsigned int nid;
	const char *node_path;

	node = GET_NODE(tree, 0U);
	node_path = S->seq;
	for (n = 1; n <= S->length; n++) {
		linktag = CHAR2LINKTAG(tree, *node_path);
		nid = transition(tree, node, node_path, linktag);
		node = GET_NODE(tree, nid);
		node_path++;
		if (IS_LEAFNODE(node))
			_match_pdict_flanks_at(NODE_P_ID(node) - 1,
				low2high, headtail, S, n,
				max_nmis, min_nmis, fixedP,
				matchpdict_buf);
	}
	return;
}

/* Does report matches */
static void walk_pdict_nonfixed_subject(ACtree *tree,
		SEXP low2high, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_nmis, int min_nmis, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	error("walk_pdict_nonfixed_subject(): implement me");
	return;
}

void _match_pdictACtree2(SEXP pptb, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_nmis, int min_nmis, int fixedP, int fixedS,
		MatchPDictBuf *matchpdict_buf)
{
	ACtree tree;
	SEXP low2high;

	tree = pptb_asACtree(pptb);
	low2high = _get_PreprocessedTB_low2high(pptb);
	if (fixedS)
		walk_pdict_subject(&tree,
			low2high, headtail, S,
			max_nmis, min_nmis, fixedP,
			matchpdict_buf);
	else
		walk_pdict_nonfixed_subject(&tree,
			low2high, headtail, S,
			max_nmis, min_nmis, fixedP,
			matchpdict_buf);
	return;
}

