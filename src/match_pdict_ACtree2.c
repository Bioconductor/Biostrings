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

#include <stdlib.h> /* for div() */
#include <limits.h> /* for UINT_MAX */

static int debug = 0;

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

SEXP ACtree2_nodebuf_max_nblock()
{
	return ScalarInteger(ACNODEBUF_MAX_NBLOCK);
}

static int ACnodeBuf_isfull(ACnodeBuf *buf)
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

	if (ACnodeBuf_isfull(buf))
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
 *   (a) ACNODEEXTBUF_MAX_NBLOCK * ACNODEEXTBUF_MAX_NELT_PER_BLOCK <= UINT_MAX + 1
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
			"As a consequence not all new links and failure links will be\n"
			"set. This might (slightly) affect speed but not the results.",
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

#define IS_ROOTNODE(tree, node) _IS_ROOTNODE(&((tree)->nodebuf), node)
#define IS_LEAFNODE(node) ((node)->attribs & ISLEAF_BIT)
#define TREE_DEPTH(tree) ((tree)->depth)
#define NODE_DEPTH(tree, node) (IS_LEAFNODE(node) ? TREE_DEPTH(tree) : _NODE_DEPTH(node))
#define GET_NODE(tree, nid) get_node_from_buf(&((tree)->nodebuf), nid)
#define CHAR2LINKTAG(tree, c) ((tree)->char2linktag[(unsigned char) (c)])
#define NEW_NODE(tree, depth) new_ACnode(tree, depth)
#define NEW_LEAFNODE(tree, P_id) new_leafACnode(tree, P_id)
#define GET_NODE_LINK(tree, node, linktag) get_ACnode_link(tree, node, linktag)
#define SET_NODE_LINK(tree, node, linktag, nid) set_ACnode_link(tree, node, linktag, nid)
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

static unsigned int get_ACnode_flink(ACtree *tree, ACnode *node)
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
	   in memory of the 2 buffers (nodebuf and nodeextbuf) is about the same */
	n = (unsigned int) (0.40 * nnodes);
	/* then we round up to the closer multiple of ACNODEBUF_MAX_NELT_PER_BLOCK
           so we don't waste space in the last block */
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
	//max_nelt = a_nice_max_nodeextbuf_nelt(get_ACnodeBuf_nelt(&(tree.nodebuf)));
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

SEXP ACtree2_print_nodes(SEXP pptb)
{
	ACtree tree;
	ACnodeBuf *nodebuf;
	ACnode *node;
	unsigned int nnodes, nid;

	tree = pptb_asACtree(pptb);
	nodebuf = &(tree.nodebuf);
	nnodes = get_ACnodeBuf_nelt(nodebuf);
	for (nid = 0U; nid < nnodes; nid++) {
		node = get_node_from_buf(nodebuf, nid);
		print_ACnode(&tree, node);
	}
	return R_NilValue;
}

SEXP ACtree2_summary(SEXP pptb)
{
	ACtree tree;
	ACnodeBuf *nodebuf;
	ACnode *node;
	unsigned int nnodes, nlink_table[MAX_CHILDREN_PER_NODE+2],
		     nid, max_nn, min_nn;
	int nleaves, nlink;

	tree = pptb_asACtree(pptb);
	nodebuf = &(tree.nodebuf);
	nnodes = get_ACnodeBuf_nelt(nodebuf);
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

/****************************************************************************
 * .Call entry point for preprocessing
 * -----------------------------------
 *
 * Arguments:
 *   tb:         the Trusted Band extracted from the original dictionary as a
 *               constant width DNAStringSet object;
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
 *                             H. MATCH FINDING                             *
 ****************************************************************************/

/*
 * We use indirect recursion for walking a sequence.
 * This indirect recursion involves the 2 core functions walk_shortseq() and
 * transition(): the former calls the latter which in turn calls the former.
 */
static unsigned int walk_shortseq(ACtree *tree, const char *seq, int seq_len);

/*
 * A trick is to have the path from the root node to the current node 'node'
 * stored at *negative* indices in 'seq_tail' i.e., if d is the depth of
 * 'node', then the path is seq_tail[-d], seq_tail[-d+1], seq_tail[-d+2], ...,
 * seq_tail[-1].
 */
static unsigned int transition(ACtree *tree, ACnode *node, int linktag, const char *seq_tail)
{
	unsigned int link, flink;
	int newpath_len;
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
		Rprintf("node_depth=%d linktag=%d path=%s\n",
			node_depth, linktag, pathbuf);
	}
#endif
*/

	if (linktag == NA_INTEGER) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%u\n", 0U);
		}
		rec_level--;
#endif
*/
		return 0U;
	}
	link = GET_NODE_LINK(tree, node, linktag);
	if (link != NOT_AN_ID) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%u\n", link);
		}
		rec_level--;
#endif
*/
		return link;
	}
	if (IS_ROOTNODE(tree, node)) {
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG]  LEAVING transition():");
			Rprintf(format, " ");
			Rprintf("link=%u\n", 0U);
		}
		rec_level--;
#endif
*/
		return 0U;
	}
	flink = GET_NODE_FLINK(tree, node);
	if (flink == NOT_AN_ID) {
		newpath_len = NODE_DEPTH(tree, node) - 1;
		newpath = seq_tail - newpath_len;
		flink = walk_shortseq(tree, newpath, newpath_len);
		SET_NODE_FLINK(tree, node, flink);
	}
	link = transition(tree, GET_NODE(tree, flink), linktag, seq_tail);
	SET_NODE_LINK(tree, node, linktag, link); /* sets a shortcut */
/*
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG]  LEAVING transition():");
		Rprintf(format, " ");
		Rprintf("link=%u\n", link);
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
static unsigned int walk_shortseq(ACtree *tree, const char *seq, int seq_len)
{
	ACnode *node;
	int n, linktag;
	unsigned int nid;
	const char *seq_tail;

	nid = 0U;
	for (n = 0, seq_tail = seq; n < seq_len; n++, seq_tail++) {
		node = GET_NODE(tree, nid);
		linktag = CHAR2LINKTAG(tree, *seq_tail);
		nid = transition(tree, node, linktag, seq_tail);
	}
	return nid;
}

/* Does report matches */
static void walk_tb_subject(ACtree *tree, const cachedCharSeq *S,
		TBMatchBuf *tb_matches)
{
	ACnode *node;
	int n, linktag;
	unsigned int nid;
	const char *S_tail;

	node = GET_NODE(tree, 0U);
	for (n = 1, S_tail = S->seq; n <= S->length; n++, S_tail++) {
		linktag = CHAR2LINKTAG(tree, *S_tail);
		nid = transition(tree, node, linktag, S_tail);
		node = GET_NODE(tree, nid);
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
	error("walk_tb_nonfixed_subject(): implement me");
	return;
}

/* Entry point for the MATCH FINDING section */
void _match_tbACtree2(SEXP pptb, const cachedCharSeq *S, int fixedS,
		TBMatchBuf *tb_matches)
{
	ACtree tree;

	tree = pptb_asACtree(pptb);
	if (fixedS)
		walk_tb_subject(&tree, S, tb_matches);
	else
		walk_tb_nonfixed_subject(&tree, S, tb_matches);
	return;
}



/****************************************************************************
 *                          I. MORE MATCH FINDING                           *
 ****************************************************************************/

/* Does report matches */
static void walk_pdict_subject(ACtree *tree,
		SEXP low2high, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_mis, int min_mis, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	ACnode *node;
	int n, linktag;
	unsigned int nid;
	const char *S_tail;

	node = GET_NODE(tree, 0U);
	for (n = 1, S_tail = S->seq; n <= S->length; n++, S_tail++) {
		linktag = CHAR2LINKTAG(tree, *S_tail);
		nid = transition(tree, node, linktag, S_tail);
		node = GET_NODE(tree, nid);
		if (IS_LEAFNODE(node))
			_match_pdict_flanks_at(NODE_P_ID(node) - 1,
				low2high, headtail, S, n,
				max_mis, min_mis, fixedP,
				matchpdict_buf);
	}
	return;
}

/* Does report matches */
static void walk_pdict_nonfixed_subject(ACtree *tree,
		SEXP low2high, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_mis, int min_mis, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	error("walk_pdict_nonfixed_subject(): implement me");
	return;
}

void _match_pdictACtree2(SEXP pptb, HeadTail *headtail,
		const cachedCharSeq *S,
		int max_mis, int min_mis, int fixedP, int fixedS,
		MatchPDictBuf *matchpdict_buf)
{
	ACtree tree;
	SEXP low2high;

	tree = pptb_asACtree(pptb);
	low2high = _get_PreprocessedTB_low2high(pptb);
	if (fixedS)
		walk_pdict_subject(&tree,
			low2high, headtail, S,
			max_mis, min_mis, fixedP,
			matchpdict_buf);
	else
		walk_pdict_nonfixed_subject(&tree,
			low2high, headtail, S,
			max_mis, min_mis, fixedP,
			matchpdict_buf);
	return;
}

