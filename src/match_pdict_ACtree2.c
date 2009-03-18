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

#define MAX_CHILDREN_PER_NODE 4  /* do NOT change this */



/****************************************************************************
 *    A. GENERAL UTILITIES FOR COMPUTING THE NB OF NODES AND EXTENSIONS     *
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

static unsigned int count_max_needed_nextensions_at_pp_time(int nleaves, int depth)
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
 * and it looked correct... pfff :-b */
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
			n2 = count_max_needed_nextensions_at_pp_time(nleaves, depth);
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

/*
 * OptMaxNN: Optimistic Max Needed nb of Nodes.
 * get_OptMaxNN() is a decreasing function of 'max_needed_nnodes' bounded by
 * 'max_nn' and 'min_nn2' (max_nn >= get_OptMaxNN() > min_nn2).
 */
static unsigned int get_OptMaxNN(unsigned int max_needed_nnodes,
		unsigned int min_needed_nnodes)
{
	double max_nn, min_nn2, x;

	max_nn = max_needed_nnodes;
	min_nn2 = min_needed_nnodes + 0.70 * (max_nn - min_needed_nnodes);
	x = max_needed_nnodes / 300000000.00;
	return (unsigned int) (min_nn2 + (max_nn - min_nn2) / (1.00 + x));
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
 * Note that when 'extbuf_length == 0.4 * nodebuf_nelt', the 'extensions'
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
 *                         B. ACnode AND ACnodeBuf                          *
 ****************************************************************************/

#define LINKTAG_BITSHIFT 28
#define MAX_DEPTH ((1 << LINKTAG_BITSHIFT) - 1)
#define ISLEAF_BIT (1 << 30)
#define ISEXTENDED_BIT (ISLEAF_BIT << 1) /* strongest bit for 32-bit integers */
#define MAX_P_ID (ISLEAF_BIT - 1)  /* P_id values are encoded on 30 bits */

#define ISEXTENDED(node) ((node)->attribs & ISEXTENDED_BIT)
#define ISLEAF(node) ((node)->attribs & ISLEAF_BIT)
#define _NODE_DEPTH(node) ((node)->attribs & MAX_DEPTH)
/* result of NODE_P_ID() is undefined on a non-leaf node */
#define NODE_P_ID(node) ((node)->attribs & MAX_P_ID)

#define NULLNODE UINT_MAX

/*
 * A node id (nid) is represented by an unsigned int. The ANSI C standard only
 * guarantees that an (unsigned) int will take at least 2 bytes in memory i.e.
 * that UINT_MAX will always be >= 2^16-1 even though, on most modern
 * platforms (including PC and Mac), an int will take at least 4 bytes in
 * memory i.e. UINT_MAX will be 2^32-1 or more.
 * The same apply to extension ids (eid).
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

#define ISROOT(node, nodebuf) ((node) == (nodebuf)->block[0])

SEXP ACtree2_nodebuf_max_nblock()
{
	return ScalarInteger(ACNODEBUF_MAX_NBLOCK);
}

static ACnode *get_node_from_buf(ACnodeBuf *nodebuf, unsigned int nid)
{
	unsigned int b, i;

	b = nid >> 22U;
	i = nid & (ACNODEBUF_MAX_NELT_PER_BLOCK - 1U);
	return nodebuf->block[b] + i;
}

static unsigned int get_ACnodeBuf_length(ACnodeBuf *nodebuf)
{
	return (unsigned int) *(nodebuf->nblock) * ACNODEBUF_MAX_NELT_PER_BLOCK;
}

static unsigned int get_ACnodeBuf_nelt(ACnodeBuf *nodebuf)
{
	int nblock;

	nblock = *(nodebuf->nblock);
	if (nblock == 0)
		return 0U;
	return (unsigned int) (nblock - 1) * ACNODEBUF_MAX_NELT_PER_BLOCK
	       + *(nodebuf->lastblock_nelt);
}

static ACnodeBuf new_ACnodeBuf(SEXP bab)
{
	ACnodeBuf nodebuf;
	SEXP bab_blocks;
	int nblock, b;

	nodebuf.bab = bab;
	nblock = *(nodebuf.nblock = _get_BAB_nblock_ptr(bab));
	nodebuf.lastblock_nelt = _get_BAB_lastblock_nelt_ptr(bab);
	bab_blocks = _get_BAB_blocks(bab);
	for (b = 0; b < nblock; b++)
		nodebuf.block[b] = (ACnode *) INTEGER(VECTOR_ELT(bab_blocks, b));
	return nodebuf;
}

static void extend_nodebuf(ACnodeBuf *nodebuf)
{
	int length;
	SEXP bab_block;

	length = ACNODEBUF_MAX_NELT_PER_BLOCK * INTS_PER_NODE;
	bab_block = _IntegerBAB_addblock(nodebuf->bab, length);
	/* sync 'nodebuf->block' with 'nodebuf->bab' */
	nodebuf->block[*(nodebuf->nblock) - 1] = (ACnode *) INTEGER(bab_block);
	return;
}



/****************************************************************************
 *                       C. ACnodeExt AND ACnodeExtBuf                      *
 ****************************************************************************/

typedef struct acnodeext {
	unsigned int link_nid[MAX_CHILDREN_PER_NODE];
	unsigned int flink_nid;
} ACnodeExt;

#define INTS_PER_EXTENSION (sizeof(ACnodeExt) / sizeof(int))

/*
 * We must have:
 *   (a) ACNODEEXTBUF_MAX_NBLOCK * ACNODEEXTBUF_MAX_NELT_PER_BLOCK <= UINT_MAX + 1
 *   (b) ACNODEEXTBUF_MAX_NELT_PER_BLOCK * INTS_PER_EXTENSION <= INT_MAX
 * The following settings are assuming UINT_MAX >= 2^32-1 and
 * INT_MAX >= 2^31-1. They result in blocks of size 80 MB.
 */
#define ACNODEEXTBUF_MAX_NBLOCK 1024  /* = 2^10 */
#define ACNODEEXTBUF_MAX_NELT_PER_BLOCK 4194304U  /* = 2^22 */

typedef struct acnodeextbuf {
	SEXP bab;  /* Big Atomic Buffer */
	int *nblock;
	int *lastblock_nelt;
	ACnodeExt *block[ACNODEEXTBUF_MAX_NBLOCK];
} ACnodeExtBuf;

SEXP ACtree2_extbuf_max_nblock()
{
	return ScalarInteger(ACNODEEXTBUF_MAX_NBLOCK);
}

static ACnodeExt *get_extension_from_buf(ACnodeExtBuf *extbuf, unsigned int eid)
{
	unsigned int b, i;

	b = eid >> 22U;
	i = eid & (ACNODEEXTBUF_MAX_NELT_PER_BLOCK - 1U);
	return extbuf->block[b] + i;
}

static unsigned int get_ACnodeExtBuf_length(ACnodeExtBuf *extbuf)
{
	return (unsigned int) *(extbuf->nblock) * ACNODEEXTBUF_MAX_NELT_PER_BLOCK;
}

static unsigned int get_ACnodeExtBuf_nelt(ACnodeExtBuf *extbuf)
{
	int nblock;

	nblock = *(extbuf->nblock);
	if (nblock == 0)
		return 0U;
	return (unsigned int) (nblock - 1) * ACNODEEXTBUF_MAX_NELT_PER_BLOCK
	       + *(extbuf->lastblock_nelt);
}

static ACnodeExtBuf new_ACnodeExtBuf(SEXP bab)
{
	ACnodeExtBuf extbuf;
	SEXP bab_blocks;
	int nblock, b;

	extbuf.bab = bab;
	nblock = *(extbuf.nblock = _get_BAB_nblock_ptr(bab));
	extbuf.lastblock_nelt = _get_BAB_lastblock_nelt_ptr(bab);
	bab_blocks = _get_BAB_blocks(bab);
	for (b = 0; b < nblock; b++)
		extbuf.block[b] = (ACnodeExt *) INTEGER(VECTOR_ELT(bab_blocks, b));
	return extbuf;
}

static void extend_extbuf(ACnodeExtBuf *extbuf)
{
	int length;
	SEXP bab_block;

	length = ACNODEEXTBUF_MAX_NELT_PER_BLOCK * INTS_PER_EXTENSION;
	bab_block = _IntegerBAB_addblock(extbuf->bab, length);
	/* sync 'extbuf->block' with 'extbuf->bab' */
	extbuf->block[*(extbuf->nblock) - 1] = (ACnodeExt *) INTEGER(bab_block);
	return;
}



/****************************************************************************
 *                                 D. ACtree                                *
 ****************************************************************************/

SEXP debug_match_pdict_ACtree2()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
	if (debug) {
		Rprintf("[DEBUG] debug_match_pdict_ACtree2():\n");
/*
		Rprintf("  INTS_PER_NODE=%d MAX_NNODES=%d\n"
			"  INTS_PER_EXTENSION=%d MAX_NEXTENSIONS=%d\n",
			INTS_PER_NODE, MAX_NNODES,
			INTS_PER_EXTENSION, MAX_NEXTENSIONS);
*/
		Rprintf("  INTS_PER_NODE=%d INTS_PER_EXTENSION=%d\n",
			INTS_PER_NODE, INTS_PER_EXTENSION);
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

typedef struct actree {
	int depth;  /* this is the depth of all leaf nodes */
	ACnodeBuf nodebuf;
	ACnodeExtBuf extbuf;
	ByteTrTable char2linktag;
} ACtree;

#define TREE_DEPTH(tree) ((tree)->depth)
#define NODE_DEPTH(tree, node) (ISLEAF(node) ? TREE_DEPTH(tree) : _NODE_DEPTH(node))

#define NODEBUF_LENGTH(tree) get_ACnodeBuf_length(&((tree)->nodebuf))
#define NODEBUF_NELT(tree) get_ACnodeBuf_nelt(&((tree)->nodebuf))
#define GET_NODE(tree, nid) get_node_from_buf(&((tree)->nodebuf), nid)

#define EXTBUF_LENGTH(tree) get_ACnodeExtBuf_length(&((tree)->extbuf))
#define EXTBUF_NELT(tree) get_ACnodeExtBuf_nelt(&((tree)->extbuf))
#define GET_EXTENSION(tree, eid) get_extension_from_buf(&((tree)->extbuf), eid)

#define CHAR2LINKTAG(tree, c) ((tree)->char2linktag[(unsigned char) (c)])



/****************************************************************************
 *                          E. LOW-LEVEL UTILITIES                          *
 ****************************************************************************/

static unsigned int new_nid(ACtree *tree)
{
	unsigned int nid;

	nid = NODEBUF_NELT(tree);
	if (nid >= NODEBUF_LENGTH(tree))
		extend_nodebuf(&(tree->nodebuf));
	(*(tree->nodebuf.lastblock_nelt))++;
	return nid;
}

static unsigned int new_eid(ACtree *tree)
{
	unsigned int eid;

	eid = EXTBUF_NELT(tree);
	if (eid >= EXTBUF_LENGTH(tree))
		extend_extbuf(&(tree->extbuf));
	(*(tree->extbuf.lastblock_nelt))++;
	return eid;
}

static void extend_ACnode(ACtree *tree, ACnode *node)
{
	ACnodeExt *extension;
	unsigned int eid;
	int i, linktag;

	eid = new_eid(tree);
	extension = GET_EXTENSION(tree, eid);
	for (i = 0; i < MAX_CHILDREN_PER_NODE; i++)
		extension->link_nid[i] = NULLNODE;
	extension->flink_nid = NULLNODE;
	if (node->nid_or_eid != NULLNODE) {
		/* this is correct because 'node' cannot be a leaf node */
		linktag = node->attribs >> LINKTAG_BITSHIFT;
		extension->link_nid[linktag] = node->nid_or_eid;
	}
	node->nid_or_eid = eid;
	/* sets the "ISEXTENDED" bit to 1 */
	node->attribs |= ISEXTENDED_BIT;
	return;
}



/****************************************************************************
 *                           F. Aho-Corasick tree API                       *
 ****************************************************************************/

static unsigned int new_ACnode(ACtree *tree, int depth)
{
	unsigned int nid;
	ACnode *node;

	if (depth >= TREE_DEPTH(tree))
		error("new_ACnode(): depth >= TREE_DEPTH(tree)");
	nid = new_nid(tree);
	node = GET_NODE(tree, nid);
	/* this sets the "ISEXTENDED" and "ISLEAF" bits to 0 */
	node->attribs = depth;
	node->nid_or_eid = NULLNODE;
	return nid;
}

static unsigned int new_leafACnode(ACtree *tree, int P_id)
{
	unsigned int nid;
	ACnode *node;

	nid = new_nid(tree);
	node = GET_NODE(tree, nid);
	/* this sets the "ISEXTENDED" bit to 0 and "ISLEAF" bit to 1 */
	node->attribs = ISLEAF_BIT | P_id;
	node->nid_or_eid = NULLNODE;
	return nid;
}

static unsigned int get_ACnode_link(ACtree *tree, ACnode *node, int linktag)
{
	ACnodeExt *extension;

	if (node->nid_or_eid == NULLNODE)
		return NULLNODE;
	if (ISEXTENDED(node)) {
		extension = GET_EXTENSION(tree, node->nid_or_eid);
		return extension->link_nid[linktag];
	}
	/* the node has no extension and is not a leaf node */
	if (linktag == (node->attribs >> LINKTAG_BITSHIFT))
		return node->nid_or_eid;
	return NULLNODE;
}

/*
 * set_ACnode_flink() should always have been called before set_ACnode_link()
 * on a leaf node. So we can assume that, if set_ACnode_link() is called on a
 * leaf node, then this node is already extended.
 */
static void set_ACnode_link(ACtree *tree, ACnode *node, int linktag, unsigned int nid)
{
	ACnodeExt *extension;

	if (node->nid_or_eid == NULLNODE) {
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

static unsigned int get_ACnode_flink(ACtree *tree, ACnode *node)
{
	ACnodeExt *extension;

	if (!ISEXTENDED(node))
		return NULLNODE;
	extension = GET_EXTENSION(tree, node->nid_or_eid);
	return extension->flink_nid;
}

static void set_ACnode_flink(ACtree *tree, ACnode *node, unsigned int nid)
{
	ACnodeExt *extension;

	if (!ISEXTENDED(node))
		extend_ACnode(tree, node);
	extension = GET_EXTENSION(tree, node->nid_or_eid);
	extension->flink_nid = nid;
	return;
}

static ACtree new_ACtree(int tb_length, int tb_width, SEXP base_codes,
		SEXP nodebuf_ptr, SEXP extbuf_ptr)
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
	tree.extbuf = new_ACnodeExtBuf(extbuf_ptr);
	_init_byte2offset_with_INTEGER(tree.char2linktag, base_codes, 1);
	new_ACnode(&tree, 0);  /* create the root node */
	return tree;
}

static ACtree pptb_asACtree(SEXP pptb)
{
	ACtree tree;
	SEXP base_codes;

	tree.depth = _get_PreprocessedTB_width(pptb);
	tree.nodebuf = new_ACnodeBuf(_get_ACtree2_nodebuf_ptr(pptb));
	tree.extbuf = new_ACnodeExtBuf(_get_ACtree2_extbuf_ptr(pptb));
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

static int get_ACnode_nlink(ACtree *tree, ACnode *node)
{
	int nlink, linktag;

	nlink = get_ACnode_flink(tree, node) != NULLNODE;
	for (linktag = 0; linktag < MAX_CHILDREN_PER_NODE; linktag++)
		if (get_ACnode_link(tree, node, linktag) != NULLNODE)
			nlink++;
	return nlink;
}



/****************************************************************************
 *                       G. STATS AND DEBUG UTILITIES                       *
 ****************************************************************************/

SEXP ACtree2_print_nodes(SEXP pptb)
{
	ACtree tree;
	ACnode *node;
	unsigned int nnodes, nid;

	tree = pptb_asACtree(pptb);
	nnodes = NODEBUF_NELT(&tree);
	for (nid = 0U; nid < nnodes; nid++) {
		node = GET_NODE(&tree, nid);
		print_ACnode(&tree, node);
	}
	return R_NilValue;
}

SEXP ACtree2_summary(SEXP pptb)
{
	ACtree tree;
	ACnode *node;
	unsigned int nnodes, nlink_table[MAX_CHILDREN_PER_NODE+2],
		     nid, max_nn, min_nn, n1;
	int nleaves, nlink;

	tree = pptb_asACtree(pptb);
	nnodes = NODEBUF_NELT(&tree);
	Rprintf("  Total nb of nodes = %u\n", nnodes);
	for (nlink = 0; nlink < MAX_CHILDREN_PER_NODE+2; nlink++)
		nlink_table[nlink] = 0U;
	nleaves = 0;
	for (nid = 0U; nid < nnodes; nid++) {
		node = GET_NODE(&tree, nid);
		nlink = get_ACnode_nlink(&tree, node);
		nlink_table[nlink]++;
		if (ISLEAF(node))
			nleaves++;
	}
	for (nlink = 0; nlink < MAX_CHILDREN_PER_NODE+2; nlink++)
		Rprintf("  - %u nodes (%.2f%) with %d links\n",
			nlink_table[nlink],
			100.00 * nlink_table[nlink] / nnodes,
			nlink);
	Rprintf("  Nb of leaf nodes (nleaves) = %d\n", nleaves);
	max_nn = count_max_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	min_nn = count_min_needed_nnodes(nleaves, TREE_DEPTH(&tree));
	n1 = get_OptMaxNN(max_nn, min_nn);
	Rprintf("  - max_needed_nnodes(nleaves, TREE_DEPTH) = %u\n", max_nn);
	Rprintf("  - min_needed_nnodes(nleaves, TREE_DEPTH) = %u\n", min_nn);
	Rprintf("  - OptMaxNN(nleaves, TREE_DEPTH) = %u\n", n1);
	return R_NilValue;
}



/****************************************************************************
 *                             H. PREPROCESSING                             *
 ****************************************************************************/

static void add_pattern(ACtree *tree, const RoSeq *P, int P_offset)
{
	int P_id, depth, dmax, linktag;
	unsigned int nid1, nid2;
	ACnode *node1, *node2;

	P_id = P_offset + 1;
	dmax = TREE_DEPTH(tree) - 1;
	for (depth = 0, nid1 = 0U; depth <= dmax; depth++, nid1 = nid2) {
		node1 = GET_NODE(tree, nid1);
		linktag = CHAR2LINKTAG(tree, P->elts[depth]);
		if (linktag == NA_INTEGER)
			error("non base DNA letter found in Trusted Band "
			      "for pattern %d", P_id);
		nid2 = get_ACnode_link(tree, node1, linktag);
		if (depth < dmax) {
			if (nid2 != NULLNODE)
				continue;
			nid2 = new_ACnode(tree, depth + 1);
			set_ACnode_link(tree, node1, linktag, nid2);
			continue;
		}
		if (nid2 != NULLNODE) {
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
 */

SEXP ACtree2_build(SEXP tb, SEXP dup2unq0, SEXP base_codes,
		SEXP nodebuf_ptr, SEXP extbuf_ptr)
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
			tree = new_ACtree(tb_length, tb_width, base_codes,
					nodebuf_ptr, extbuf_ptr);
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
	SET_ELEMENT(ans, 0, R_NilValue);

	/* set the "dup2unq" element */
	PROTECT(ans_elt = _dup2unq_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}



/****************************************************************************
 *                             I. MATCH FINDING                             *
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
		Rprintf("nid=%u node_depth=%d linktag=%d path=%s\n",
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
			Rprintf("link=%u\n", 0U);
		}
		rec_level--;
#endif
*/
		return 0U;
	}
	link = get_ACnode_link(tree, node, linktag);
	if (link != NULLNODE) {
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
	if (ISROOT(node, &(tree->nodebuf))) {
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
	flink = get_ACnode_flink(tree, node);
	if (flink == NULLNODE) {
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
static void walk_subject(ACtree *tree, const RoSeq *S)
{
	ACnode *node;
	int n, linktag;
	unsigned int nid;
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

