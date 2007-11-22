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


/****************************************************************************/
static int debug = 0;

SEXP match_ACuldna_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_ACuldna.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_ACuldna.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Manipulation of the buffer of duplicates
 */

typedef struct dupsbuf_line {
	int *vals;
	int countmax;
	int count;
} DupsBufLine;

static DupsBufLine *dupsbuf;
static int dupsbuf_countmax, dupsbuf_count;

static void dupsbuf_reset()
{
	/* No memory leak here, because we use transient storage allocation */
	dupsbuf = NULL;
	dupsbuf_countmax = dupsbuf_count = 0;
	return;
}

static int dupsbuf_appendtoline(DupsBufLine *line, int P_id)
{
	long new_countmax;

	if (line->count >= line->countmax) {
		if (line->countmax == 0)
			new_countmax = 1000;
		else
			new_countmax = 2 * line->countmax;
		line->vals = Srealloc((char *) line->vals, new_countmax,
				(long) line->countmax, int);
		line->countmax = new_countmax;
	}
	line->vals[line->count++] = P_id;
	return line->count;
}

static int dupsbuf_append(int P_id1, int P_id2)
{
	int i;
	DupsBufLine *line;
	long new_countmax;

	for (i = 0, line = dupsbuf; i < dupsbuf_count; i++, line++) {
		if (line->vals[0] == P_id1)
			return dupsbuf_appendtoline(line, P_id2);
	}
	if (dupsbuf_count >= dupsbuf_countmax) {
		/* Buffer is full */
		if (dupsbuf_countmax == 0)
			new_countmax = 1000;
		else
			new_countmax = 2 * dupsbuf_countmax;
		dupsbuf = Srealloc((char *) dupsbuf, new_countmax,
					(long) dupsbuf_countmax, DupsBufLine);
		dupsbuf_countmax = new_countmax;
	}
	line = dupsbuf + dupsbuf_count++;
	line->countmax = line->count = 0;
	line->vals = NULL;
	dupsbuf_appendtoline(line, P_id1);
	return dupsbuf_appendtoline(line, P_id2);
}


/****************************************************************************
 * Initialization of the Aho-Corasick 4-ary tree
 * =============================================
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

typedef struct ac_node {
        int ac_node1_id;
        int ac_node2_id;
        int ac_node3_id;
        int ac_node4_id;
	int flink;
        int P_id;
} ACNode;

static ACNode *ACnodebuf;
static int ACnodebuf_countmax, ACnodebuf_count;
static int ACnodebuf_pattern_length;
static int ACnodebuf_base_codes[4];

void ACnodebuf_reset()
{
	int i;

	/* No memory leak here, because we use transient storage allocation */
	ACnodebuf = NULL;
	ACnodebuf_countmax = ACnodebuf_count = 0;
	ACnodebuf_pattern_length = -1;
	for (i = 0; i < 4; i++)
		ACnodebuf_base_codes[i] = -1;
	return;
}

static int ACnodebuf_newNode()
{
	long new_countmax;
	ACNode *node;

	if (ACnodebuf_count >= ACnodebuf_countmax) {
		/* Buffer is full */
		if (ACnodebuf_countmax == 0)
			new_countmax = 50000;
		else
			new_countmax = 2 * ACnodebuf_countmax;
		ACnodebuf = Srealloc((char *) ACnodebuf, new_countmax,
					(long) ACnodebuf_countmax, ACNode);
		ACnodebuf_countmax = new_countmax;
	}
	node = ACnodebuf + ACnodebuf_count;
	node->ac_node1_id = -1;
	node->ac_node2_id = -1;
	node->ac_node3_id = -1;
	node->ac_node4_id = -1;
	node->flink = -1;
	node->P_id = -1;
	return ACnodebuf_count++;
}

static ACNode *ACnodebuf_tryMovingToChild(char c, int childslot, int *ac_node_id)
{
	if (ACnodebuf_base_codes[childslot] == -1) {
		ACnodebuf_base_codes[childslot] = c;
	} else {
		if (c != ACnodebuf_base_codes[childslot])
			return NULL;
	}
	if (*ac_node_id == -1)
		*ac_node_id = ACnodebuf_newNode();
	return ACnodebuf + *ac_node_id;
}

static void ACnodebuf_addPattern(const char *P, int nP, int P_id)
{
	ACNode *node, *child_node;
	int n;

	if (ACnodebuf_pattern_length == -1) {
		if (nP == 0)
			error("dictionary contains empty patterns");
		ACnodebuf_pattern_length = nP;
	} else {
		if (nP != ACnodebuf_pattern_length)
			error("all patterns in dictionary must have the same length");
	}
	node = ACnodebuf;
	for (n = 0; n < nP; n++, node = child_node) {
		child_node = ACnodebuf_tryMovingToChild(P[n], 0, &(node->ac_node1_id));
		if (child_node != NULL)
			continue;
		child_node = ACnodebuf_tryMovingToChild(P[n], 1, &(node->ac_node2_id));
		if (child_node != NULL)
			continue;
		child_node = ACnodebuf_tryMovingToChild(P[n], 2, &(node->ac_node3_id));
		if (child_node != NULL)
			continue;
		child_node = ACnodebuf_tryMovingToChild(P[n], 3, &(node->ac_node4_id));
		if (child_node != NULL)
			continue;
		error("dictionary contains more than 4 distinct letters");
	}
	if (node->P_id == -1)
		node->P_id = P_id;
	else
		dupsbuf_append(node->P_id, P_id);
	return;
}

static void ACuldna_init()
{
	dupsbuf_reset();
	ACnodebuf_reset();
	ACnodebuf_newNode();
	ACnodebuf_addPattern("aabaz", 5, 1);
	ACnodebuf_addPattern("abzab", 5, 2);
	ACnodebuf_addPattern("aabaz", 5, 3);
	return;
}


/****************************************************************************
 * Exact matching
 * ======================================================
 */

static int ACuldna_exact_search()
{
	return 0;
}


/****************************************************************************
 * .Call entry points: "ACuldna_init_with_StrVect"
 *                 and "ACuldna_init_with_BStringList"
 *
 * Arguments:
 *   'dict': a string vector (aka character vector) containing the
 *           uniform-length dictionary for ACuldna_init_with_StrVect.
 *           A list of (pattern@data@xp, pattern@offset, pattern@length)
 *           triplets containing the uniform-length dictionary for
 *           ACuldna_init_with_BStringList.
 *
 * Return an R list with the following elements:
 *   - AC_tree: XInteger object containing the Aho-Corasick 4-ary tree built
 *         from 'dict'.
 *   - AC_base_codes: integer vector containing the 4 character codes (ASCII)
 *         attached to the 4 child slots of any node in the AC_tree object.
 *   - dups: an unnamed (and eventually empty) list of integer vectors
 *         containing the indices of the duplicated words found in 'dict'.
 *
 ****************************************************************************/

SEXP ACuldna_init_with_StrVect(SEXP dict)
{
	int subj_offset, subj_length, pat_length, c1, c2, c3, c4;
	const Rbyte *subj;
	SEXP buf, ans, ans_names, ans_elt;

	ACuldna_init();
	error("Not ready yet!\n");

	PROTECT(ans = NEW_LIST(3));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("AC_tree"));
	SET_STRING_ELT(ans_names, 1, mkChar("AC_base_codes"));
	SET_STRING_ELT(ans_names, 2, mkChar("dups"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "AC_tree" element */
	PROTECT(ans_elt = NEW_NUMERIC(4));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "AC_base_codes" element */
	PROTECT(ans_elt = NEW_INTEGER(4));
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "dups" element */
	PROTECT(ans_elt = NEW_INTEGER(4));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

SEXP ACuldna_init_with_BStringList(SEXP dict)
{
	SEXP ans;

	error("Not ready yet!\n");
	return ans;
}

