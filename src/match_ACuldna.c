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

typedef struct pattern {
	const char *P;
	int nP;
} Pattern;


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

#define INTS_PER_ACNODE (sizeof(ACNode) / sizeof(int))

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
	for (n = 0, node = ACnodebuf; n < nP; n++, node = child_node) {
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

static void ACuldna_init(const Pattern *patterns, int dict_length)
{
	const Pattern *p;
	int i;

	dupsbuf_reset();
	ACnodebuf_reset();
	ACnodebuf_newNode();
	for (i = 0, p = patterns; i < dict_length; i++, p++)
		ACnodebuf_addPattern(p->P, p->nP, i + 1);
	return;
}

static SEXP ACuldna_init_mkSEXP()
{
	SEXP ans, ans_names, ans_elt, tag, ans_elt_elt;
	DupsBufLine *line;
	int tag_length, i;

	PROTECT(ans = NEW_LIST(3));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("AC_tree_xp"));
	SET_STRING_ELT(ans_names, 1, mkChar("AC_base_codes"));
	SET_STRING_ELT(ans_names, 2, mkChar("dups"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "AC_tree_xp" element */
	PROTECT(ans_elt = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	tag_length = ACnodebuf_count * INTS_PER_ACNODE;
	PROTECT(tag = NEW_INTEGER(tag_length));
	memcpy((char *) INTEGER(tag), ACnodebuf, tag_length * sizeof(int));
	R_SetExternalPtrTag(ans_elt, tag);
	UNPROTECT(1);
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "AC_base_codes" element */
	PROTECT(ans_elt = NEW_INTEGER(4));
	for (i = 0; i < 4; i++)
		INTEGER(ans_elt)[i] = ACnodebuf_base_codes[i];
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	/* set the "dups" element */
	PROTECT(ans_elt = NEW_LIST(dupsbuf_count));
	for (i = 0, line = dupsbuf; i < dupsbuf_count; i++, line++) {
		PROTECT(ans_elt_elt = NEW_INTEGER(line->count));
		memcpy((char *) INTEGER(ans_elt_elt), line->vals, line->count * sizeof(int));
		SET_ELEMENT(ans_elt, i, ans_elt_elt);
		UNPROTECT(1);
	}
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
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
 *   - AC_tree_xp: "externalptr" object pointing to the Aho-Corasick 4-ary
 *         tree built from 'dict'.
 *   - AC_base_codes: integer vector containing the 4 character codes (ASCII)
 *         attached to the 4 child slots of any node in the tree pointed by
 *         AC_tree_xp.
 *   - dups: an unnamed (and eventually empty) list of integer vectors
 *         containing the indices of the duplicated words found in 'dict'.
 *
 ****************************************************************************/

SEXP ACuldna_init_with_StrVect(SEXP dict)
{
	Pattern *patterns, *p;
	SEXP dict_elt;
	int i;

	patterns = Salloc((long) LENGTH(dict), Pattern);
	for (i = 0, p = patterns; i < LENGTH(dict); i++, p++) {
		dict_elt = STRING_ELT(dict, i);
		p->P = CHAR(dict_elt);
		p->nP = LENGTH(dict_elt);
	}
	ACuldna_init(patterns, LENGTH(dict));
	return ACuldna_init_mkSEXP();
}

SEXP ACuldna_init_with_BStringList(SEXP dict)
{
	SEXP ans;

	error("Not ready yet!\n");
	return ans;
}

