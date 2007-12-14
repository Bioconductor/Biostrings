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

SEXP match_ULdna_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_ULdna.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_ULdna.c'\n");
#endif
	return R_NilValue;
}

typedef struct pattern {
	const char *P;
	int nP;
} Pattern;


/****************************************************************************
 * Buffer of duplicates
 */

IBBuf dupsbuf;

static void dupsbuf_append(int P_id1, int P_id2)
{
	int i;
	IBuf *line, new_line;

	for (i = 0, line = dupsbuf.ibufs; i < dupsbuf.count; i++, line++) {
		if (line->vals[0] == P_id1) {
			_IBuf_insert_at(line, line->count, P_id2);
			return;
		}
	}
	_IBuf_init(&new_line);
	_IBuf_insert_at(&new_line, new_line.count, P_id1);
	_IBuf_insert_at(&new_line, new_line.count, P_id2);
	_IBBuf_insert_at(&dupsbuf, dupsbuf.count, new_line);
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
static int ACnodebuf_maxcount, ACnodebuf_count;
static int ACnodebuf_pattern_length;
static int ACnodebuf_base_codes[4];

void ACnodebuf_reset()
{
	int i;

	/* No memory leak here, because we use transient storage allocation */
	ACnodebuf = NULL;
	ACnodebuf_maxcount = ACnodebuf_count = 0;
	ACnodebuf_pattern_length = -1;
	for (i = 0; i < 4; i++)
		ACnodebuf_base_codes[i] = -1;
	return;
}

static int ACnodebuf_newNode()
{
	long new_maxcount;
	ACNode *node;

	if (ACnodebuf_count >= ACnodebuf_maxcount) {
		/* Buffer is full */
		if (ACnodebuf_maxcount == 0)
			new_maxcount = 50000;
		else
			new_maxcount = 2 * ACnodebuf_maxcount;
		ACnodebuf = Srealloc((char *) ACnodebuf, new_maxcount,
					(long) ACnodebuf_maxcount, ACNode);
		ACnodebuf_maxcount = new_maxcount;
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

static void ULdna_init(const Pattern *patterns, int dict_length)
{
	const Pattern *p;
	int i;

	_IBBuf_init(&dupsbuf);
	ACnodebuf_reset();
	ACnodebuf_newNode();
	for (i = 0, p = patterns; i < dict_length; i++, p++)
		ACnodebuf_addPattern(p->P, p->nP, i + 1);
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

IBBuf startsbuf;

static void ULdna_exact_search(int uldna_len,
		ACNode *ACtree, int ACtree_length,
		const char *S, int nS)
{
	int i;
	IBuf *ibuf_p;

	startsbuf.ibufs = Salloc((long) uldna_len, IBuf);
	for (i = 0, ibuf_p = startsbuf.ibufs; i < uldna_len; i++, ibuf_p++)
		_IBuf_init(ibuf_p);
	startsbuf.maxcount = startsbuf.count = uldna_len;
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
 * Turn the ACnodebuf_base_codes array into an R vector of 4 integers
 */
static SEXP base_codes_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(4));
	memcpy(INTEGER(ans), ACnodebuf_base_codes, 4 * sizeof(int));
	UNPROTECT(1);
	return ans;
}

/*
 * Return an R list with the following elements:
 *   - ACtree_xp: "externalptr" object pointing to the Aho-Corasick 4-ary
 *         tree built from 'dict'.
 *   - AC_base_codes: integer vector containing the 4 character codes (ASCII)
 *         attached to the 4 child slots of any node in the tree pointed by
 *         ACtree_xp.
 *   - dups: an unnamed (and eventually empty) list of integer vectors
 *         containing the indices of the duplicated words found in 'dict'.
 */
static SEXP ULdna_asLIST()
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
	PROTECT(ans_elt = _IBBuf_asLIST(&dupsbuf));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry points: "ULdna_init_with_StrVect"
 *                 and "ULdna_init_with_BStringList"
 *
 * Arguments:
 *   'dict': a string vector (aka character vector) containing the
 *           uniform-length dictionary for ULdna_init_with_StrVect.
 *           A list of (pattern@data@xp, pattern@offset, pattern@length)
 *           triplets containing the uniform-length dictionary for
 *           ULdna_init_with_BStringList.
 *
 * See ULdna_asLIST() for a description of the returned SEXP.
 *
 ****************************************************************************/

SEXP ULdna_init_with_StrVect(SEXP dict)
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
	ULdna_init(patterns, LENGTH(dict));
	return ULdna_asLIST();
}

SEXP ULdna_init_with_BStringList(SEXP dict)
{
	SEXP ans;

	error("Not ready yet!\n");
	return ans;
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
	dupsbuf = _LIST_asIBBuf(uldna_dups);
	tag = R_ExternalPtrTag(ACtree_xp);
	ACtree = (ACNode *) INTEGER(tag);
	ACtree_length = LENGTH(tag) / INTS_PER_ACNODE;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;

	ULdna_exact_search(uldna_len, ACtree, ACtree_length, (char *) subj, subj_length);

	return _IBBuf_asLIST(&startsbuf);
}

