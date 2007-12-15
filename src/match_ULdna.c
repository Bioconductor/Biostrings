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

IBBuf dups_bbuf;

static void dups_bbuf_append(int P_id1, int P_id2)
{
	int i;
	IBuf *line, new_line;

	for (i = 0, line = dups_bbuf.ibufs; i < dups_bbuf.count; i++, line++) {
		if (line->vals[0] == P_id1) {
			_IBuf_insert_at(line, line->count, P_id2);
			return;
		}
	}
	_IBuf_init(&new_line);
	_IBuf_insert_at(&new_line, new_line.count, P_id1);
	_IBuf_insert_at(&new_line, new_line.count, P_id2);
	_IBBuf_insert_at(&dups_bbuf, dups_bbuf.count, new_line);
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

static ACNode *ACnodebuf;
static int ACnodebuf_maxcount, ACnodebuf_count;
static int ACnodebuf_pattern_length;
static int ACnodebuf_base_codes[ALPHABET_LENGTH];

static void ACnodebuf_init()
{
	int i;

	/* No memory leak here, because we use transient storage allocation */
	ACnodebuf = NULL;
	ACnodebuf_maxcount = ACnodebuf_count = 0;
	ACnodebuf_pattern_length = -1;
	for (i = 0; i < ALPHABET_LENGTH; i++)
		ACnodebuf_base_codes[i] = -1;
	return;
}

static void ACnodebuf_init_node(ACNode *node, int parent_id)
{
	int i;

	node->parent_id = parent_id;
	for (i = 0; i < ALPHABET_LENGTH; i++)
		node->child_id[i] = -1;
	node->flink = -1;
	node->P_id = -1;
	return;
}

static void ACnodebuf_get_more_room()
{
	long new_maxcount;

	if (ACnodebuf_maxcount == 0)
		new_maxcount = 50000;
	else
		new_maxcount = 2 * ACnodebuf_maxcount;
	ACnodebuf = Srealloc((char *) ACnodebuf, new_maxcount,
				(long) ACnodebuf_maxcount, ACNode);
	ACnodebuf_maxcount = new_maxcount;
	return;
}

static int ACnodebuf_append_node(int parent_id)
{
	ACNode *node;

	if (ACnodebuf_count >= ACnodebuf_maxcount)
		ACnodebuf_get_more_room();
	node = ACnodebuf + ACnodebuf_count;
	ACnodebuf_init_node(node, parent_id);
	return ACnodebuf_count++;
}

static int ACnodebuf_tryMovingToChild(int node_id, int childslot, char c)
{
	ACNode *node;
	int *child_id;

	if (ACnodebuf_base_codes[childslot] == -1) {
		ACnodebuf_base_codes[childslot] = c;
	} else {
		if (c != ACnodebuf_base_codes[childslot])
			return -1;
	}
	node = ACnodebuf + node_id;
	child_id = node->child_id + childslot;
	if (*child_id == -1)
		*child_id = ACnodebuf_append_node(node_id);
	return *child_id;
}

static void ACnodebuf_addPattern(const char *P, int nP, int P_id)
{
	int n, node_id, child_id, i;
	char c;
	ACNode *node;

	if (ACnodebuf_pattern_length == -1) {
		if (nP == 0)
			error("dictionary contains empty patterns");
		ACnodebuf_pattern_length = nP;
	} else {
		if (nP != ACnodebuf_pattern_length)
			error("all patterns in dictionary must have the same length");
	}
	for (n = 0, node_id = 0; n < nP; n++, node_id = child_id) {
		c = P[n];
		for (i = 0; i < ALPHABET_LENGTH; i++) {
			child_id = ACnodebuf_tryMovingToChild(node_id, i, c);
			if (child_id != -1)
				break;
		}
		if (child_id == -1)
			error("dictionary contains more than %d distinct letters", ALPHABET_LENGTH);
	}
	node = ACnodebuf + node_id;
	if (node->P_id == -1)
		node->P_id = P_id;
	else
		dups_bbuf_append(node->P_id, P_id);
	return;
}

static void ULdna_init(const Pattern *patterns, int dict_length)
{
	const Pattern *p;
	int i;

	_IBBuf_init(&dups_bbuf);
	ACnodebuf_init();
	ACnodebuf_append_node(0);
	for (i = 0, p = patterns; i < dict_length; i++, p++)
		ACnodebuf_addPattern(p->P, p->nP, i + 1);
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

IBBuf starts_bbuf;

static int get_child_id(ACNode *node, char c)
{
	int i;

	for (i = 0; i < ALPHABET_LENGTH; i++)
		if (c == ACnodebuf_base_codes[i])
			return node->child_id[i];
	return -1;
}

static ACNode *follow_flink(ACNode *node, ACNode *node0)
{
	/* not ready yet */
	return node0;
}

static void ULdna_exact_search(int uldna_len,
		ACNode *ACtree, int ACtree_length,
		const char *S, int nS)
{
	int i, n, child_id;
	IBuf *starts_buf;
	ACNode *node;
	char c;

	starts_bbuf.ibufs = Salloc((long) uldna_len, IBuf);
	for (i = 0, starts_buf = starts_bbuf.ibufs; i < uldna_len; i++, starts_buf++)
		_IBuf_init(starts_buf);
	starts_bbuf.maxcount = starts_bbuf.count = uldna_len;
	node = ACtree;
	for (n = 0; n < nS; n++) {
		c = S[n];
		child_id = get_child_id(node, c);
		while (child_id == -1) {
			if (node == ACtree)
				goto continue0;
			node = follow_flink(node, ACtree);
			child_id = get_child_id(node, c);
		}
		node = ACtree + child_id;
		if (node->P_id != -1) {
			starts_buf = starts_bbuf.ibufs + node->P_id - 1;
			_IBuf_insert_at(starts_buf, starts_buf->count, n + 1);
		}
		continue0: ;
	}
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
 * Turn the ACnodebuf_base_codes array into an R vector of ALPHABET_LENGTH integers
 */
static SEXP base_codes_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ALPHABET_LENGTH));
	memcpy(INTEGER(ans), ACnodebuf_base_codes, ALPHABET_LENGTH * sizeof(int));
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
	PROTECT(ans_elt = _IBBuf_asLIST(&dups_bbuf));
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
	dups_bbuf = _LIST_asIBBuf(uldna_dups);
	tag = R_ExternalPtrTag(ACtree_xp);
	ACtree = (ACNode *) INTEGER(tag);
	ACtree_length = LENGTH(tag) / INTS_PER_ACNODE;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;

	ULdna_exact_search(uldna_len, ACtree, ACtree_length, (char *) subj, subj_length);

	return _IBBuf_asLIST(&starts_bbuf);
}

