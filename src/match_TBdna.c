/****************************************************************************
 *                    Fast matching of a DNA dictionary                     *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_match_TBdna()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_TBdna.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_TBdna.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Match reporting
 * ===============
 */

static int match_reporting_mode; // 0, 1 or 2
static IntBuf match_count; // used when mode == 0 and initialized when mode == 2
static IntBBuf ends_bbuf;  // used when mode >= 1

static void init_match_reporting(int not_tail, int is_count_only, int length)
{
	match_reporting_mode = is_count_only ? (not_tail ? 0 : 2) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		match_count = _new_IntBuf(length, length);
	if (match_reporting_mode >= 1)
		ends_bbuf = _new_IntBBuf(length, length);
	return;
}

static void report_match(int poffset, int end)
{
	IntBuf *ends_buf;

	if (match_reporting_mode == 0) {
		match_count.elts[poffset]++;
		return;
	}
	ends_buf = ends_bbuf.elts + poffset;
	_IntBuf_insert_at(ends_buf, ends_buf->nelt, end);
	return;
}

static void report_matches_for_dups(const int *dups, int length)
{
	int poffset, *val;
	IntBuf *ends_buf;

	if (match_reporting_mode == 0) {
		for (poffset = 0, val = match_count.elts;
		     poffset < length;
		     poffset++, val++, dups++) {
			if (*dups == 0)
				continue;
			*val = *(match_count.elts + *dups - 1);
		}
		return;
	}
	for (poffset = 0, ends_buf = ends_bbuf.elts;
	     poffset < length;
	     poffset++, ends_buf++, dups++) {
		if (*dups == 0)
			continue;
		*ends_buf = *(ends_bbuf.elts + *dups - 1);
	}
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

static int code2childoffset_chrtrtable[CHRTRTABLE_LENGTH];

static int get_child_id(const ACNode *node, char c)
{
	int offset;

	offset = code2childoffset_chrtrtable[(unsigned char) c];
	if (offset == -1)
		return -1;
	return node->child_id[offset];
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
 * follow_string() function.
 */
static void set_shortcut(ACNode *basenode, char c, int node_id)
{
	int offset, *slot;

	offset = code2childoffset_chrtrtable[(unsigned char) c];
	if (offset == -1)
		return;
	slot = basenode->child_id + offset;
	if (*slot == -1)
		*slot = node_id;
	return;
}

static int follow_string(ACNode *node0, const int *base_codes, const char *S, int nS)
{
	int basenode_id, node_id, child_id, n, subcall_nS;
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
			Rprintf("[DEBUG] follow_string():");
			sprintf(format, "%%%ds", 1 + 2*rec_level);
			Rprintf(format, " ");
			snprintf(pathbuf, basenode->depth + 1, "%s", S - basenode->depth);
			Rprintf("On basenode_id=%d (basepath=%s), reading S[%d]=%c\n", basenode_id, pathbuf, n, *S);
		}
#endif
		node_id = basenode_id;
		node = basenode;
		while (1) {
			child_id = get_child_id(node, *S);
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
				subcall_nS = node->depth - 1;
				node->flink = follow_string(node0, base_codes, S - subcall_nS, subcall_nS);
				rec_level--;
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					Rprintf("[DEBUG] follow_string():");
					Rprintf(format, " ");
					Rprintf("setting failure link %d -> %d\n", node_id, node->flink);
				}
#endif
			}
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] follow_string():");
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
			Rprintf("[DEBUG] follow_string():");
			Rprintf(format, " ");
			Rprintf("moving to basenode %d\n", basenode_id);
		}
#endif
		// Finding a match cannot happen during a nested call to follow_string()
		// so there is no need to check that rec_level is 0
		if (basenode->P_id != -1)
			report_match(basenode->P_id - 1, n + 1);
	}
	return basenode_id;
}

static void CWdna_exact_search(ACNode *node0, const int *base_codes, const char *S, int nS)
{
	_init_chrtrtable(base_codes, MAX_CHILDREN_PER_ACNODE, code2childoffset_chrtrtable);
	follow_string(node0, base_codes, S, nS);
	return;
}


/****************************************************************************
 * Inexact matching on the tails of the TBdna_PDict object
 * =======================================================
 */

static void TBdna_match_pattern_tail(RoSeq pattern_tail, RoSeq S,
		int poffset, const int *dups, int dups_len,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int dup0, i, end;
	IntBuf *ends_buf, *ends_buf0;

	ends_buf = ends_buf0 = ends_bbuf.elts + poffset;
	dup0 = dups[poffset];
	if (dup0 != 0)
		ends_buf0 = ends_bbuf.elts + dup0 - 1;
	for (i = 0; i < ends_buf0->nelt; i++) {
		end = ends_buf0->elts[i];
		if (_is_matching(pattern_tail, S, end, max_mm, fixedP, fixedS)) {
			/* Match */
			if (is_count_only) {
				match_count.elts[poffset]++;
				continue;
			}
			if (dup0 == 0) {
				ends_buf0->elts[i] += pattern_tail.nelt;
				continue;
			}
			end += pattern_tail.nelt;
			_IntBuf_insert_at(ends_buf, ends_buf->nelt, end);
			continue;
		}
		/* Mismatch */
		if (is_count_only)
			continue;
		if (dup0 != 0)
			continue;
		/*
		 * We need to shrink the buffer we are walking on! This is safe
		 * because shrinking an IntBuf object should never trigger reallocation.
		 */
		_IntBuf_delete_at(ends_buf0, i--);
	}
	return;
}

static void TBdna_match_tail(SEXP tail, RoSeq S,
		const int *dups, int dups_len,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int poffset;
	CachedXStringSet cached_tail;
	RoSeq pattern_tail;

	cached_tail = _new_CachedXStringSet(tail);
	// The duplicated must be treated BEFORE the first pattern they
	// duplicate, hence we must walk from last to first.
	for (poffset = dups_len - 1; poffset >= 0; poffset--) {
		pattern_tail = _get_CachedXStringSet_elt_asRoSeq(&cached_tail, poffset);
		TBdna_match_pattern_tail(pattern_tail, S,
			poffset, dups, dups_len,
			max_mm, fixedP, fixedS, is_count_only);
	}
	return;
}


/****************************************************************************
 * .Call entry point: "match_TBdna"
 *
 * Arguments:
 *   'actree_nodes_xp': uldna_pdict@actree@nodes@xp
 *   'actree_base_codes': uldna_pdict@actree@base_codes
 *   'pdict_dups': pdict@dups
 *   'pdict_head_XStringSet': pdict@head (can be NULL)
 *   'pdict_tail_XStringSet': pdict@tail (can be NULL)
 *   'subject_XString': subject
 *   'max_mismatch': max.mismatch (max nb of mismatches in the tail)
 *   'fixed': logical vector of length 2
 *   'count_only': TRUE or FALSE
 *   'envir': environment to be populated with the matches
 *
 * Return an R object that will be assigned to the 'ends' slot of the
 * MIndex object returned by matchPDict(). Refer to the description
 * of this slot in the matchPDict.R file for the details.
 *
 ****************************************************************************/

SEXP match_TBdna(SEXP actree_nodes_xp, SEXP actree_base_codes,
		SEXP pdict_dups, SEXP pdict_head_XStringSet, SEXP pdict_tail_XStringSet,
		SEXP subject_XString,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	ACNode *actree_nodes;
	RoSeq S;
	int is_count_only, no_head, no_tail;

	actree_nodes = (ACNode *) INTEGER(R_ExternalPtrTag(actree_nodes_xp));
	S = _get_XString_asRoSeq(subject_XString);
	is_count_only = LOGICAL(count_only)[0];
	no_head = pdict_head_XStringSet == R_NilValue;
	no_tail = pdict_tail_XStringSet == R_NilValue;

	if (!no_head)
		error("matchPDict() doesn't support PDict objects with a head yet, sorry!");
	init_match_reporting(no_tail, is_count_only, LENGTH(pdict_dups));
	CWdna_exact_search(actree_nodes, INTEGER(actree_base_codes), S.elts, S.nelt);
	if (no_tail) {
		report_matches_for_dups(INTEGER(pdict_dups), LENGTH(pdict_dups));
	} else {
		TBdna_match_tail(pdict_tail_XStringSet, S,
			INTEGER(pdict_dups), LENGTH(pdict_dups),
			INTEGER(max_mismatch)[0], LOGICAL(fixed)[0], LOGICAL(fixed)[1],
			is_count_only);
	}
	if (is_count_only)
		return _IntBuf_asINTEGER(&match_count);
	if (envir == R_NilValue)
		return _IntBBuf_asLIST(&ends_bbuf, 1);
	return _IntBBuf_toEnvir(&ends_bbuf, envir, 1);
}

