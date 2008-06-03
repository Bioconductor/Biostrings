/****************************************************************************
 *                                                                          *
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
static IntBBuf match_ends;  // used when mode >= 1
static IntBuf matching_poffsets;

static void init_match_reporting(int no_tail, int is_count_only, int length)
{
	match_reporting_mode = is_count_only ? (no_tail ? 0 : 2) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		match_count = _new_IntBuf(length, length, 0);
	if (match_reporting_mode >= 1)
		match_ends = _new_IntBBuf(length, length);
	matching_poffsets = _new_IntBuf(0, 0, 0);
	return;
}

static void report_match(int poffset, int end)
{
	int is_new_matching_poffset;
	IntBuf *ends_buf;

	if (match_reporting_mode == 0) {
		is_new_matching_poffset = match_count.elts[poffset]++ == 0;
	} else {
		ends_buf = match_ends.elts + poffset;
		is_new_matching_poffset = ends_buf->nelt == 0;
		_IntBuf_insert_at(ends_buf, ends_buf->nelt, end);
	}
	if (is_new_matching_poffset)
		_IntBuf_insert_at(&matching_poffsets,
				  matching_poffsets.nelt, poffset);
	return;
}

static void report_matches_for_dups(SEXP unq2dup)
{
	int n1, i, *poffset, j, *dup, k1, k2;
	SEXP dups;

	/* The value of matching_poffsets.nelt can increase during the for loop! */
	n1 = matching_poffsets.nelt;
	for (i = 0, poffset = matching_poffsets.elts; i < n1; i++, poffset++) {
		k1 = *poffset;
		dups = VECTOR_ELT(unq2dup, k1);
		if (dups == R_NilValue)
			continue;
		for (j = 0, dup = INTEGER(dups); j < LENGTH(dups); j++, dup++) {
			k2 = *dup - 1;
			if (match_reporting_mode == 0)
				match_count.elts[k2] = match_count.elts[k1];
			else
				match_ends.elts[k2] = match_ends.elts[k1];
			_IntBuf_insert_at(&matching_poffsets,
					  matching_poffsets.nelt, k2);
		}
	}
	return;
}

static void merge_matches(IntBuf *global_match_count,
		IntBBuf *global_match_ends, int view_offset)
{
	int i, *poffset;
	IntBuf *ends_buf, *global_ends_buf;

	for (i = 0, poffset = matching_poffsets.elts;
	     i < matching_poffsets.nelt;
	     i++, poffset++)
	{
		if (match_reporting_mode == 0 || match_reporting_mode == 2) {
			global_match_count->elts[*poffset] += match_count.elts[*poffset];
			match_count.elts[*poffset] = 0;
		} else {
			ends_buf = match_ends.elts + *poffset;
			global_ends_buf = global_match_ends->elts + *poffset;
			_IntBuf_append_shifted_vals(global_ends_buf,
				ends_buf->elts, ends_buf->nelt, view_offset);
		}
		if (match_reporting_mode >= 1)
			match_ends.elts[*poffset].nelt = 0;
	}
	matching_poffsets.nelt = 0;
	return;
}


/****************************************************************************
 * Exact matching
 * ==============
 */

static int slotno_chrtrtable[CHRTRTABLE_LENGTH];

static int get_child_node_id(const ACNode *node, char c)
{
	int slotno;

	slotno = slotno_chrtrtable[(unsigned char) c];
	if (slotno == -1)
		return -1;
	return node->child_id[slotno];
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
 * walk_string() function.
 */
static void set_shortcut(ACNode *node, char c, int next_node_id)
{
	int slotno, *slot;

	slotno = slotno_chrtrtable[(unsigned char) c];
	if (slotno == -1)
		return;
	slot = node->child_id + slotno;
	if (*slot == -1)
		*slot = next_node_id;
	return;
}

/*
 * We use indirect recursion for walking the Aho-Corasick tree.
 * This indirect recursion involves the 2 core functions path_to_node_id() and
 * get_next_node_id(): the latter calls the former which in turn calls the
 * latter.
 */
static int get_next_node_id(ACNode *node0, const int *base_codes,
		int node_id, const char *Stail, char c);

static int path_to_node_id(ACNode *node0, const int *base_codes,
		const char *path, int path_len)
{
	int node_id, n;
	ACNode *node;

	node_id = 0;
	for (n = 0; n < path_len; n++, path++) {
		node = node0 + node_id;
		node_id = get_next_node_id(node0, base_codes,
				node_id, path, *path);
	}
	return node_id;
}

/*
 * An important trick here is that the chars located _before_ 'Stail' will
 * always describe the path that goes from the root node to 'node_id'.
 * More precisely, if d is the depth of 'node_id', then this path is made
 * of the path[i] chars where path is 'Stail' - d and 0 <= i < d.
 */
static int get_next_node_id(ACNode *node0, const int *base_codes,
		int node_id, const char *Stail, char c)
{
	ACNode *node, *next_node;
	int next_node_id, child_node_id, subpath_len;
	const char *subpath;
#ifdef DEBUG_BIOSTRINGS
	static int rec_level = 0;
	char format[20], pathbuf[2000];
#endif

	node = node0 + node_id;
	next_node_id = node_id;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] ENTERING get_next_node_id():");
		sprintf(format, "%%%ds", 1 + 2*rec_level);
		Rprintf(format, " ");
		snprintf(pathbuf, node->depth + 1, "%s", Stail - node->depth);
		Rprintf("node_id=%d path=%s c=%c\n",
			node_id, pathbuf, c);
	}
#endif
	while (1) {
		next_node = node0 + next_node_id;
		child_node_id = get_child_node_id(next_node, c);
		if (child_node_id != -1) {
			next_node_id = child_node_id;
			break;
		}
		if (next_node_id == 0) {
			//next_node->flink = 0;
			break;
		}
		if (next_node->flink == -1) {
			rec_level++;
			subpath_len = next_node->depth - 1;
			subpath = Stail - subpath_len;
			next_node->flink = path_to_node_id(node0, base_codes,
						subpath, subpath_len);
			rec_level--;
		}
		next_node_id = next_node->flink;
	}
	set_shortcut(node, c, next_node_id);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] LEAVING get_next_node_id(): ");
		Rprintf(format, " ");
		Rprintf("next_node_id=%d\n", next_node_id);
	}
#endif
	return next_node_id;
}

static int walk_string(ACNode *node0, const int *base_codes,
		const char *S, int nS)
{
	int basenode_id, node_id, child_id, n, subwalk_nS;
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
			Rprintf("[DEBUG] walk_string():");
			sprintf(format, "%%%ds", 1 + 2*rec_level);
			Rprintf(format, " ");
			snprintf(pathbuf, basenode->depth + 1, "%s", S - basenode->depth);
			Rprintf("On basenode_id=%d (basepath=%s), reading S[%d]=%c\n", basenode_id, pathbuf, n, *S);
		}
#endif
		node_id = basenode_id;
		node = basenode;
		while (1) {
			child_id = get_child_node_id(node, *S);
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
				subwalk_nS = node->depth - 1;
				node->flink = walk_string(node0, base_codes, S - subwalk_nS, subwalk_nS);
				rec_level--;
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					Rprintf("[DEBUG] walk_string():");
					Rprintf(format, " ");
					Rprintf("setting failure link %d -> %d\n", node_id, node->flink);
				}
#endif
			}
#ifdef DEBUG_BIOSTRINGS
			if (debug) {
				Rprintf("[DEBUG] walk_string():");
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
			Rprintf("[DEBUG] walk_string():");
			Rprintf(format, " ");
			Rprintf("moving to basenode %d\n", basenode_id);
		}
#endif
		// Finding a match cannot happen during a nested call to walk_string()
		// so there is no need to check that rec_level is 0
		if (basenode->P_id != -1)
			report_match(basenode->P_id - 1, n + 1);
	}
	return basenode_id;
}

static void CWdna_exact_search(ACNode *node0,
		const int *base_codes, RoSeq S)
{
	_init_chrtrtable(base_codes, MAX_CHILDREN_PER_ACNODE, slotno_chrtrtable);
	walk_string(node0, base_codes, S.elts, S.nelt);
	return;
}

static void CWdna_exact_search_on_nonfixedS(ACNode *node0,
		const int *base_codes, RoSeq S)
{
	IntBuf cnode_ids; // buffer of current node ids
	int n, npointers, i, node_id, next_node_id, is_first, j, base, P_id;
	const char *Stail;
	char c;

	_init_chrtrtable(base_codes, MAX_CHILDREN_PER_ACNODE, slotno_chrtrtable);
	cnode_ids = _new_IntBuf(256, 0, 0);
	_IntBuf_insert_at(&cnode_ids, 0, 0);
	for (n = 1, Stail = S.elts; n <= S.nelt; n++, Stail++) {
		c = *Stail;
		npointers = cnode_ids.nelt;
		// move and split pointers
		for (i = 0; i < npointers; i++) {
			node_id = cnode_ids.elts[i];
			is_first = 1;
			for (j = 0, base = 1; j < 4; j++, base *= 2) {
				if ((((unsigned char) c) & base) != 0) {
					next_node_id = get_next_node_id(node0,
							base_codes,
							node_id, Stail, base);
					if (is_first) {
						cnode_ids.elts[i] = next_node_id;
						is_first = 0;
					} else {
						_IntBuf_insert_at(&cnode_ids,
							cnode_ids.nelt, next_node_id);
					}
				}
			}
		}
		// merge pointers and report matches
		for (i = 0; i < cnode_ids.nelt; i++) {
			node_id = cnode_ids.elts[i];
			// FIXME: This merging algo is dumb and inefficient!
			// There must be a way to do something better.
			for (j = i + 1; j < cnode_ids.nelt; j++) {
				if (cnode_ids.elts[j] == node_id)
					_IntBuf_delete_at(&cnode_ids, j--);
			}
			P_id = node0[node_id].P_id;
			if (P_id != -1)
				report_match(P_id - 1, n);
		}
		// error if too many remaining pointers
		if (cnode_ids.nelt > 4096)
			error("too many IUPAC ambiguity letters in 'subject'");
	}
	return;
}


/****************************************************************************
 * Inexact matching on the tails of the TBdna_PDict object
 * =======================================================
 */

static int match_TBdna_Ptail(RoSeq Ptail, RoSeq S, int k1, int k2,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int i, end1;
	IntBuf *ends_buf1, *ends_buf2;

	ends_buf1 = match_ends.elts + k1;
	ends_buf2 = match_ends.elts + k2;
	for (i = 0; i < ends_buf1->nelt; i++) {
		end1 = ends_buf1->elts[i];
		if (!_is_matching(Ptail, S, end1, max_mm, fixedP, fixedS)) {
			/* Mismatch */
			if (match_reporting_mode == 0)
				continue;
			if (k1 != k2)
				continue;
			/*
			 * We need to shrink the buffer we are walking on! This is safe
			 * because shrinking an IntBuf object should never trigger
			 * reallocation.
			 */
			_IntBuf_delete_at(ends_buf1, i--);
			continue;
		}
		/* Match */
		if (is_count_only) {
			match_count.elts[k2]++;
			if (match_reporting_mode == 0)
				continue;
		}
		if (k1 != k2) {
			_IntBuf_insert_at(ends_buf2, ends_buf2->nelt, end1 + Ptail.nelt);
		} else {
			ends_buf2->elts[i] += Ptail.nelt;
		}
	}
	return is_count_only ? match_count.elts[k2] : ends_buf2->nelt;
}

static void match_TBdna_tail(CachedXStringSet *cached_tail, RoSeq S,
		int pdict_len, SEXP unq2dup,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int n1, i, j, *dup, k1, k2, nmatches;
	SEXP dups;
	RoSeq Ptail;

	/* The number of elements in matching_poffsets can increase or decrease
	   during the for loop! */
	n1 = matching_poffsets.nelt;
	for (i = 0; i < n1; i++) {
		k1 = matching_poffsets.elts[i];
		dups = VECTOR_ELT(unq2dup, k1);
		if (dups != R_NilValue) {
			for (j = 0, dup = INTEGER(dups); j < LENGTH(dups); j++, dup++) {
				k2 = *dup - 1;
				Ptail = _get_CachedXStringSet_elt_asRoSeq(cached_tail, k2);
				nmatches = match_TBdna_Ptail(Ptail, S, k1, k2,
							     max_mm, fixedP, fixedS, is_count_only);
				if (nmatches != 0)
					_IntBuf_insert_at(&matching_poffsets,
							  matching_poffsets.nelt, k2);
			}
		}
		Ptail = _get_CachedXStringSet_elt_asRoSeq(cached_tail, k1);
		nmatches = match_TBdna_Ptail(Ptail, S, k1, k1,
					     max_mm, fixedP, fixedS, is_count_only);
		if (nmatches == 0) {
			_IntBuf_delete_at(&matching_poffsets, i--);
			n1--;
		}
	}
	return;
}

static void match_TBdna(ACNode *actree_nodes, const int *base_codes,
		int pdict_len, SEXP unq2dup,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		RoSeq S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	if (fixedS)
		CWdna_exact_search(actree_nodes, base_codes, S);
	else
		CWdna_exact_search_on_nonfixedS(actree_nodes, base_codes, S);
	if (cached_tail == NULL)
		report_matches_for_dups(unq2dup);
	else
		match_TBdna_tail(cached_tail, S,
			pdict_len, unq2dup,
			max_mm, fixedP, fixedS, is_count_only);
	return;
}


/****************************************************************************
 * .Call entry point: "XString_match_TBdna" and "XStringViews_match_TBdna"
 *
 * Arguments:
 *   'actree_nodes_xp': pdict@actree@nodes@xp
 *   'actree_base_codes': pdict@actree@base_codes
 *   'pdict_length': length(pdict)
 *   'pdict_unq2dup': pdict@dups@unq2dup
 *   'pdict_head': pdict@head (XStringSet or NULL)
 *   'pdict_tail': pdict@tail (XStringSet or NULL)
 *   'subject': subject
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

SEXP XString_match_TBdna(SEXP actree_nodes_xp, SEXP actree_base_codes,
		SEXP pdict_length, SEXP pdict_unq2dup,
		SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	ACNode *actree_nodes;
	CachedXStringSet cached_head, *pcached_head, cached_tail, *pcached_tail;
	RoSeq S;
	int pdict_len, no_head, no_tail,
	    fixedP, fixedS, is_count_only;

	actree_nodes = (ACNode *) INTEGER(R_ExternalPtrTag(actree_nodes_xp));
	pdict_len = INTEGER(pdict_length)[0];
	pcached_head = pcached_tail = NULL;
	no_head = pdict_head == R_NilValue;
	if (!no_head) {
		error("matchPDict() doesn't support PDict objects with "
		      "a head yet, sorry!");
		cached_head = _new_CachedXStringSet(pdict_head);
		pcached_head = &cached_head;
	}
	no_tail = pdict_tail == R_NilValue;
	if (!no_tail) {
		cached_tail = _new_CachedXStringSet(pdict_tail);
		pcached_tail = &cached_tail;
	}
	S = _get_XString_asRoSeq(subject);
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	init_match_reporting(no_tail, is_count_only, pdict_len);
	match_TBdna(actree_nodes, INTEGER(actree_base_codes),
		pdict_len, pdict_unq2dup,
		pcached_head, pcached_tail,
		S,
		INTEGER(max_mismatch)[0], fixedP, fixedS, is_count_only);

	if (is_count_only)
		return _IntBuf_asINTEGER(&match_count);
	if (envir == R_NilValue)
		return _IntBBuf_asLIST(&match_ends, 1);
	return _IntBBuf_toEnvir(&match_ends, envir, 1);
}

SEXP XStringViews_match_TBdna(SEXP actree_nodes_xp, SEXP actree_base_codes,
		SEXP pdict_length, SEXP pdict_unq2dup,
		SEXP pdict_head, SEXP pdict_tail,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	ACNode *actree_nodes;
	CachedXStringSet cached_head, *pcached_head, cached_tail, *pcached_tail;
	RoSeq S, V;
	int pdict_len, no_head, no_tail,
	    fixedP, fixedS, is_count_only,
	    nviews, i, *view_start, *view_width, view_offset;
	IntBuf global_match_count;
	IntBBuf global_match_ends;

	actree_nodes = (ACNode *) INTEGER(R_ExternalPtrTag(actree_nodes_xp));
	pdict_len = INTEGER(pdict_length)[0];
	pcached_head = pcached_tail = NULL;
	no_head = pdict_head == R_NilValue;
	if (!no_head) {
		error("matchPDict() doesn't support PDict objects with "
		      "a head yet, sorry!");
		cached_head = _new_CachedXStringSet(pdict_head);
		pcached_head = &cached_head;
	}
	no_tail = pdict_tail == R_NilValue;
	if (!no_tail) {
		cached_tail = _new_CachedXStringSet(pdict_tail);
		pcached_tail = &cached_tail;
	}
	S = _get_XString_asRoSeq(subject);
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	if (is_count_only)
		global_match_count = _new_IntBuf(pdict_len, pdict_len, 0);
	else
		global_match_ends = _new_IntBBuf(pdict_len, pdict_len);
	init_match_reporting(no_tail, is_count_only, pdict_len);
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.nelt)
			error("'subject' has out of limits views");
		V.elts = S.elts + view_offset;
		V.nelt = *view_width;
		match_TBdna(actree_nodes, INTEGER(actree_base_codes),
			pdict_len, pdict_unq2dup,
			pcached_head, pcached_tail,
			V,
			INTEGER(max_mismatch)[0], fixedP, fixedS,
			is_count_only);
		merge_matches(&global_match_count, &global_match_ends, view_offset);
	}

	if (is_count_only)
		return _IntBuf_asINTEGER(&global_match_count);
	if (envir == R_NilValue)
		return _IntBBuf_asLIST(&global_match_ends, 1);
	return _IntBBuf_toEnvir(&global_match_ends, envir, 1);
}

