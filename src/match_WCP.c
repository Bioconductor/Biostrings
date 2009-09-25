/****************************************************************************
 *            WEIGHTED CLUSTERED POSITIONS MATCHING                         *
 *		             Author: Patrick Aboyoun                                *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

typedef struct roseqslist {
	RoSeqs *elts;
	int nelt;
} RoSeqsList;


RoSeqsList _alloc_RoSeqsList(int nelt)
{
	RoSeqsList seqsList;

	seqsList.elts = Salloc((long) nelt, RoSeqs);
	seqsList.nelt = nelt;
	return seqsList;
}


/****************************************************************************
 * INPUTS:
 * key_seqs_list   - the list of dictionary keys
 * key_scores_list - the list of scores for the keys
 * key_total_nchar - the total number of characters across key sets
 * key_nchars      - the number of characters for each key set
 * pos_indices     - the position indices for the clusters
 * ndict           - the number of dictionary clusters
 * S_buffer        - an R list containing XString objects
 * S_char          - the subject character vector
 * S_nchar         - the number of characters in subject
 * S_shift         - location on subject to compute score
 */
static double compute_wcp_score(RoSeqsList *key_seqs_list,
		                        double **key_scores_list,
		                        int **key_order_list,
		                        int key_total_nchar,
		                        const int *key_nchars,
		                        const int *pos_indices,
		                        int ndict,
		                        int *match_buffer,
		                        RoSeqs *S_buffer,
		                        const char *S_char,
		                        int S_nchar, int S_shift)
{
	int i, match_loc, S_order, *key_nchar, *subscript, *key_order;
	RoSeqs *key_seqs;
	double *key_scores, score;

	if (S_shift < 0 || S_shift > (S_nchar - key_total_nchar))
		error("trying to compute the score from an invalid starting position");
	S_char += S_shift;

	score = 0;
	subscript = (int *) pos_indices;
	S_order = 0;
	for (i = 0, key_nchar = (int *) key_nchars; i < ndict; i++, key_nchar++) {
		S_buffer->elts[0].length = *key_nchar;
		key_seqs = &key_seqs_list->elts[i];
		key_scores = key_scores_list[i];
		key_order = key_order_list[i];

		Ocopy_byteblocks_from_subscript(subscript, *key_nchar,
				                (char *) S_buffer->elts[0].seq, *key_nchar,
				                S_char, S_nchar, sizeof(Rbyte));
		_get_RoSeqs_match(S_buffer, key_seqs, 0, &S_order, key_order, match_buffer, &match_loc);
		if (match_loc > 0) {
			score += key_scores[match_loc-1];
		}
		subscript += *key_nchar;
	}

	return score;
}


/*
 * --- .Call ENTRY POINT ---
 * WCP_score_starting_at() arguments are assumed to be:
 *   wcp: the weighted clustered positions (WCP) object
 *   subject: an XString object containing the subject sequence
 *   starting_at: an integer vector of arbitrary length (NAs accepted)
 */
SEXP WCP_score_starting_at(SEXP wcp, SEXP subject, SEXP starting_at)
{
	int i, j, ndict, *pos_indices, *start_elt;
	int *clust_ends, *end, curr_end, prev_end;
	int nkey, max_nkey, max_key_nchar, key_total_nchar;
	int *key_nchars, *nchar, *key_order;
	int *match_buffer, **key_order_list;
	cachedCharSeq S;
	RoSeqs S_buffer;
	RoSeqsList key_seqs_list;
	SEXP dict_list, clust_bins, dict_elt, key, ans;
	double **key_scores_list, *ans_elt;

	dict_list = GET_SLOT(GET_SLOT(wcp, install("dictList")), install("listData"));
	clust_bins = GET_SLOT(GET_SLOT(wcp, install("clusters")), install("bins"));
	pos_indices = INTEGER(GET_SLOT(clust_bins, install("unlistData")));
	clust_ends =
		INTEGER(GET_SLOT(GET_SLOT(clust_bins, install("partitioning")),
				         install("end")));

	prev_end = 0;
	max_nkey = 0;
	max_key_nchar = 0;
	key_total_nchar = 0;
	ndict = LENGTH(dict_list);
	key_nchars = (int *) R_alloc((long) ndict, sizeof(int));
	key_seqs_list = _alloc_RoSeqsList(ndict);
	key_scores_list = (double **) R_alloc((long) ndict, sizeof(double *));
	key_order_list = (int **) R_alloc((long) ndict, sizeof(int *));
	for (i = 0, nchar = key_nchars, end = clust_ends; i < ndict;
	     i++, nchar++, end++) {
		curr_end = *end;
		*nchar = curr_end - prev_end;
		key_total_nchar += *nchar;
		max_key_nchar = *nchar > max_key_nchar ? *nchar : max_key_nchar;
		prev_end = curr_end;

		dict_elt = VECTOR_ELT(dict_list, i);
		key = GET_SLOT(dict_elt, install("key"));
		nkey = _get_XStringSet_length(key);
		max_nkey = nkey > max_nkey ? nkey : max_nkey;
		key_seqs_list.elts[i] = _new_RoSeqs_from_XStringSet(nkey, key);
		key_scores_list[i] =
			REAL(VECTOR_ELT(GET_SLOT(GET_SLOT(dict_elt, install("table")),
					                 install("listData")), 0));
		key_order_list[i] = (int *) R_alloc((long) nkey, sizeof(int));
		for (j = 0, key_order = key_order_list[i]; j < nkey;
		     j++, key_order++)
			*key_order = j;
	}
	S_buffer = _alloc_RoSeqs(1);
	S_buffer.elts[0].seq = (char *) R_alloc((long) max_key_nchar, sizeof(char));
	S_buffer.elts[0].length = max_key_nchar;
	match_buffer = (int *) R_alloc((long) max_nkey, sizeof(int));

	S = cache_XRaw(subject);
	PROTECT(ans = NEW_NUMERIC(LENGTH(starting_at)));
	for (i = 0, start_elt = INTEGER(starting_at), ans_elt = REAL(ans);
	     i < LENGTH(starting_at); i++, start_elt++, ans_elt++) {
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_REAL;
			continue;
		}
		*ans_elt = compute_wcp_score(&key_seqs_list, key_scores_list,
				                     key_order_list,
				                     key_total_nchar, key_nchars,
				                     pos_indices, ndict, match_buffer,
				                     &S_buffer, S.seq, S.length, *start_elt - 1);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * XString_match_WCP() arguments are assumed to be:
 *   wcp: the weighted clustered positions (WCP) object
 *   subject: an XString object containing the subject sequence
 *   min_score: a single double (not NA)
 *   count_only: a single logical (not NA)
 */
SEXP XString_match_WCP(SEXP wcp, SEXP subject, SEXP min_score, SEXP count_only)
{
	int n1, n2, ndict, is_count_only, *pos_indices;
	int *clust_ends, *end, curr_end, prev_end;
	int nkey, max_nkey, max_key_nchar, key_total_nchar;
	int *key_nchars, *nchar, *key_order;
	int *match_buffer, **key_order_list;
	cachedCharSeq S;
	RoSeqs S_buffer;
	RoSeqsList key_seqs_list;
	SEXP dict_list, clust_bins, dict_elt, key;
	double **key_scores_list, minscore;

	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];

	dict_list = GET_SLOT(GET_SLOT(wcp, install("dictList")), install("listData"));
	clust_bins = GET_SLOT(GET_SLOT(wcp, install("clusters")), install("bins"));
	pos_indices = INTEGER(GET_SLOT(clust_bins, install("unlistData")));
	clust_ends =
		INTEGER(GET_SLOT(GET_SLOT(clust_bins, install("partitioning")),
				         install("end")));

	prev_end = 0;
	max_nkey = 0;
	max_key_nchar = 0;
	key_total_nchar = 0;
	ndict = LENGTH(dict_list);
	key_nchars = (int *) R_alloc((long) ndict, sizeof(int));
	key_seqs_list = _alloc_RoSeqsList(ndict);
	key_scores_list = (double **) R_alloc((long) ndict, sizeof(double *));
	key_order_list = (int **) R_alloc((long) ndict, sizeof(int *));
	for (n1 = 0, nchar = key_nchars, end = clust_ends; n1 < ndict;
	     n1++, nchar++, end++) {
		curr_end = *end;
		*nchar = curr_end - prev_end;
		key_total_nchar += *nchar;
		max_key_nchar = *nchar > max_key_nchar ? *nchar : max_key_nchar;
		prev_end = curr_end;

		dict_elt = VECTOR_ELT(dict_list, n1);
		key = GET_SLOT(dict_elt, install("key"));
		nkey = _get_XStringSet_length(key);
		max_nkey = nkey > max_nkey ? nkey : max_nkey;
		key_seqs_list.elts[n1] = _new_RoSeqs_from_XStringSet(nkey, key);
		key_scores_list[n1] =
			REAL(VECTOR_ELT(GET_SLOT(GET_SLOT(dict_elt, install("table")),
					                 install("listData")), 0));
		key_order_list[n1] = (int *) R_alloc((long) nkey, sizeof(int));
		for (n2 = 0, key_order = key_order_list[n1]; n2 < nkey;
		     n2++, key_order++)
			*key_order = n2;
	}
	S_buffer = _alloc_RoSeqs(1);
	S_buffer.elts[0].seq = (char *) R_alloc((long) max_key_nchar, sizeof(char));
	S_buffer.elts[0].length = max_key_nchar;
	match_buffer = (int *) R_alloc((long) max_nkey, sizeof(int));

	S = cache_XRaw(subject);
	_init_match_reporting(is_count_only ? mkString("COUNTONLY") : mkString("ASIRANGES"));
	for (n1 = 0, n2 = key_total_nchar; n2 <= S.length; n1++, n2++) {
		if (compute_wcp_score(&key_seqs_list, key_scores_list,
				              key_order_list,
				              key_total_nchar, key_nchars,
				              pos_indices, ndict, match_buffer,
                              &S_buffer, S.seq, S.length, n1) >= minscore)
			_report_match(n1 + 1, key_total_nchar);
	}
	// The SEXP returned by reported_matches_asSEXP() is UNPROTECTED
	// but you don't have to PROTECT it here since you are returning it
	// right away.
	return _reported_matches_asSEXP();
}

/*
 * --- .Call ENTRY POINT ---
 * XStringViews_match_WCP() arguments are assumed to be:
 *   wcp: the weighted clustered positions (WCP) object
 *   subject: an XStringViews object containing the subject sequence
 *   views_start, views_width: 2 integer vectors describing views on 'subject'.
 *   base_codes: named integer vector of length 4 obtained with
 *       xscodes(subject, baseOnly=TRUE)
 *   min_score: a single double (not NA)
 *   count_only: a single logical (not NA)
 */
SEXP XStringViews_match_WCP(SEXP wcp,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP min_score, SEXP count_only)
{
	int n1, n2, ndict, is_count_only, *pos_indices;
	int nviews, i, *view_start, *view_width, view_offset;
	int *clust_ends, *end, curr_end, prev_end;
	int nkey, max_nkey, max_key_nchar, key_total_nchar;
	int *key_nchars, *nchar, *key_order;
	int *match_buffer, **key_order_list;
	cachedCharSeq S, S_view;
	RoSeqs S_buffer;
	RoSeqsList key_seqs_list;
	SEXP dict_list, clust_bins, dict_elt, key;
	double **key_scores_list, minscore;

	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];

	dict_list = GET_SLOT(GET_SLOT(wcp, install("dictList")), install("listData"));
	clust_bins = GET_SLOT(GET_SLOT(wcp, install("clusters")), install("bins"));
	pos_indices = INTEGER(GET_SLOT(clust_bins, install("unlistData")));
	clust_ends =
		INTEGER(GET_SLOT(GET_SLOT(clust_bins, install("partitioning")),
				         install("end")));

	prev_end = 0;
	max_nkey = 0;
	max_key_nchar = 0;
	key_total_nchar = 0;
	ndict = LENGTH(dict_list);
	key_nchars = (int *) R_alloc((long) ndict, sizeof(int));
	key_seqs_list = _alloc_RoSeqsList(ndict);
	key_scores_list = (double **) R_alloc((long) ndict, sizeof(double *));
	key_order_list = (int **) R_alloc((long) ndict, sizeof(int *));
	for (n1 = 0, nchar = key_nchars, end = clust_ends; n1 < ndict;
	     n1++, nchar++, end++) {
		curr_end = *end;
		*nchar = curr_end - prev_end;
		key_total_nchar += *nchar;
		max_key_nchar = *nchar > max_key_nchar ? *nchar : max_key_nchar;
		prev_end = curr_end;

		dict_elt = VECTOR_ELT(dict_list, n1);
		key = GET_SLOT(dict_elt, install("key"));
		nkey = _get_XStringSet_length(key);
		max_nkey = nkey > max_nkey ? nkey : max_nkey;
		key_seqs_list.elts[n1] = _new_RoSeqs_from_XStringSet(nkey, key);
		key_scores_list[n1] =
			REAL(VECTOR_ELT(GET_SLOT(GET_SLOT(dict_elt, install("table")),
					                 install("listData")), 0));
		key_order_list[n1] = (int *) R_alloc((long) nkey, sizeof(int));
		for (n2 = 0, key_order = key_order_list[n1]; n2 < nkey;
		     n2++, key_order++)
			*key_order = n2;
	}
	S_buffer = _alloc_RoSeqs(1);
	S_buffer.elts[0].seq = (char *) R_alloc((long) max_key_nchar, sizeof(char));
	S_buffer.elts[0].length = max_key_nchar;
	match_buffer = (int *) R_alloc((long) max_nkey, sizeof(int));

	S = cache_XRaw(subject);
	nviews = LENGTH(views_start);
	_init_match_reporting(is_count_only ? mkString("COUNTONLY") : mkString("ASIRANGES"));
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S.seq + view_offset;
		S_view.length = *view_width;
		_shift_match_on_reporting(view_offset);
		for (n1 = 0, n2 = key_total_nchar; n2 <= S_view.length; n1++, n2++) {
			if (compute_wcp_score(&key_seqs_list, key_scores_list,
					              key_order_list,
					              key_total_nchar, key_nchars,
					              pos_indices, ndict, match_buffer,
					              &S_buffer,
					              S_view.seq, S_view.length, n1) >= minscore)
				_report_match(n1 + 1, key_total_nchar);
		}
	}
	// The SEXP returned by reported_matches_asSEXP() is UNPROTECTED
	// but you don't have to PROTECT it here since you are returning it
	// right away.
	return _reported_matches_asSEXP();
}
