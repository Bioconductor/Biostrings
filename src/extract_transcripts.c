#include "Biostrings.h"
#include "IRanges_interface.h"

static SEXP make_ans_width(const cachedCharSeq *X, SEXP exonStarts, SEXP exonEnds)
{
	int ntranscripts, i, nstarts, nends, cumwidth, j, start, end, width;
	SEXP ans_width, starts, ends;

	ntranscripts = LENGTH(exonStarts);
	PROTECT(ans_width = NEW_INTEGER(ntranscripts));
	for (i = 0; i < ntranscripts; i++) {
		starts = VECTOR_ELT(exonStarts, i);
		if (starts == R_NilValue) {
			nstarts = 0;
		} else if (IS_INTEGER(starts)) {
			nstarts = LENGTH(starts);
		} else {
			UNPROTECT(1);
			error("'exonStarts' has invalid elements");
		}
		ends = VECTOR_ELT(exonEnds, i);
		if (ends == R_NilValue) {
			nends = 0;
		} else if (IS_INTEGER(ends)) {
			nends = LENGTH(ends);
		} else {
			UNPROTECT(1);
			error("'exonEnds' has invalid elements");
		}
		if (nstarts != nends) {
			UNPROTECT(1);
			error("'exonStarts' and 'exonEnds' are not compatible");
		}
		cumwidth = 0;
		for (j = 0; j < nstarts; j++) {
			start = INTEGER(starts)[j];
			if (start < 1 || start > X->length + 1) {
				UNPROTECT(1);
				error("'exonStarts' contains \"out of limits\" values");
			}
			end = INTEGER(ends)[j];
			if (end < 0 || end > X->length) {
				UNPROTECT(1);
				error("'exonEnds' contains \"out of limits\" values");
			}
			width = end - start + 1;
			if (width < 0) {
				UNPROTECT(1);
				error("'exonStarts/exonEnds' define exons "
				      "of negative length");
			}
			cumwidth += width;
		}
		INTEGER(ans_width)[i] = cumwidth;
	}
	UNPROTECT(1);
	return ans_width;
}


/*
 * --- .Call ENTRY POINT ---
 */
SEXP extract_transcripts(SEXP x, SEXP exonStarts, SEXP exonEnds, SEXP strand,
		SEXP reorder_exons_on_minus_strand, SEXP lkup)
{
	cachedCharSeq X, Y;
	SEXP ans_width, ans, starts, ends, strand_elt;
	cachedXVectorList cached_ans;
	int ans_length, i, nexons, cumwidth, j, start, end, width;

	X = cache_XRaw(x);
	PROTECT(ans_width = make_ans_width(&X, exonStarts, exonEnds));
	PROTECT(ans = alloc_XRawList("DNAStringSet", "DNAString", ans_width));
	cached_ans = cache_XVectorList(ans);
	ans_length = get_cachedXVectorList_length(&cached_ans);
	for (i = 0; i < ans_length; i++) {
		starts = VECTOR_ELT(exonStarts, i);
		if (starts == R_NilValue || (nexons = LENGTH(starts)) == 0)
			continue;
		ends = VECTOR_ELT(exonEnds, i);
		strand_elt = STRING_ELT(strand, i);
		Y = get_cachedXRawList_elt(&cached_ans, i);
		cumwidth = 0;
		for (j = 0; j < nexons; j++) {
			if (CHAR(strand_elt)[0] != '+'
			 && LOGICAL(reorder_exons_on_minus_strand)[0]) {
				start = INTEGER(starts)[nexons - 1 - j];
				end = INTEGER(ends)[nexons - 1 - j];
			} else {
				start = INTEGER(starts)[j];
				end = INTEGER(ends)[j];
			}
			start--; end--;
			width = end - start + 1;
			/* Y.seq is a const char * so we need to cast it to char *
                           before we can write to it */
			if (CHAR(strand_elt)[0] == '+') {
				Ocopy_bytes_from_i1i2_with_lkup(start, end,
					(char *) Y.seq + cumwidth, width,
					X.seq, X.length,
					NULL, 0);
			} else {
				Orevcopy_bytes_from_i1i2_with_lkup(start, end,
					(char *) Y.seq + cumwidth, width,
					X.seq, X.length,
					INTEGER(lkup), LENGTH(lkup));
			}
			cumwidth += width;
		}
	}
	UNPROTECT(2);
	return ans;
}

