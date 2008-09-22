/****************************************************************************
 *                          Read/write fasta files                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <ctype.h> /* for isspace() */

static int debug = 0;

SEXP debug_fasta_io()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

#define LINEBUF_SIZE 20001
static char errmsg_buf[200];

/*
 * Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be able to eventually support longer sympols.
 */
static const char *comment_markup = ";", *desc_markup = ">";

/*
 * Return the number of chars that remain in the buffer after we've removed
 * the right spaces ('\n', '\r', '\t', ' ', etc...).
 */
static int rtrim(char *linebuf)
{
	int i;

	i = strlen(linebuf) - 1;
	while (i >= 0 && isspace(linebuf[i])) i--;
	linebuf[++i] = 0;
	return i;
}

/*
 * The LENGTHONLY storage handlers keep only the lengths of the description
 * lines and the sequences.
 */

static IntAE desc_lengths_buf, seq_lengths_buf;

static void add_desc_LENGTHONLY(int seqno, const RoSeq *dataline)
{
	IntAE_insert_at(&desc_lengths_buf, desc_lengths_buf.nelt, dataline->nelt);
	return;
}

static void add_empty_seq_LENGTHONLY(int seqno)
{
	IntAE_insert_at(&seq_lengths_buf, seq_lengths_buf.nelt, 0);
	return;
}

static void append_to_last_seq_LENGTHONLY(const RoSeq *dataline)
{
	seq_lengths_buf.elts[seq_lengths_buf.nelt - 1] += dataline->nelt;
	return;
}


/*
 * The CHARAEAE storage handlers use 2 Auto-Extending buffers to keep all the
 * data coming from the FASTA file in a single pass. The downside is that
 * the content of the 2 Auto-Extending buffers must then be turned into an
 * SEXP before it can be returned to the user.
 */

static CharAEAE descs_buf, seqs_buf;

static void add_desc_CHARAEAE(int seqno, const RoSeq *dataline)
{
	// This works only because dataline->elts is nul-terminated!
	append_string_to_CharAEAE(&descs_buf, dataline->elts);
	return;
}

static void add_empty_seq_CHARAEAE(int seqno)
{
	append_string_to_CharAEAE(&seqs_buf, "");
	return;
}

static void append_to_last_seq_CHARAEAE(const RoSeq *dataline)
{
	// This works only because dataline->elts is nul-terminated!
	append_string_to_CharAE(seqs_buf.elts + seqs_buf.nelt - 1, dataline->elts);
	return;
}


/*
 * Other storage handlers.
 */

static SEXP ans_names;

static void add_desc1(int seqno, const RoSeq *dataline)
{
	// This works only because dataline->elts is nul-terminated!
	SET_STRING_ELT(ans_names, seqno, mkChar(dataline->elts));
	return;
}

/*
 * Return the number of sequences in the file (>= 0) or -1 if an error occured.
 * Ignore empty lines and lines starting with 'comment_markup' like
 * in the original Pearson FASTA format.
 * This function is agnostic about how the data that are read are stored in
 * memory, how they will be returned to the user and how they will look to him.
 * This is delegated to 3 storage handlers: add_desc(), add_empty_seq() and
 * append_to_last_seq().
 */
static int read_fasta_file(FILE *stream,
		void (*add_desc)(int seqno, const RoSeq *dataline),
		void (*add_empty_seq)(int seqno),
		void (*append_to_last_seq)(const RoSeq *dataline))
{
	int comment_markup_length, desc_markup_length, lineno, seqno;
	char linebuf[LINEBUF_SIZE];
	RoSeq dataline;

	comment_markup_length = strlen(comment_markup);
	desc_markup_length = strlen(desc_markup);
	lineno = seqno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.nelt = rtrim(linebuf);
		// dataline.nelt > LINEBUF_SIZE - 1 should never happen
		if (dataline.nelt >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long", lineno);
			return -1;
		}
		if (dataline.nelt == 0)
			continue; // we ignore empty lines
		if (strncmp(linebuf, comment_markup, comment_markup_length) == 0)
			continue; // we ignore comment lines
		dataline.elts = linebuf;
		if (strncmp(linebuf, desc_markup, desc_markup_length) == 0) {
			if (add_desc != NULL) {
				dataline.elts += desc_markup_length;
				dataline.nelt -= desc_markup_length;
				add_desc(seqno, &dataline);
			}
			if (add_empty_seq != NULL)
				add_empty_seq(seqno);
			seqno++;
		} else {
			if (seqno == 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 desc_markup, lineno);
				return -1;
			}
			if (append_to_last_seq != NULL)
				append_to_last_seq(&dataline);
		}
	}
	return seqno;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP fasta_info(SEXP filepath, SEXP use_descs)
{
	const char *path;
	FILE *stream;
	void (*add_desc)(int seqno, const RoSeq *dataline);
	RoSeqs descs;
	SEXP ans, ans_names;

	path = translateChar(STRING_ELT(filepath, 0));
	if ((stream = fopen(path, "r")) == NULL)
		error("cannot open file '%s'", path);
	if (LOGICAL(use_descs)[0]) {
		add_desc = &add_desc_CHARAEAE;
		descs_buf = new_CharAEAE(0, 0);
	} else {
		add_desc = NULL;
	}
	seq_lengths_buf = new_IntAE(0, 0, 0);
	if (read_fasta_file(stream, add_desc,
			&add_empty_seq_LENGTHONLY,
			&append_to_last_seq_LENGTHONLY) < 0)
	{
		error("%s", errmsg_buf);
	}
	PROTECT(ans = IntAE_asINTEGER(&seq_lengths_buf));
	if (LOGICAL(use_descs)[0]) {
		descs = _new_RoSeqs_from_CharAEAE(&descs_buf);
		PROTECT(ans_names = _new_STRSXP_from_RoSeqs(&descs, R_NilValue));
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


#define FASTALINE_MAX 20000

/*
 RawPtr_loadFASTA().
 Load a FASTA file into an RawPtr object.

 Return a named list of 4 elements:

   - "nbyte" (single integer): number of bytes that were written to the RawPtr
     object. RawPtr_loadFASTA() starts to write at position 1 (the first byte)
     in the RawPtr object. If nbyte_max is the length of the RawPtr object, we
     should always have 1 <= nbyte <= nbyte_max. An error is raised if the
     FASTA file contains no data or if the size of the data exceeds the
     capacity (i.e. the length) of the RawPtr object i.e. if RawPtr_loadFASTA()
     tries to write more than nbyte_max bytes to the RawPtr object (it does not
     try to resize it). 

   - "start" and "end": 2 integer vectors of same length L (L should always be
     >= 1) containing the start/end positions of the FASTA sequences relatively
     to the first byte of the RawPtr object. The start/end positions define a set
     of views such that:
       (a) there is at least one view,
       (b) all views have a width >= 1 ('all(start <= end)' is TRUE),
       (c) the views are not overlapping and are sorted from left to right,
       (d) the leftmost view starts at position 1 (start[1] == 1) and the
           rightmost view ends at position nbyte (end[L] == nbyte).

   - "desc" (character vector): descriptions of the FASTA sequences. The length
     of 'desc' is L too.

 RawPtr_loadFASTA() is designed to be very fast:

     library(Biostrings)
     filepath <- "~hpages/BSgenome_srcdata/UCSC/hg18/chr1.fa"
     filesize <- file.info(filepath)$size
     if (is.na(filesize))
         stop(filepath, ": file not found")
     filesize <- as.integer(filesize)
     if (is.na(filesize))
         stop(filepath, ": file is too big")
     rawptr <- RawPtr(filesize)

     # Takes < 1 second on lamb1!
     .Call("RawPtr_loadFASTA", rawptr@xp, filepath, "", NULL, PACKAGE="Biostrings")

     # or use safe wrapper:
     Biostrings:::RawPtr.loadFASTA(rawptr, filepath)

   Compare to the time it takes to load Hsapiens$chr1 from the
   BSgenome.Hsapiens.UCSC.hg18 package: 1.340 second on lamb1!
   Conclusion: we could put and load directly the FASTA files in the
   BSgenome.* packages. The DNAString and XStringViews instances would be
   created on the fly. No need to store them in .rda files anymore!

*/

SEXP RawPtr_loadFASTA(SEXP rawptr_xp, SEXP filepath, SEXP collapse, SEXP lkup)
{
	SEXP ans, ans_names, ans_elt, dest;
	const char *path, *coll;
	FILE *infile;
	long int lineno;
	char line[FASTALINE_MAX+1], desc[FASTALINE_MAX+1];
	int nbyte_max, gaplen, line_len, status, view_start, i1, i2;
	char c0;

	dest = R_ExternalPtrTag(rawptr_xp);
	nbyte_max = LENGTH(dest);
	path = translateChar(STRING_ELT(filepath, 0));
	coll = CHAR(STRING_ELT(collapse, 0));
	gaplen = strlen(coll);

	if ((infile = fopen(path, "r")) == NULL)
		error("cannot open file");
	lineno = i1 = 0;
	status = 0; /* 0: expecting desc; 1: expecting seq; 2: no expectation */
	_init_match_reporting(0);
	while ((line_len = fgets_rtrimmed(line, FASTALINE_MAX+1, infile)) != -1) {
	/* while (fgets(line, FASTALINE_MAX+1, infile) != NULL) { */
		lineno++;
		if (line_len >= FASTALINE_MAX) {
			fclose(infile);
			error("file contains lines that are too long");
		}
		if (line_len == 0)
			continue;
		c0 = line[0];
		if (c0 == ';')
			continue;
		if (c0 != '>') {
			if (status == 0) {
				fclose(infile);
				error("file does not seem to be FASTA");
			}
			i2 = i1 + line_len - 1;
			if (lkup == R_NilValue)
				IRanges_memcpy_to_i1i2(i1, i2,
					(char *) RAW(dest), nbyte_max,
					line, line_len, sizeof(char));
			else
				IRanges_charcpy_to_i1i2_with_lkup(i1, i2,
					(char *) RAW(dest), nbyte_max,
					line, line_len,
					INTEGER(lkup), LENGTH(lkup));
			i1 = i2 + 1;
			status = 2;
			continue;
		}
		if (status == 1) {
			fclose(infile);
			error("file does not seem to be FASTA");
		}
		if (status == 2) {
			_report_view(view_start, i1, desc);
			if (gaplen != 0) {
				i2 = i1 + gaplen - 1;
				IRanges_memcpy_to_i1i2(i1, i2,
					(char *) RAW(dest), nbyte_max,
					coll, gaplen, sizeof(char));
				i1 = i2 + 1;
			}
		}
		view_start = i1 + 1;
		strcpy(desc, line + 1);
		status = 1;
	}
	fclose(infile);
	if (status != 2)
		error("file does not seem to be FASTA");
	_report_view(view_start, i1, desc);

	PROTECT(ans = NEW_LIST(4));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(4));
	SET_STRING_ELT(ans_names, 0, mkChar("nbyte"));
	SET_STRING_ELT(ans_names, 1, mkChar("start"));
	SET_STRING_ELT(ans_names, 2, mkChar("end"));
	SET_STRING_ELT(ans_names, 3, mkChar("desc"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "nbyte" element */
	PROTECT(ans_elt = NEW_INTEGER(1));
	INTEGER(ans_elt)[0] = i1;
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = _reported_match_starts_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = _reported_match_ends_asINTEGER());
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);
	/* set the "desc" element */
	PROTECT(ans_elt = _reported_view_names_asCHARACTER());
	SET_ELEMENT(ans, 3, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

