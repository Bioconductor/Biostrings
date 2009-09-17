/****************************************************************************
 *                          Read/write fasta files                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <ctype.h> /* for isspace() */
#include <stdlib.h> /* for abs() */
#include <math.h> /* for log() */
#include <limits.h> /* for INT_MAX */

static int debug = 0;

SEXP debug_XStringSet_io()
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

/* This function should go somewhere else */
static int overflow_mult_int(int a, int b)
{
	a = abs(a);
	b = abs(b);
	return a != 0 && b != 0 && (log((double) a) + log((double) b) >= log((double) INT_MAX));
}

static int nfile;
static FILE **files;

static void open_files(SEXP paths)
{
	const char *path_str;

	nfile = 0;
	files = (FILE **) malloc(LENGTH(paths) * sizeof(FILE *));
	if (files == NULL)
		error("malloc() failed");
	for (nfile = 0; nfile < LENGTH(paths); nfile++) {
		path_str = translateChar(STRING_ELT(paths, nfile));
		if ((files[nfile] = fopen(path_str, "r")) == NULL)
			error("cannot open file '%s'", path_str);
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP io_cleanup()
{
	int i;

	for (i = 0; i < nfile; i++)
		fclose(files[i]);
	free(files);
	return R_NilValue;
}

#define LINEBUF_SIZE 20001
static char errmsg_buf[200];

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


/****************************************************************************
 * Read FASTA.
 */

/*
 * Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be able to eventually support longer sympols.
 */
static const char *FASTA_comment_markup = ";", *FASTA_desc_markup = ">";

/*
 * The LENGTHONLY storage handlers keep only the lengths of the description
 * lines and the sequences.
 */

static IntAE desc_lengths_buf, seq_lengths_buf;

static void add_desc_LENGTHONLY(int recno, const cachedCharSeq *dataline)
{
	IntAE_insert_at(&desc_lengths_buf, desc_lengths_buf.nelt, dataline->length);
	return;
}

static void add_empty_seq_LENGTHONLY(int recno)
{
	IntAE_insert_at(&seq_lengths_buf, seq_lengths_buf.nelt, 0);
	return;
}

static void append_to_last_seq_LENGTHONLY(const cachedCharSeq *dataline)
{
	seq_lengths_buf.elts[seq_lengths_buf.nelt - 1] += dataline->length;
	return;
}


/*
 * The CHARAEAE storage handlers use 2 Auto-Extending buffers to keep all the
 * data coming from the FASTA file in a single pass. The downside is that
 * the content of the 2 Auto-Extending buffers must then be turned into an
 * SEXP before it can be returned to the user.
 */

static CharAEAE descs_buf, seqs_buf;

static void add_desc_CHARAEAE(int recno, const cachedCharSeq *dataline)
{
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAEAE(&descs_buf, dataline->seq);
	return;
}

static void add_empty_seq_CHARAEAE(int recno)
{
	append_string_to_CharAEAE(&seqs_buf, "");
	return;
}

static void append_to_last_seq_CHARAEAE(const cachedCharSeq *dataline)
{
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAE(seqs_buf.elts + seqs_buf.nelt - 1, dataline->seq);
	return;
}


/*
 * Other storage handlers.
 */

static SEXP ans_names;

static void add_desc1(int recno, const cachedCharSeq *dataline)
{
	// This works only because dataline->seq is nul-terminated!
	SET_STRING_ELT(ans_names, recno, mkChar(dataline->seq));
	return;
}

/*
 * Return the number of sequences (1 sequence per FASTA record) in the file
 * (>= 0) or -1 if an error occured.
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like
 * in the original Pearson FASTA format.
 * This function is agnostic about how the data that are read are stored in
 * memory, how they will be returned to the user and how they will look to him.
 * This is delegated to 3 storage handlers: add_desc(), add_empty_seq() and
 * append_to_last_seq().
 */
static int parse_FASTA_file(FILE *stream,
		void (*add_desc)(int recno, const cachedCharSeq *dataline),
		void (*add_empty_seq)(int recno),
		void (*append_to_last_seq)(const cachedCharSeq *dataline))
{
	int FASTA_comment_markup_length, FASTA_desc_markup_length,
	    lineno, recno;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTA_comment_markup_length = strlen(FASTA_comment_markup);
	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	lineno = recno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrim(linebuf);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long", lineno);
			return -1;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		if (strncmp(linebuf, FASTA_comment_markup, FASTA_comment_markup_length) == 0)
			continue; // we ignore comment lines
		dataline.seq = linebuf;
		if (strncmp(linebuf, FASTA_desc_markup, FASTA_desc_markup_length) == 0) {
			if (add_desc != NULL) {
				dataline.seq += FASTA_desc_markup_length;
				dataline.length -= FASTA_desc_markup_length;
				add_desc(recno, &dataline);
			}
			if (add_empty_seq != NULL)
				add_empty_seq(recno);
			recno++;
		} else {
			if (recno == 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTA_desc_markup, lineno);
				return -1;
			}
			if (append_to_last_seq != NULL)
				append_to_last_seq(&dataline);
		}
	}
	return recno;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP fasta_info(SEXP filepath, SEXP use_descs)
{
	const char *path;
	FILE *stream;
	void (*add_desc)(int recno, const cachedCharSeq *dataline);
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
	if (parse_FASTA_file(stream, add_desc,
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
 SharedRaw_loadFASTA().
 Load a FASTA file into an SharedRaw object.

 Return a named list of 4 elements:

   - "nbyte" (single integer): number of bytes that were written to the SharedRaw
     object. SharedRaw_loadFASTA() starts to write at position 1 (the first byte)
     in the SharedRaw object. If nbyte_max is the length of the SharedRaw object, we
     should always have 1 <= nbyte <= nbyte_max. An error is raised if the
     FASTA file contains no data or if the size of the data exceeds the
     capacity (i.e. the length) of the SharedRaw object i.e. if SharedRaw_loadFASTA()
     tries to write more than nbyte_max bytes to the SharedRaw object (it does not
     try to resize it). 

   - "start" and "end": 2 integer vectors of same length L (L should always be
     >= 1) containing the start/end positions of the FASTA sequences relatively
     to the first byte of the SharedRaw object. The start/end positions define a set
     of views such that:
       (a) there is at least one view,
       (b) all views have a width >= 1 ('all(start <= end)' is TRUE),
       (c) the views are not overlapping and are sorted from left to right,
       (d) the leftmost view starts at position 1 (start[1] == 1) and the
           rightmost view ends at position nbyte (end[L] == nbyte).

   - "desc" (character vector): descriptions of the FASTA sequences. The length
     of 'desc' is L too.

 SharedRaw_loadFASTA() is designed to be very fast:

     library(Biostrings)
     filepath <- "~hpages/BSgenome_srcdata/UCSC/hg18/chr1.fa"
     filesize <- file.info(filepath)$size
     if (is.na(filesize))
         stop(filepath, ": file not found")
     filesize <- as.integer(filesize)
     if (is.na(filesize))
         stop(filepath, ": file is too big")
     rawptr <- SharedRaw(filesize)

     # Takes < 1 second on lamb1!
     .Call("SharedRaw_loadFASTA", rawptr@xp, filepath, "", NULL, PACKAGE="Biostrings")

     # or use safe wrapper:
     Biostrings:::SharedRaw.loadFASTA(rawptr, filepath)

   Compare to the time it takes to load Hsapiens$chr1 from the
   BSgenome.Hsapiens.UCSC.hg18 package: 1.340 second on lamb1!
   Conclusion: we could put and load directly the FASTA files in the
   BSgenome.* packages. The DNAString and XStringViews instances would be
   created on the fly. No need to store them in .rda files anymore!

*/

SEXP SharedRaw_loadFASTA(SEXP rawptr_xp, SEXP filepath, SEXP collapse, SEXP lkup)
{
	SEXP ans, ans_names, ans_elt, dest;
	const char *path, *coll;
	FILE *infile;
	long int lineno;
	char line[FASTALINE_MAX+1], desc[FASTALINE_MAX+1];
	int nbyte_max, gaplen, line_len, status, view_start, i1, i2;
	char c0;

	error("SharedRaw_loadFASTA() is not ready yet");
	dest = R_ExternalPtrTag(rawptr_xp);
	nbyte_max = LENGTH(dest);
	path = translateChar(STRING_ELT(filepath, 0));
	coll = CHAR(STRING_ELT(collapse, 0));
	gaplen = strlen(coll);

	if ((infile = fopen(path, "r")) == NULL)
		error("cannot open file");
	lineno = i1 = 0;
	status = 0; /* 0: expecting desc; 1: expecting seq; 2: no expectation */
	//_init_match_reporting(0);
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
			//_report_view(view_start, i1, desc);
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
	//_report_view(view_start, i1, desc);

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
	//PROTECT(ans_elt = _reported_match_starts_asINTEGER());
	//SET_ELEMENT(ans, 1, ans_elt);
	//UNPROTECT(1);
	/* set the "end" element */
	//PROTECT(ans_elt = _reported_match_ends_asINTEGER());
	//SET_ELEMENT(ans, 2, ans_elt);
	//UNPROTECT(1);
	/* set the "desc" element */
	//PROTECT(ans_elt = _reported_view_names_asCHARACTER());
	//SET_ELEMENT(ans, 3, ans_elt);
	//UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Read FASTQ.
 */

static const char *FASTQ_line1_markup = "@", *FASTQ_line3_markup = "+";

/*
 * The WIDTHONLY storage handlers keep only the common width of the sequences
 * (NA if the set of sequences is not rectangular).
 */

static int FASTQ_width;

static void add_seq_WIDTHONLY(int recno, const cachedCharSeq *dataline)
{
	if (recno == 0) {
		FASTQ_width = dataline->length;
		return;
	}
	if (FASTQ_width == NA_INTEGER || dataline->length == FASTQ_width)
		return;
	FASTQ_width = NA_INTEGER;
	return;
}

static SEXP FASTQ_seqbuf, FASTQ_seqbuf_shared, FASTQ_qualbuf, FASTQ_qualbuf_shared;

static void append_seq_to_FASTQ_seqbuf(int recno, const cachedCharSeq *dataline)
{
	const ByteTrTable *byte2code;

	byte2code = get_enc_byte2code("DNAString");
	_write_RoSeq_to_SharedRaw(FASTQ_seqbuf_shared, recno * FASTQ_width, dataline, byte2code);
	return;
}

static void append_qual_to_FASTQ_qualbuf(int recno, const cachedCharSeq *dataline)
{
	_write_RoSeq_to_SharedRaw(FASTQ_qualbuf_shared, recno * FASTQ_width, dataline, NULL);
	return;
}

/*
 * Return the number of sequences (1 sequence per FASTQ record) in the file
 * (>= 0) or -1 if an error occured. Ignore empty lines.
 * This function is agnostic about how the data that are read are stored in
 * memory, how they will be returned to the user and how they will look to him.
 * This is delegated to 3 storage handlers: add_desc(), add_empty_seq() and
 * append_to_last_seq().
 */
static int parse_FASTQ_file(FILE *stream,
		void (*add_seqid)(int recno, const cachedCharSeq *dataline),
		void (*add_seq)(int recno, const cachedCharSeq *dataline),
		void (*add_qualid)(int recno, const cachedCharSeq *dataline),
		void (*add_qual)(int recno, const cachedCharSeq *dataline))
{
	int FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineno, recno, lineinrecno;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineno = recno = lineinrecno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrim(linebuf);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long", lineno);
			return -1;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		dataline.seq = linebuf;
		lineinrecno++;
		if (lineinrecno > 4)
			lineinrecno = 1;
		switch (lineinrecno) {
		    case 1:
			if (strncmp(linebuf, FASTQ_line1_markup, FASTQ_line1_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line1_markup, lineno);
				return -1;
			}
			if (add_seqid != NULL) {
				dataline.seq += FASTQ_line1_markup_length;
				dataline.length -= FASTQ_line1_markup_length;
				add_seqid(recno, &dataline);
			}
		    break;
		    case 2:
			if (add_seq != NULL)
				add_seq(recno, &dataline);
		    break;
		    case 3:
			if (strncmp(linebuf, FASTQ_line3_markup, FASTQ_line3_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line3_markup, lineno);
				return -1;
			}
			if (add_qualid != NULL) {
				dataline.seq += FASTQ_line3_markup_length;
				dataline.length -= FASTQ_line3_markup_length;
				add_qualid(recno, &dataline);
			}
		    break;
		    case 4:
			if (add_qual != NULL)
				add_qual(recno, &dataline);
			recno++;
		    break;
		}
	}
	return recno;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP fastq_geometry(SEXP filepath)
{
	int fileno;
	int totalnseq, nseq;
	SEXP ans;

	open_files(filepath);
	totalnseq = 0;
	for (fileno = 0; fileno < nfile; fileno++) {
		nseq = parse_FASTQ_file(files[fileno], NULL, add_seq_WIDTHONLY, NULL, NULL);
		if (nseq < 0)
			error("reading FASTQ file %s: %s", STRING_ELT(filepath, fileno), errmsg_buf);
		totalnseq += nseq;
	}
	PROTECT(ans = NEW_INTEGER(2));
	INTEGER(ans)[0] = totalnseq;
	INTEGER(ans)[1] = FASTQ_width;
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP read_fastq(SEXP filepath, SEXP drop_quality)
{
	SEXP ans_geom, ans;
	int fileno, buf_length;

	PROTECT(ans_geom = fastq_geometry(filepath));
	if (INTEGER(ans_geom)[1] == NA_INTEGER)
		error("read_fastq(): FASTQ files with variable sequence lengths are not supported yet");
	if (overflow_mult_int(INTEGER(ans_geom)[0], INTEGER(ans_geom)[1]))
		error("read_fastq(): FASTQ files contain more data an XStringSet object can hold, sorry!");
	buf_length = INTEGER(ans_geom)[0] * INTEGER(ans_geom)[1];
	PROTECT(FASTQ_seqbuf = _alloc_XString("DNAString", buf_length));
	FASTQ_seqbuf_shared = get_XVector_shared(FASTQ_seqbuf);
	if (!LOGICAL(drop_quality)[0]) {
		PROTECT(FASTQ_qualbuf = _alloc_XString("BString", buf_length));
		FASTQ_qualbuf_shared = get_XVector_shared(FASTQ_qualbuf);
	}
	for (fileno = 0; fileno < nfile; fileno++) {
		rewind(files[fileno]);
		parse_FASTQ_file(files[fileno],
			NULL, append_seq_to_FASTQ_seqbuf,
			NULL, LOGICAL(drop_quality)[0] ? NULL : append_qual_to_FASTQ_qualbuf);
	}
	if (!LOGICAL(drop_quality)[0]) {
		PROTECT(ans = NEW_LIST(3));
		SET_ELEMENT(ans, 0, ans_geom);
		SET_ELEMENT(ans, 1, FASTQ_seqbuf);
		SET_ELEMENT(ans, 2, FASTQ_qualbuf);
		UNPROTECT(4);
	} else {
		PROTECT(ans = NEW_LIST(2));
		SET_ELEMENT(ans, 0, ans_geom);
		SET_ELEMENT(ans, 1, FASTQ_seqbuf);
		UNPROTECT(3);
	}
	return ans;
}

