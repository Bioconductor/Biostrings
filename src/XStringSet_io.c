/****************************************************************************
 *                       Read/write FASTA/FASTQ files                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#include <ctype.h> /* for isspace() */
#include <stdlib.h> /* for abs() */
#include <math.h> /* for log() */
#include <limits.h> /* for INT_MAX */

#include <sys/types.h>
#include <sys/stat.h> /* for fstat() */

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

/* Code taken from do_url() in R/src/main/connections.c */
static void get_file_ztype(const char *path, int *ztype, int *subtype)
{
	FILE *fp;
	char buf[7];
	int res;

	*ztype = -1;
	*subtype = 0;
	if ((fp = fopen(path, "rb")) == NULL)
		return;
	memset(buf, 0, 7);
	res = fread(buf, 5, 1, fp);
	fclose(fp);
	if (res != 1)
		return;
	if (buf[0] == '\x1f' && buf[1] == '\x8b')
		*ztype = 0;
	else if (strncmp(buf, "BZh", 3) == 0)
		*ztype = 1;
	else if (buf[0] == '\xFD' && strncmp(buf+1, "7zXZ", 4) == 0)
		*ztype = 2;
	else if ((buf[0] == '\xFF') && strncmp(buf+1, "LZMA", 4) == 0) {
		*ztype = 2;
		*subtype = 1;
	} else if (memcmp(buf, "]\0\0\200\0", 5) == 0) {
		*ztype = 2;
		*subtype = 1;
	}
	return;
}

static int ninputfiles;

static struct inputfile {
	FILE *fp;
	//gzFile gfp;
} *inputfiles;

static void open_inputfiles(SEXP filepath)
{
	SEXP filepath_elt;
	const char *path;
	int ret, ztype, subtype;
	struct stat buf;

	if (!IS_CHARACTER(filepath))
		error("'filepath' must be a character vector");
	ninputfiles = 0;
	inputfiles = (struct inputfile *) malloc(LENGTH(filepath) * sizeof(struct inputfile));
	if (inputfiles == NULL)
		error("malloc() failed in open_inputfiles()");
	while (ninputfiles < LENGTH(filepath)) {
		filepath_elt = STRING_ELT(filepath, ninputfiles);
		if (filepath_elt == NA_STRING)
			error("'filepath' contains NAs");
		path = R_ExpandFileName(translateChar(filepath_elt));
		get_file_ztype(path, &ztype, &subtype);
		switch (ztype) {
		/* No compression */
		    case -1:
			inputfiles[ninputfiles].fp = fopen(path, "r");
			if (inputfiles[ninputfiles].fp == NULL)
				error("cannot open file '%s'", path);
			ninputfiles++;
			ret = fstat(fileno(inputfiles[ninputfiles].fp), &buf);
			if (ret != 0)
				error("Biostrings internal error in open_inputfiles(): ",
				      "cannot stat file '%s'", path);
			if (S_ISDIR(buf.st_mode))
				error("file '%s' is a directory", path);
		    break;
		/* gzfile */
		    case 0:
			error("cannot open file '%s' (gzip-compressed files "
			      "are not supported yet, sorry!)", path);
			//open the file with gzFile gfp = gzopen(path, "r");
			//close it with gzclose();
		    break;
		/* bzfile */
		    case 1:
			error("cannot open file '%s' (bzip2-compressed files "
			      "are not supported yet, sorry!)", path);
			//requires #include <bzlib.h>
			//open the file with BZFILE* bfp = BZ2_bzReadOpen(...);
			//close it with BZ2_bzReadClose();
		    break;
		/* xzfile */
		    case 2:
			error("cannot open file '%s' (LZMA-compressed files "
			      "are not supported yet, sorry!)", path);
			//requires #include <lzma.h>
			//opening/closing this type of file seems quite complicated
		    break;
		    default:
			error("Biostrings internal error in open_files(): ",
			      "invalid ztype value %d", ztype);
		}
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP io_cleanup()
{
	int fn;

	for (fn = 0; fn < ninputfiles; fn++)
		fclose(inputfiles[fn].fp);
	free(inputfiles);
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
 * The LENGTHONLY storage handlers only store the lengths of the sequences.
 */

static IntAE seq_lengths_buf;

static void add_empty_seq_LENGTHONLY(int recno)
{
	IntAE_insert_at(&seq_lengths_buf, seq_lengths_buf.nelt, 0);
	return;
}

static void append_to_last_seq_LENGTHONLY(int recno, const cachedCharSeq *dataline)
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

static CharAEAE descs_buf;

static void add_desc_CHARAEAE(int recno, const cachedCharSeq *dataline)
{
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAEAE(&descs_buf, dataline->seq);
	return;
}


/*
 * The FASTA_seqbuf storage handlers.
 */

static cachedXVectorList FASTA_seqbuf;
static cachedCharSeq FASTA_seqbuf_elt;
static const int *FASTA_lkup;
static int FASTA_lkup_length;

static void add_empty_seq_to_FASTA_seqbuf(int recno)
{
	FASTA_seqbuf_elt = get_cachedXRawList_elt(&FASTA_seqbuf, recno);
	FASTA_seqbuf_elt.length = 0;
	return;
}

static void append_to_last_seq_in_FASTA_seqbuf(int recno, const cachedCharSeq *dataline)
{
	int i1;

	i1 = FASTA_seqbuf_elt.length;
	FASTA_seqbuf_elt.length += dataline->length;
	/* FASTA_seqbuf_elt.seq is a const char * so we need to cast it to
           char * before we can write to it */
	Ocopy_bytes_to_i1i2_with_lkup(i1, FASTA_seqbuf_elt.length - 1,
		(char *) FASTA_seqbuf_elt.seq, FASTA_seqbuf_elt.length,
		dataline->seq, dataline->length,
		FASTA_lkup, FASTA_lkup_length);
	return;
}


/*
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
 * the original Pearson FASTA format.
 * This function is agnostic about how the data that are read are stored in
 * memory, how they will be returned to the user and how they will look to him.
 * This is delegated to 3 storage handlers: add_desc(), add_empty_seq() and
 * append_to_last_seq().
 */
static const char *parse_FASTA_file(FILE *stream, int *recno,
		void (*add_desc)(int recno, const cachedCharSeq *dataline),
		void (*add_empty_seq)(int recno),
		void (*append_to_last_seq)(int recno, const cachedCharSeq *dataline))
{
	int FASTA_comment_markup_length, FASTA_desc_markup_length, lineno;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTA_comment_markup_length = strlen(FASTA_comment_markup);
	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	lineno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrim(linebuf);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long", lineno);
			return errmsg_buf;
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
				add_desc(*recno, &dataline);
			}
			if (add_empty_seq != NULL)
				add_empty_seq(*recno);
			(*recno)++;
		} else {
			if (*recno == 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTA_desc_markup, lineno);
				return errmsg_buf;
			}
			if (append_to_last_seq != NULL)
				append_to_last_seq(*recno, &dataline);
		}
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fasta_info(SEXP filepath, SEXP use_descs)
{
	void (*add_desc)(int recno, const cachedCharSeq *dataline);
	int fn, recno;
	SEXP ans, ans_names;
	RoSeqs descs;
	const char *errmsg;

	seq_lengths_buf = new_IntAE(0, 0, 0);
	if (LOGICAL(use_descs)[0]) {
		add_desc = &add_desc_CHARAEAE;
		descs_buf = new_CharAEAE(0, 0);
	} else {
		add_desc = NULL;
	}
	open_inputfiles(filepath);
	for (fn = recno = 0; fn < ninputfiles; fn++) {
		errmsg = parse_FASTA_file(inputfiles[fn].fp, &recno,
				add_desc,
				&add_empty_seq_LENGTHONLY,
				&append_to_last_seq_LENGTHONLY);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      STRING_ELT(filepath, fn), errmsg_buf);
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

/* --- .Call ENTRY POINT --- */
SEXP read_fasta_in_XStringSet(SEXP filepath, SEXP set_names,
		SEXP elementType, SEXP lkup)
{
	SEXP ans, ans_width, ans_names;
	const char *element_type;
	char classname[40]; // longest string will be "DNAStringSet"
	int fn, recno;

	PROTECT(ans_width = fasta_info(filepath, set_names));
	PROTECT(ans_names = GET_NAMES(ans_width));
	SET_NAMES(ans_width, R_NilValue);
	element_type = CHAR(STRING_ELT(elementType, 0));
	if (snprintf(classname, sizeof(classname), "%sSet", element_type)
	    >= sizeof(classname))
		error("Biostrings internal error in "
		      "read_fasta_in_XStringSet(): 'elementType' too long");
	PROTECT(ans = alloc_XRawList(classname, element_type, ans_width));
	_set_XStringSet_names(ans, ans_names);
	FASTA_seqbuf = cache_XVectorList(ans);
	if (lkup == R_NilValue) {
		FASTA_lkup = NULL;
	} else {
		FASTA_lkup = INTEGER(lkup);
		FASTA_lkup_length = LENGTH(lkup);
	}
	for (fn = recno = 0; fn < ninputfiles; fn++) {
		rewind(inputfiles[fn].fp);
		parse_FASTA_file(inputfiles[fn].fp, &recno,
			NULL,
			&add_empty_seq_to_FASTA_seqbuf,
			&append_to_last_seq_in_FASTA_seqbuf);
	}
	UNPROTECT(3);
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
	Ocopy_cachedCharSeq_to_SharedRaw_offset(
		FASTQ_seqbuf_shared, recno * FASTQ_width,
		dataline,
		*byte2code, BYTETRTABLE_LENGTH);
	return;
}

static void append_qual_to_FASTQ_qualbuf(int recno, const cachedCharSeq *dataline)
{
	Ocopy_cachedCharSeq_to_SharedRaw_offset(
		FASTQ_qualbuf_shared, recno * FASTQ_width,
		dataline,
		NULL, 0);
	return;
}

/*
 * Ignore empty lines.
 * This function is agnostic about how the data that are read are stored in
 * memory, how they will be returned to the user and how they will look to him.
 * This is delegated to 4 storage handlers: add_seqid(), add_seq(), add_qualid()
 * and add_qual().
 */
static const char *parse_FASTQ_file(FILE *stream, int *recno,
		void (*add_seqid)(int recno, const cachedCharSeq *dataline),
		void (*add_seq)(int recno, const cachedCharSeq *dataline),
		void (*add_qualid)(int recno, const cachedCharSeq *dataline),
		void (*add_qual)(int recno, const cachedCharSeq *dataline))
{
	int FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineno, lineinrecno;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineno = lineinrecno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrim(linebuf);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long", lineno);
			return errmsg_buf;
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
				return errmsg_buf;
			}
			if (add_seqid != NULL) {
				dataline.seq += FASTQ_line1_markup_length;
				dataline.length -= FASTQ_line1_markup_length;
				add_seqid(*recno, &dataline);
			}
		    break;
		    case 2:
			if (add_seq != NULL)
				add_seq(*recno, &dataline);
		    break;
		    case 3:
			if (strncmp(linebuf, FASTQ_line3_markup, FASTQ_line3_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (add_qualid != NULL) {
				dataline.seq += FASTQ_line3_markup_length;
				dataline.length -= FASTQ_line3_markup_length;
				add_qualid(*recno, &dataline);
			}
		    break;
		    case 4:
			if (add_qual != NULL)
				add_qual(*recno, &dataline);
			(*recno)++;
		    break;
		}
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fastq_geometry(SEXP filepath)
{
	SEXP ans;
	int fn, recno;
	const char *errmsg;

	FASTQ_width = NA_INTEGER;
	open_inputfiles(filepath);
	for (fn = recno = 0; fn < ninputfiles; fn++) {
		errmsg = parse_FASTQ_file(inputfiles[fn].fp, &recno,
				NULL, add_seq_WIDTHONLY,
				NULL, NULL);
		if (errmsg != NULL)
			error("reading FASTQ file %s: %s",
			      STRING_ELT(filepath, fn), errmsg_buf);
	}
	PROTECT(ans = NEW_INTEGER(2));
	INTEGER(ans)[0] = recno;
	INTEGER(ans)[1] = FASTQ_width;
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP read_fastq(SEXP filepath, SEXP drop_quality)
{
	SEXP ans_geom, ans;
	int fn, recno, buf_length;

	PROTECT(ans_geom = fastq_geometry(filepath));
	if (INTEGER(ans_geom)[0] == 0) {
		buf_length = 0;
	} else {
		if (INTEGER(ans_geom)[1] == NA_INTEGER)
			error("read_fastq(): FASTQ files with variable sequence "
			      "lengths are not supported yet");
		if (overflow_mult_int(INTEGER(ans_geom)[0], INTEGER(ans_geom)[1]))
			error("read_fastq(): FASTQ files contain more data an "
			      "XStringSet object can hold, sorry!");
		buf_length = INTEGER(ans_geom)[0] * INTEGER(ans_geom)[1];
	}
	PROTECT(FASTQ_seqbuf = alloc_XRaw("DNAString", buf_length));
	FASTQ_seqbuf_shared = get_XVector_shared(FASTQ_seqbuf);
	if (!LOGICAL(drop_quality)[0]) {
		PROTECT(FASTQ_qualbuf = alloc_XRaw("BString", buf_length));
		FASTQ_qualbuf_shared = get_XVector_shared(FASTQ_qualbuf);
	}
	for (fn = recno = 0; fn < ninputfiles; fn++) {
		rewind(inputfiles[fn].fp);
		parse_FASTQ_file(inputfiles[fn].fp, &recno,
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

