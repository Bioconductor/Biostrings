/****************************************************************************
 *                       Read/write FASTA/FASTQ files                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

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

#define LINEBUF_SIZE 20001
static char errmsg_buf[200];



/****************************************************************************
 *  A. FASTA FORMAT                                                         *
 ****************************************************************************/

/*
 * Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be able to eventually support longer sympols.
 */
static const char *FASTA_comment_markup = ";", *FASTA_desc_markup = ">";


/****************************************************************************
 * Reading FASTA files.
 */

typedef struct fasta_loader {
	const int *lkup;
	int lkup_length;
	void (*load_desc_line)(struct fasta_loader *loader,
			       const cachedCharSeq *dataline);
	void (*load_empty_seq)(struct fasta_loader *loader);
	void (*load_seq_line)(struct fasta_loader *loader,
			      const cachedCharSeq *dataline);
	int nrec;
	void *ext;  /* loader extension (optional) */
} FASTAloader;

/*
 * The FASTAINFO loader only loads the lengths (and optionally the names)
 * of the sequences.
 */

typedef struct fastainfo_loader_ext {
	CharAEAE ans_names_buf;
	IntAE seqlengths_buf;
} FASTAINFO_loaderExt;

static void FASTAINFO_load_desc_line(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	FASTAINFO_loaderExt *loader_ext;
	CharAEAE *ans_names_buf;

	loader_ext = loader->ext;
	ans_names_buf = &(loader_ext->ans_names_buf);
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAEAE(ans_names_buf, dataline->seq);
	return;
}

static void FASTAINFO_load_empty_seq(FASTAloader *loader)
{
	FASTAINFO_loaderExt *loader_ext;
	IntAE *seqlengths_buf;

	loader_ext = loader->ext;
	seqlengths_buf = &(loader_ext->seqlengths_buf);
	IntAE_insert_at(seqlengths_buf, IntAE_get_nelt(seqlengths_buf), 0);
	return;
}

static void FASTAINFO_load_seq_line(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	FASTAINFO_loaderExt *loader_ext;
	IntAE *seqlengths_buf;

	loader_ext = loader->ext;
	seqlengths_buf = &(loader_ext->seqlengths_buf);
	seqlengths_buf->elts[IntAE_get_nelt(seqlengths_buf) - 1] += dataline->length;
	return;
}

static FASTAINFO_loaderExt new_FASTAINFO_loaderExt()
{
	FASTAINFO_loaderExt loader_ext;

	loader_ext.ans_names_buf = new_CharAEAE(0, 0);
	loader_ext.seqlengths_buf = new_IntAE(0, 0, 0);
	return loader_ext;
}

static FASTAloader new_FASTAINFO_loader(SEXP lkup, int load_descs,
		FASTAINFO_loaderExt *loader_ext)
{
	FASTAloader loader;

	if (lkup == R_NilValue) {
		loader.lkup = NULL;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_length = LENGTH(lkup);
	}
	loader.load_desc_line = load_descs ? &FASTAINFO_load_desc_line : NULL;
	loader.load_empty_seq = &FASTAINFO_load_empty_seq;
	loader.load_seq_line = &FASTAINFO_load_seq_line;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

/*
 * The FASTA loader.
 */

typedef struct fasta_loader_ext {
	cachedXVectorList cached_ans;
	cachedCharSeq cached_ans_elt;
} FASTA_loaderExt;

static void FASTA_load_empty_seq(FASTAloader *loader)
{
	FASTA_loaderExt *loader_ext;
	cachedCharSeq *cached_ans_elt;

	loader_ext = loader->ext;
	cached_ans_elt = &(loader_ext->cached_ans_elt);
	*cached_ans_elt = get_cachedXRawList_elt(&(loader_ext->cached_ans),
						 loader->nrec);
	cached_ans_elt->length = 0;
	return;
}

static void FASTA_load_seq_line(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	FASTA_loaderExt *loader_ext;
	cachedCharSeq *cached_ans_elt;

	loader_ext = loader->ext;
	cached_ans_elt = &(loader_ext->cached_ans_elt);
	/* cached_ans_elt->seq is a const char * so we need to cast it to
	   char * before we can write to it */
	memcpy((char *) cached_ans_elt->seq + cached_ans_elt->length,
	       dataline->seq, dataline->length * sizeof(char));
	cached_ans_elt->length += dataline->length;
	return;
}

static FASTA_loaderExt new_FASTA_loaderExt(SEXP ans)
{
	FASTA_loaderExt loader_ext;

	loader_ext.cached_ans = cache_XVectorList(ans);
	return loader_ext;
}

static FASTAloader new_FASTA_loader(SEXP lkup, FASTA_loaderExt *loader_ext)
{
	FASTAloader loader;

	if (lkup == R_NilValue) {
		loader.lkup = NULL;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_length = LENGTH(lkup);
	}
	loader.load_desc_line = NULL;
	loader.load_empty_seq = &FASTA_load_empty_seq;
	loader.load_seq_line = &FASTA_load_seq_line;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

static int translate(cachedCharSeq *dataline, const int *lkup, int lkup_length)
{
	char *dest;
	int nbinvalid, i, j, key, val;

	/* dataline->seq is a const char * so we need to cast it to
	   char * before we can write to it */
	dest = (char *) dataline->seq;
	nbinvalid = j = 0;
	for (i = 0; i < dataline->length; i++) {
		key = (unsigned char) dataline->seq[i];
		if (key >= lkup_length || (val = lkup[key]) == NA_INTEGER) {
			nbinvalid++;
			continue;
		}
		dest[j++] = val;
	}
	dataline->length = j;
	return nbinvalid;
}

/*
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
 * the original Pearson FASTA format.
 */
static const char *parse_FASTA_file(FILE *stream, int *recno,
		int nrec, int skip, FASTAloader *loader)
{
	int FASTA_comment_markup_length, FASTA_desc_markup_length,
	    lineno, load_record, new_record;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTA_comment_markup_length = strlen(FASTA_comment_markup);
	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	lineno = 0;
	load_record = -1;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = delete_trailing_LF_or_CRLF(linebuf, -1);
		// dataline.length >= LINEBUF_SIZE should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		linebuf[dataline.length] = '\0';
		if (strncmp(linebuf, FASTA_comment_markup,
				FASTA_comment_markup_length) == 0)
			continue; // we ignore comment lines
		dataline.seq = linebuf;
		new_record = strncmp(linebuf, FASTA_desc_markup,
				FASTA_desc_markup_length) == 0;
		if (new_record) {
			load_record = *recno >= skip;
			if (load_record && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_record = load_record && loader != NULL;
			if (load_record && loader->load_desc_line != NULL) {
				dataline.seq += FASTA_desc_markup_length;
				dataline.length -= FASTA_desc_markup_length;
				loader->load_desc_line(loader, &dataline);
			}
			if (load_record && loader->load_empty_seq != NULL)
				loader->load_empty_seq(loader);
			if (load_record)
				loader->nrec++;
			(*recno)++;
			continue;
		}
		if (load_record == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "\"%s\" expected at beginning of line %d",
				 FASTA_desc_markup, lineno);
			return errmsg_buf;
		}
		if (load_record && loader->load_seq_line != NULL) {
			if (loader->lkup != NULL)
				translate(&dataline,
					  loader->lkup, loader->lkup_length);
			loader->load_seq_line(loader, &dataline);
		}
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fasta_info(SEXP efp_list, SEXP nrec, SEXP skip, SEXP use_names, SEXP lkup)
{
	int nrec0, skip0, load_descs, i, recno;
	FASTAINFO_loaderExt loader_ext;
	FASTAloader loader;
	FILE *stream;
	SEXP ans, ans_names;
	const char *errmsg;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	load_descs = LOGICAL(use_names)[0];
	loader_ext = new_FASTAINFO_loaderExt();
	loader = new_FASTAINFO_loader(lkup, load_descs, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		errmsg = parse_FASTA_file(stream, &recno,
				nrec0, skip0, &loader);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(efp_list), i)),
			      errmsg_buf);
	}
	PROTECT(ans = new_INTEGER_from_IntAE(&(loader_ext.seqlengths_buf)));
	if (load_descs) {
		PROTECT(ans_names =
		    new_CHARACTER_from_CharAEAE(&(loader_ext.ans_names_buf)));
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP read_fasta_in_XStringSet(SEXP efp_list, SEXP nrec, SEXP skip,
		SEXP use_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, i, recno;
	SEXP ans_width, ans_names, ans;
	const char *element_type;
	char classname[40];  /* longest string should be "DNAStringSet" */
	FASTA_loaderExt loader_ext;
	FASTAloader loader;
	FILE *stream;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	PROTECT(ans_width = fasta_info(efp_list, nrec, skip, use_names, lkup));
	PROTECT(ans_names = GET_NAMES(ans_width));
	SET_NAMES(ans_width, R_NilValue);
	element_type = CHAR(STRING_ELT(elementType, 0));
	if (snprintf(classname, sizeof(classname), "%sSet", element_type)
	    >= sizeof(classname))
	{
		UNPROTECT(2);
		error("Biostrings internal error in "
		      "read_fasta_in_XStringSet(): "
		      "'classname' buffer too small");
	}
	PROTECT(ans = alloc_XRawList(classname, element_type, ans_width));
	_set_XStringSet_names(ans, ans_names);
	loader_ext = new_FASTA_loaderExt(ans);
	loader = new_FASTA_loader(lkup, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		rewind(stream);
		parse_FASTA_file(stream, &recno, nrec0, skip0, &loader);
	}
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Writing FASTA files.
 */

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fasta(SEXP x, SEXP efp_list, SEXP width, SEXP lkup)
{
	cachedXStringSet X;
	int x_length, width0, lkup_length, i, j1, j2, dest_nbytes;
	FILE *stream;
	const int *lkup0;
	SEXP x_names, desc;
	cachedCharSeq X_elt;
	char linebuf[LINEBUF_SIZE];

	X = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&X);
	stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, 0));
	width0 = INTEGER(width)[0];
	if (width0 >= LINEBUF_SIZE)
		error("'width' must be <= %d", LINEBUF_SIZE - 1);
	linebuf[width0] = 0;
	if (lkup == R_NilValue) {
		lkup0 = NULL;
		lkup_length = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_length = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		if (fputs(FASTA_desc_markup, stream) == EOF)
			error("write error");
		if (x_names != R_NilValue) {
			desc = STRING_ELT(x_names, i);
			if (desc == NA_STRING)
				error("'names(x)' contains NAs");
			if (fputs(CHAR(desc), stream) == EOF)
				error("write error");
		}
		if (fputs("\n", stream) == EOF)
			error("write error");
		X_elt = _get_cachedXStringSet_elt(&X, i);
		for (j1 = 0; j1 < X_elt.length; j1 += width0) {
			j2 = j1 + width0;
			if (j2 > X_elt.length)
				j2 = X_elt.length;
			dest_nbytes = j2 - j1;
			j2--;
			Ocopy_bytes_from_i1i2_with_lkup(j1, j2,
				linebuf, dest_nbytes,
				X_elt.seq, X_elt.length,
				lkup0, lkup_length);
			linebuf[dest_nbytes] = 0;
			if (fputs(linebuf, stream) == EOF
			 || fputs("\n", stream) == EOF)
				error("write error");
		}
	}
	return R_NilValue;
}



/****************************************************************************
 *  B. FASTQ FORMAT                                                         *
 ****************************************************************************/

static const char *FASTQ_line1_markup = "@", *FASTQ_line3_markup = "+";


/****************************************************************************
 * Reading FASTQ files.
 */

typedef struct fastq_loader {
	void (*load_seqid)(struct fastq_loader *loader,
			   const cachedCharSeq *dataline);
	void (*load_seq)(struct fastq_loader *loader,
			 const cachedCharSeq *dataline);
	void (*load_qualid)(struct fastq_loader *loader,
			    const cachedCharSeq *dataline);
	void (*load_qual)(struct fastq_loader *loader,
			  const cachedCharSeq *dataline);
	int nrec;
	void *ext;  /* loader extension (optional) */
} FASTQloader;

/*
 * The FASTQGEOM loader keeps only the common width of the sequences
 * (NA if the set of sequences is not rectangular).
 */

typedef struct fastqgeom_loader_ext {
	int width;
} FASTQGEOM_loaderExt;

static void FASTQGEOM_load_seq(FASTQloader *loader,
		const cachedCharSeq *dataline)
{
	FASTQGEOM_loaderExt *loader_ext;

	loader_ext = loader->ext;
	if (loader->nrec == 0) {
		loader_ext->width = dataline->length;
		return;
	}
	if (loader_ext->width == NA_INTEGER
	 || dataline->length == loader_ext->width)
		return;
	loader_ext->width = NA_INTEGER;
	return;
}

static FASTQGEOM_loaderExt new_FASTQGEOM_loaderExt()
{
	FASTQGEOM_loaderExt loader_ext;

	loader_ext.width = NA_INTEGER;
	return loader_ext;
}

static FASTQloader new_FASTQGEOM_loader(FASTQGEOM_loaderExt *loader_ext)
{
	FASTQloader loader;

	loader.load_seqid = NULL;
	loader.load_seq = FASTQGEOM_load_seq;
	loader.load_qualid = NULL;
	loader.load_qual = NULL;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

/*
 * The FASTQ loader.
 */

typedef struct fastq_loader_ext {
	CharAEAE ans_names_buf;
	cachedXVectorList cached_ans;
	const int *lkup;
	int lkup_length;
} FASTQ_loaderExt;

static void FASTQ_load_seqid(FASTQloader *loader,
		const cachedCharSeq *dataline)
{
	FASTQ_loaderExt *loader_ext;
	CharAEAE *ans_names_buf;

	loader_ext = loader->ext;
	ans_names_buf = &(loader_ext->ans_names_buf);
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAEAE(ans_names_buf, dataline->seq);
	return;
}

static void FASTQ_load_seq(FASTQloader *loader, const cachedCharSeq *dataline)
{
	FASTQ_loaderExt *loader_ext;
	cachedCharSeq cached_ans_elt;

	loader_ext = loader->ext;
	cached_ans_elt = get_cachedXRawList_elt(&(loader_ext->cached_ans),
						loader->nrec);
	/* cached_ans_elt.seq is a const char * so we need to cast it to
	   char * before we can write to it */
	Ocopy_bytes_to_i1i2_with_lkup(0, cached_ans_elt.length - 1,
		(char *) cached_ans_elt.seq, cached_ans_elt.length,
		dataline->seq, dataline->length,
		loader_ext->lkup, loader_ext->lkup_length);
	return;
}

static FASTQ_loaderExt new_FASTQ_loaderExt(SEXP ans, SEXP lkup)
{
	FASTQ_loaderExt loader_ext;

	loader_ext.ans_names_buf =
		new_CharAEAE(_get_XStringSet_length(ans), 0);
	loader_ext.cached_ans = cache_XVectorList(ans);
	if (lkup == R_NilValue) {
		loader_ext.lkup = NULL;
	} else {
		loader_ext.lkup = INTEGER(lkup);
		loader_ext.lkup_length = LENGTH(lkup);
	}
	return loader_ext;
}

static FASTQloader new_FASTQ_loader(int load_seqids,
		FASTQ_loaderExt *loader_ext)
{
	FASTQloader loader;

	loader.load_seqid = load_seqids ? &FASTQ_load_seqid : NULL;
	loader.load_seq = FASTQ_load_seq;
	loader.load_qualid = NULL;
	loader.load_qual = NULL;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}


/*
 * Ignore empty lines.
 */
static const char *parse_FASTQ_file(FILE *stream, int *recno,
		int nrec, int skip, FASTQloader *loader)
{
	int FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineno, lineinrecno, load_record;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineno = lineinrecno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = delete_trailing_LF_or_CRLF(linebuf, -1);
		// dataline.length >= LINEBUF_SIZE should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		linebuf[dataline.length] = '\0';
		dataline.seq = linebuf;
		lineinrecno++;
		if (lineinrecno > 4)
			lineinrecno = 1;
		switch (lineinrecno) {
		    case 1:
			if (strncmp(linebuf, FASTQ_line1_markup,
					FASTQ_line1_markup_length) != 0) {
			    snprintf(errmsg_buf, sizeof(errmsg_buf),
				     "\"%s\" expected at beginning of line %d",
				     FASTQ_line1_markup, lineno);
				return errmsg_buf;
			}
			load_record = *recno >= skip;
			if (load_record && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_record = load_record && loader != NULL;
			if (load_record && nrec >= 0)
				load_record = *recno < skip + nrec;
			if (load_record && loader->load_seqid != NULL) {
				dataline.seq += FASTQ_line1_markup_length;
				dataline.length -= FASTQ_line1_markup_length;
				loader->load_seqid(loader, &dataline);
			}
		    break;
		    case 2:
			if (load_record && loader->load_seq != NULL)
				loader->load_seq(loader, &dataline);
		    break;
		    case 3:
			if (strncmp(linebuf, FASTQ_line3_markup,
					FASTQ_line3_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (load_record && loader->load_qualid != NULL) {
				dataline.seq += FASTQ_line3_markup_length;
				dataline.length -= FASTQ_line3_markup_length;
				loader->load_qualid(loader, &dataline);
			}
		    break;
		    case 4:
			if (load_record && loader->load_qual != NULL)
				loader->load_qual(loader, &dataline);
			if (load_record)
				loader->nrec++;
			(*recno)++;
		    break;
		}
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fastq_geometry(SEXP efp_list, SEXP nrec, SEXP skip)
{
	int nrec0, skip0, i, recno;
	FASTQGEOM_loaderExt loader_ext;
	FASTQloader loader;
	FILE *stream;
	const char *errmsg;
	SEXP ans;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	loader_ext = new_FASTQGEOM_loaderExt();
	loader = new_FASTQGEOM_loader(&loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		errmsg = parse_FASTQ_file(stream, &recno,
				nrec0, skip0, &loader);
		if (errmsg != NULL)
			error("reading FASTQ file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(efp_list), i)),
			      errmsg_buf);
	}
	PROTECT(ans = NEW_INTEGER(2));
	INTEGER(ans)[0] = loader.nrec;
	INTEGER(ans)[1] = loader_ext.width;
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP read_fastq_in_XStringSet(SEXP efp_list, SEXP nrec, SEXP skip,
		SEXP use_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, load_seqids, ans_length, i, recno;
	SEXP ans_geom, ans_width, ans, ans_names;
	const char *element_type;
	char classname[40];  /* longest string should be "DNAStringSet" */
	FASTQ_loaderExt loader_ext;
	FASTQloader loader;
	FILE *stream;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	load_seqids = LOGICAL(use_names)[0];
	PROTECT(ans_geom = fastq_geometry(efp_list, nrec, skip));
	ans_length = INTEGER(ans_geom)[0];
	PROTECT(ans_width = NEW_INTEGER(ans_length));
	if (ans_length != 0) {
		if (INTEGER(ans_geom)[1] == NA_INTEGER) {
			UNPROTECT(2);
			error("read_fastq_in_XStringSet(): FASTQ files with "
			      "variable sequence lengths are not supported yet");
		}
		for (recno = 0; recno < ans_length; recno++)
			INTEGER(ans_width)[recno] = INTEGER(ans_geom)[1];
	}
	element_type = CHAR(STRING_ELT(elementType, 0));
	if (snprintf(classname, sizeof(classname), "%sSet", element_type)
	    >= sizeof(classname))
	{
		UNPROTECT(2);
		error("Biostrings internal error in "
		      "read_fastq_in_XStringSet(): "
		      "'classname' buffer too small");
	}
	PROTECT(ans = alloc_XRawList(classname, element_type, ans_width));
	loader_ext = new_FASTQ_loaderExt(ans, lkup);
	loader = new_FASTQ_loader(load_seqids, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		rewind(stream);
		parse_FASTQ_file(stream, &recno, nrec0, skip0, &loader);
	}
	if (load_seqids) {
		PROTECT(ans_names =
		    new_CHARACTER_from_CharAEAE(&(loader_ext.ans_names_buf)));
		_set_XStringSet_names(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Writing FASTQ files.
 */

static const char *get_FASTQ_rec_id(SEXP x_names, SEXP q_names, int i)
{
	SEXP seqid, qualid;

	seqid = NA_STRING;
	if (x_names != R_NilValue) {
		seqid = STRING_ELT(x_names, i);
		if (seqid == NA_STRING)
			error("'names(x)' contains NAs");
	}
	if (q_names != R_NilValue) {
		qualid = STRING_ELT(q_names, i);
		if (qualid == NA_STRING)
			error("'names(qualities)' contains NAs");
		if (seqid == NA_STRING)
			seqid = qualid;
		/* Comparing the *content* of 2 CHARSXP's with
		   'qualid != seqid' only works because of the global CHARSXP
		   cache. Otherwise we would need to do:
		     LENGTH(qualid) != LENGTH(seqid) ||
		       memcmp(CHAR(qualid), CHAR(seqid), LENGTH(seqid))
		 */
		else if (qualid != seqid)
			error("when 'x' and 'qualities' both have "
			      "names, they must be identical");
	}
	if (seqid == NA_STRING)
		error("either 'x' or 'qualities' must have names");
	return CHAR(seqid);
}

static int write_FASTQ_id(FILE *stream, const char *markup, const char *id)
{
	return fputs(markup, stream) == EOF
	    || fputs(id, stream) == EOF
	    || fputs("\n", stream) == EOF;
}

static int write_FASTQ_seq(FILE *stream, const char *linebuf)
{
	return fputs(linebuf, stream) == EOF
	    || fputs("\n", stream) == EOF;
}

static int write_FASTQ_qual(FILE *stream, int seqlen,
		const cachedXStringSet *Q, int i)
{
	cachedCharSeq Q_elt;
	int j;

	Q_elt = _get_cachedXStringSet_elt(Q, i);
	if (Q_elt.length != seqlen)
		error("'x' and 'quality' must have the same width");
	for (j = 0; j < seqlen; j++) {
		if (fputc((int) *(Q_elt.seq++), stream) == EOF)
			return 1;
	}
	return fputs("\n", stream) == EOF;
}

static int write_FASTQ_fakequal(FILE *stream, int seqlen)
{
	int j;

	for (j = 0; j < seqlen; j++) {
		if (fputc((int) ';', stream) == EOF)
			return 1;
	}
	return fputs("\n", stream) == EOF;
}

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fastq(SEXP x, SEXP efp_list, SEXP qualities, SEXP lkup)
{
	cachedXStringSet X, Q;
	int x_length, lkup_length, i;
	FILE *stream;
	const int *lkup0;
	SEXP x_names, q_names;
	const char *id;
	cachedCharSeq X_elt;
	char linebuf[LINEBUF_SIZE];

	X = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&X);
	if (qualities != R_NilValue) {
		Q = _cache_XStringSet(qualities);
		if (_get_cachedXStringSet_length(&Q) != x_length)
			error("'x' and 'qualities' must have the same length");
		q_names = get_XVectorList_names(qualities);
	} else {
		q_names = R_NilValue;
	}
	stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, 0));
	if (lkup == R_NilValue) {
		lkup0 = NULL;
		lkup_length = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_length = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		id = get_FASTQ_rec_id(x_names, q_names, i);
		X_elt = _get_cachedXStringSet_elt(&X, i);
		Ocopy_bytes_from_i1i2_with_lkup(0, X_elt.length - 1,
			linebuf, X_elt.length,
			X_elt.seq, X_elt.length,
			lkup0, lkup_length);
		linebuf[X_elt.length] = 0;
		if (write_FASTQ_id(stream, FASTQ_line1_markup, id)
		 || write_FASTQ_seq(stream, linebuf)
		 || write_FASTQ_id(stream, FASTQ_line3_markup, id))
			error("write error");
		if (qualities != R_NilValue) {
			if (write_FASTQ_qual(stream, X_elt.length, &Q, i))
				error("write error");
		} else {
			if (write_FASTQ_fakequal(stream, X_elt.length))
				error("write error");
		}
	}
	return R_NilValue;
}

