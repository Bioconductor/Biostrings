/****************************************************************************
 *                       Read/write FASTA/FASTQ files                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

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

#define IOBUF_SIZE 20002
static char errmsg_buf[200];

static int has_prefix(const char *s, const char *prefix)
{
	int i = 0;
	char c;

	while ((c = prefix[i]) != '\0') {
		if (s[i] != c)
			return 0;
		i++;
	}
	return 1;
}



/****************************************************************************
 *  A. FASTA FORMAT                                                         *
 ****************************************************************************/

/*
 * Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be ready to support longer sympols just in case.
 */
static const char *FASTA_comment_markup = ";", *FASTA_desc_markup = ">";


/****************************************************************************
 * Reading FASTA files.
 */

typedef struct fasta_loader {
	const int *lkup;
	int lkup_length;
	void (*load_desc_line)(struct fasta_loader *loader,
			       const Chars_holder *desc_line);
	void (*load_empty_seq)(struct fasta_loader *loader);
	void (*load_seq_data)(struct fasta_loader *loader,
			      const Chars_holder *seq_data);
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
		const Chars_holder *desc_line)
{
	FASTAINFO_loaderExt *loader_ext;
	CharAEAE *ans_names_buf;

	loader_ext = loader->ext;
	ans_names_buf = &(loader_ext->ans_names_buf);
	// This works only because desc_line->seq is nul-terminated!
	append_string_to_CharAEAE(ans_names_buf, desc_line->seq);
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

static void FASTAINFO_load_seq_data(FASTAloader *loader,
		const Chars_holder *seq_data)
{
	FASTAINFO_loaderExt *loader_ext;
	IntAE *seqlengths_buf;

	loader_ext = loader->ext;
	seqlengths_buf = &(loader_ext->seqlengths_buf);
	seqlengths_buf->elts[IntAE_get_nelt(seqlengths_buf) - 1] +=
		seq_data->length;
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
		loader.lkup_length = 0;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_length = LENGTH(lkup);
	}
	loader.load_desc_line = load_descs ? &FASTAINFO_load_desc_line : NULL;
	loader.load_empty_seq = &FASTAINFO_load_empty_seq;
	loader.load_seq_data = &FASTAINFO_load_seq_data;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

/*
 * The FASTA loader.
 */

typedef struct fasta_loader_ext {
	XVectorList_holder ans_holder;
	Chars_holder ans_elt_holder;
} FASTA_loaderExt;

static void FASTA_load_empty_seq(FASTAloader *loader)
{
	FASTA_loaderExt *loader_ext;
	Chars_holder *ans_elt_holder;

	loader_ext = loader->ext;
	ans_elt_holder = &(loader_ext->ans_elt_holder);
	*ans_elt_holder = get_elt_from_XRawList_holder(&(loader_ext->ans_holder),
						 loader->nrec);
	ans_elt_holder->length = 0;
	return;
}

static void FASTA_load_seq_data(FASTAloader *loader,
		const Chars_holder *seq_data)
{
	FASTA_loaderExt *loader_ext;
	Chars_holder *ans_elt_holder;

	loader_ext = loader->ext;
	ans_elt_holder = &(loader_ext->ans_elt_holder);
	/* ans_elt_holder->seq is a const char * so we need to cast it to
	   char * in order to write to it */
	memcpy((char *) ans_elt_holder->seq + ans_elt_holder->length,
	       seq_data->seq, seq_data->length * sizeof(char));
	ans_elt_holder->length += seq_data->length;
	return;
}

static FASTA_loaderExt new_FASTA_loaderExt(SEXP ans)
{
	FASTA_loaderExt loader_ext;

	loader_ext.ans_holder = hold_XVectorList(ans);
	return loader_ext;
}

static FASTAloader new_FASTA_loader(SEXP lkup, FASTA_loaderExt *loader_ext)
{
	FASTAloader loader;

	if (lkup == R_NilValue) {
		loader.lkup = NULL;
		loader.lkup_length = 0;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_length = LENGTH(lkup);
	}
	loader.load_desc_line = NULL;
	loader.load_empty_seq = &FASTA_load_empty_seq;
	loader.load_seq_data = &FASTA_load_seq_data;
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

static int translate(Chars_holder *seq_data, const int *lkup, int lkup_length)
{
	char *dest;
	int nbinvalid, i, j, key, val;

	/* seq_data->seq is a const char * so we need to cast it to
	   char * before we can write to it */
	dest = (char *) seq_data->seq;
	nbinvalid = j = 0;
	for (i = 0; i < seq_data->length; i++) {
		key = (unsigned char) seq_data->seq[i];
		if (key >= lkup_length || (val = lkup[key]) == NA_INTEGER) {
			nbinvalid++;
			continue;
		}
		dest[j++] = val;
	}
	seq_data->length = j;
	return nbinvalid;
}

/*
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
 * the original Pearson FASTA format.
 */
static const char *parse_FASTA_file(SEXP efp, int *recno, int *ninvalid,
		int nrec, int skip, int seek_first_rec, FASTAloader *loader)
{
	int lineno, EOL_in_buf, EOL_in_prev_buf, ret_code,
	    FASTA_desc_markup_length, load_rec, is_new_rec;
	char buf[IOBUF_SIZE];
	Chars_holder data;

	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	load_rec = -1;
	for (lineno = EOL_in_prev_buf = 1;
	     (ret_code = ExternalFilePtr_gets(efp, buf, IOBUF_SIZE,
					      &EOL_in_buf));
	     lineno += EOL_in_prev_buf = EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (seek_first_rec) {
			if (EOL_in_prev_buf
			 && has_prefix(buf, FASTA_desc_markup)) {
				seek_first_rec = 0;
			} else {
				continue;
			}
		}
		if (EOL_in_buf) {
			data.length = delete_trailing_LF_or_CRLF(buf, -1);
		} else {
			data.length = IOBUF_SIZE - 1;
		}
		data.seq = buf;
		if (!EOL_in_prev_buf)
			goto parse_seq_data;
		if (data.length == 0)
			continue; // we ignore empty lines
		if (has_prefix(data.seq, FASTA_comment_markup)) {
			if (!EOL_in_buf) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "cannot read line %d, "
					 "line is too long", lineno);
				return errmsg_buf;
			}
			continue; // we ignore comment lines
		}
		buf[data.length] = '\0';
		is_new_rec = has_prefix(data.seq, FASTA_desc_markup);
		if (is_new_rec) {
			if (!EOL_in_buf) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "cannot read line %d, "
					 "line is too long", lineno);
				return errmsg_buf;
			}
			load_rec = *recno >= skip;
			if (load_rec && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_rec = load_rec && loader != NULL;
			if (load_rec && loader->load_desc_line != NULL) {
				data.seq += FASTA_desc_markup_length;
				data.length -= FASTA_desc_markup_length;
				loader->load_desc_line(loader, &data);
			}
			if (load_rec && loader->load_empty_seq != NULL)
				loader->load_empty_seq(loader);
			if (load_rec)
				loader->nrec++;
			(*recno)++;
			continue;
		}
		if (load_rec == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "\"%s\" expected at beginning of line %d",
				 FASTA_desc_markup, lineno);
			return errmsg_buf;
		}
		parse_seq_data:
		if (load_rec && loader->load_seq_data != NULL) {
			if (loader->lkup != NULL)
				*ninvalid += translate(&data,
						       loader->lkup,
						       loader->lkup_length);
			loader->load_seq_data(loader, &data);
		}
	}
	if (seek_first_rec) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "no FASTA record found");
		return errmsg_buf;
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fasta_info(SEXP efp_list, SEXP nrec, SEXP skip, SEXP seek_first_rec,
		SEXP use_names, SEXP lkup)
{
	int nrec0, skip0, seek_rec0, load_descs, i, recno, ninvalid;
	FASTAINFO_loaderExt loader_ext;
	FASTAloader loader;
	SEXP efp, ans, ans_names;
	const char *errmsg;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	load_descs = LOGICAL(use_names)[0];
	loader_ext = new_FASTAINFO_loaderExt();
	loader = new_FASTAINFO_loader(lkup, load_descs, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		efp = VECTOR_ELT(efp_list, i);
		ninvalid = 0;
		errmsg = parse_FASTA_file(efp, &recno, &ninvalid,
					  nrec0, skip0, seek_rec0, &loader);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(efp_list), i)),
			      errmsg_buf);
		if (ninvalid != 0)
			warning("reading FASTA file %s: ignored %d "
				"invalid one-letter sequence codes",
				CHAR(STRING_ELT(GET_NAMES(efp_list), i)),
				ninvalid);
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
		SEXP seek_first_rec,
		SEXP use_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, seek_rec0, i, recno, ninvalid;
	SEXP efp, ans_width, ans_names, ans;
	const char *element_type;
	char classname[40];  /* longest string should be "DNAStringSet" */
	FASTA_loaderExt loader_ext;
	FASTAloader loader;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	PROTECT(ans_width = fasta_info(efp_list, nrec, skip, seek_first_rec,
				       use_names, lkup));
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
	recno = ninvalid = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		efp = VECTOR_ELT(efp_list, i);
		ExternalFilePtr_rewind(efp);
		parse_FASTA_file(efp, &recno, &ninvalid,
				 nrec0, skip0, seek_rec0, &loader);
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
	XStringSet_holder X;
	int x_length, width0, lkup_length, i, j1, j2, dest_nbytes;
	const int *lkup0;
	SEXP efp, x_names, desc;
	Chars_holder X_elt;
	char buf[IOBUF_SIZE];

	X = _hold_XStringSet(x);
	x_length = _get_length_from_XStringSet_holder(&X);
	efp = VECTOR_ELT(efp_list, 0);
	width0 = INTEGER(width)[0];
	if (width0 >= IOBUF_SIZE)
		error("'width' must be <= %d", IOBUF_SIZE - 1);
	buf[width0] = 0;
	if (lkup == R_NilValue) {
		lkup0 = NULL;
		lkup_length = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_length = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		ExternalFilePtr_puts(efp, FASTA_desc_markup);
		if (x_names != R_NilValue) {
			desc = STRING_ELT(x_names, i);
			if (desc == NA_STRING)
				error("'names(x)' contains NAs");
			ExternalFilePtr_puts(efp, CHAR(desc));
		}
		ExternalFilePtr_puts(efp, "\n");
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		for (j1 = 0; j1 < X_elt.length; j1 += width0) {
			j2 = j1 + width0;
			if (j2 > X_elt.length)
				j2 = X_elt.length;
			dest_nbytes = j2 - j1;
			j2--;
			Ocopy_bytes_from_i1i2_with_lkup(j1, j2,
				buf, dest_nbytes,
				X_elt.seq, X_elt.length,
				lkup0, lkup_length);
			buf[dest_nbytes] = 0;
			ExternalFilePtr_puts(efp, buf);
			ExternalFilePtr_puts(efp, "\n");
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
			   const Chars_holder *seqid);
	void (*load_seq)(struct fastq_loader *loader,
			 const Chars_holder *seq);
	void (*load_qualid)(struct fastq_loader *loader,
			    const Chars_holder *qualid);
	void (*load_qual)(struct fastq_loader *loader,
			  const Chars_holder *qual);
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

static void FASTQGEOM_load_seq(FASTQloader *loader, const Chars_holder *seq)
{
	FASTQGEOM_loaderExt *loader_ext;

	loader_ext = loader->ext;
	if (loader->nrec == 0) {
		loader_ext->width = seq->length;
		return;
	}
	if (loader_ext->width == NA_INTEGER
	 || seq->length == loader_ext->width)
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
	XVectorList_holder ans_holder;
	const int *lkup;
	int lkup_length;
} FASTQ_loaderExt;

static void FASTQ_load_seqid(FASTQloader *loader, const Chars_holder *seqid)
{
	FASTQ_loaderExt *loader_ext;
	CharAEAE *ans_names_buf;

	loader_ext = loader->ext;
	ans_names_buf = &(loader_ext->ans_names_buf);
	// This works only because seqid->seq is nul-terminated!
	append_string_to_CharAEAE(ans_names_buf, seqid->seq);
	return;
}

static void FASTQ_load_seq(FASTQloader *loader, const Chars_holder *seq)
{
	FASTQ_loaderExt *loader_ext;
	Chars_holder ans_elt_holder;

	loader_ext = loader->ext;
	ans_elt_holder = get_elt_from_XRawList_holder(&(loader_ext->ans_holder),
						loader->nrec);
	/* ans_elt_holder.seq is a const char * so we need to cast it to
	   char * before we can write to it */
	Ocopy_bytes_to_i1i2_with_lkup(0, ans_elt_holder.length - 1,
		(char *) ans_elt_holder.seq, ans_elt_holder.length,
		seq->seq, seq->length,
		loader_ext->lkup, loader_ext->lkup_length);
	return;
}

static FASTQ_loaderExt new_FASTQ_loaderExt(SEXP ans, SEXP lkup)
{
	FASTQ_loaderExt loader_ext;

	loader_ext.ans_names_buf =
		new_CharAEAE(_get_XStringSet_length(ans), 0);
	loader_ext.ans_holder = hold_XVectorList(ans);
	if (lkup == R_NilValue) {
		loader_ext.lkup = NULL;
		loader_ext.lkup_length = 0;
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
static const char *parse_FASTQ_file(SEXP efp, int *recno,
		int nrec, int skip, int seek_first_rec, FASTQloader *loader)
{
	int lineno, EOL_in_buf, EOL_in_prev_buf, ret_code,
	    FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineinrecno, load_rec;
	char buf[IOBUF_SIZE];
	Chars_holder data;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineinrecno = 0;
	for (lineno = EOL_in_prev_buf = 1;
	     (ret_code = ExternalFilePtr_gets(efp, buf, IOBUF_SIZE,
					      &EOL_in_buf));
	     lineno += EOL_in_prev_buf = EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (seek_first_rec) {
			if (EOL_in_prev_buf
			 && has_prefix(buf, FASTQ_line1_markup)) {
				seek_first_rec = 0;
			} else {
				continue;
			}
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		data.length = delete_trailing_LF_or_CRLF(buf, -1);
		if (data.length == 0)
			continue; // we ignore empty lines
		buf[data.length] = '\0';
		data.seq = buf;
		lineinrecno++;
		if (lineinrecno > 4)
			lineinrecno = 1;
		switch (lineinrecno) {
		    case 1:
			if (!has_prefix(buf, FASTQ_line1_markup)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
				     "\"%s\" expected at beginning of line %d",
				     FASTQ_line1_markup, lineno);
				return errmsg_buf;
			}
			load_rec = *recno >= skip;
			if (load_rec && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_rec = load_rec && loader != NULL;
			if (load_rec && nrec >= 0)
				load_rec = *recno < skip + nrec;
			if (load_rec && loader->load_seqid != NULL) {
				data.seq += FASTQ_line1_markup_length;
				data.length -= FASTQ_line1_markup_length;
				loader->load_seqid(loader, &data);
			}
		    break;
		    case 2:
			if (load_rec && loader->load_seq != NULL)
				loader->load_seq(loader, &data);
		    break;
		    case 3:
			if (!has_prefix(buf, FASTQ_line3_markup)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of "
					 "line %d", FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (load_rec && loader->load_qualid != NULL) {
				data.seq += FASTQ_line3_markup_length;
				data.length -= FASTQ_line3_markup_length;
				loader->load_qualid(loader, &data);
			}
		    break;
		    case 4:
			if (load_rec && loader->load_qual != NULL)
				loader->load_qual(loader, &data);
			if (load_rec)
				loader->nrec++;
			(*recno)++;
		    break;
		}
	}
	if (seek_first_rec) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "no FASTQ record found");
		return errmsg_buf;
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fastq_geometry(SEXP efp_list, SEXP nrec, SEXP skip, SEXP seek_first_rec)
{
	int nrec0, skip0, seek_rec0, i, recno;
	FASTQGEOM_loaderExt loader_ext;
	FASTQloader loader;
	const char *errmsg;
	SEXP efp, ans;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	loader_ext = new_FASTQGEOM_loaderExt();
	loader = new_FASTQGEOM_loader(&loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		efp = VECTOR_ELT(efp_list, i);
		errmsg = parse_FASTQ_file(efp, &recno,
					  nrec0, skip0, seek_rec0, &loader);
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
		SEXP seek_first_rec,
		SEXP use_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, seek_rec0, load_seqids, ans_length, i, recno;
	SEXP efp, ans_geom, ans_width, ans, ans_names;
	const char *element_type;
	char classname[40];  /* longest string should be "DNAStringSet" */
	FASTQ_loaderExt loader_ext;
	FASTQloader loader;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	load_seqids = LOGICAL(use_names)[0];
	PROTECT(ans_geom = fastq_geometry(efp_list, nrec, skip,
					  seek_first_rec));
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
		efp = VECTOR_ELT(efp_list, i);
		ExternalFilePtr_rewind(efp);
		parse_FASTQ_file(efp, &recno,
				 nrec0, skip0, seek_rec0, &loader);
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

static void write_FASTQ_id(SEXP efp, const char *markup, const char *id)
{
	ExternalFilePtr_puts(efp, markup);
	ExternalFilePtr_puts(efp, id);
	ExternalFilePtr_puts(efp, "\n");
}

static void write_FASTQ_seq(SEXP efp, const char *buf)
{
	ExternalFilePtr_puts(efp, buf);
	ExternalFilePtr_puts(efp, "\n");
}

static void write_FASTQ_qual(SEXP efp, int seqlen,
		const XStringSet_holder *Q, int i)
{
	Chars_holder Q_elt;
	int j;

	Q_elt = _get_elt_from_XStringSet_holder(Q, i);
	if (Q_elt.length != seqlen)
		error("'x' and 'quality' must have the same width");
	for (j = 0; j < seqlen; j++)
		ExternalFilePtr_putc(efp, (int) *(Q_elt.seq++));
	ExternalFilePtr_puts(efp, "\n");
}

static void write_FASTQ_fakequal(SEXP efp, int seqlen)
{
	int j;

	for (j = 0; j < seqlen; j++)
		ExternalFilePtr_putc(efp, (int) ';');
	ExternalFilePtr_puts(efp, "\n");
}

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fastq(SEXP x, SEXP efp_list, SEXP qualities, SEXP lkup)
{
	XStringSet_holder X, Q;
	int x_length, lkup_length, i;
	const int *lkup0;
	SEXP efp, x_names, q_names;
	const char *id;
	Chars_holder X_elt;
	char buf[IOBUF_SIZE];

	X = _hold_XStringSet(x);
	x_length = _get_length_from_XStringSet_holder(&X);
	if (qualities != R_NilValue) {
		Q = _hold_XStringSet(qualities);
		if (_get_length_from_XStringSet_holder(&Q) != x_length)
			error("'x' and 'qualities' must have the same length");
		q_names = get_XVectorList_names(qualities);
	} else {
		q_names = R_NilValue;
	}
	efp = VECTOR_ELT(efp_list, 0);
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
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		Ocopy_bytes_from_i1i2_with_lkup(0, X_elt.length - 1,
			buf, X_elt.length,
			X_elt.seq, X_elt.length,
			lkup0, lkup_length);
		buf[X_elt.length] = 0;
		write_FASTQ_id(efp, FASTQ_line1_markup, id);
		write_FASTQ_seq(efp, buf);
		write_FASTQ_id(efp, FASTQ_line3_markup, id);
		if (qualities != R_NilValue) {
			write_FASTQ_qual(efp, X_elt.length, &Q, i);
		} else {
			write_FASTQ_fakequal(efp, X_elt.length);
		}
	}
	return R_NilValue;
}

