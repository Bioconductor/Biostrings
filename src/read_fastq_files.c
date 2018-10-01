/****************************************************************************
 *                          Read/write FASTQ files                          *
 *                                 --------                                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"


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

/* TODO: Move this to XVector (together with _copy_CHARSXP_to_Chars_holder). */
static void copy_Chars_holder(Chars_holder *dest, const Chars_holder *src,
			      const int *lkup, int lkup_length)
{
	char *dest_ptr;

	/* dest->ptr is a (const char *) so we need to cast it to (char *)
	   before we can write to it */
	dest_ptr = (char *) dest->ptr;
	if (lkup == NULL) {
		memcpy(dest_ptr, src->ptr, dest->length);
	} else {
		Ocopy_bytes_to_i1i2_with_lkup(0, dest->length - 1,
			dest_ptr, dest->length,
			src->ptr, src->length,
			lkup, lkup_length);
	}
	return;
}


/****************************************************************************
 * FASTQ parser
 */

static const char *FASTQ_line1_markup = "@", *FASTQ_line3_markup = "+";

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
 * The FASTQ SEQLEN loader.
 * Used in parse_FASTQ_file() to load the read lengths only.
 */

typedef struct seqlen_fastq_loader_ext {
	IntAE *seqlength_buf;
} SEQLEN_FASTQloaderExt;

static SEQLEN_FASTQloaderExt new_SEQLEN_FASTQloaderExt()
{
	SEQLEN_FASTQloaderExt loader_ext;

	loader_ext.seqlength_buf = new_IntAE(0, 0, 0);
	return loader_ext;
}

static void FASTQ_SEQLEN_load_seq(FASTQloader *loader, const Chars_holder *seq)
{
	SEQLEN_FASTQloaderExt *loader_ext;
	IntAE *seqlength_buf;

	loader_ext = loader->ext;
	seqlength_buf = loader_ext->seqlength_buf;
	IntAE_insert_at(seqlength_buf, IntAE_get_nelt(seqlength_buf),
			seq->length);
	return;
}

static FASTQloader new_FASTQloader_with_SEQLEN_ext(
		SEQLEN_FASTQloaderExt *loader_ext)
{
	FASTQloader loader;

	loader.load_seqid = NULL;
	loader.load_seq = FASTQ_SEQLEN_load_seq;
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
	CharAEAE *seqid_buf;
	XVectorList_holder seq_holder;
	const int *lkup;
	int lkup_length;
	XVectorList_holder qual_holder;
} FASTQloaderExt;

static FASTQloaderExt new_FASTQloaderExt(SEXP sequences, SEXP lkup,
		SEXP qualities)
{
	FASTQloaderExt loader_ext;

	loader_ext.seqid_buf =
		new_CharAEAE(_get_XStringSet_length(sequences), 0);
	loader_ext.seq_holder = hold_XVectorList(sequences);
	if (lkup == R_NilValue) {
		loader_ext.lkup = NULL;
		loader_ext.lkup_length = 0;
	} else {
		loader_ext.lkup = INTEGER(lkup);
		loader_ext.lkup_length = LENGTH(lkup);
	}
	if (qualities != R_NilValue)
		loader_ext.qual_holder = hold_XVectorList(qualities);
	return loader_ext;
}

static void FASTQ_load_seqid(FASTQloader *loader, const Chars_holder *seqid)
{
	FASTQloaderExt *loader_ext;
	CharAEAE *seqid_buf;

	loader_ext = loader->ext;
	seqid_buf = loader_ext->seqid_buf;
	// This works only because seqid->ptr is nul-terminated!
	CharAEAE_append_string(seqid_buf, seqid->ptr);
	return;
}

static void FASTQ_load_seq(FASTQloader *loader, const Chars_holder *seq)
{
	FASTQloaderExt *loader_ext;
	Chars_holder seq_elt_holder;

	loader_ext = loader->ext;
	seq_elt_holder = get_elt_from_XRawList_holder(
					&(loader_ext->seq_holder),
					loader->nrec);
	copy_Chars_holder(&seq_elt_holder, seq,
			  loader_ext->lkup, loader_ext->lkup_length);
	return;
}

static void FASTQ_load_qual(FASTQloader *loader, const Chars_holder *qual)
{
	FASTQloaderExt *loader_ext;
	Chars_holder qual_elt_holder;

	loader_ext = loader->ext;
	qual_elt_holder = get_elt_from_XRawList_holder(
					&(loader_ext->qual_holder),
					loader->nrec);
	copy_Chars_holder(&qual_elt_holder, qual, NULL, 0);
	return;
}

static FASTQloader new_FASTQloader(int load_seqids, int load_quals,
				   FASTQloaderExt *loader_ext)
{
	FASTQloader loader;

	loader.load_seqid = load_seqids ? &FASTQ_load_seqid : NULL;
	loader.load_seq = FASTQ_load_seq;
	if (load_quals) {
		/* Quality ids are always ignored for now. */
		//loader.load_qualid = &FASTQ_load_qualid;
		loader.load_qualid = NULL;
		loader.load_qual = &FASTQ_load_qual;
	} else {
		loader.load_qualid = NULL;
		loader.load_qual = NULL;
	}
	loader.nrec = 0;
	loader.ext = loader_ext;
	return loader;
}

/*
 * Ignore empty lines.
 */
static const char *parse_FASTQ_file(SEXP filexp,
		int nrec, int skip, int seek_first_rec,
		FASTQloader *loader,
		int *recno, long long int *offset)
{
	int lineno, EOL_in_buf, EOL_in_prev_buf, ret_code, nbyte_in,
	    FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineinrecno, load_rec, seq_len;
	char buf[IOBUF_SIZE];
	Chars_holder data;
	long long int prev_offset;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineinrecno = 0;
	lineno = 0;
	EOL_in_buf = 1;
	while (1) {
		if (EOL_in_buf)
			lineno++;
		EOL_in_prev_buf = EOL_in_buf;
		ret_code = filexp_gets(filexp, buf, IOBUF_SIZE, &EOL_in_buf);
		if (ret_code == 0)
			break;
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		nbyte_in = strlen(buf);
		data.length = delete_trailing_LF_or_CRLF(buf, nbyte_in);
		prev_offset = *offset;
		*offset += nbyte_in;
		if (seek_first_rec) {
			if (EOL_in_prev_buf
			 && has_prefix(buf, FASTQ_line1_markup)) {
				seek_first_rec = 0;
			} else {
				continue;
			}
		}
		if (data.length == 0)
			continue; // we ignore empty lines
		buf[data.length] = '\0';
		data.ptr = buf;
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
			if (load_rec && nrec >= 0 && *recno >= skip + nrec) {
				/* Calls to filexp_seek() are costly on
				   compressed files and the cost increases as
				   we advance in the file. This is not a
				   problem when reading the entire file but
				   becomes one when reading a compressed file
				   by chunk. */
				filexp_seek(filexp, prev_offset, SEEK_SET);
				return NULL;
			}
			load_rec = load_rec && loader != NULL;
			if (load_rec && nrec >= 0)
				load_rec = *recno < skip + nrec;
			if (load_rec && loader->load_seqid != NULL) {
				data.ptr += FASTQ_line1_markup_length;
				data.length -= FASTQ_line1_markup_length;
				loader->load_seqid(loader, &data);
			}
		    break;
		    case 2:
			if (load_rec) {
				seq_len = data.length;
				if  (loader->load_seq != NULL)
					loader->load_seq(loader, &data);
			}
		    break;
		    case 3:
			if (!has_prefix(buf, FASTQ_line3_markup)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of "
					 "line %d", FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (load_rec && loader->load_qualid != NULL) {
				data.ptr += FASTQ_line3_markup_length;
				data.length -= FASTQ_line3_markup_length;
				loader->load_qualid(loader, &data);
			}
		    break;
		    case 4:
			if (load_rec) {
				if (data.length != seq_len) {
					snprintf(errmsg_buf, sizeof(errmsg_buf),
						 "length of quality string "
						 "at line %d\n  differs from "
						 "length of corresponding "
						 "sequence", lineno);
					return errmsg_buf;
				}
				if (loader->load_qual != NULL)
					loader->load_qual(loader, &data);
			}
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


/****************************************************************************
 * read_fastq_files()
 */

static SEXP get_fastq_seqlengths(SEXP filexp_list,
		int nrec, int skip, int seek_first_rec)
{
	SEQLEN_FASTQloaderExt loader_ext;
	FASTQloader loader;
	int recno, i;
	SEXP filexp;
	long long int offset0, offset;
	const char *errmsg;

	loader_ext = new_SEQLEN_FASTQloaderExt();
	loader = new_FASTQloader_with_SEQLEN_ext(&loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		/* Calls to filexp_tell() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		offset0 = offset = filexp_tell(filexp);
		errmsg = parse_FASTQ_file(filexp, nrec, skip, seek_first_rec,
					  &loader,
					  &recno, &offset);
		/* Calls to filexp_seek() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		filexp_seek(filexp, offset0, SEEK_SET);
		if (errmsg != NULL)
			error("reading FASTQ file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(filexp_list), i)),
			      errmsg_buf);
	}
	return new_INTEGER_from_IntAE(loader_ext.seqlength_buf);
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   filexp_list:     A list of external pointers.
 *   nrec:            An integer vector of length 1.
 *   skip:            An integer vector of length 1.
 *   seek_first_rec:  A logical vector of length 1.
 */
SEXP fastq_seqlengths(SEXP filexp_list,
		SEXP nrec, SEXP skip, SEXP seek_first_rec)
{
	int nrec0, skip0, seek_rec0;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	return get_fastq_seqlengths(filexp_list, nrec0, skip0, seek_rec0);
}

/* --- .Call ENTRY POINT ---
 * Return an XStringSet object if 'with_qualities' is FALSE, or a list of 2
 * parallel XStringSet objects of the same shape if 'with_qualities' is TRUE.
 * We use a 2-pass algo so we can pre-alloc the exact amount of memory for the
 * XStringSet object that we load with the string data during the 2nd pass.
 * Because a XStringSet object cannot be grown in-place, a 1-pass algo would
 * need to load the string data in a growing object first (e.g. CharAEAE)
 * then turn it into an XStringSet object (this would require a copy of the
 * entire string data).
 */
SEXP read_fastq_files(SEXP filexp_list, SEXP nrec, SEXP skip,
		SEXP seek_first_rec,
		SEXP use_names, SEXP elementType, SEXP lkup,
		SEXP with_qualities)
{
	int nrec0, skip0, seek_rec0, load_seqids, load_quals, recno, i;
	SEXP filexp, seqlengths, sequences, seqids, qualities, ans;
	FASTQloaderExt loader_ext;
	FASTQloader loader;
	long long int offset;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	load_seqids = LOGICAL(use_names)[0];
	load_quals = LOGICAL(with_qualities)[0];
	/* 1st pass */
	PROTECT(seqlengths = get_fastq_seqlengths(filexp_list,
						  nrec0, skip0, seek_rec0));
	/* Allocation */
	PROTECT(sequences = _alloc_XStringSet(CHAR(STRING_ELT(elementType, 0)),
					      seqlengths));
	if (load_quals) {
		PROTECT(qualities = _alloc_XStringSet("BString", seqlengths));
	} else {
		qualities = R_NilValue;
	}
	/* 2nd pass */
	loader_ext = new_FASTQloaderExt(sequences, lkup, qualities);
	loader = new_FASTQloader(load_seqids, load_quals, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		/* Calls to filexp_tell() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		offset = filexp_tell(filexp);
		parse_FASTQ_file(filexp, nrec0, skip0, seek_rec0,
				 &loader,
				 &recno, &offset);
	}
	if (load_seqids) {
		PROTECT(seqids =
			new_CHARACTER_from_CharAEAE(loader_ext.seqid_buf));
		_set_XStringSet_names(sequences, seqids);
		UNPROTECT(1);
	}
	if (!load_quals) {
		UNPROTECT(2);
		return sequences;
	}
	PROTECT(ans = NEW_LIST(2));
	SET_ELEMENT(ans, 0, sequences);
	SET_ELEMENT(ans, 1, qualities);
	UNPROTECT(4);
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

static void write_FASTQ_id(SEXP filexp, const char *markup, const char *id)
{
	filexp_puts(filexp, markup);
	filexp_puts(filexp, id);
	filexp_puts(filexp, "\n");
}

static void write_FASTQ_seq(SEXP filexp, const char *buf)
{
	filexp_puts(filexp, buf);
	filexp_puts(filexp, "\n");
}

static void write_FASTQ_qual(SEXP filexp, int seqlen,
		const XStringSet_holder *Q, int i)
{
	Chars_holder Q_elt;
	int j;

	Q_elt = _get_elt_from_XStringSet_holder(Q, i);
	if (Q_elt.length != seqlen)
		error("'x' and 'quality' must have the same width");
	for (j = 0; j < seqlen; j++)
		filexp_putc(filexp, (int) *(Q_elt.ptr++));
	filexp_puts(filexp, "\n");
}

static void write_FASTQ_fakequal(SEXP filexp, int seqlen)
{
	int j;

	for (j = 0; j < seqlen; j++)
		filexp_putc(filexp, (int) ';');
	filexp_puts(filexp, "\n");
}

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fastq(SEXP x, SEXP filexp_list,
		SEXP qualities, SEXP lkup)
{
	XStringSet_holder X, Q;
	int x_length, lkup_length, i;
	const int *lkup0;
	SEXP filexp, x_names, q_names;
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
	filexp = VECTOR_ELT(filexp_list, 0);
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
			X_elt.ptr, X_elt.length,
			lkup0, lkup_length);
		buf[X_elt.length] = 0;
		write_FASTQ_id(filexp, FASTQ_line1_markup, id);
		write_FASTQ_seq(filexp, buf);
		write_FASTQ_id(filexp, FASTQ_line3_markup, id);
		if (qualities != R_NilValue) {
			write_FASTQ_qual(filexp, X_elt.length, &Q, i);
		} else {
			write_FASTQ_fakequal(filexp, X_elt.length);
		}
	}
	return R_NilValue;
}

