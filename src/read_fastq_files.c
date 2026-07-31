/****************************************************************************
 *                          Read/write FASTQ files                          *
 *                                 --------                                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"


#define IOBUF_SIZE 200001
static char iobuf[IOBUF_SIZE];

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

static int translate(Chars_holder *x, const int *lkup, int lkup_len)
{
	char *dest;
	int nbinvalid, i, j, c;

	/* x->ptr is a (const char *) so we need to cast it to (char *)
	   before we can write to it */
	dest = (char *) x->ptr;
	nbinvalid = j = 0;
	for (i = 0; i < x->length; i++) {
		c = translate_byte(x->ptr[i], lkup, lkup_len);
		if (c == NA_INTEGER) {
			nbinvalid++;
			continue;
		}
		dest[j++] = (char) c;
	}
	x->length = j;
	return nbinvalid;
}

static void append_Chars_holder(Chars_holder *dest, const Chars_holder *src)
{
	/* dest->ptr is a (const char *) so we need to cast it to (char *)
	   in order to write to it */
	memcpy((char *) dest->ptr + dest->length,
	       src->ptr, src->length * sizeof(char));
	dest->length += src->length;
	return;
}


/****************************************************************************
 * FASTQ parser
 */

static const char *FASTQ_line1_markup = "@", *FASTQ_line3_markup = "+";

typedef struct fastq_loader {
	void (*new_seqid_hook)(struct fastq_loader *loader,
			       const Chars_holder *seqid);
	void (*new_empty_seq_hook)(struct fastq_loader *loader);
	const char *(*append_seq_hook)(struct fastq_loader *loader,
				       Chars_holder *seq_data);
	void (*new_qualid_hook)(struct fastq_loader *loader,
				const Chars_holder *qualid);
        void (*new_empty_qual_hook)(struct fastq_loader *loader);
	const char *(*append_qual_hook)(struct fastq_loader *loader,
					const Chars_holder *qual_data);
	const int *lkup;
	int lkup_len;
	void *ext;  /* loader extension (optional) */
} FASTQloader;

/*
 * The FASTQ SEQLEN loader.
 * Used in parse_FASTQ_file() to load the read lengths only.
 * Does NOT check that the read sequences are valid and completely ignores
 * the quality sequences.
 */

typedef struct seqlen_fastq_loader_ext {
	IntAE *seqlength_buf;
} SEQLEN_FASTQloaderExt;

static SEQLEN_FASTQloaderExt new_SEQLEN_FASTQloaderExt(void)
{
	SEQLEN_FASTQloaderExt loader_ext;

	loader_ext.seqlength_buf = new_IntAE(0, 0, 0);
	return loader_ext;
}

static void FASTQ_SEQLEN_new_empty_seq_hook(FASTQloader *loader)
{
	SEQLEN_FASTQloaderExt *loader_ext;
	IntAE *seqlength_buf;

	loader_ext = loader->ext;
	seqlength_buf = loader_ext->seqlength_buf;
	IntAE_insert_at(seqlength_buf, IntAE_get_nelt(seqlength_buf), 0);
	return;
}

/* Unlike FASTQ_append_seq_hook(), FASTQ_SEQLEN_append_seq_hook() does NOT
   check that the read sequence is valid. */
static const char *FASTQ_SEQLEN_append_seq_hook(FASTQloader *loader,
		Chars_holder *seq_data)
{
	SEQLEN_FASTQloaderExt *loader_ext;
	IntAE *seqlength_buf;

	loader_ext = loader->ext;
	seqlength_buf = loader_ext->seqlength_buf;
	seqlength_buf->elts[IntAE_get_nelt(seqlength_buf) - 1] +=
		seq_data->length;
	return NULL;
}

static FASTQloader new_FASTQloader_with_SEQLEN_ext(SEXP lkup,
		SEQLEN_FASTQloaderExt *loader_ext)
{
	FASTQloader loader;

	loader.new_seqid_hook = NULL;
	loader.new_empty_seq_hook = FASTQ_SEQLEN_new_empty_seq_hook;
	loader.append_seq_hook = FASTQ_SEQLEN_append_seq_hook;
	loader.new_qualid_hook = NULL;
	loader.new_empty_qual_hook = NULL;
	loader.append_qual_hook = NULL;
	if (lkup == R_NilValue) {
		loader.lkup = NULL;
		loader.lkup_len = 0;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_len = LENGTH(lkup);
	}
	loader.ext = loader_ext;
	return loader;
}

/*
 * The FASTQ loader.
 * Unlike the FASTQ SEQLEN loader, the FASTQ loader returns an error if a
 * read sequence is invalid or if a quality sequence is longer than the
 * corresponding read sequence. Note that a quality sequence shorter than
 * the corresponding read sequence is not considered an error. It's padded
 * to the length of the read with whatever bytes were present in the
 * BStringSet object where the quality sequences are copied (this object
 * is pre-allocated and its geometry always reflects the read lengths).
 * Unfortunately this is likely to cause problems downstream.
 */

typedef struct fastq_loader_ext {
	CharAEAE *seqid_buf;
	XVectorList_holder seq_holder;
	Chars_holder seq_elt_holder;
	int nseq;
	XVectorList_holder qual_holder;
	Chars_holder qual_elt_holder;
	int nqual;
} FASTQloaderExt;

static FASTQloaderExt new_FASTQloaderExt(SEXP sequences, SEXP qualities)
{
	FASTQloaderExt loader_ext;

	loader_ext.seqid_buf =
		new_CharAEAE(_get_XStringSet_length(sequences), 0);
	loader_ext.seq_holder = hold_XVectorList(sequences);
	loader_ext.nseq = 0;
	if (qualities != R_NilValue) {
		loader_ext.qual_holder = hold_XVectorList(qualities);
		loader_ext.nqual = 0;
	}
	return loader_ext;
}

static void FASTQ_new_seqid_hook(FASTQloader *loader,
				 const Chars_holder *seqid)
{
	FASTQloaderExt *loader_ext;
	CharAEAE *seqid_buf;

	loader_ext = loader->ext;
	seqid_buf = loader_ext->seqid_buf;
	// This works only because seqid->ptr is nul-terminated!
	CharAEAE_append_string(seqid_buf, seqid->ptr);
	return;
}

static void FASTQ_new_empty_seq_hook(FASTQloader *loader)
{
	FASTQloaderExt *loader_ext;
	Chars_holder *seq_elt_holder;

	loader_ext = loader->ext;
	seq_elt_holder = &(loader_ext->seq_elt_holder);
	*seq_elt_holder = get_elt_from_XRawList_holder(
					&(loader_ext->seq_holder),
					loader_ext->nseq);
	seq_elt_holder->length = 0;
	loader_ext->nseq++;
	return;
}

/* Unlike FASTQ_SEQLEN_append_seq_hook(), FASTQ_append_seq_hook() does check
   that the read sequence is valid. */
static const char *FASTQ_append_seq_hook(FASTQloader *loader,
		Chars_holder *seq_data)
{
	FASTQloaderExt *loader_ext;
	Chars_holder *seq_elt_holder;
	int ninvalid;

	loader_ext = loader->ext;
	seq_elt_holder = &(loader_ext->seq_elt_holder);
	if (loader->lkup != NULL) {
		ninvalid = translate(seq_data,
				     loader->lkup,
				     loader->lkup_len);
		if (ninvalid != 0)
			return "read sequence contains invalid letters";
	}
	append_Chars_holder(seq_elt_holder, seq_data);
	return NULL;
}

static void FASTQ_new_empty_qual_hook(FASTQloader *loader)
{
	FASTQloaderExt *loader_ext;
	Chars_holder *qual_elt_holder;

	loader_ext = loader->ext;
	qual_elt_holder = &(loader_ext->qual_elt_holder);
	*qual_elt_holder = get_elt_from_XRawList_holder(
					&(loader_ext->qual_holder),
					loader_ext->nqual);
	qual_elt_holder->length = 0;
	loader_ext->nqual++;
	return;
}

/* Check that the quality sequence is not longer than the read sequence.
   This prevents writing beyond the pre-allocated memory for the quality
   sequence. */
static const char *FASTQ_append_qual_hook(FASTQloader *loader,
		const Chars_holder *qual_data)
{
	FASTQloaderExt *loader_ext;
	Chars_holder *seq_elt_holder, *qual_elt_holder;

	loader_ext = loader->ext;
	seq_elt_holder = &(loader_ext->seq_elt_holder);
	qual_elt_holder = &(loader_ext->qual_elt_holder);
	if (qual_elt_holder->length + qual_data->length >
	    seq_elt_holder->length)
	{
		return "quality sequence is longer than read sequence";
	}
	append_Chars_holder(qual_elt_holder, qual_data);
	return NULL;
}

static FASTQloader new_FASTQloader(int load_seqids, int load_quals,
		SEXP lkup, FASTQloaderExt *loader_ext)
{
	FASTQloader loader;

	loader.new_seqid_hook = load_seqids ? &FASTQ_new_seqid_hook : NULL;
	loader.new_empty_seq_hook = FASTQ_new_empty_seq_hook;
	loader.append_seq_hook = FASTQ_append_seq_hook;
	if (load_quals) {
		/* Quality ids are always ignored for now. */
		//loader.new_qualid_hook = &FASTQ_new_qualid_hook;
		loader.new_qualid_hook = NULL;
		loader.new_empty_qual_hook = &FASTQ_new_empty_qual_hook;
		loader.append_qual_hook = &FASTQ_append_qual_hook;
	} else {
		loader.new_qualid_hook = NULL;
		loader.new_empty_qual_hook = NULL;
		loader.append_qual_hook = NULL;
	}
	if (lkup == R_NilValue) {
		loader.lkup = NULL;
		loader.lkup_len = 0;
	} else {
		loader.lkup = INTEGER(lkup);
		loader.lkup_len = LENGTH(lkup);
	}
	loader.ext = loader_ext;
	return loader;
}

/* Ignore empty lines. */
static const char *parse_FASTQ_file(SEXP filexp,
		int nrec, int skip, int seek_first_rec,
		FASTQloader *loader,
		int *recno, long long int *offset)
{
	int lineno, EOL_in_buf, EOL_in_prev_buf, ret_code, nbyte_in,
	    FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineinrecno, dont_load;
	Chars_holder data;
	long long int prev_offset;
	const char *errmsg;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineinrecno = 0;
	lineno = 0;
	EOL_in_buf = 1;
	while (1) {
		if (EOL_in_buf)
			lineno++;
		EOL_in_prev_buf = EOL_in_buf;
		ret_code = filexp_gets(filexp, iobuf, IOBUF_SIZE, &EOL_in_buf);
		if (ret_code == 0)
			break;
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (EOL_in_buf) {
			nbyte_in = strlen(iobuf);
			data.length =
				delete_trailing_LF_or_CRLF(iobuf, nbyte_in);
		} else {
			data.length = nbyte_in = IOBUF_SIZE - 1;
		}
		prev_offset = *offset;
		*offset += nbyte_in;
		if (seek_first_rec) {
			if (EOL_in_prev_buf
			 && has_prefix(iobuf, FASTQ_line1_markup)) {
				seek_first_rec = 0;
			} else {
				continue;
			}
		}
		data.ptr = iobuf;
		if (EOL_in_prev_buf) {
			if (data.length == 0)
				continue;  // we ignore empty lines
			lineinrecno++;
			if (lineinrecno > 4)
				lineinrecno = 1;
		}
		if (!EOL_in_buf && (lineinrecno == 1 || lineinrecno == 3)) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, "
				 "line is too long", lineno);
			return errmsg_buf;
		}
		iobuf[data.length] = '\0';
		errmsg = NULL;
		switch (lineinrecno) {
		    case 1:
			if (nrec >= 0 && *recno >= skip + nrec) {
				/* Calls to filexp_seek() are costly on
				   compressed files and the cost increases as
				   we advance in the file. This is not a
				   problem when reading the entire file but
				   becomes one when reading a compressed file
				   by chunk. */
				filexp_seek(filexp, prev_offset, SEEK_SET);
				*offset = prev_offset;
				return NULL;
			}
			if (!has_prefix(iobuf, FASTQ_line1_markup)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
				    "\"%s\" expected at beginning of line %d",
				    FASTQ_line1_markup, lineno);
				return errmsg_buf;
			}
			dont_load = *recno < skip || loader == NULL;
			if (dont_load || loader->new_seqid_hook == NULL)
				continue;
			data.ptr += FASTQ_line1_markup_length;
			data.length -= FASTQ_line1_markup_length;
			loader->new_seqid_hook(loader, &data);
			continue;
		    case 2:
			if (dont_load || loader->new_empty_seq_hook == NULL)
				continue;
			if (EOL_in_prev_buf)
				loader->new_empty_seq_hook(loader);
			if (loader->append_seq_hook != NULL)
				errmsg = loader->append_seq_hook(loader, &data);
			break;
		    case 3:
			if (!has_prefix(iobuf, FASTQ_line3_markup)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
				    "\"%s\" expected at beginning of line %d",
				    FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (dont_load || loader->new_qualid_hook == NULL)
				continue;
			data.ptr += FASTQ_line3_markup_length;
			data.length -= FASTQ_line3_markup_length;
			loader->new_qualid_hook(loader, &data);
			continue;
		    case 4:
			if (EOL_in_buf)
				(*recno)++;
			if (dont_load || loader->new_empty_qual_hook == NULL)
				continue;
			if (EOL_in_prev_buf)
				loader->new_empty_qual_hook(loader);
			if (loader->append_qual_hook != NULL)
				errmsg = loader->append_qual_hook(loader,
								  &data);
			break;
		}
		if (errmsg != NULL) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "line %d: %s", lineno, errmsg);
			return errmsg_buf;
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
	loader = new_FASTQloader_with_SEQLEN_ext(R_NilValue, &loader_ext);
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
	const char *errmsg;

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
	loader_ext = new_FASTQloaderExt(sequences, qualities);
	loader = new_FASTQloader(load_seqids, load_quals, lkup, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		/* Calls to filexp_tell() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		offset = filexp_tell(filexp);
		errmsg = parse_FASTQ_file(filexp, nrec0, skip0, seek_rec0,
					  &loader,
					  &recno, &offset);
		if (errmsg != NULL) {
			UNPROTECT(load_quals ? 3 : 2);
			error("reading FASTQ file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(filexp_list), i)),
			      errmsg_buf);
		}
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

static void write_FASTQ_seq(SEXP filexp, const Chars_holder X,
		const int *lkup0, int lkup_len)
{
	// write sequences in chunks in case sequence is longer than I/O buffer
	int i1, i2, bytes_to_write, bytes_remaining;
	bytes_remaining = X.length;
	i1 = 0;
	while(bytes_remaining){
		bytes_to_write = bytes_remaining >= IOBUF_SIZE ? IOBUF_SIZE-1 : bytes_remaining;
		i2 = i1+bytes_to_write-1;
		Ocopy_bytes_from_i1i2_with_lkup(i1, i2,
			iobuf, bytes_to_write,
			X.ptr, X.length,
			lkup0, lkup_len);
		iobuf[bytes_to_write] = 0;
		filexp_puts(filexp, iobuf);
		i1 = i2 + 1;
		bytes_remaining -= bytes_to_write;
	}
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
	int x_length, lkup_len, i;
	const int *lkup0;
	SEXP filexp, x_names, q_names;
	const char *id;
	Chars_holder X_elt;

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
		lkup_len = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_len = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		id = get_FASTQ_rec_id(x_names, q_names, i);
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		write_FASTQ_id(filexp, FASTQ_line1_markup, id);
		write_FASTQ_seq(filexp, X_elt, lkup0, lkup_len);
		write_FASTQ_id(filexp, FASTQ_line3_markup, id);
		if (qualities != R_NilValue) {
			write_FASTQ_qual(filexp, X_elt.length, &Q, i);
		} else {
			write_FASTQ_fakequal(filexp, X_elt.length);
		}
	}
	return R_NilValue;
}

