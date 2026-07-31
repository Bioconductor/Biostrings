/****************************************************************************
 *                          Read/write FASTA files                          *
 *                                 --------                                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

#include <math.h>  /* for llround */


#define IOBUF_SIZE 200002
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
 * FASTA parser
 */

/* Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be ready to support longer sympols just in case. */
static const char *FASTA_comment_markup = ";", *FASTA_desc_markup = ">";

typedef struct fasta_loader {
	void (*new_desc_hook)(struct fasta_loader *loader,
			      int recno, long long int offset,
			      const Chars_holder *desc_line);
	void (*new_empty_seq_hook)(struct fasta_loader *loader);
	void (*append_seq_hook)(struct fasta_loader *loader,
				const Chars_holder *seq_data);
	const int *lkup;
	int lkup_len;
	void *ext;  /* loader extension (optional) */
} FASTAloader;

/*
 * The FASTA INDEX loader.
 * Used in parse_FASTA_file() to build the FASTA index only.
 */

typedef struct index_fasta_loader_ext {
	IntAE *recno_buf;
	LLongAE *offset_buf;
	CharAEAE *desc_buf;
	IntAE *seqlength_buf;
} INDEX_FASTAloaderExt;

static INDEX_FASTAloaderExt new_INDEX_FASTAloaderExt(void)
{
	INDEX_FASTAloaderExt loader_ext;

	loader_ext.recno_buf = new_IntAE(0, 0, 0);
	loader_ext.offset_buf = new_LLongAE(0, 0, 0);
	loader_ext.desc_buf = new_CharAEAE(0, 0);
	loader_ext.seqlength_buf = new_IntAE(0, 0, 0);
	return loader_ext;
}

static void FASTA_INDEX_new_desc_hook(FASTAloader *loader,
		int recno, long long int offset,
		const Chars_holder *desc_line)
{
	INDEX_FASTAloaderExt *loader_ext;
	IntAE *recno_buf;
	LLongAE *offset_buf;
	CharAEAE *desc_buf;

	loader_ext = loader->ext;

	recno_buf = loader_ext->recno_buf;
	IntAE_insert_at(recno_buf, IntAE_get_nelt(recno_buf), recno + 1);

	offset_buf = loader_ext->offset_buf;
	LLongAE_insert_at(offset_buf, LLongAE_get_nelt(offset_buf), offset);

	desc_buf = loader_ext->desc_buf;
	// This works only because desc_line->ptr is nul-terminated!
	CharAEAE_append_string(desc_buf, desc_line->ptr);
	return;
}

static void FASTA_INDEX_new_empty_seq_hook(FASTAloader *loader)
{
	INDEX_FASTAloaderExt *loader_ext;
	IntAE *seqlength_buf;

	loader_ext = loader->ext;
	seqlength_buf = loader_ext->seqlength_buf;
	IntAE_insert_at(seqlength_buf, IntAE_get_nelt(seqlength_buf), 0);
	return;
}

static void FASTA_INDEX_append_seq_hook(FASTAloader *loader,
		const Chars_holder *seq_data)
{
	INDEX_FASTAloaderExt *loader_ext;
	IntAE *seqlength_buf;

	loader_ext = loader->ext;
	seqlength_buf = loader_ext->seqlength_buf;
	seqlength_buf->elts[IntAE_get_nelt(seqlength_buf) - 1] +=
		seq_data->length;
	return;
}

static FASTAloader new_FASTAloader_with_INDEX_ext(int load_descs, SEXP lkup,
		INDEX_FASTAloaderExt *loader_ext)
{
	FASTAloader loader;

	loader.new_desc_hook = load_descs ? &FASTA_INDEX_new_desc_hook : NULL;
	loader.new_empty_seq_hook = &FASTA_INDEX_new_empty_seq_hook;
	loader.append_seq_hook = &FASTA_INDEX_append_seq_hook;
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
 * The FASTA loader.
 */

typedef struct fasta_loader_ext {
	XVectorList_holder seq_holder;
	Chars_holder seq_elt_holder;
	int nseq;
} FASTAloaderExt;

static FASTAloaderExt new_FASTAloaderExt(SEXP sequences)
{
	FASTAloaderExt loader_ext;

	loader_ext.seq_holder = hold_XVectorList(sequences);
	loader_ext.seq_elt_holder.ptr = NULL;  // only for -Wuninitialized
	loader_ext.seq_elt_holder.length = 0;  // only for -Wuninitialized
	loader_ext.nseq = 0;
	return loader_ext;
}

static void FASTA_new_empty_seq_hook(FASTAloader *loader)
{
	FASTAloaderExt *loader_ext;
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

static void FASTA_append_seq_hook(FASTAloader *loader,
		const Chars_holder *seq_data)
{
	FASTAloaderExt *loader_ext;
	Chars_holder *seq_elt_holder;

	loader_ext = loader->ext;
	seq_elt_holder = &(loader_ext->seq_elt_holder);
	append_Chars_holder(seq_elt_holder, seq_data);
	return;
}

static FASTAloader new_FASTAloader(SEXP lkup, FASTAloaderExt *loader_ext)
{
	FASTAloader loader;

	loader.new_desc_hook = NULL;
	loader.new_empty_seq_hook = &FASTA_new_empty_seq_hook;
	loader.append_seq_hook = &FASTA_append_seq_hook;
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

/* Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
   the original Pearson FASTA format. */
static const char *parse_FASTA_file(SEXP filexp,
		int nrec, int skip, int seek_first_rec,
		FASTAloader *loader,
		int *recno, long long int *offset, long long int *ninvalid)
{
	int lineno, EOL_in_buf, EOL_in_prev_buf, ret_code, nbyte_in,
	    FASTA_desc_markup_length, dont_load, is_comment, is_desc;
	Chars_holder data;
	long long int prev_offset;

	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	lineno = 0;
	EOL_in_buf = 1;
	dont_load = -1;
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
			 && has_prefix(iobuf, FASTA_desc_markup)) {
				seek_first_rec = 0;
			} else {
				continue;
			}
		}
		data.ptr = iobuf;
		if (EOL_in_prev_buf) {
			if (data.length == 0)
				continue;  // we ignore empty lines
			is_comment = has_prefix(iobuf, FASTA_comment_markup);
			is_desc = has_prefix(iobuf, FASTA_desc_markup);
			if (!EOL_in_buf && (is_comment || is_desc)) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "cannot read line %d, "
					 "line is too long", lineno);
				return errmsg_buf;
			}
			if (is_comment)
				continue;  // we ignore comment lines
		}
		iobuf[data.length] = '\0';
		if (EOL_in_prev_buf && is_desc) {
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
			dont_load = *recno < skip || loader == NULL;
			if (!dont_load && loader->new_desc_hook != NULL) {
				data.ptr += FASTA_desc_markup_length;
				data.length -= FASTA_desc_markup_length;
				loader->new_desc_hook(loader,
						      *recno, prev_offset,
						      &data);
			}
			if (!dont_load && loader->new_empty_seq_hook != NULL)
				loader->new_empty_seq_hook(loader);
			(*recno)++;
			continue;
		}
		if (dont_load == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "\"%s\" expected at beginning of line %d",
				 FASTA_desc_markup, lineno);
			return errmsg_buf;
		}
		if (dont_load || loader->new_empty_seq_hook == NULL)
			continue;
		if (loader->append_seq_hook != NULL) {
			if (loader->lkup != NULL)
				*ninvalid += translate(&data,
						       loader->lkup,
						       loader->lkup_len);
			loader->append_seq_hook(loader, &data);
		}
	}
	if (seek_first_rec) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "no FASTA record found");
		return errmsg_buf;
	}
	return NULL;
}


/****************************************************************************
 * read_fasta_files()
 */

static SEXP get_fasta_seqlengths(SEXP filexp_list,
		int nrec, int skip, int seek_first_rec,
		int use_names, SEXP lkup)
{
	INDEX_FASTAloaderExt loader_ext;
	FASTAloader loader;
	int recno, i;
	SEXP filexp, ans, ans_names;
	const char *filename, *errmsg;
	long long int offset0, offset, ninvalid;

	loader_ext = new_INDEX_FASTAloaderExt();
	loader = new_FASTAloader_with_INDEX_ext(use_names, lkup, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		filename = CHAR(STRING_ELT(GET_NAMES(filexp_list), i));
		/* Calls to filexp_tell() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		offset0 = offset = filexp_tell(filexp);
		ninvalid = 0LL;
		errmsg = parse_FASTA_file(filexp, nrec, skip, seek_first_rec,
					  &loader,
					  &recno, &offset, &ninvalid);
		/* Calls to filexp_seek() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		filexp_seek(filexp, offset0, SEEK_SET);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      filename, errmsg_buf);
		if (ninvalid != 0LL)
			warning("reading FASTA file %s: ignored %lld "
				"invalid one-letter sequence codes",
				filename, ninvalid);
	}
	PROTECT(ans = new_INTEGER_from_IntAE(loader_ext.seqlength_buf));
	if (use_names) {
		PROTECT(ans_names =
			new_CHARACTER_from_CharAEAE(loader_ext.desc_buf));
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * We use a 2-pass algo so we can pre-alloc the exact amount of memory for the
 * XStringSet object that we load with the string data during the 2nd pass.
 * Because a XStringSet object cannot be grown in-place, a 1-pass algo would
 * need to load the string data in a growing object first (e.g. CharAEAE)
 * then turn it into an XStringSet object (this would require a copy of the
 * entire string data).
 */
SEXP read_fasta_files(SEXP filexp_list,
		      SEXP nrec, SEXP skip, SEXP seek_first_rec,
		      SEXP use_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, seek_rec0, use_names0, recno, i;
	SEXP seqlengths, ans, filexp;
	FASTAloaderExt loader_ext;
	FASTAloader loader;
	long long int offset, ninvalid;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	use_names0 = LOGICAL(use_names)[0];
	/* 1st pass */
	PROTECT(seqlengths = get_fasta_seqlengths(filexp_list,
						  nrec0, skip0, seek_rec0,
						  use_names0, lkup));
	/* Allocation */
	PROTECT(ans = _alloc_XStringSet(CHAR(STRING_ELT(elementType, 0)),
					seqlengths));
	/* 2nd pass */
	loader_ext = new_FASTAloaderExt(ans);
	loader = new_FASTAloader(lkup, &loader_ext);
	recno = 0;
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		/* Calls to filexp_tell() are costly on compressed files
		   and the cost increases as we advance in the file.
		   This is not a problem when reading the entire file but
		   becomes one when reading a compressed file by chunk. */
		offset = filexp_tell(filexp);
		parse_FASTA_file(filexp, nrec0, skip0, seek_rec0,
				 &loader,
				 &recno, &offset, &ninvalid);
	}
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * read_fasta_blocks()
 */

static SEXP make_fasta_index_data_frame(const IntAE *recno_buf,
					const IntAE *fileno_buf,
					const LLongAE *offset_buf,
					const CharAEAE *desc_buf,
					const IntAE *seqlength_buf)
{
	SEXP df, colnames, tmp;
	int i;

	PROTECT(df = NEW_LIST(5));

	PROTECT(colnames = NEW_CHARACTER(5));
	PROTECT(tmp = mkChar("recno"));
	SET_STRING_ELT(colnames, 0, tmp);
	UNPROTECT(1);
	PROTECT(tmp = mkChar("fileno"));
	SET_STRING_ELT(colnames, 1, tmp);
	UNPROTECT(1);
	PROTECT(tmp = mkChar("offset"));
	SET_STRING_ELT(colnames, 2, tmp);
	UNPROTECT(1);
	PROTECT(tmp = mkChar("desc"));
	SET_STRING_ELT(colnames, 3, tmp);
	UNPROTECT(1);
	PROTECT(tmp = mkChar("seqlength"));
	SET_STRING_ELT(colnames, 4, tmp);
	UNPROTECT(1);
	SET_NAMES(df, colnames);
	UNPROTECT(1);

	PROTECT(tmp = new_INTEGER_from_IntAE(recno_buf));
	SET_ELEMENT(df, 0, tmp);
	UNPROTECT(1);

	PROTECT(tmp = new_INTEGER_from_IntAE(fileno_buf));
	SET_ELEMENT(df, 1, tmp);
	UNPROTECT(1);

	PROTECT(tmp = NEW_NUMERIC(LLongAE_get_nelt(offset_buf)));
	for (i = 0; i < LENGTH(tmp); i++)
		REAL(tmp)[i] = (double) offset_buf->elts[i];
	SET_ELEMENT(df, 2, tmp);
	UNPROTECT(1);

	PROTECT(tmp = new_CHARACTER_from_CharAEAE(desc_buf));
	SET_ELEMENT(df, 3, tmp);
	UNPROTECT(1);

	PROTECT(tmp = new_INTEGER_from_IntAE(seqlength_buf));
	SET_ELEMENT(df, 4, tmp);
	UNPROTECT(1);

	/* list_as_data_frame() performs IN-PLACE coercion */
	list_as_data_frame(df, IntAE_get_nelt(recno_buf));
	UNPROTECT(1);
	return df;
}

/* --- .Call ENTRY POINT --- */
SEXP fasta_index(SEXP filexp_list,
		 SEXP nrec, SEXP skip, SEXP seek_first_rec, SEXP lkup)
{
	int nrec0, skip0, seek_rec0, i, recno, old_nrec, new_nrec, k;
	INDEX_FASTAloaderExt loader_ext;
	FASTAloader loader;
	IntAE *seqlength_buf, *fileno_buf;
	SEXP filexp;
	long long int offset, ninvalid;
	const char *errmsg;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	seek_rec0 = LOGICAL(seek_first_rec)[0];
	loader_ext = new_INDEX_FASTAloaderExt();
	loader = new_FASTAloader_with_INDEX_ext(1, lkup, &loader_ext);
	seqlength_buf = loader_ext.seqlength_buf;
	fileno_buf = new_IntAE(0, 0, 0);
	for (i = recno = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		offset = filexp_tell(filexp);
		ninvalid = 0LL;
		errmsg = parse_FASTA_file(filexp, nrec0, skip0, seek_rec0,
					  &loader,
					  &recno, &offset, &ninvalid);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      CHAR(STRING_ELT(GET_NAMES(filexp_list), i)),
			      errmsg_buf);
		if (ninvalid != 0LL)
			warning("reading FASTA file %s: ignored %lld "
				"invalid one-letter sequence codes",
				CHAR(STRING_ELT(GET_NAMES(filexp_list), i)),
				ninvalid);
		old_nrec = IntAE_get_nelt(fileno_buf);
		new_nrec = IntAE_get_nelt(seqlength_buf);
		for (k = old_nrec; k < new_nrec; k++)
			IntAE_insert_at(fileno_buf, k, i + 1);
	}
	return make_fasta_index_data_frame(loader_ext.recno_buf,
					   fileno_buf,
					   loader_ext.offset_buf,
					   loader_ext.desc_buf,
					   seqlength_buf);
}

/* --- .Call ENTRY POINT ---
 * "FASTA blocks" are groups of consecutive FASTA records.
 * Args:
 *   seqlengths:  Integer vector with 1 element per record to read. The
 *                elements must be placed in the order that the records are
 *                going to be read.
 *   filexp_list: A list of N "File External Pointers" (see src/io_utils.c in
 *                the XVector package) with 1 element per input file. Files are
 *                going to be accessed from first to last in the list.
 *   nrec_list:   A list of N integer vectors (1 element per input file).
 *                Each integer vector has 1 value per FASTA block, which is the
 *                number of records in the block. IMPORTANT: Even if not
 *                required, the blocks in each integer vector should preferably
 *                be placed in the same order as in the file. This ensures that
 *                the calls to filexp_seek() are always moving forward
 *                in the file and are therefore efficient (moving backward on a
 *                compressed file can be *extremely* slow).
 *   offset_list: A list of N numeric vectors (1 element per input file) with
 *                the same shape as 'nrec_list', i.e. each numeric vector has 1
 *                value per FASTA block. This value is the offset of the block
 *                (i.e. the offset of its first record) relative to the start
 *                of the file. Measured in bytes.
 *   elementType: The elementType of the XStringSet to return (its class is
 *                inferred from this).
 *   lkup:        Lookup table for encoding the incoming sequence bytes.
 */
SEXP read_fasta_blocks(SEXP seqlengths,
		SEXP filexp_list, SEXP nrec_list, SEXP offset_list,
		SEXP elementType, SEXP lkup)
{
	SEXP ans, filexp, nrec, offset;
	FASTAloaderExt loader_ext;
	FASTAloader loader;
	int i, j, nrec_j, recno;
	long long int offset_j, ninvalid;

	PROTECT(ans = _alloc_XStringSet(CHAR(STRING_ELT(elementType, 0)),
					seqlengths));
	loader_ext = new_FASTAloaderExt(ans);
	loader = new_FASTAloader(lkup, &loader_ext);
	for (i = 0; i < LENGTH(filexp_list); i++) {
		filexp = VECTOR_ELT(filexp_list, i);
		nrec = VECTOR_ELT(nrec_list, i);
		offset = VECTOR_ELT(offset_list, i);
		for (j = 0; j < LENGTH(nrec); j++) {
			nrec_j = INTEGER(nrec)[j];
			offset_j = llround(REAL(offset)[j]);
			filexp_seek(filexp, offset_j, SEEK_SET);
			recno = 0;
			ninvalid = 0LL;
			parse_FASTA_file(filexp, nrec_j, 0, 0,
					 &loader,
					 &recno, &offset_j, &ninvalid);
		}
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Writing FASTA files.
 */

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fasta(SEXP x, SEXP filexp_list, SEXP width, SEXP lkup)
{
	XStringSet_holder X;
	int x_length, width0, lkup_len, i, j1, j2, dest_nbytes;
	int nmax_write_bytes, nwrote_bytes;
	const int *lkup0;
	SEXP filexp, x_names, desc;
	Chars_holder X_elt;

	X = _hold_XStringSet(x);
	x_length = _get_length_from_XStringSet_holder(&X);
	filexp = VECTOR_ELT(filexp_list, 0);
	width0 = INTEGER(width)[0];
	nmax_write_bytes = width0 >= IOBUF_SIZE ? IOBUF_SIZE - 1 : width0;
	iobuf[width0] = 0;
	if (lkup == R_NilValue) {
		lkup0 = NULL;
		lkup_len = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_len = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		filexp_puts(filexp, FASTA_desc_markup);
		if (x_names != R_NilValue) {
			desc = STRING_ELT(x_names, i);
			if (desc == NA_STRING)
				error("'names(x)' contains NAs");
			filexp_puts(filexp, CHAR(desc));
		}
		filexp_puts(filexp, "\n");
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		j1 = 0;
		nwrote_bytes = 0;
		while(j1 < X_elt.length){
			j2 = j1 + nmax_write_bytes;
			if (j2 > X_elt.length)
				j2 = X_elt.length;
			dest_nbytes = j2 - j1;
			if(nwrote_bytes + dest_nbytes > width0){
				// align write to size width0
				dest_nbytes = width0 - nwrote_bytes;
				j2 = j1 + dest_nbytes;
			}
			j2--;
			Ocopy_bytes_from_i1i2_with_lkup(j1, j2,
				iobuf, dest_nbytes,
				X_elt.ptr, X_elt.length,
				lkup0, lkup_len);
			iobuf[dest_nbytes] = 0;
			nwrote_bytes += dest_nbytes;
			j1 = j2+1;
			filexp_puts(filexp, iobuf);
			if(j1 == X_elt.length || nwrote_bytes == width0){
				// add newline if at end of line (width0) or end of sequence
				filexp_puts(filexp, "\n");
				nwrote_bytes = 0;
			}
		}
	}
	return R_NilValue;
}
