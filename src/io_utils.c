/****************************************************************************
 *                           I/O low-level utils                            *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdlib.h> /* for malloc(), free() */
#include <ctype.h> /* for isspace() */

#include <sys/types.h>
#include <sys/stat.h> /* for fstat() */


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

static FILE *open_file(const char *path, const char *mode)
{
	FILE *fp;
	int ret;
	struct stat buf;

	fp = fopen(path, mode);
	if (fp == NULL)
		error("cannot open file '%s'", path);
	ret = fstat(fileno(fp), &buf);
	if (ret != 0) {
		fclose(fp);
		error("Biostrings internal error in open_file(): "
		      "cannot stat file '%s'", path);
	}
	if (S_ISDIR(buf.st_mode)) {
		fclose(fp);
		error("file '%s' is a directory", path);
	}
	return fp;
}


/****************************************************************************
 * File_holder structs
 */

typedef struct file_holder {
	const char *path;
	const char *expath;
	const char *mode;
	int ztype;
	int subtype;
	void *file;
} File_holder;

static File_holder new_File_holder(const char *path, const char *expath,
		const char *mode, const char *compress, int level)
{
	File_holder file_holder;
	int ztype, subtype;
	void *file;

	file_holder.path = path;
	file_holder.expath = expath;
	file_holder.mode = mode;
	if (strcmp(mode, "r") == 0) {
		get_file_ztype(expath, &ztype, &subtype);
		file_holder.ztype = ztype;
		file_holder.subtype = subtype;
		switch (ztype) {
		/* No compression */
		    case -1:
		/* gzfile */
		    case 0:
			file = gzopen(expath, "r");
			if (file == NULL)
				error("cannot open file '%s'", expath);
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
			//opening/closing this type of file seems quite
			//complicated
			break;
		    default:
			error("Biostrings internal error in "
			      "new_File_holder(): "
			      "invalid ztype value %d", ztype);
		}
	} else {
		file_holder.ztype = -1;
		file = open_file(expath, mode);
	}
	file_holder.file = file;
	return file_holder;
}

static void File_holder_close(const File_holder *file_holder)
{
	const char *mode;
	int ztype;
	void *file;

	mode = file_holder->mode;
	ztype = file_holder->ztype;
	file = file_holder->file;
	if (strcmp(mode, "r") == 0) {
		switch (ztype) {
		/* No compression */
		    case -1:
		/* gzfile */
		    case 0:
			gzclose((gzFile) file);
			break;
		    default:
			error("Biostrings internal error "
			      "in File_holder_close(): "
			      "invalid ztype value %d", ztype);
		}
	} else {
		fclose((FILE *) file);
	}
	return;
}

/*
 * Similar to fgets()/gzgets(), except that it returns a code instead of
 * NULL/Z_NULL or a pointer to the buffer. The returned code can be:
 *    2: if reading stopped after an EOF or a newline,
 *    1: if reading stopped because buffer was full and no EOF or newline was
 *       read in,
 *    0: if end of file occurred while no character was read,
 *   -1: on read error.
 */
static int File_holder_gets(const File_holder *file_holder,
		char *buf, int buf_size, int *EOL_in_buf)
{
	int ztype;
	void *file;

	buf[buf_size - 1] = 'N'; // any non '\0' would do

	ztype = file_holder->ztype;
	file = file_holder->file;
	switch (ztype) {
	/* No compression */
	    case -1:
	/* gzfile */
	    case 0:
		if (gzgets((gzFile) file, buf, buf_size) == Z_NULL) {
			//if (ferror(file) != 0 || feof(file) == 0)
			//	return -1;
			return 0;
		}
		break;
	    default:
		error("Biostrings internal error in File_holder_gets(): "
		      "invalid ztype value %d", ztype);
	}
	*EOL_in_buf = buf[buf_size - 1] == 'N' || buf[buf_size - 2] == '\n';
	return *EOL_in_buf ? 2 : 1;
}

static void File_holder_rewind(const File_holder *file_holder)
{
	int ztype;
	void *file;

	ztype = file_holder->ztype;
	file = file_holder->file;
	switch (ztype) {
	/* No compression */
	    case -1:
	/* gzfile */
	    case 0:
		gzrewind((gzFile) file);
		break;
	    default:
		error("Biostrings internal error in File_holder_rewind(): "
		      "invalid ztype value %d", ztype);
	}
	return;
}

static int File_holder_puts(const File_holder *file_holder, const char *s)
{
	int ztype, n;
	void *file;

	ztype = file_holder->ztype;
	file = file_holder->file;
	switch (ztype) {
	/* No compression */
	    case -1:
		if ((n = fputs(s, (FILE *) file)) >= 0)
			return n;
		break;
	/* gzfile */
	    case 0:
		if ((n = gzputs((gzFile) file, s)) >= 0)
			return n;
		break;
	    default:
		error("Biostrings internal error in File_holder_puts(): "
		      "invalid ztype value %d", ztype);
	}
	error("write error");
}

static void File_holder_putc(const File_holder *file_holder, int c)
{
	int ztype;
	void *file;

	ztype = file_holder->ztype;
	file = file_holder->file;
	switch (ztype) {
	/* No compression */
	    case -1:
		if (fputc(c, (FILE *) file) != EOF)
			return;
		break;
	/* gzfile */
	    case 0:
		if (gzputc((gzFile) file, c) != -1)
			return;
		break;
	    default:
		error("Biostrings internal error in File_holder_putc(): "
		      "invalid ztype value %d", ztype);
	}
	error("write error");
}


/****************************************************************************
 * Low-level manipulation of External File Pointers
 */

int ExternalFilePtr_gets(SEXP efp, char *buf, int buf_size, int *EOL_in_buf)
{
	return File_holder_gets(R_ExternalPtrAddr(efp),
				buf, buf_size, EOL_in_buf);
}

void ExternalFilePtr_rewind(SEXP efp)
{
	return File_holder_rewind(R_ExternalPtrAddr(efp));
}

int ExternalFilePtr_puts(SEXP efp, const char *s)
{
	return File_holder_puts(R_ExternalPtrAddr(efp), s);
}

void ExternalFilePtr_putc(SEXP efp, int c)
{
	return File_holder_putc(R_ExternalPtrAddr(efp), c);
}

static SEXP new_ExternalFilePtr(SEXP filepath,
		const char *mode, const char *compress, int level)
{
	SEXP filepath0, ans, string;
	const char *expath;
	File_holder file_holder, *ans_p;

	if (!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' is NA");
	expath = R_ExpandFileName(translateChar(filepath0));
	file_holder = new_File_holder(CHAR(filepath0), expath,
				      mode, compress, level);
	ans_p = (File_holder *) malloc(sizeof(File_holder));
	if (ans_p == NULL)
		error("Biostrings internal error in new_ExternalFilePtr(): "
		      "malloc() failed");
	*ans_p = file_holder;
	PROTECT(ans = R_MakeExternalPtr(ans_p, R_NilValue, R_NilValue));
	PROTECT(string = mkString(expath));
	setAttrib(ans, install("expath"), string);
	UNPROTECT(2);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Returns an external pointer.
 * From R:
 *   x <- .Call("new_input_ExternalFilePtr", "path/to/some/file",
 *              PACKAGE="Biostrings")
 *   reg.finalizer(x,
 *       function(e) .Call("finalize_ExternalFilePtr", e, PACKAGE="Biostrings"),
 *       onexit=TRUE)
 */
SEXP new_input_ExternalFilePtr(SEXP filepath)
{
	return new_ExternalFilePtr(filepath, "r", NULL, 0);
}

/* --- .Call ENTRY POINT ---
 * Returns an external pointer.
 */
SEXP new_output_ExternalFilePtr(SEXP filepath, SEXP append,
		SEXP compress, SEXP compression_level)
{
	const char *mode, *compress0;
	int level;

	mode = LOGICAL(append)[0] ? "a" : "w";
	compress0 = CHAR(STRING_ELT(compress, 0));
	level = INTEGER(compression_level)[0];
	return new_ExternalFilePtr(filepath, mode, compress0, level);
}

/* --- .Call ENTRY POINT ---
 * Closes the file pointed to by e.
 */
SEXP finalize_ExternalFilePtr(SEXP efp)
{
	File_holder *file_holder;

	file_holder = R_ExternalPtrAddr(efp);
	if (file_holder != NULL) {
		File_holder_close(file_holder);
		free(file_holder);
		R_SetExternalPtrAddr(efp, NULL);
	}
	return R_NilValue;
}


/****************************************************************************
 * Other stuff.
 */

/*
 * Doesn't actually delete anything but returns the length the 'buf' char
 * array would have after deletion of the LF ("\n") or CRLF ("\r\n") chars
 * located at its end.
 * If 'buf_len' is -1, then 'buf' must be a C string (i.e. null-terminated).
 */
int delete_trailing_LF_or_CRLF(const char *buf, int buf_len)
{
	if (buf_len == -1)
		buf_len = strlen(buf);
	if (buf_len == 0)
		return 0;
	if (buf[--buf_len] != '\n')
		return ++buf_len;
	if (buf_len == 0)
		return 0;
	if (buf[--buf_len] != '\r')
		return ++buf_len;
	return buf_len;
}

