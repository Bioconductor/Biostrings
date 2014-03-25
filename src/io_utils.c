/****************************************************************************
 *                           I/O low-level utils                            *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdlib.h> /* for malloc(), free() */
#include <sys/stat.h> /* for fstat() */

#include <zlib.h>
#ifndef _WIN32
#include <bzlib.h>
#endif

#define UNCOMPRESSED		0
#define GZ_TYPE			1
#define BZ_TYPE			2 /* detected but not supported yet */
#define XZ_TYPE			3 /* detected but not supported yet */

/* Code taken from do_url() in R/src/main/connections.c */
static void get_file_ztype(const char *path, int *ztype, int *subtype)
{
	FILE *fp;
	char buf[7];
	int res;

	*ztype = UNCOMPRESSED;
	*subtype = 0;
	if ((fp = fopen(path, "rb")) == NULL)
		return;
	memset(buf, 0, 7);
	res = fread(buf, 5, 1, fp);
	fclose(fp);
	if (res != 1)
		return;
	if (buf[0] == '\x1f' && buf[1] == '\x8b')
		*ztype = GZ_TYPE;
	else if (strncmp(buf, "BZh", 3) == 0)
		*ztype = BZ_TYPE;
	else if (buf[0] == '\xFD' && strncmp(buf+1, "7zXZ", 4) == 0)
		*ztype = XZ_TYPE;
	else if ((buf[0] == '\xFF') && strncmp(buf+1, "LZMA", 4) == 0) {
		*ztype = XZ_TYPE;
		*subtype = 1;
	} else if (memcmp(buf, "]\0\0\200\0", 5) == 0) {
		*ztype = XZ_TYPE;
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
		error(INTERNAL_ERR_IN "open_file(): "
		      "cannot stat file '%s'", path);
	}
	if (S_ISDIR(buf.st_mode)) {
		fclose(fp);
		error("file '%s' is a directory", path);
	}
	return fp;
}


/****************************************************************************
 * ZFile structs
 */

typedef struct zfile {
	const char *path;
	const char *expath;
	const char *mode;
	int ztype;
	int subtype;
	void *file;
} ZFile;

static ZFile new_ZFile(const char *path, const char *expath,
		const char *mode, const char *compress, int level)
{
	ZFile zfile;
	int ztype, subtype;
	void *file;

	zfile.path = path;
	zfile.expath = expath;
	zfile.mode = mode;
	if (strcmp(mode, "r") == 0) {
		/* Open file for reading. */
		get_file_ztype(expath, &ztype, &subtype);
		zfile.ztype = ztype;
		zfile.subtype = subtype;
		switch (ztype) {
		    case UNCOMPRESSED:
		    case GZ_TYPE:
			file = gzopen(expath, "r");
			break;
		    case BZ_TYPE:
#ifndef _WIN32
			file = BZ2_bzopen(expath, "rb");
#else
			error("cannot open file '%s'\n"
                              "  bzip2-compressed files are not supported "
                              "on Windows, sorry!", expath);
#endif
			break;
		    case XZ_TYPE:
#ifndef _WIN32
			error("cannot open file '%s'\n"
                              "  LZMA-compressed files are not supported "
                              "yet, sorry!", expath);
			//requires #include <lzma.h>
			//opening/closing this type of file seems quite
			//complicated
#else
			error("cannot open file '%s'\n"
                              "  LZMA-compressed files are not supported "
                              "on Windows, sorry!", expath);
#endif
			break;
		    default:
			error(INTERNAL_ERR_IN "new_ZFile(): "
			      "invalid ztype value %d", ztype);
		}
		if (file == NULL)
			error("cannot open file '%s'", expath);
	} else {
		/* Open file for writing or appending. */
		zfile.ztype = UNCOMPRESSED;
		file = open_file(expath, mode);
	}
	zfile.file = file;
	return zfile;
}

static void ZFile_close(const ZFile *zfile)
{
	const char *mode;
	int ztype;
	void *file;

	mode = zfile->mode;
	ztype = zfile->ztype;
	file = zfile->file;
	if (strcmp(mode, "r") == 0) {
		switch (ztype) {
		    case UNCOMPRESSED:
		    case GZ_TYPE:
			gzclose((gzFile) file);
			break;
#ifndef _WIN32
		    case BZ_TYPE:
			BZ2_bzclose((BZFILE *) file);
			break;
#endif
		    default:
			error(INTERNAL_ERR_IN "ZFile_close(): "
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
static int ZFile_gets(const ZFile *zfile,
		char *buf, int buf_size, int *EOL_in_buf)
{
	int max_buf_len, ztype, i;
	void *file;
	char *buf_p;

	max_buf_len = buf_size - 1;
	buf[max_buf_len] = 'N'; // any non '\0' would do

	ztype = zfile->ztype;
	file = zfile->file;
	switch (ztype) {
	    case UNCOMPRESSED:
	    case GZ_TYPE:
		if (gzgets((gzFile) file, buf, buf_size) == Z_NULL) {
			//if (ferror(file) != 0 || feof(file) == 0)
			//	return -1;
			return 0;
		}
		break;
#ifndef _WIN32
	    case BZ_TYPE:
		for (i = 0, buf_p = buf; i < max_buf_len; i++, buf_p++) {
			if (BZ2_bzread((BZFILE *) file, buf_p, 1) == 0) {
				if (i == 0)
					return 0;
				break;
			}
			if (*buf_p == '\n') {
				buf_p++;
				break;
			}
		}
		*buf_p = '\0';
		break;
#endif
	    default:
		error(INTERNAL_ERR_IN "ZFile_gets(): "
		      "invalid ztype value %d", ztype);
	}

	*EOL_in_buf = buf[max_buf_len] == 'N' || buf[max_buf_len - 1] == '\n';
	return *EOL_in_buf ? 2 : 1;
}

static void ZFile_rewind(ZFile *zfile)
{
	int ztype;
	void *file;

	ztype = zfile->ztype;
	file = zfile->file;
	switch (ztype) {
	    case UNCOMPRESSED:
	    case GZ_TYPE:
		gzrewind((gzFile) file);
		break;
#ifndef _WIN32
	    case BZ_TYPE:
		/* No rewind in the bzip2 library. So we close and re-open
		   the file. */
		BZ2_bzclose((BZFILE *) file);
		file = BZ2_bzopen(zfile->expath, "rb");
		if (file == NULL)
			error(INTERNAL_ERR_IN "ZFile_rewind(): "
			      "cannot re-open file '%s'", zfile->expath);
		zfile->file = file;
		break;
#endif
	    default:
		error(INTERNAL_ERR_IN "ZFile_rewind(): "
		      "invalid ztype value %d", ztype);
	}
	return;
}

static int ZFile_puts(const ZFile *zfile, const char *s)
{
	int ztype, n;
	void *file;

	ztype = zfile->ztype;
	file = zfile->file;
	switch (ztype) {
	    case UNCOMPRESSED:
		if ((n = fputs(s, (FILE *) file)) >= 0)
			return n;
		break;
	    case GZ_TYPE:
		if ((n = gzputs((gzFile) file, s)) >= 0)
			return n;
		break;
	    default:
		error(INTERNAL_ERR_IN "ZFile_puts(): "
		      "invalid ztype value %d", ztype);
	}
	error("write error");
}

static void ZFile_putc(const ZFile *zfile, int c)
{
	int ztype;
	void *file;

	ztype = zfile->ztype;
	file = zfile->file;
	switch (ztype) {
	    case UNCOMPRESSED:
		if (fputc(c, (FILE *) file) != EOF)
			return;
		break;
	    case GZ_TYPE:
		if (gzputc((gzFile) file, c) != -1)
			return;
		break;
	    default:
		error(INTERNAL_ERR_IN "ZFile_putc(): "
		      "invalid ztype value %d", ztype);
	}
	error("write error");
}


/****************************************************************************
 * Low-level manipulation of External File Pointers
 */

int ExternalFilePtr_gets(SEXP efp, char *buf, int buf_size, int *EOL_in_buf)
{
	return ZFile_gets(R_ExternalPtrAddr(efp),
				buf, buf_size, EOL_in_buf);
}

void ExternalFilePtr_rewind(SEXP efp)
{
	return ZFile_rewind(R_ExternalPtrAddr(efp));
}

int ExternalFilePtr_puts(SEXP efp, const char *s)
{
	return ZFile_puts(R_ExternalPtrAddr(efp), s);
}

void ExternalFilePtr_putc(SEXP efp, int c)
{
	return ZFile_putc(R_ExternalPtrAddr(efp), c);
}

static SEXP new_ExternalFilePtr(SEXP filepath,
		const char *mode, const char *compress, int level)
{
	SEXP filepath0, ans, string;
	const char *expath;
	ZFile zfile, *extptraddr;

	if (!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' is NA");
	expath = R_ExpandFileName(translateChar(filepath0));
	zfile = new_ZFile(CHAR(filepath0), expath, mode, compress, level);
	extptraddr = (ZFile *) malloc(sizeof(ZFile));
	if (extptraddr == NULL) {
		ZFile_close(&zfile);
		error(INTERNAL_ERR_IN "new_ExternalFilePtr(): "
		      "malloc() failed");
	}
	*extptraddr = zfile;
	PROTECT(ans = R_MakeExternalPtr(extptraddr, R_NilValue, R_NilValue));
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
	ZFile *zfile;

	zfile = R_ExternalPtrAddr(efp);
	if (zfile != NULL) {
		ZFile_close(zfile);
		free(zfile);
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

