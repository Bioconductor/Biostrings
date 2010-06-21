/****************************************************************************
 *                           I/O low-level utils                            *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

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

static FILE *open_input_file(const char *path)
{
	FILE *stream;
	int ret, ztype, subtype;
	struct stat buf;

	get_file_ztype(path, &ztype, &subtype);
	switch (ztype) {
	/* No compression */
	    case -1:
		stream = fopen(path, "r");
		if (stream == NULL)
			error("cannot open file '%s'", path);
		ret = fstat(fileno(stream), &buf);
		if (ret != 0)
			error("Biostrings internal error in open_input_file(): "
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
		error("Biostrings internal error in open_input_file(): ",
		      "invalid ztype value %d", ztype);
	}
	return stream;
}

/* --- .Call ENTRY POINT ---
 * Returns an external pointer.
 * From R:
 *   x <- .Call("new_ExternalFilePtr", "path/to/some/file",
 *              PACKAGE="Biostrings")
 *   reg.finalizer(x,
 *       function(e) .Call("ExternalFilePtr_close", e, PACKAGE="Biostrings"),
 *       onexit=TRUE)
 */
SEXP new_ExternalFilePtr(SEXP filepath)
{
	SEXP filepath_elt;
	const char *path;
	FILE *stream;

	if (!IS_CHARACTER(filepath) || LENGTH(filepath) == 0)
		error("'filepath' must be a non-empty character vector");
	filepath_elt = STRING_ELT(filepath, 0);
	if (filepath_elt == NA_STRING)
		error("'filepath' is NA");
	// Maybe do this in R, not here.
	path = R_ExpandFileName(translateChar(filepath_elt));
	stream = open_input_file(path);
	return R_MakeExternalPtr(stream, R_NilValue, R_NilValue);
}

/* --- .Call ENTRY POINT ---
 * Closes the file pointed to by e.
 */
SEXP ExternalFilePtr_close(SEXP e)
{
	if (R_ExternalPtrAddr(e) != NULL) {
		//Rprintf("closing file...\n");
		fclose(R_ExternalPtrAddr(e));
		R_SetExternalPtrAddr(e, NULL);
	}
	return R_NilValue;
}


/****************************************************************************
 * Other stuff.
 */

/*
 * Right-trims the string stored into the buffer pointed to by linebuf.
 * If size is -1, then the size of the buffer is determined by a call to
 * strlen().
 * Return the number of chars that remain in the buffer after we've removed
 * the right spaces ('\n', '\r', '\t', ' ', etc...).
 */
int rtrimline(char *linebuf, int size)
{
	if (size == -1)
		size = strlen(linebuf);
	size--;
	while (size >= 0 && isspace(linebuf[size])) size--;
	linebuf[++size] = 0;
	return size;
}

/* Like fgets() except that:
 *   - The string stored into the buffer pointed to by s is right-trimmed i.e.
 *     all the rightmost white-space characters are removed.
 *   - It returns the length of s on success and -1 on error or when end of
 *     file occurs while no characters have been read.
 */
int fgets_rtrimmed(char *s, int size, FILE *stream)
{
	long pos0;
	char *s1;
	int line_len;

	pos0 = ftell(stream);
	s1 = fgets(s, size, stream);
	if (s1 == NULL)
		return -1;
	/* 2 almost equivalent ways to get the length of the current line,
	   "almost" because they will differ if a line contains embedded
	   NUL characters */
	line_len = ftell(stream) - pos0; /* should be faster than strlen() */
	/* line_len = strlen(s); */
	return rtrimline(s, line_len);
}

