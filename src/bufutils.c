/*
 * Functions for low-level manipulation of the "temporary buffers".
 *
 * Except for Biostrings_debug_bufutils(), the functions defined in
 * this file are NOT .Call methods (but they are used by .Call methods
 * defined in other .c files) so THEY DON'T NEED TO BE REGISTERED in
 * R_init_Biostrings.c. They are prefixed with a "_" (underscore) to
 * emphasize the fact that they are used internally within the Biostrings
 * shared lib.
 */
#include "Biostrings.h"
#include <S.h> /* for Srealloc() */

#define MAX_BUFSIZEINC (128 * 1024 * 1024)
#define MAX_BUFSIZE (8 * MAX_BUFSIZEINC)


static int debug = 0;

SEXP Biostrings_debug_bufutils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'bufutils.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'bufutils.c'\n");
#endif
	return R_NilValue;
}

static int get_new_maxcount(int maxcount)
{
	if (maxcount >= MAX_BUFSIZE)
		error("get_new_maxcount(): MAX_BUFSIZE reached");
	if (maxcount == 0)
		return 1024;
	if (maxcount <= MAX_BUFSIZEINC)
		return 2 * maxcount;
	return maxcount + MAX_BUFSIZEINC;
}


/****************************************************************************
 * IBuf functions
 */

void _IBuf_init(IBuf *ibuf)
{
	/* No memory leak here, because we use transient storage allocation */
	ibuf->vals = NULL;
	ibuf->maxcount = ibuf->count = 0;
	return;
}

void _IBuf_get_more_room(IBuf *ibuf)
{
	long new_maxcount;

	new_maxcount = get_new_maxcount(ibuf->maxcount);
	ibuf->vals = Srealloc((char *) ibuf->vals, new_maxcount,
					(long) ibuf->maxcount, int);
	ibuf->maxcount = new_maxcount;
	return;
}

void _IBuf_insert_at(IBuf *ibuf, int at, int val)
{
	int i, j;

	if (ibuf->count >= ibuf->maxcount)
		_IBuf_get_more_room(ibuf);
	j = ibuf->count++;
	for (i = j - 1; i >= at; i--, j--)
		ibuf->vals[j] = ibuf->vals[i];
	ibuf->vals[j] = val;
	return;
}

SEXP _IBuf_asINTEGER(IBuf *ibuf)
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ibuf->count));
	memcpy(INTEGER(ans), ibuf->vals, sizeof(int) * ibuf->count);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * IBBuf functions
 */

void _IBBuf_init(IBBuf *ibbuf)
{
	/* No memory leak here, because we use transient storage allocation */
	ibbuf->ibufs = NULL;
	ibbuf->maxcount = ibbuf->count = 0;
	return;
}

void _IBBuf_get_more_room(IBBuf *ibbuf)
{
	long new_maxcount;

	new_maxcount = get_new_maxcount(ibbuf->maxcount);
	ibbuf->ibufs = Srealloc((char *) ibbuf->ibufs, new_maxcount,
					(long) ibbuf->maxcount, IBuf);
	ibbuf->maxcount = new_maxcount;
	return;
}

void _IBBuf_insert_at(IBBuf *ibbuf, int at, IBuf ibuf)
{
	int i, j;

	if (ibbuf->count >= ibbuf->maxcount)
		_IBBuf_get_more_room(ibbuf);
	j = ibbuf->count++;
	for (i = j - 1; i >= at; i--, j--)
		ibbuf->ibufs[j] = ibbuf->ibufs[i];
	ibbuf->ibufs[j] = ibuf;
	return;
}

SEXP _IBBuf_asLIST(IBBuf *ibbuf)
{
	SEXP ans, ans_elt;
	int i;
	IBuf *ibuf;

	PROTECT(ans = NEW_LIST(ibbuf->count));
	for (i = 0, ibuf = ibbuf->ibufs; i < ibbuf->count; i++, ibuf++) {
		PROTECT(ans_elt = _IBuf_asINTEGER(ibuf));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

