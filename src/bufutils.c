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
#include <S.h> /* for Salloc() and Srealloc() */

#define MAX_BUF_LENGTHINC (128 * 1024 * 1024)
#define MAX_BUF_LENGTH (8 * MAX_BUF_LENGTHINC)


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
	if (maxcount == 0)
		return 256;
	if (maxcount <= 256 * 1024)
		return 4 * maxcount;
	if (maxcount <= MAX_BUF_LENGTHINC)
		return 2 * maxcount;
	if (maxcount == MAX_BUF_LENGTH)
		error("get_new_maxcount(): MAX_BUF_LENGTH reached");
	maxcount += MAX_BUF_LENGTHINC;
	if (maxcount <= MAX_BUF_LENGTH)
		return maxcount;
	return MAX_BUF_LENGTH;
}


/****************************************************************************
 * IBuf functions
 */

void _IBuf_init(IBuf *ibuf, int maxcount, int count)
{
	int i, *val;

	/* No memory leak here, because we use transient storage allocation */
	if (maxcount == 0)
		ibuf->vals = NULL;
	else
		ibuf->vals = Salloc((long) maxcount, int);
	ibuf->maxcount = maxcount;
	for (i = 0, val = ibuf->vals; i < count; i++, val++)
		*val = 0;
	ibuf->count = count;
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

IBuf _INTEGER_asIBuf(SEXP x)
{
	IBuf ibuf;
	int x_length;

	x_length = LENGTH(x);
	ibuf.vals = Salloc((long) x_length, int);
	memcpy(ibuf.vals, INTEGER(x), sizeof(int) * x_length);
	ibuf.maxcount = ibuf.count = x_length;
	return ibuf;
}


/****************************************************************************
 * IBBuf functions
 */

void _IBBuf_init(IBBuf *ibbuf, int maxcount, int count)
{
	int i;
	IBuf *ibuf;

	/* No memory leak here, because we use transient storage allocation */
	if (maxcount == 0)
		ibbuf->ibufs = NULL;
	else
		ibbuf->ibufs = Salloc((long) maxcount, IBuf);
	ibbuf->maxcount = maxcount;
	for (i = 0, ibuf = ibbuf->ibufs; i < count; i++, ibuf++)
		_IBuf_init(ibuf, 0, 0);
	ibbuf->count = count;
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

IBBuf _LIST_asIBBuf(SEXP x)
{
	IBBuf ibbuf;
	int x_length, i;
	IBuf *ibuf_p;

	x_length = LENGTH(x);
	ibbuf.ibufs = Salloc((long) x_length, IBuf);
	for (i = 0, ibuf_p = ibbuf.ibufs; i < x_length; i++, ibuf_p++)
		*ibuf_p = _INTEGER_asIBuf(VECTOR_ELT(x, i));
	ibbuf.maxcount = ibbuf.count = x_length;
	return ibbuf;
}


/****************************************************************************
 * CBuf functions
 */

void _CBuf_init(CBuf *cbuf, int maxcount)
{
	/* No memory leak here, because we use transient storage allocation */
	if (maxcount == 0)
		cbuf->vals = NULL;
	else
		cbuf->vals = Salloc((long) maxcount, char);
	cbuf->maxcount = maxcount;
	cbuf->count = 0;
	return;
}

void _CBuf_get_more_room(CBuf *cbuf)
{
	long new_maxcount;

	new_maxcount = get_new_maxcount(cbuf->maxcount);
	cbuf->vals = Srealloc((char *) cbuf->vals, new_maxcount,
					(long) cbuf->maxcount, char);
	cbuf->maxcount = new_maxcount;
	return;
}

void _CBuf_insert_at(CBuf *cbuf, int at, char val)
{
	int i, j;

	if (cbuf->count >= cbuf->maxcount)
		_CBuf_get_more_room(cbuf);
	j = cbuf->count++;
	for (i = j - 1; i >= at; i--, j--)
		cbuf->vals[j] = cbuf->vals[i];
	cbuf->vals[j] = val;
	return;
}

SEXP _CBuf_asRAW(CBuf *cbuf)
{
	SEXP ans;

	if (sizeof(Rbyte) != sizeof(char)) // should never happen!
		error("_CBuf_asRAW(): sizeof(Rbyte) != sizeof(char)");
	PROTECT(ans = NEW_RAW(cbuf->count));
	memcpy(RAW(ans), cbuf->vals, sizeof(char) * cbuf->count);
	UNPROTECT(1);
	return ans;
}

