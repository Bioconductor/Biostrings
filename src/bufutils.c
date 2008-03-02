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

static void _IBuf_extend(IBuf *ibuf)
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
	int *val1, *val2, i1;

	if (ibuf->count >= ibuf->maxcount)
		_IBuf_extend(ibuf);
	val1 = ibuf->vals + ibuf->count;
	val2 = val1 + 1;
	for (i1 = ibuf->count++; i1 >= at; i1--)
		*(val2--) = *(val1--);
	*val2 = val;
	return;
}

void _IBuf_delete_at(IBuf *ibuf, int at)
{
	int *val1, *val2, i2;

	val1 = ibuf->vals + at;
	val2 = val1 + 1;
	for (i2 = at + 1; i2 < ibuf->count; i2++)
		*(val1++) = *(val2++);
	ibuf->count--;
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

	_IBuf_init(&ibuf, LENGTH(x), 0);
	memcpy(ibuf.vals, INTEGER(x), sizeof(int) * LENGTH(x));
	ibuf.count = ibuf.maxcount;
	return ibuf;
}

IBuf _CHARACTER_asIBuf(SEXP x, int keyshift)
{
	IBuf ibuf;
	int *val;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CHARACTER_asIBuf(): BEGIN ... LENGTH(x)=%d keyshift=%d\n",
			LENGTH(x), keyshift);
	}
#endif
	_IBuf_init(&ibuf, LENGTH(x), 0);
	for (ibuf.count = 0, val = ibuf.vals;
	     ibuf.count < ibuf.maxcount;
	     ibuf.count++, val++) {
		sscanf(CHAR(STRING_ELT(x, ibuf.count)), "%d", val);
		*val += keyshift;
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (ibuf.count < 100 || ibuf.count >= ibuf.maxcount - 100)
				Rprintf("[DEBUG] _CHARACTER_asIBuf(): ibuf.count=%d key=%s *val=%d\n",
					ibuf.count, CHAR(STRING_ELT(x, ibuf.count)), *val);
		}
#endif
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CHARACTER_asIBuf(): END\n");
	}
#endif
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

static void _IBBuf_extend(IBBuf *ibbuf)
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
		_IBBuf_extend(ibbuf);
	j = ibbuf->count++;
	for (i = j - 1; i >= at; i--, j--)
		ibbuf->ibufs[j] = ibbuf->ibufs[i];
	ibbuf->ibufs[j] = ibuf;
	return;
}

/*
 * mode: 0 -> integer(0), 1-> NULL, 2 -> NA
 */
SEXP _IBBuf_asLIST(IBBuf *ibbuf, int mode)
{
	SEXP ans, ans_elt;
	int i;
	IBuf *ibuf;

	PROTECT(ans = NEW_LIST(ibbuf->count));
	for (i = 0, ibuf = ibbuf->ibufs; i < ibbuf->count; i++, ibuf++) {
		if (ibuf->count == 0 && mode != 0) {
			if (mode == 1) {
				PROTECT(ans_elt = R_NilValue);
			} else {
				// Not sure new LOGICALs are initialized with NAs,
				// need to check! If not, then LOGICAL(ans_elt)[0]
				// must be set to NA but I don't know how to do this :-/
				PROTECT(ans_elt = NEW_LOGICAL(1));
			}
		} else {
			PROTECT(ans_elt = _IBuf_asINTEGER(ibuf));
		}
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

IBBuf _LIST_asIBBuf(SEXP x)
{
	IBBuf ibbuf;
	IBuf *ibuf_p;

	_IBBuf_init(&ibbuf, LENGTH(x), 0);
	for (ibbuf.count = 0, ibuf_p = ibbuf.ibufs;
	     ibbuf.count < ibbuf.maxcount;
	     ibbuf.count++, ibuf_p++) {
		*ibuf_p = _INTEGER_asIBuf(VECTOR_ELT(x, ibbuf.count));
	}
	return ibbuf;
}

SEXP _IBBuf_toEnvir(IBBuf *ibbuf, SEXP envir, int keyshift)
{
	int i;
	IBuf *ibuf;
	char key[11];
	SEXP value;

#ifdef DEBUG_BIOSTRINGS
	int nkey = 0, cum_length = 0;
	if (debug) {
		Rprintf("[DEBUG] _IBBuf_toEnvir(): BEGIN ... ibbuf->count=%d keyshift=%d\n",
			ibbuf->count, keyshift);
	}
#endif
	for (i = 0, ibuf = ibbuf->ibufs; i < ibbuf->count; i++, ibuf++) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (i < 100 || i >= ibbuf->count - 100)
				Rprintf("[DEBUG] _IBBuf_toEnvir(): nkey=%d ibbuf->ibufs[%d].count=%d\n",
					nkey, i, ibuf->count);
		}
#endif
		if (ibuf->count == 0)
			continue;
		//snprintf(key, sizeof(key), "%d", i + keyshift);
		snprintf(key, sizeof(key), "%010d", i + keyshift);
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (i < 100 || i >= ibbuf->count - 100)
				Rprintf("[DEBUG] _IBBuf_toEnvir(): installing key=%s ... ", key);
		}
#endif
		PROTECT(value = _IBuf_asINTEGER(ibuf));
		defineVar(install(key), value, envir);
		UNPROTECT(1);
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			nkey++;
			cum_length += ibuf->count;
			if (i < 100 || i >= ibbuf->count - 100)
				Rprintf("OK (nkey=%d cum_length=%d)\n", nkey, cum_length);
		}
#endif
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _IBBuf_toEnvir(): END (nkey=%d cum_length=%d)\n", nkey, cum_length);
	}
#endif
	return envir;
}


/****************************************************************************
 * InterBuf functions
 */

void _InterBuf_init(InterBuf *interbuf, int maxcount, int count)
{
	_IBuf_init(&(interbuf->start), maxcount, count);
	_IBuf_init(&(interbuf->width), maxcount, count);
	return;
}

void _InterBuf_insert_at(InterBuf *interbuf, int at, int start, int width)
{
	_IBuf_insert_at(&(interbuf->start), at, start);
	_IBuf_insert_at(&(interbuf->width), at, width);
	return;
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

static void _CBuf_extend(CBuf *cbuf)
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
		_CBuf_extend(cbuf);
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

