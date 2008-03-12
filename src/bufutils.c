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

#define MAX_BUFLENGTH_INC (128 * 1024 * 1024)
#define MAX_BUFLENGTH (8 * MAX_BUFLENGTH_INC)


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

static int get_new_buflength(int buflength)
{
	if (buflength >= MAX_BUFLENGTH)
		error("get_new_buflength(): MAX_BUFLENGTH reached");
	if (buflength == 0)
		return 256;
	if (buflength <= 256 * 1024)
		return 4 * buflength;
	if (buflength <= MAX_BUFLENGTH_INC)
		return 2 * buflength;
	buflength += MAX_BUFLENGTH_INC;
	if (buflength <= MAX_BUFLENGTH)
		return buflength;
	return MAX_BUFLENGTH;
}


/****************************************************************************
 * IBuf functions
 */

IBuf _new_IBuf(int buflength, int nelt)
{
	IBuf ibuf;
	int *elt;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		ibuf.elts = NULL;
	else
		ibuf.elts = Salloc((long) buflength, int);
	ibuf.buflength = buflength;
	for (ibuf.nelt = 0, elt = ibuf.elts;
	     ibuf.nelt < nelt;
	     ibuf.nelt++, elt++)
		*elt = 0;
	return ibuf;
}

static void _IBuf_extend(IBuf *ibuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(ibuf->buflength);
	ibuf->elts = Srealloc((char *) ibuf->elts, new_buflength,
					(long) ibuf->buflength, int);
	ibuf->buflength = new_buflength;
	return;
}

void _IBuf_insert_at(IBuf *ibuf, int at, int val)
{
	int *elt1, *elt2;
	int i1;

	if (ibuf->nelt >= ibuf->buflength)
		_IBuf_extend(ibuf);
	elt1 = ibuf->elts + ibuf->nelt;
	elt2 = elt1 + 1;
	for (i1 = ibuf->nelt++; i1 >= at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = val;
	return;
}

void _IBuf_delete_at(IBuf *ibuf, int at)
{
	int *elt1, *elt2;
	int i2;

	elt1 = ibuf->elts + at;
	elt2 = elt1 + 1;
	for (i2 = at + 1; i2 < ibuf->nelt; i2++)
		*(elt1++) = *(elt2++);
	ibuf->nelt--;
	return;
}

SEXP _IBuf_asINTEGER(IBuf *ibuf)
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ibuf->nelt));
	memcpy(INTEGER(ans), ibuf->elts, sizeof(int) * ibuf->nelt);
	UNPROTECT(1);
	return ans;
}

IBuf _INTEGER_asIBuf(SEXP x)
{
	IBuf ibuf;

	ibuf = _new_IBuf(LENGTH(x), 0);
	memcpy(ibuf.elts, INTEGER(x), sizeof(int) * LENGTH(x));
	ibuf.nelt = ibuf.buflength;
	return ibuf;
}

IBuf _CHARACTER_asIBuf(SEXP x, int keyshift)
{
	IBuf ibuf;
	int *elt;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CHARACTER_asIBuf(): BEGIN ... "
			"LENGTH(x)=%d keyshift=%d\n",
			LENGTH(x), keyshift);
	}
#endif
	ibuf = _new_IBuf(LENGTH(x), 0);
	for (ibuf.nelt = 0, elt = ibuf.elts;
	     ibuf.nelt < ibuf.buflength;
	     ibuf.nelt++, elt++) {
		sscanf(CHAR(STRING_ELT(x, ibuf.nelt)), "%d", elt);
		*elt += keyshift;
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (ibuf.nelt < 100
			 || ibuf.nelt >= ibuf.buflength - 100)
				Rprintf("[DEBUG] _CHARACTER_asIBuf(): "
					"ibuf.nelt=%d key=%s *elt=%d\n",
					ibuf.nelt,
					CHAR(STRING_ELT(x, ibuf.nelt)), *elt);
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

IBBuf _new_IBBuf(int buflength, int nelt)
{
	IBBuf ibbuf;
	IBuf *elt;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		ibbuf.elts = NULL;
	else
		ibbuf.elts = Salloc((long) buflength, IBuf);
	ibbuf.buflength = buflength;
	for (ibbuf.nelt = 0, elt = ibbuf.elts;
	     ibbuf.nelt < nelt;
	     ibbuf.nelt++, elt++)
		*elt = _new_IBuf(0, 0);
	return ibbuf;
}

static void _IBBuf_extend(IBBuf *ibbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(ibbuf->buflength);
	ibbuf->elts = Srealloc((char *) ibbuf->elts, new_buflength,
					(long) ibbuf->buflength, IBuf);
	ibbuf->buflength = new_buflength;
	return;
}

void _IBBuf_insert_at(IBBuf *ibbuf, int at, IBuf ibuf)
{
	IBuf *elt1, *elt2;
	int i1;

	if (ibbuf->nelt >= ibbuf->buflength)
		_IBBuf_extend(ibbuf);
	elt1 = ibbuf->elts + ibbuf->nelt;
	elt2 = elt1 + 1;
	for (i1 = ibbuf->nelt++; i1 >= at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = ibuf;
	return;
}

/*
 * mode: 0 -> integer(0), 1-> NULL, 2 -> NA
 */
SEXP _IBBuf_asLIST(IBBuf *ibbuf, int mode)
{
	SEXP ans, ans_elt;
	int i;
	IBuf *elt;

	PROTECT(ans = NEW_LIST(ibbuf->nelt));
	for (i = 0, elt = ibbuf->elts; i < ibbuf->nelt; i++, elt++) {
		if (elt->nelt == 0 && mode != 0) {
			if (mode == 1) {
				PROTECT(ans_elt = R_NilValue);
			} else {
				// Not sure new LOGICALs are initialized with NAs,
				// need to check! If not, then LOGICAL(ans_elt)[0]
				// must be set to NA but I don't know how to do this :-/
				PROTECT(ans_elt = NEW_LOGICAL(1));
			}
		} else {
			PROTECT(ans_elt = _IBuf_asINTEGER(elt));
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
	IBuf *elt;

	ibbuf = _new_IBBuf(LENGTH(x), 0);
	for (ibbuf.nelt = 0, elt = ibbuf.elts;
	     ibbuf.nelt < ibbuf.buflength;
	     ibbuf.nelt++, elt++) {
		*elt = _INTEGER_asIBuf(VECTOR_ELT(x, ibbuf.nelt));
	}
	return ibbuf;
}

SEXP _IBBuf_toEnvir(IBBuf *ibbuf, SEXP envir, int keyshift)
{
	int i;
	IBuf *elt;
	char key[11];
	SEXP value;

#ifdef DEBUG_BIOSTRINGS
	int nkey = 0, cum_length = 0;
	if (debug) {
		Rprintf("[DEBUG] _IBBuf_toEnvir(): BEGIN ... "
			"ibbuf->nelt=%d keyshift=%d\n",
			ibbuf->nelt, keyshift);
	}
#endif
	for (i = 0, elt = ibbuf->elts; i < ibbuf->nelt; i++, elt++) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (i < 100 || i >= ibbuf->nelt - 100)
				Rprintf("[DEBUG] _IBBuf_toEnvir(): "
					"nkey=%d ibbuf->elts[%d].nelt=%d\n",
					nkey, i, elt->nelt);
		}
#endif
		if (elt->nelt == 0)
			continue;
		//snprintf(key, sizeof(key), "%d", i + keyshift);
		snprintf(key, sizeof(key), "%010d", i + keyshift);
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (i < 100 || i >= ibbuf->nelt - 100)
				Rprintf("[DEBUG] _IBBuf_toEnvir(): "
					"installing key=%s ... ", key);
		}
#endif
		PROTECT(value = _IBuf_asINTEGER(elt));
		defineVar(install(key), value, envir);
		UNPROTECT(1);
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			nkey++;
			cum_length += elt->nelt;
			if (i < 100 || i >= ibbuf->nelt - 100)
				Rprintf("OK (nkey=%d cum_length=%d)\n",
					nkey, cum_length);
		}
#endif
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _IBBuf_toEnvir(): END "
			"(nkey=%d cum_length=%d)\n", nkey, cum_length);
	}
#endif
	return envir;
}


/****************************************************************************
 * RangesBuf functions
 */

RangesBuf _new_RangesBuf(int buflength, int nelt)
{
	RangesBuf rangesbuf;

	rangesbuf.start = _new_IBuf(buflength, nelt);
	rangesbuf.width = _new_IBuf(buflength, nelt);
	return rangesbuf;
}

void _RangesBuf_insert_at(RangesBuf *rangesbuf, int at, int start, int width)
{
	_IBuf_insert_at(&(rangesbuf->start), at, start);
	_IBuf_insert_at(&(rangesbuf->width), at, width);
	return;
}


/****************************************************************************
 * CBuf functions
 */

CBuf _new_CBuf(int buflength)
{
	CBuf cbuf;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		cbuf.elts = NULL;
	else
		cbuf.elts = Salloc((long) buflength, char);
	cbuf.buflength = buflength;
	cbuf.nelt = 0;
	return cbuf;
}

CBuf _new_CBuf_from_string(const char *string)
{
	CBuf cbuf;
	int buflength;

	buflength = strlen(string) + 1;
	cbuf = _new_CBuf(buflength);
	memcpy(cbuf.elts, string, buflength);
        return cbuf;
}

static void _CBuf_extend(CBuf *cbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(cbuf->buflength);
	cbuf->elts = Srealloc((char *) cbuf->elts, new_buflength,
					(long) cbuf->buflength, char);
	cbuf->buflength = new_buflength;
	return;
}

void _CBuf_insert_at(CBuf *cbuf, int at, char c)
{
	char *elt1, *elt2;
	int i1;

	if (cbuf->nelt >= cbuf->buflength)
		_CBuf_extend(cbuf);
	elt1 = cbuf->elts + cbuf->nelt;
	elt2 = elt1 + 1;
	for (i1 = cbuf->nelt++; i1 >= at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = c;
	return;
}

SEXP _CBuf_asRAW(CBuf *cbuf)
{
	SEXP ans;

	if (sizeof(Rbyte) != sizeof(char)) // should never happen!
		error("_CBuf_asRAW(): sizeof(Rbyte) != sizeof(char)");
	PROTECT(ans = NEW_RAW(cbuf->nelt));
	memcpy(RAW(ans), cbuf->elts, sizeof(char) * cbuf->nelt);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * CBBuf functions
 */

CBBuf _new_CBBuf(int buflength, int nelt)
{
	CBBuf cbbuf;
	CBuf *elt;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		cbbuf.elts = NULL;
	else
		cbbuf.elts = Salloc((long) buflength, CBuf);
	cbbuf.buflength = buflength;
	for (cbbuf.nelt = 0, elt = cbbuf.elts;
	     cbbuf.nelt < nelt;
	     cbbuf.nelt++, elt++)
		*elt = _new_CBuf(0);
	return cbbuf;
}

static void _CBBuf_extend(CBBuf *cbbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(cbbuf->buflength);
	cbbuf->elts = Srealloc((char *) cbbuf->elts, new_buflength,
				(long) cbbuf->buflength, CBuf);
	cbbuf->buflength = new_buflength;
	return;
}

void _CBBuf_insert_at(CBBuf *cbbuf, int at, CBuf cbuf)
{
	CBuf *elt1, *elt2;
	int i1;

	if (cbbuf->nelt >= cbbuf->buflength)
		_CBBuf_extend(cbbuf);
	elt1 = cbbuf->elts + cbbuf->nelt;
	elt2 = elt1 + 1;
	for (i1 = cbbuf->nelt++; i1 >= at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = cbuf;
	return;
}

void _append_string_to_CBBuf(CBBuf *cbbuf, const char *string)
{
	CBuf cbuf;

	cbuf = _new_CBuf_from_string(string);
	_CBBuf_insert_at(cbbuf, cbbuf->nelt, cbuf);
	return;
}

