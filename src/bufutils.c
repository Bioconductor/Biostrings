/*
 * Low-level manipulation of the extendable buffers.
 *
 * Except for debug_bufutils(), the functions defined in
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

SEXP debug_bufutils()
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
 * IntBuf functions
 */

void _IntBuf_set_val(IntBuf *ibuf, int val)
{
	int i, *elt;

	for (i = 0, elt = ibuf->elts; i < ibuf->nelt; i++, elt++)
		*elt = val;
	return;
}

IntBuf _new_IntBuf(int buflength, int nelt)
{
	IntBuf ibuf;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		ibuf.elts = NULL;
	else
		ibuf.elts = Salloc((long) buflength, int);
	ibuf.buflength = buflength;
	ibuf.nelt = nelt;
	_IntBuf_set_val(&ibuf, 0);
	return ibuf;
}

static void IntBuf_extend(IntBuf *ibuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(ibuf->buflength);
	ibuf->elts = Srealloc((char *) ibuf->elts, new_buflength,
					(long) ibuf->buflength, int);
	ibuf->buflength = new_buflength;
	return;
}

void _IntBuf_insert_at(IntBuf *ibuf, int at, int val)
{
	int *elt1, *elt2;
	int i1;

	if (ibuf->nelt >= ibuf->buflength)
		IntBuf_extend(ibuf);
	elt2 = ibuf->elts + ibuf->nelt;
	elt1 = elt2 - 1;
	for (i1 = ibuf->nelt++; i1 > at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = val;
	return;
}

void _IntBuf_append(IntBuf *ibuf, int *vals, int nval)
{
	int new_nelt, *dest;

	new_nelt = ibuf->nelt + nval;
	while (new_nelt > ibuf->buflength)
		IntBuf_extend(ibuf);
	dest = ibuf->elts + ibuf->nelt;
	memcpy(dest, vals, nval * sizeof(int));
	ibuf->nelt = new_nelt;
	return;
}

void _IntBuf_delete_at(IntBuf *ibuf, int at)
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

void _IntBuf_sum_val(IntBuf *ibuf, int val)
{
	int i, *elt;

	for (i = 0, elt = ibuf->elts; i < ibuf->nelt; i++, elt++)
		*elt += val;
	return;
}

/*
 * Left and right IntBuf objects must have the same length. This is
 * NOT checked!
 */
void _IntBuf_sum_IntBuf(IntBuf *ibuf1, IntBuf *ibuf2)
{
	int i, *elt1, *elt2;

	for (i = 0, elt1 = ibuf1->elts, elt2 = ibuf2->elts;
	     i < ibuf1->nelt;
	     i++, elt1++, elt2++)
		*elt1 += *elt2;
	return;
}

SEXP _IntBuf_asINTEGER(IntBuf *ibuf)
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(ibuf->nelt));
	memcpy(INTEGER(ans), ibuf->elts, sizeof(int) * ibuf->nelt);
	UNPROTECT(1);
	return ans;
}

IntBuf _INTEGER_asIntBuf(SEXP x)
{
	IntBuf ibuf;

	ibuf = _new_IntBuf(LENGTH(x), 0);
	memcpy(ibuf.elts, INTEGER(x), sizeof(int) * LENGTH(x));
	ibuf.nelt = ibuf.buflength;
	return ibuf;
}

IntBuf _CHARACTER_asIntBuf(SEXP x, int keyshift)
{
	IntBuf ibuf;
	int *elt;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CHARACTER_asIntBuf(): BEGIN ... "
			"LENGTH(x)=%d keyshift=%d\n",
			LENGTH(x), keyshift);
	}
#endif
	ibuf = _new_IntBuf(LENGTH(x), 0);
	for (ibuf.nelt = 0, elt = ibuf.elts;
	     ibuf.nelt < ibuf.buflength;
	     ibuf.nelt++, elt++) {
		sscanf(CHAR(STRING_ELT(x, ibuf.nelt)), "%d", elt);
		*elt += keyshift;
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (ibuf.nelt < 100
			 || ibuf.nelt >= ibuf.buflength - 100)
				Rprintf("[DEBUG] _CHARACTER_asIntBuf(): "
					"ibuf.nelt=%d key=%s *elt=%d\n",
					ibuf.nelt,
					CHAR(STRING_ELT(x, ibuf.nelt)), *elt);
		}
#endif
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CHARACTER_asIntBuf(): END\n");
	}
#endif
	return ibuf;
}


/****************************************************************************
 * IntBBuf functions
 */

IntBBuf _new_IntBBuf(int buflength, int nelt)
{
	IntBBuf ibbuf;
	IntBuf *elt;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		ibbuf.elts = NULL;
	else
		ibbuf.elts = Salloc((long) buflength, IntBuf);
	ibbuf.buflength = buflength;
	for (ibbuf.nelt = 0, elt = ibbuf.elts;
	     ibbuf.nelt < nelt;
	     ibbuf.nelt++, elt++)
		*elt = _new_IntBuf(0, 0);
	return ibbuf;
}

static void IntBBuf_extend(IntBBuf *ibbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(ibbuf->buflength);
	ibbuf->elts = Srealloc((char *) ibbuf->elts, new_buflength,
					(long) ibbuf->buflength, IntBuf);
	ibbuf->buflength = new_buflength;
	return;
}

void _IntBBuf_insert_at(IntBBuf *ibbuf, int at, IntBuf ibuf)
{
	IntBuf *elt1, *elt2;
	int i1;

	if (ibbuf->nelt >= ibbuf->buflength)
		IntBBuf_extend(ibbuf);
	elt2 = ibbuf->elts + ibbuf->nelt;
	elt1 = elt2 - 1;
	for (i1 = ibbuf->nelt++; i1 > at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = ibuf;
	return;
}

/*
 * Left and right IntBBuf objects must have the same length. This is
 * NOT checked!
 */
void _IntBBuf_eltwise_append(IntBBuf *ibbuf1, IntBBuf *ibbuf2)
{
	int i;
	IntBuf *elt1, *elt2;

	for (i = 0, elt1 = ibbuf1->elts, elt2 = ibbuf2->elts;
	     i < ibbuf1->nelt;
	     i++, elt1++, elt2++)
		_IntBuf_append(elt1, elt2->elts, elt2->nelt);
	return;
}

void _IntBBuf_sum_val(IntBBuf *ibbuf, int val)
{
	int i;
	IntBuf *elt;

	for (i = 0, elt = ibbuf->elts; i < ibbuf->nelt; i++, elt++)
		_IntBuf_sum_val(elt, val);
	return;
}

/*
 * mode: 0 -> integer(0), 1-> NULL, 2 -> NA
 */
SEXP _IntBBuf_asLIST(IntBBuf *ibbuf, int mode)
{
	SEXP ans, ans_elt;
	int i;
	IntBuf *elt;

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
			PROTECT(ans_elt = _IntBuf_asINTEGER(elt));
		}
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

IntBBuf _LIST_asIntBBuf(SEXP x)
{
	IntBBuf ibbuf;
	IntBuf *elt;

	ibbuf = _new_IntBBuf(LENGTH(x), 0);
	for (ibbuf.nelt = 0, elt = ibbuf.elts;
	     ibbuf.nelt < ibbuf.buflength;
	     ibbuf.nelt++, elt++) {
		*elt = _INTEGER_asIntBuf(VECTOR_ELT(x, ibbuf.nelt));
	}
	return ibbuf;
}

SEXP _IntBBuf_toEnvir(IntBBuf *ibbuf, SEXP envir, int keyshift)
{
	int i;
	IntBuf *elt;
	char key[11];
	SEXP value;

#ifdef DEBUG_BIOSTRINGS
	int nkey = 0, cum_length = 0;
	if (debug) {
		Rprintf("[DEBUG] _IntBBuf_toEnvir(): BEGIN ... "
			"ibbuf->nelt=%d keyshift=%d\n",
			ibbuf->nelt, keyshift);
	}
#endif
	for (i = 0, elt = ibbuf->elts; i < ibbuf->nelt; i++, elt++) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			if (i < 100 || i >= ibbuf->nelt - 100)
				Rprintf("[DEBUG] _IntBBuf_toEnvir(): "
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
				Rprintf("[DEBUG] _IntBBuf_toEnvir(): "
					"installing key=%s ... ", key);
		}
#endif
		PROTECT(value = _IntBuf_asINTEGER(elt));
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
		Rprintf("[DEBUG] _IntBBuf_toEnvir(): END "
			"(nkey=%d cum_length=%d)\n", nkey, cum_length);
	}
#endif
	return envir;
}


/****************************************************************************
 * RangeBuf functions
 */

RangeBuf _new_RangeBuf(int buflength, int nelt)
{
	RangeBuf rangebuf;

	rangebuf.start = _new_IntBuf(buflength, nelt);
	rangebuf.width = _new_IntBuf(buflength, nelt);
	return rangebuf;
}

void _RangeBuf_insert_at(RangeBuf *rangebuf, int at, int start, int width)
{
	_IntBuf_insert_at(&(rangebuf->start), at, start);
	_IntBuf_insert_at(&(rangebuf->width), at, width);
	return;
}


/****************************************************************************
 * CharBuf functions
 */

CharBuf _new_CharBuf(int buflength)
{
	CharBuf cbuf;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		cbuf.elts = NULL;
	else
		cbuf.elts = Salloc((long) buflength, char);
	cbuf.buflength = buflength;
	cbuf.nelt = 0;
	return cbuf;
}

CharBuf _new_CharBuf_from_string(const char *string)
{
	CharBuf cbuf;
	int buflength;

	buflength = strlen(string);
	cbuf = _new_CharBuf(buflength);
	memcpy(cbuf.elts, string, buflength);
	cbuf.nelt = buflength;
	return cbuf;
}

static void CharBuf_extend(CharBuf *cbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(cbuf->buflength);
	cbuf->elts = Srealloc((char *) cbuf->elts, new_buflength,
					(long) cbuf->buflength, char);
	cbuf->buflength = new_buflength;
	return;
}

void _CharBuf_insert_at(CharBuf *cbuf, int at, char c)
{
	char *elt1, *elt2;
	int i1;

	if (cbuf->nelt >= cbuf->buflength)
		CharBuf_extend(cbuf);
	elt2 = cbuf->elts + cbuf->nelt;
	elt1 = elt2 - 1;
	for (i1 = cbuf->nelt++; i1 > at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = c;
	return;
}

SEXP _CharBuf_asRAW(CharBuf *cbuf)
{
	SEXP ans;

	if (sizeof(Rbyte) != sizeof(char)) // should never happen!
		error("_CharBuf_asRAW(): sizeof(Rbyte) != sizeof(char)");
	PROTECT(ans = NEW_RAW(cbuf->nelt));
	memcpy(RAW(ans), cbuf->elts, sizeof(char) * cbuf->nelt);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * CharBBuf functions
 */

CharBBuf _new_CharBBuf(int buflength, int nelt)
{
	CharBBuf cbbuf;
	CharBuf *elt;

	/* No memory leak here, because we use transient storage allocation */
	if (buflength == 0)
		cbbuf.elts = NULL;
	else
		cbbuf.elts = Salloc((long) buflength, CharBuf);
	cbbuf.buflength = buflength;
	for (cbbuf.nelt = 0, elt = cbbuf.elts;
	     cbbuf.nelt < nelt;
	     cbbuf.nelt++, elt++)
		*elt = _new_CharBuf(0);
	return cbbuf;
}

static void CharBBuf_extend(CharBBuf *cbbuf)
{
	long new_buflength;

	new_buflength = get_new_buflength(cbbuf->buflength);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] CharBBuf_extend(): BEGIN\n");
		Rprintf("[DEBUG] CharBBuf_extend(): "
			"cbbuf->elts=%p buflength=%d new_buflength=%d\n",
			cbbuf->elts, cbbuf->buflength, new_buflength);
	}
#endif
	cbbuf->elts = Srealloc((char *) cbbuf->elts, new_buflength,
				(long) cbbuf->buflength, CharBuf);
	cbbuf->buflength = new_buflength;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] CharBBuf_extend(): END (cbbuf->elts=%p)\n",
			cbbuf->elts);
	}
#endif
	return;
}

void _CharBBuf_insert_at(CharBBuf *cbbuf, int at, CharBuf cbuf)
{
	CharBuf *elt1, *elt2;
	int i1;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CharBBuf_insert_at(): BEGIN\n");
	}
#endif
	if (cbbuf->nelt >= cbbuf->buflength)
		CharBBuf_extend(cbbuf);
	elt2 = cbbuf->elts + cbbuf->nelt;
	elt1 = elt2 - 1;
	for (i1 = cbbuf->nelt++; i1 > at; i1--)
		*(elt2--) = *(elt1--);
	*elt2 = cbuf;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _CharBBuf_insert_at(): END\n");
	}
#endif
	return;
}

void _append_string_to_CharBBuf(CharBBuf *cbbuf, const char *string)
{
	CharBuf cbuf;

	cbuf = _new_CharBuf_from_string(string);
	_CharBBuf_insert_at(cbbuf, cbbuf->nelt, cbuf);
	return;
}

