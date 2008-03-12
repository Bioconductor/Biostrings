#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/*
 * Extendable buffers used for temporary storage of incoming data whose size
 * is not known in advance.
 * They are NOT an attempt to reinvent an SEXP subsystem. Some notable
 * differences are: (a) they are extendable (i.e. they are automatically
 * reallocated when more room is needed to add a new element), (b) they are
 * much faster, and (c) they don't require any PROTECT/UNPROTECT mechanism.
 */

typedef struct ibuf {
	int *elts;
	int buflength;
	int nelt;
} IntBuf; // Extendable buffer of integers

typedef struct ibbuf {
	IntBuf *elts;
	int buflength;
	int nelt;
} IntBBuf; // Extendable buffer of extendable buffers of integers

typedef struct rangesbuf {
	IntBuf start;
	IntBuf width;
} RangesBuf; // Extendable buffer of integer ranges

typedef struct cbuf {
	char *elts;
	int buflength;
	int nelt;
} CharBuf; // Extendable buffer of chars

typedef struct cbbuf {
        CharBuf *elts;
        int buflength;
        int nelt;
} CharBBuf; // Extendable buffer of extendable buffers of chars

typedef struct ccarr {
	const char *elts;
	int nelt;
} ConstCharArr;

typedef struct ccaarr {
	ConstCharArr *elts;
	int nelt;
} ConstCharAArr;


/*
 * Match reporting modes (more modes will be added soon...)
 */
#define COUNT_MRMODE	1
#define START_MRMODE	2

