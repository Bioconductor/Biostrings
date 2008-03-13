#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/*
 * Extendable buffers used for temporary storage of incoming data whose size
 * is not known in advance:
 *
 *   o IntBuf:   extendable buffer of ints
 *   o IntBBuf:  extendable buffer of extendable buffers of ints
 *   o RangeBuf: extendable buffer of integer ranges
 *   o CharBuf:  extendable buffer of chars
 *   o CharBBuf: extendable buffer of extendable buffers of chars
 *
 * They are NOT an attempt to reinvent an SEXP subsystem. Some notable
 * differences are: (a) they are extendable (i.e. they are automatically
 * reallocated when more room is needed to add a new element), (b) they are
 * much faster, and (c) they don't require any PROTECT/UNPROTECT mechanism.
 */

typedef struct ibuf {
	int buflength;
	int *elts;
	int nelt;
} IntBuf;

typedef struct ibbuf {
	int buflength;
	IntBuf *elts;
	int nelt;
} IntBBuf; 

typedef struct rangebuf {
	IntBuf start;
	IntBuf width;
} RangeBuf;

typedef struct cbuf {
	int buflength;
	char *elts;
	int nelt;
} CharBuf; 

typedef struct cbbuf {
        int buflength;
        CharBuf *elts;
        int nelt;
} CharBBuf; 


/*
 * Two additional types:
 *
 *   o CharArr:  array of const chars (think of this as a pointer to a non
 *                  null-terminated sequence of chars)
 *   o CharAArr: array of arrays of const chars
 */
typedef struct carr {
	const char *elts;
	int nelt;
} CharArr;

typedef struct caarr {
	CharArr *elts;
	int nelt;
} CharAArr;


/*
 * Match reporting modes (more modes will be added soon...)
 */
#define COUNT_MRMODE	1
#define START_MRMODE	2

