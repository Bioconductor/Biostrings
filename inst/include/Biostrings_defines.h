#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/*
 * Extendable buffers used for temporary storage of incoming data whose size
 * is not known in advance.
 * They are NOT an attempt to reinvent an SEXP subsystem. Some notable
 * differences are: they are extendable (i.e. they are automatically
 * reallocated when more room is needed to add a new element), they are much
 * faster, and they don't need any PROTECT/UNPROTECT mechanism.
 */

typedef struct ibuf {
	int *elts;
	int buflength;
	int nelt;
} IBuf; // Extendable buffer of integers

typedef struct ibbuf {
	IBuf *elts;
	int buflength;
	int nelt;
} IBBuf; // Extendable buffer of extendable buffers of integers

typedef struct rangesbuf {
	IBuf start;
	IBuf width;
} RangesBuf; // Extendable buffer of integer ranges

typedef struct cbuf {
	char *elts;
	int buflength;
	int nelt;
} CBuf; // Extendable buffer of chars

typedef struct cbbuf {
        CBuf *elts;
        int buflength;
        int nelt;
} CBBuf; // Extendable buffer of extendable buffers of chars


/*
 * Match reporting modes
 */
#define COUNT_MRMODE	1
#define START_MRMODE	2

// more modes to come soon...


