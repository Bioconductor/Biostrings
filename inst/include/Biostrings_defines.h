#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/*
 * Extendable buffers used for temporary storage of data whose size is not
 * known in advance.
 */

typedef struct ibuf {
	int *vals;
	int maxcount;
	int count;
} IBuf; // Extendable buffer of integers

typedef struct ibbuf {
	IBuf *ibufs;
	int maxcount;
	int count;
} IBBuf; // Extendable buffer of arrays of integers

typedef struct rangesbuf {
	IBuf start;
	IBuf width;
} RangesBuf; // Extendable buffer of integer ranges

typedef struct cbuf {
	char *vals;
	int maxcount;
	int count;
} CBuf; // Extendable buffer of chars

typedef struct sbuf {
	char **strings;
	int maxcount;
	int count;
} SBuf; // Extendable buffer of strings

typedef struct namedsbuf {
	SBuf sbuf;
	SBuf names;
} NamedSBuf; // Extendable buffer of named strings


/*
 * Match reporting modes
 */
#define COUNT_MRMODE	1
#define START_MRMODE	2

// more modes to come soon...


