#include "Biostrings_interface.h"

typedef void (*NamedSBuf_init_FUNTYPE)
	(NamedSBuf *, int, int);

void NamedSBuf_init(NamedSBuf *namedsbuf, int maxcount, int count)
{
	static NamedSBuf_init_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (NamedSBuf_init_FUNTYPE)
			R_GetCCallable("Biostrings", "_NamedSBuf_init");
	return fun(namedsbuf, maxcount, count);
}

typedef void (*NamedSBuf_insert_at_FUNTYPE)
	(NamedSBuf *, int, const char *, const char *);

void NamedSBuf_insert_at(NamedSBuf *namedsbuf, int at,
		const char *string, const char *name)
{
	static NamedSBuf_insert_at_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (NamedSBuf_insert_at_FUNTYPE)
			R_GetCCallable("Biostrings", "_NamedSBuf_insert_at");
	return fun(namedsbuf, at, string, name);
}

char DNAencode(char c)
{
	static char (*fun)(char) = NULL;

	if (fun == NULL)
		fun = (char (*)(char)) R_GetCCallable("Biostrings", "_DNAencode");
	return fun(c);
}

char DNAdecode(char code)
{
	static char (*fun)(char) = NULL;

	if (fun == NULL)
		fun = (char (*)(char)) R_GetCCallable("Biostrings", "_DNAdecode");
	return fun(code);
}

char RNAencode(char c)
{
	static char (*fun)(char) = NULL;

	if (fun == NULL)
		fun = (char (*)(char)) R_GetCCallable("Biostrings", "_RNAencode");
	return fun(c);
}

char RNAdecode(char code)
{
	static char (*fun)(char) = NULL;

	if (fun == NULL)
		fun = (char (*)(char)) R_GetCCallable("Biostrings", "_RNAdecode");
	return fun(code);
}

const char *get_XString_charseq(SEXP x, int *length)
{
	static const char *(*fun)(SEXP, int *) = NULL;

	if (fun == NULL)
		fun = (const char *(*)(SEXP, int *)) R_GetCCallable("Biostrings", "_get_XString_charseq");
	return fun(x, length);
}

int get_XStringList_length(SEXP x)
{
	static int (*fun)(SEXP) = NULL;

	if (fun == NULL)
		fun = (int (*)(SEXP)) R_GetCCallable("Biostrings", "_get_XStringList_length");
	return fun(x);
}

const char *get_XStringList_charseq(SEXP x, int i, int *nchar)
{
	static const char *(*fun)(SEXP, int, int *) = NULL;

	if (fun == NULL)
		fun = (const char *(*)(SEXP, int, int *)) R_GetCCallable("Biostrings", "_get_XStringList_charseq");
	return fun(x, i, nchar);
}

int get_XStringSet_length(SEXP x)
{
	static int (*fun)(SEXP) = NULL;

	if (fun == NULL)
		fun = (int (*)(SEXP)) R_GetCCallable("Biostrings", "_get_XStringSet_length");
	return fun(x);
}

const char *get_XStringSet_charseq(SEXP x, int i, int *nchar)
{
	static const char *(*fun)(SEXP, int, int *) = NULL;

	if (fun == NULL)
		fun = (const char *(*)(SEXP, int, int *)) R_GetCCallable("Biostrings", "_get_XStringSet_charseq");
	return fun(x, i, nchar);
}

void init_match_reporting(int mrmode)
{
	static void (*fun)(int) = NULL;

	if (fun == NULL)
		fun = (void (*)(int)) R_GetCCallable("Biostrings", "_init_match_reporting");
	return fun(mrmode);
}

int report_match(int start, int end)
{
	static int (*fun)(int, int) = NULL;

	if (fun == NULL)
		fun = (int (*)(int, int)) R_GetCCallable("Biostrings", "_report_match");
	return fun(start, end);
}

SEXP reported_matches_asSEXP()
{
	static SEXP (*fun)() = NULL;

	if (fun == NULL)
		fun = (SEXP (*)()) R_GetCCallable("Biostrings", "_reported_matches_asSEXP");
	return fun();
}

