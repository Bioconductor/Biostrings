#include "Biostrings_interface.h"

typedef CharBBuf (*new_CharBBuf_FUNTYPE)(int, int);
CharBBuf new_CharBBuf(int buflength, int nelt)
{
	static new_CharBBuf_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (new_CharBBuf_FUNTYPE)
			R_GetCCallable("Biostrings", "_new_CharBBuf");
	return fun(buflength, nelt);
}

typedef void (*append_string_to_CharBBuf_FUNTYPE)(CharBBuf *, const char *);
void append_string_to_CharBBuf(CharBBuf *cbbuf, const char *string)
{
	static append_string_to_CharBBuf_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (append_string_to_CharBBuf_FUNTYPE)
			R_GetCCallable("Biostrings", "_append_string_to_CharBBuf");
	return fun(cbbuf, string);
}

typedef ConstCharAArr (*new_ConstCharAArr_from_CharBBuf_FUNTYPE)(CharBBuf);
ConstCharAArr new_ConstCharAArr_from_CharBBuf(CharBBuf cbbuf)
{
	static new_ConstCharAArr_from_CharBBuf_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (new_ConstCharAArr_from_CharBBuf_FUNTYPE)
			R_GetCCallable("Biostrings", "_new_ConstCharAArr_from_CharBBuf");
	return fun(cbbuf);
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

typedef SEXP (*new_XStringSet_from_seqsnames_FUNTYPE)(const char *, ConstCharAArr, ConstCharAArr);
SEXP new_XStringSet_from_seqsnames(const char *baseClass, ConstCharAArr seqs, ConstCharAArr names)
{
	static new_XStringSet_from_seqsnames_FUNTYPE fun = NULL;

	if (fun == NULL)
		fun = (new_XStringSet_from_seqsnames_FUNTYPE)
			R_GetCCallable("Biostrings", "_new_XStringSet_from_seqsnames");
	return fun(baseClass, seqs, names);
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

