#include "Biostrings_interface.h"

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

const char *get_BString_charseq(SEXP x, int *length)
{
	static const char *(*fun)(SEXP, int *) = NULL;

	if (fun == NULL)
		fun = (const char *(*)(SEXP, int *)) R_GetCCallable("Biostrings", "_get_BString_charseq");
	return fun(x, length);
}

void init_match_reporting(int mode)
{
	static void (*fun)(int) = NULL;

	if (fun == NULL)
		fun = (void (*)(int)) R_GetCCallable("Biostrings", "_init_match_reporting");
	return fun(mode);
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

