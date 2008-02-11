#include "Biostrings_defines.h"

int DNAencode(char c);

char DNAdecode(char code);

const char *get_BString_charseq(SEXP x, int *length);

void init_match_reporting(int mode);

int report_match(int start, int end);

SEXP reported_matches_asSEXP();

