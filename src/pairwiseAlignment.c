#include "Biostrings.h"

#define  GLOBAL_ALIGNMENT 1
#define   LOCAL_ALIGNMENT 2
#define OVERLAP_ALIGNMENT 3

#define REPLACEMENT 0
#define INSERTION   1
#define DELETION    2
#define NOTHING     3

#define F_MATRIX(i, j) (fMatrix[(nCharString2 + 1) * i + j])
#define TRACE_MATRIX(i, j) (traceMatrix[(nCharString2 + 1) * i + j])

#define SET_LOOKUP_VALUE(lookupTable, length, key) \
{ \
	unsigned char lookupKey = (unsigned char) (key); \
	if (lookupKey >= (length) || (lookupValue = (lookupTable)[lookupKey]) == NA_INTEGER) { \
		error("key %d not in lookup table", (int) lookupKey); \
	} \
}

static int nCharAligned = 0;
static char *align1Buffer, *align2Buffer, *align1, *align2;

/* Returns the score of the alignment */
static int pairwiseAlignment(
		RoSeq stringElements1,
		RoSeq stringElements2,
		const int *matchScores,
		const int *matchScoresDim,
		const int *lookupTable,
		int lookupTableLength,
		int gapExtension,
		char gapCode,
		int typeCode)
{
	int i, j, iMinus1, jMinus1;

	/* Step 1:  Get information on input strings */
	int nCharString1, nCharString2;
	nCharString1 = stringElements1.nelt;
	nCharString2 = stringElements2.nelt;

	/* Step 2:  Create objects for scores and traceback values */
	int *fMatrix = (int *) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(int));
	if (typeCode == GLOBAL_ALIGNMENT) {
		for (i = 0; i <= nCharString1; i++)
			F_MATRIX(i, 0) = i * gapExtension;
		for (j = 0; j <= nCharString2; j++)
			F_MATRIX(0, j) = j * gapExtension;
	} else if (typeCode == LOCAL_ALIGNMENT || typeCode == OVERLAP_ALIGNMENT) {
		for (i = 0; i <= nCharString1; i++)
			F_MATRIX(i, 0) = 0;
		for (j = 0; j <= nCharString2; j++)
			F_MATRIX(0, j) = 0;
	}

	int traceValue;
	int **traceMatrix = (int **) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(int*));
	int *possibleTraceValues = (int *) R_alloc((long) 4, sizeof(int));
	possibleTraceValues[REPLACEMENT] = REPLACEMENT;
	possibleTraceValues[INSERTION] = INSERTION;
	possibleTraceValues[DELETION] = DELETION;
	possibleTraceValues[NOTHING] = NOTHING;

	/* Step 3:  Generate scores and traceback values */
	int score;
	int startIndex = -1;
	int startScore = -2147483646;
	for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
		for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
			int lookupValue;
			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, stringElements1.elts[iMinus1]);
			int element1 = lookupValue;
			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, stringElements2.elts[jMinus1]);
			int element2 = lookupValue;
			int scoreReplacement = F_MATRIX(iMinus1, jMinus1) + matchScores[matchScoresDim[0] * element1 + element2];
			int scoreDeletion    = F_MATRIX(iMinus1, j) + gapExtension;
			int scoreInsertion   = F_MATRIX(i, jMinus1) + gapExtension;
			if (scoreDeletion >= scoreInsertion) {
				score = scoreDeletion;
				traceValue = DELETION;
			} else {
				score = scoreInsertion;
				traceValue = INSERTION;
			}
			if (scoreReplacement >= score) {
				score = scoreReplacement;
				traceValue = REPLACEMENT;
			}
			if (typeCode == LOCAL_ALIGNMENT) {
				if (score < 0) {
					score = 0;
					traceValue = NOTHING;
				}
				if (score > startScore) {
					startIndex = (nCharString2 + 1) * i + j;
					startScore = score;
				}
			}
			F_MATRIX(i, j) = score;
			TRACE_MATRIX(i, j) = &possibleTraceValues[traceValue];
		}
	}
	if (typeCode == GLOBAL_ALIGNMENT) {
		startIndex = (nCharString2 + 1) * nCharString1 + nCharString2;
		startScore = F_MATRIX(nCharString1, nCharString2);
	} else if (typeCode == OVERLAP_ALIGNMENT) {
		for (i = 1; i <= nCharString1; i++) {
			score = F_MATRIX(i, nCharString2);
			if (score > startScore) {
				startIndex = (nCharString2 + 1) * i + nCharString2;
				startScore = score;
			}
	    }
		for (j = 1; j <= nCharString2; j++) {
			score = F_MATRIX(nCharString1, j);
			if (score > startScore) {
				startIndex = (nCharString2 + 1) * nCharString1 + j;
				startScore = score;
			}
	    }
	}

	/* Step 4:  Get a starting location for the traceback */
	nCharAligned = 0;
	int alignmentBufferSize = nCharString1 + nCharString2;
	align1Buffer = (char *) R_alloc((long) alignmentBufferSize, sizeof(char));
	align2Buffer = (char *) R_alloc((long) alignmentBufferSize, sizeof(char));
	align1 = align1Buffer + alignmentBufferSize;
	align2 = align2Buffer + alignmentBufferSize;
	int startRow = startIndex / (nCharString2 + 1);
	int startCol = startIndex - (nCharString2 + 1) * startRow;
	if (typeCode == OVERLAP_ALIGNMENT && (startRow < nCharString1 || startCol < nCharString2)) {
		if (startRow == nCharString1) {
			nCharAligned += nCharString2 - startCol;
			for (j = 1; j <= nCharString2 - startCol; j++) {
				align1--;
				align2--;
				*align1 = gapCode;
				*align2 = stringElements2.elts[nCharString2 - j];
			}
	    } else {
			nCharAligned += nCharString1 - startRow;
			for (i = 1; i <= nCharString1 - startRow; i++) {
				align1--;
				align2--;
				*align1 = stringElements1.elts[nCharString1 - i];
				*align2 = gapCode;
			}
	    }
	}

	/* Step 5:  Traceback through the score matrix */
	i = startRow;
	j = startCol;
	while (((typeCode == GLOBAL_ALIGNMENT || typeCode == OVERLAP_ALIGNMENT) && (i >= 1 || j >= 1)) ||
			(typeCode == LOCAL_ALIGNMENT && F_MATRIX(i, j) > 0)) {
		nCharAligned++;
		align1--;
		align2--;
		iMinus1 = i - 1;
		jMinus1 = j - 1;
		if (j == 0)
			traceValue = DELETION;
		else if (i == 0)
			traceValue = INSERTION;
		else
			traceValue = *TRACE_MATRIX(i, j);
		switch (traceValue) {
		    case DELETION:
			*align1 = stringElements1.elts[iMinus1];
			*align2 = gapCode;
			i--;
			break;
		    case INSERTION:
			*align1 = gapCode;
			*align2 = stringElements2.elts[jMinus1];
			j--;
			break;
		    case REPLACEMENT:
			*align1 = stringElements1.elts[iMinus1];
			*align2 = stringElements2.elts[jMinus1];
			i--;
			j--;
			break;
		    case NOTHING:
			break;
		    default:
			error("unknown traceback code %d", traceValue);
			break;
		}
	}

	return startScore;
}

/*
 * INPUTS
 * 'string1', 'string2':  left and right XString objects
 * 'matchScores':  scoring matrix for matches/mismatches (integer matrix)
 * 'matchScoresDim':  dimension of 'matchScores' (integer vector of length 2
 * 'lookupTable':  lookup table for translating XString bytes to scoring matrix
 *                 indices (integer vector)
 * 'gapExtension':  gap cost or penalty (integer vector of length 1)
 * 'gapCode':  encoded value of the '-' letter (raw vector of length 1)
 * 'typeCode':  type of pairwise alignment
 *          (integer vector of length 1; 1 = 'global', 2 = 'local', 3 = 'overlap')
 * 
 * OUTPUT
 * Return a named list with 3 elements: 2 "externalptr" objects describing
 * the alignments and 1 integer vector containing the alignment score.
 * Note that the 2 XString objects to align should contain no gaps.
 */
SEXP align_pairwiseAlignment(
		SEXP string1,
		SEXP string2,
		SEXP matchScores,
		SEXP matchScoresDim,
		SEXP lookupTable,
		SEXP gapExtension,
		SEXP gapCode,
		SEXP typeCode)
{
	int score;
	RoSeq stringElements1, stringElements2;
	SEXP answer, answerNames, answerElements, tag;

	stringElements1 = _get_XString_asRoSeq(string1);
	stringElements2 = _get_XString_asRoSeq(string2);
	score = pairwiseAlignment(
			stringElements1,
			stringElements2,
			INTEGER(matchScores),
			INTEGER(matchScoresDim),
			INTEGER(lookupTable),
			LENGTH(lookupTable),
			INTEGER(gapExtension)[0],
			(char) RAW(gapCode)[0],
			INTEGER(typeCode)[0]);

	PROTECT(answer = NEW_LIST(3));
	/* set the names */
	PROTECT(answerNames = NEW_CHARACTER(3));
	SET_STRING_ELT(answerNames, 0, mkChar("align1"));
	SET_STRING_ELT(answerNames, 1, mkChar("align2"));
	SET_STRING_ELT(answerNames, 2, mkChar("score"));
	SET_NAMES(answer, answerNames);
	UNPROTECT(1);
	/* set the "align1" element */
	PROTECT(tag = NEW_RAW(nCharAligned));
	memcpy((char *) RAW(tag), align1, nCharAligned * sizeof(char));
	PROTECT(answerElements = _new_XRaw(tag));
	SET_ELEMENT(answer, 0, answerElements);
	UNPROTECT(2);
	/* set the "align2" element */
	PROTECT(tag = NEW_RAW(nCharAligned));
	memcpy((char *) RAW(tag), align2, nCharAligned * sizeof(char));
	PROTECT(answerElements = _new_XRaw(tag));
	SET_ELEMENT(answer, 1, answerElements);
	UNPROTECT(2);
	/* set the "score" element */
	PROTECT(answerElements = NEW_INTEGER(1));
	INTEGER(answerElements)[0] = score;
	SET_ELEMENT(answer, 2, answerElements);
	UNPROTECT(1);
	/* answer is ready */
	UNPROTECT(1);
	return answer;
}
