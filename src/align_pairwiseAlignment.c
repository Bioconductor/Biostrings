#include <float.h>
#include "Biostrings.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

#define NEGATIVE_INFINITY (- FLT_MAX)

#define  GLOBAL_ALIGNMENT 1
#define   LOCAL_ALIGNMENT 2
#define OVERLAP_ALIGNMENT 3

#define SUBSTITUTION 'S'
#define INSERTION    'I'
#define DELETION     'D'

#define F_MATRIX(i, j) (fMatrix[(nCharString2 + 1) * i + j])
#define H_MATRIX(i, j) (hMatrix[(nCharString2 + 1) * i + j])
#define V_MATRIX(i, j) (vMatrix[(nCharString2 + 1) * i + j])
#define TRACE_MATRIX(i, j) (traceMatrix[(nCharString2 + 1) * i + j])

#define SAFE_SUM(x, y) (x == NEGATIVE_INFINITY ? x : (x + y))

#define SET_LOOKUP_VALUE(lookupTable, length, key) \
{ \
	unsigned char lookupKey = (unsigned char) (key); \
	if (lookupKey >= (length) || (lookupValue = (lookupTable)[lookupKey]) == NA_INTEGER) { \
		error("key %d not in lookup table", (int) lookupKey); \
	} \
}


/* Structure to hold alignment information */
struct AlignInfo {
	int endMatch;
	int widthMatch;
	int* endInserts;
	int* widthInserts;
	int lengthInserts;
};
void function(struct AlignInfo *);


/* Returns the score of the optimal pairwise alignment */
static float pairwiseAlignment(
		RoSeq stringElements1,
		RoSeq stringElements2,
		RoSeq qualityElements1,
		RoSeq qualityElements2,
		int typeCode,
		int scoreOnly,
		int endGap,
		float gapOpening,
		float gapExtension,
		const int *qualityLookupTable,
		int qualityLookupTableLength,
		const double *qualityMatchMatrix,
		const double *qualityMismatchMatrix,
		const int *qualityMatrixDim,
		const int *constantLookupTable,
		int constantLookupTableLength,
		const double *constantMatrix,
		const int *constantMatrixDim,
		struct AlignInfo *align1InfoPtr,
		struct AlignInfo *align2InfoPtr)
{
	int i, j, iMinus1, jMinus1;

	/* Step 1:  Get information on input XString objects */
	int nCharString1 = stringElements1.nelt;
	int nCharString2 = stringElements2.nelt;
	int nQuality1 = qualityElements1.nelt;
	int nQuality2 = qualityElements2.nelt;

	/* Step 2:  Create objects for scores and traceback values */
	float *fMatrix, *hMatrix, *vMatrix;
	fMatrix = (float *) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(float));
	if (gapOpening == 0.0) {
		hMatrix = fMatrix;
		vMatrix = fMatrix;
	} else {
		hMatrix = (float *) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(float));
		vMatrix = (float *) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(float));		
	}
	if (typeCode == GLOBAL_ALIGNMENT && endGap) {
		for (i = 0; i <= nCharString1; i++)
			H_MATRIX(i, 0) = gapOpening + i * gapExtension;
		for (j = 0; j <= nCharString2; j++)
			V_MATRIX(0, j) = gapOpening + j * gapExtension;
	} else {
		for (i = 0; i <= nCharString1; i++)
			H_MATRIX(i, 0) = 0.0;
		for (j = 0; j <= nCharString2; j++)
			V_MATRIX(0, j) = 0.0;
	}
	if (gapOpening < 0) {
		F_MATRIX(0, 0) = 0;
		for (i = 1; i <= nCharString1; i++) {
			F_MATRIX(i, 0) = NEGATIVE_INFINITY;
			V_MATRIX(i, 0) = NEGATIVE_INFINITY;
		}
		for (j = 1; j <= nCharString2; j++) {
			F_MATRIX(0, j) = NEGATIVE_INFINITY;
			H_MATRIX(0, j) = NEGATIVE_INFINITY;
		}
	}
	char *traceMatrix = (char *) R_alloc((long) (nCharString1 + 1) * (nCharString2 + 1), sizeof(char));
	for (i = 0; i <= nCharString1; i++)
		TRACE_MATRIX(i, 0) = DELETION;
	for (j = 0; j <= nCharString2; j++)
		TRACE_MATRIX(0, j) = INSERTION;

	/* Step 3:  Generate scores and traceback values */
	RoSeq sequence1, sequence2;
	int scalar1, scalar2;
	const int *lookupTable;
	int lookupTableLength;
	const double *matchMatrix, *mismatchMatrix;
	const int *matrixDim;
	if (nQuality1 == 0) {
		sequence1 = stringElements1;
		sequence2 = stringElements2;
		scalar1 = (nCharString1 == 1);
		scalar2 = (nCharString2 == 1);
		lookupTable = constantLookupTable;
		lookupTableLength = constantLookupTableLength;
		matchMatrix = constantMatrix;
		mismatchMatrix = constantMatrix;
		matrixDim = constantMatrixDim;
	} else {
		sequence1 = qualityElements1;
		sequence2 = qualityElements2;
		scalar1 = (nQuality1 == 1);
		scalar2 = (nQuality2 == 1);
		lookupTable = qualityLookupTable;
		lookupTableLength = qualityLookupTableLength;
		matchMatrix = qualityMatchMatrix;
		mismatchMatrix = qualityMismatchMatrix;
		matrixDim = qualityMatrixDim;
	}
	align1InfoPtr->endMatch = -1;
	align2InfoPtr->endMatch = -1;
	align1InfoPtr->widthMatch = 0;
	align2InfoPtr->widthMatch = 0;
	int lookupValue, element1, element2;
	float substitutionValue, score, maxScore = NEGATIVE_INFINITY;
	if (gapOpening == 0) {
		for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iMinus1]);
			element1 = lookupValue;
			for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
				SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jMinus1]);
				element2 = lookupValue;
				if (stringElements1.elts[iMinus1] == stringElements2.elts[jMinus1])
					substitutionValue = (float) matchMatrix[matrixDim[0] * element1 + element2];
				else
					substitutionValue = (float) mismatchMatrix[matrixDim[0] * element1 + element2];

				float scoreSubstitution = F_MATRIX(iMinus1, jMinus1) + substitutionValue;
				float scoreInsertion    = V_MATRIX(i, jMinus1) + gapExtension;
				float scoreDeletion     = H_MATRIX(iMinus1, j) + gapExtension;
				if (typeCode == LOCAL_ALIGNMENT) {
					if (scoreSubstitution < 0.0)
						scoreSubstitution = 0.0;
					if (scoreSubstitution > maxScore) {
						align1InfoPtr->endMatch = i;
						align2InfoPtr->endMatch = j;
						maxScore = scoreSubstitution;
					}
				}
				if (scoreSubstitution >= MAX(scoreInsertion, scoreDeletion)) {
					F_MATRIX(i, j) = scoreSubstitution;
					TRACE_MATRIX(i, j) = SUBSTITUTION;
				} else if (scoreInsertion >= scoreDeletion) {
					F_MATRIX(i, j) = scoreInsertion;
					TRACE_MATRIX(i, j) = INSERTION;
				} else {
					F_MATRIX(i, j) = scoreDeletion;
					TRACE_MATRIX(i, j) = DELETION;
				}
			}
		}
	} else {
		for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iMinus1]);
			element1 = lookupValue;
			for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
				SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jMinus1]);
				element2 = lookupValue;
				if (stringElements1.elts[iMinus1] == stringElements2.elts[jMinus1])
					substitutionValue = (float) matchMatrix[matrixDim[0] * element1 + element2];
				else
					substitutionValue = (float) mismatchMatrix[matrixDim[0] * element1 + element2];

				F_MATRIX(i, j) =
					SAFE_SUM(MAX(F_MATRIX(iMinus1, jMinus1), MAX(H_MATRIX(iMinus1, jMinus1), V_MATRIX(iMinus1, jMinus1))),
					         substitutionValue);
				if (typeCode == LOCAL_ALIGNMENT) {
					if (F_MATRIX(i, j) < 0.0)
						F_MATRIX(i, j) = 0.0;
					if (F_MATRIX(i, j) > maxScore) {
						align1InfoPtr->endMatch = i;
						align2InfoPtr->endMatch = j;
						maxScore = F_MATRIX(i, j);
					}
				}
				H_MATRIX(i, j) = 
					MAX(SAFE_SUM(F_MATRIX(iMinus1, j), gapOpening + gapExtension),
					    SAFE_SUM(H_MATRIX(iMinus1, j), gapExtension));
				V_MATRIX(i, j) =
					MAX(SAFE_SUM(F_MATRIX(i, jMinus1), gapOpening + gapExtension),
					    SAFE_SUM(V_MATRIX(i, jMinus1), gapExtension));

				float scoreSubstitution = F_MATRIX(i, j);
				float scoreInsertion    = V_MATRIX(i, j);
				float scoreDeletion     = H_MATRIX(i, j);

				if (scoreSubstitution >= MAX(scoreInsertion, scoreDeletion))
					TRACE_MATRIX(i, j) = SUBSTITUTION;
				else if (scoreInsertion >= scoreDeletion)
					TRACE_MATRIX(i, j) = INSERTION;
				else
					TRACE_MATRIX(i, j) = DELETION;
			}
		}
	}
	if (typeCode == GLOBAL_ALIGNMENT) {
		if (typeCode == GLOBAL_ALIGNMENT && endGap == 0) {
			for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
				if (H_MATRIX(i, nCharString2) < H_MATRIX(iMinus1, nCharString2)) {
					H_MATRIX(i, nCharString2) = H_MATRIX(iMinus1, nCharString2);
					if (H_MATRIX(i, nCharString2) >= MAX(F_MATRIX(i, nCharString2), V_MATRIX(i, nCharString2)))
						TRACE_MATRIX(i, nCharString2) = DELETION;
				}
			}
			for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
				if (V_MATRIX(nCharString1, j) < V_MATRIX(nCharString1, jMinus1)) {
					V_MATRIX(nCharString1, j) = V_MATRIX(nCharString1, jMinus1);
					if (V_MATRIX(nCharString1, j) >= MAX(F_MATRIX(nCharString1, j), H_MATRIX(nCharString1, j)))
						TRACE_MATRIX(nCharString1, j) = INSERTION;
				}
			}
		}
		align1InfoPtr->endMatch = nCharString1;
		align2InfoPtr->endMatch = nCharString2;
		maxScore =
			MAX(F_MATRIX(nCharString1, nCharString2),
			MAX(H_MATRIX(nCharString1, nCharString2),
			    V_MATRIX(nCharString1, nCharString2)));
	} else if (typeCode == OVERLAP_ALIGNMENT) {
		for (i = nCharString1; i >= 1; i--) {
			score =
				MAX(F_MATRIX(i, nCharString2),
				MAX(H_MATRIX(i, nCharString2),
				    V_MATRIX(i, nCharString2)));
			if (score > maxScore) {
				align1InfoPtr->endMatch = i;
				align2InfoPtr->endMatch = nCharString2;
				maxScore = score;
			}
	    }
		for (j = nCharString2; j >= 1; j--) {
			score =
				MAX(F_MATRIX(nCharString1, j),
				MAX(H_MATRIX(nCharString1, j),
				    V_MATRIX(nCharString1, j)));
			if (score > maxScore) {
				align1InfoPtr->endMatch = nCharString1;
				align2InfoPtr->endMatch = j;
				maxScore = score;
			}
	    }
	}

	if (scoreOnly == 0) {
		/* Step 4:  Traceback through the score matrix */
		char previousAction = '?';
		i = align1InfoPtr->endMatch;
		j = align2InfoPtr->endMatch;
		while ((i > 0 && j > 0) && !(typeCode == LOCAL_ALIGNMENT && F_MATRIX(i, j) == 0)) {
			iMinus1 = i - 1;
			jMinus1 = j - 1;
			char action = TRACE_MATRIX(i, j);
			switch (action) {
		    	case DELETION:
		    		if (j == nCharString2) {
		    			align1InfoPtr->endMatch--;
		    		} else {
		    			align1InfoPtr->widthMatch++;
			    		if (previousAction != DELETION) {
							align2InfoPtr->endInserts--;
							align2InfoPtr->widthInserts--;
							align2InfoPtr->lengthInserts++;
							*align2InfoPtr->endInserts = j;
			    		}
		    		}
		    		*align2InfoPtr->widthInserts += 1;
		    		i--;
		    		break;
		    	case INSERTION:
		    		if (i == nCharString1) {
		    			align2InfoPtr->endMatch--;
		    		} else {
		    			align2InfoPtr->widthMatch++;
			    		if (previousAction != INSERTION) {
							align1InfoPtr->endInserts--;
							align1InfoPtr->widthInserts--;
							align1InfoPtr->lengthInserts++;
							*align1InfoPtr->endInserts = i;
			    		}
		    		}
		    		*align1InfoPtr->widthInserts += 1;
		    		j--;
		    		break;
		    	case SUBSTITUTION:
		    		align1InfoPtr->widthMatch++;
		    		align2InfoPtr->widthMatch++;
		    		i--;
		    		j--;
		    		break;
		    	default:
		    		error("unknown traceback code %d", action);
		    		break;
			}
			previousAction = action;
		}
		int offset1 = align1InfoPtr->endMatch - align1InfoPtr->widthMatch;
		if (offset1 > 0 && align1InfoPtr->lengthInserts > 0) {
			for (i = 0; i < align1InfoPtr->lengthInserts; i++)
				align1InfoPtr->endInserts[i] -= offset1;
		}
		int offset2 = align2InfoPtr->endMatch - align2InfoPtr->widthMatch;
		if (offset2 > 0 && align2InfoPtr->lengthInserts > 0) {
			for (j = 0; j < align2InfoPtr->lengthInserts; j++)
				align2InfoPtr->endInserts[j] -= offset2;
		}
	}

	return maxScore;
}

/*
 * INPUTS
 * 'string1', 'string2':     left and right XString objects for reads
 * 'quality1', 'quality2':   left and right BString objects for quality scores
 * 'typeCode':               type of pairwise alignment
 *                           (integer vector of length 1;
 *                            1 = 'global', 2 = 'local', 3 = 'overlap')
 * 'scoreOnly':              denotes whether or not to only return the scores
 *                           of the optimal pairwise alignment
 *                           (logical vector of length 1)
 * 'endGap':                 denote whether or not to penalize end gaps like
 *                           internal gaps.
 *                           (logical vector of length 1)
 * 'gapOpening':             gap opening cost or penalty
 *                           (double vector of length 1)
 * 'gapExtension':           gap extension cost or penalty
 *                           (double vector of length 1)
 * 'qualityLookupTable':     lookup table for translating BString bytes to
 *                           quality-based scoring matrix indices
 *                           (integer vector)
 * 'qualityMatchMatrix':     quality-based substitution matrix for matches
 * 'qualityMismatchMatrix':  quality-based substitution matrix for matches
 * 'qualityMatrixDim':       dimension of 'qualityMatchMatrix' and
 *                           'qualityMismatchMatrix'
 *                           (integer vector of lenth 2)
 * 'constantLookupTable':    lookup table for translating XString bytes to scoring
 *                           matrix indices
 *                           (integer vector)
 * 'constantMatrix':         constant substitution matrix for matches/mismatches
 *                           (double matrix)
 * 'constantMatrixDim':      dimension of 'constantMatrix'
 *                           (integer vector of length 2)
 * 
 * OUTPUT
 * Return a named list with 3 elements: 2 "externalptr" objects describing
 * the alignments and 1 integer vector containing the alignment score.
 * Note that the 2 XString objects to align should contain no gaps.
 */
SEXP align_pairwiseAlignment(
		SEXP string1,
		SEXP string2,
		SEXP quality1,
		SEXP quality2,
		SEXP typeCode,
		SEXP scoreOnly,
		SEXP endGap,
		SEXP gapOpening,
		SEXP gapExtension,
		SEXP qualityLookupTable,
		SEXP qualityMatchMatrix,
		SEXP qualityMismatchMatrix,
		SEXP qualityMatrixDim,
		SEXP constantLookupTable,
		SEXP constantMatrix,
		SEXP constantMatrixDim)
{
	RoSeq stringElements1 = _get_XString_asRoSeq(string1);
	RoSeq stringElements2 = _get_XString_asRoSeq(string2);
	RoSeq qualityElements1 = _get_XString_asRoSeq(quality1);
	RoSeq qualityElements2 = _get_XString_asRoSeq(quality2);
	struct AlignInfo align1Info, align2Info;
	int alignmentBufferSize = MIN(stringElements1.nelt, stringElements2.nelt) + 1;
	int* endInserts1Buffer   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* endInserts2Buffer   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* widthInserts1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* widthInserts2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	for(int i = 0; i < alignmentBufferSize; i++) {
		endInserts1Buffer[i] = 0;
		endInserts2Buffer[i] = 0;
		widthInserts1Buffer[i] = 0;
		widthInserts2Buffer[i] = 0;
	}
	align1Info.lengthInserts = 0;
	align2Info.lengthInserts = 0;
	align1Info.endInserts = endInserts1Buffer + alignmentBufferSize;
	align2Info.endInserts = endInserts2Buffer + alignmentBufferSize;
	align1Info.widthInserts = widthInserts1Buffer + alignmentBufferSize;
	align2Info.widthInserts = widthInserts2Buffer + alignmentBufferSize;
	float score = pairwiseAlignment(
			stringElements1,
			stringElements2,
			qualityElements1,
			qualityElements2,
			INTEGER(typeCode)[0],
			LOGICAL(scoreOnly)[0],
			LOGICAL(endGap)[0],
			(float) REAL(gapOpening)[0],
			(float) REAL(gapExtension)[0],
			INTEGER(qualityLookupTable),
			LENGTH(qualityLookupTable),
			REAL(qualityMatchMatrix),
			REAL(qualityMismatchMatrix),
			INTEGER(qualityMatrixDim),
			INTEGER(constantLookupTable),
			LENGTH(constantLookupTable),
			REAL(constantMatrix),
			INTEGER(constantMatrixDim),
			&align1Info,
			&align2Info);

	SEXP answer, answerNames, answerElement1, answerElement2;
	PROTECT(answer = NEW_LIST(9));
	/* set the names */
	PROTECT(answerNames = NEW_CHARACTER(9));
	SET_STRING_ELT(answerNames, 0, mkChar("endMatch1"));
	SET_STRING_ELT(answerNames, 1, mkChar("widthMatch1"));
	SET_STRING_ELT(answerNames, 2, mkChar("endInserts1"));
	SET_STRING_ELT(answerNames, 3, mkChar("widthInserts1"));
	SET_STRING_ELT(answerNames, 4, mkChar("endMatch2"));
	SET_STRING_ELT(answerNames, 5, mkChar("widthMatch2"));
	SET_STRING_ELT(answerNames, 6, mkChar("endInserts2"));
	SET_STRING_ELT(answerNames, 7, mkChar("widthInserts2"));
	SET_STRING_ELT(answerNames, 8, mkChar("score"));
	SET_NAMES(answer, answerNames);
	UNPROTECT(1);
	/* set the "endMatch1" element */
	PROTECT(answerElement1 = NEW_INTEGER(1));
	INTEGER(answerElement1)[0] = align1Info.endMatch;
	SET_ELEMENT(answer, 0, answerElement1);
	UNPROTECT(1);
	/* set the "widthMatch1" element */
	PROTECT(answerElement1 = NEW_INTEGER(1));
	INTEGER(answerElement1)[0] = align1Info.widthMatch;
	SET_ELEMENT(answer, 1, answerElement1);
	UNPROTECT(1);
	/* set the "endInserts1" and "widthInserts1" elements */
	PROTECT(answerElement1 = NEW_INTEGER(align1Info.lengthInserts));
	PROTECT(answerElement2 = NEW_INTEGER(align1Info.lengthInserts));
	for(int i = 0; i < align1Info.lengthInserts; i++) {
		INTEGER(answerElement1)[i] = align1Info.endInserts[i];
		INTEGER(answerElement2)[i] = align1Info.widthInserts[i];
	}
	SET_ELEMENT(answer, 2, answerElement1);
	SET_ELEMENT(answer, 3, answerElement2);
	UNPROTECT(2);
	/* set the "endMatch2" element */
	PROTECT(answerElement1 = NEW_INTEGER(1));
	INTEGER(answerElement1)[0] = align2Info.endMatch;
	SET_ELEMENT(answer, 4, answerElement1);
	UNPROTECT(1);
	/* set the "widthMatch2" element */
	PROTECT(answerElement1 = NEW_INTEGER(1));
	INTEGER(answerElement1)[0] = align2Info.widthMatch;
	SET_ELEMENT(answer, 5, answerElement1);
	UNPROTECT(1);
	/* set the "endInserts2" and "widthInserts2" elements */
	PROTECT(answerElement1 = NEW_INTEGER(align2Info.lengthInserts));
	PROTECT(answerElement2 = NEW_INTEGER(align2Info.lengthInserts));
	for(int i = 0; i < align2Info.lengthInserts; i++) {
		INTEGER(answerElement1)[i] = align2Info.endInserts[i];
		INTEGER(answerElement2)[i] = align2Info.widthInserts[i];
	}
	SET_ELEMENT(answer, 6, answerElement1);
	SET_ELEMENT(answer, 7, answerElement2);
	UNPROTECT(2);
	/* set the "score" element */
	PROTECT(answerElement1 = NEW_NUMERIC(1));
	REAL(answerElement1)[0] = score;
	SET_ELEMENT(answer, 8, answerElement1);
	UNPROTECT(1);
	/* answer is ready */
	UNPROTECT(1);
	return answer;
}
