#include <float.h>
#include "Biostrings.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

#define NEGATIVE_INFINITY (- FLT_MAX)

#define   GLOBAL_ALIGNMENT 1
#define    LOCAL_ALIGNMENT 2
#define  OVERLAP_ALIGNMENT 3
#define OVERLAP1_ALIGNMENT 4
#define OVERLAP2_ALIGNMENT 5

#define SUBSTITUTION 'S'
#define INSERTION    'I'
#define DELETION     'D'

#define F_MATRIX(i, j) (fMatrix[nCharString2Plus1 * i + j])
#define H_MATRIX(i, j) (hMatrix[nCharString2Plus1 * i + j])
#define V_MATRIX(i, j) (vMatrix[nCharString2Plus1 * i + j])

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
	int startMatch;
	int widthMatch;
	int* startInserts;
	int* widthInserts;
	int lengthInserts;
	float* profile;
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
	int endGap1 = (typeCode == GLOBAL_ALIGNMENT || typeCode == OVERLAP2_ALIGNMENT);
	int endGap2 = (typeCode == GLOBAL_ALIGNMENT || typeCode == OVERLAP1_ALIGNMENT);

	/* Step 1:  Get information on input XString objects */
	int nCharString1 = stringElements1.nelt;
	int nCharString2 = stringElements2.nelt;
	int nCharString1Plus1 = nCharString1 + 1;
	int nCharString2Plus1 = nCharString2 + 1;
	int nCharString1Minus1 = nCharString1 - 1;
	int nCharString2Minus1 = nCharString2 - 1;
	int nQuality1 = qualityElements1.nelt;
	int nQuality2 = qualityElements2.nelt;

	/* Step 2:  Create objects for scores and traceback values */
	float *fMatrix = (float *) R_alloc((long) nCharString1Plus1 * nCharString2Plus1, sizeof(float));
	float *hMatrix = (float *) R_alloc((long) nCharString1Plus1 * nCharString2Plus1, sizeof(float));
	float *vMatrix = (float *) R_alloc((long) nCharString1Plus1 * nCharString2Plus1, sizeof(float));
	if (endGap1) {
		for (i = 0; i <= nCharString1; i++)
			H_MATRIX(i, 0) = gapOpening + i * gapExtension;
	} else {
		for (i = 0; i <= nCharString1; i++)
			H_MATRIX(i, 0) = 0.0;
	}
	if (endGap2) {
		for (j = 0; j <= nCharString2; j++)
			V_MATRIX(0, j) = gapOpening + j * gapExtension;
	} else {
		for (j = 0; j <= nCharString2; j++)
			V_MATRIX(0, j) = 0.0;
	}
	F_MATRIX(0, 0) = 0.0;
	for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
		F_MATRIX(i, 0) = NEGATIVE_INFINITY;
		V_MATRIX(i, 0) = NEGATIVE_INFINITY;
		align1InfoPtr->profile[iMinus1] = NEGATIVE_INFINITY;
	}
	for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
		F_MATRIX(0, j) = NEGATIVE_INFINITY;
		H_MATRIX(0, j) = NEGATIVE_INFINITY;
		align2InfoPtr->profile[jMinus1] = NEGATIVE_INFINITY;
	}

	/* Step 3:  Perform main alignment operations */
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
	align1InfoPtr->startMatch = -1;
	align2InfoPtr->startMatch = -1;
	align1InfoPtr->widthMatch = 0;
	align2InfoPtr->widthMatch = 0;
	int lookupValue, element1, element2, iElt, jElt;
	float substitutionValue, maxScore = NEGATIVE_INFINITY;
	for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
		SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iElt]);
		element1 = lookupValue;
		for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
			element2 = lookupValue;
			if (stringElements1.elts[iElt] == stringElements2.elts[jElt])
				substitutionValue = (float) matchMatrix[matrixDim[0] * element1 + element2];
			else
				substitutionValue = (float) mismatchMatrix[matrixDim[0] * element1 + element2];

			/* Step 3a:  Generate substitution, insertion, and deletion scores */
			F_MATRIX(i, j) =
				SAFE_SUM(MAX(F_MATRIX(iMinus1, jMinus1), MAX(H_MATRIX(iMinus1, jMinus1), V_MATRIX(iMinus1, jMinus1))),
				         substitutionValue);
			H_MATRIX(i, j) = 
				MAX(SAFE_SUM(MAX(F_MATRIX(iMinus1, j), V_MATRIX(iMinus1, j)), gapOpening + gapExtension),
				    SAFE_SUM(H_MATRIX(iMinus1, j), gapExtension));
			V_MATRIX(i, j) =
				MAX(SAFE_SUM(MAX(F_MATRIX(i, jMinus1), H_MATRIX(i, jMinus1)), gapOpening + gapExtension),
				    SAFE_SUM(V_MATRIX(i, jMinus1), gapExtension));

			if (typeCode == LOCAL_ALIGNMENT) {
				F_MATRIX(i, j) = MAX(F_MATRIX(i, j), 0.0);

				/* Step 3b:  Generate profile scores for local alignments */
				align1InfoPtr->profile[iElt] = MAX(F_MATRIX(i, j), align1InfoPtr->profile[iElt]);
				align2InfoPtr->profile[jElt] = MAX(F_MATRIX(i, j), align2InfoPtr->profile[jElt]);

				/* Step 3c:  Get the optimal score for local alignments */
				if (F_MATRIX(i, j) > maxScore) {
					align1InfoPtr->startMatch = iElt + 1;
					align2InfoPtr->startMatch = jElt + 1;
					maxScore = F_MATRIX(i, j);
				}
			}
		}
	}

	if (typeCode != LOCAL_ALIGNMENT) {
		/* Step 3d:  Adjust insertion and deletion scores when no end gaps */
		if (!endGap1) {
			for (i = 2, iMinus1 = 1; i <= nCharString1; i++, iMinus1++) {
				H_MATRIX(i, nCharString2) =
					MAX(F_MATRIX(iMinus1, nCharString2),
					MAX(H_MATRIX(iMinus1, nCharString2),
						V_MATRIX(iMinus1, nCharString2)));
			}
		}

		if (!endGap2) {
			for (j = 2, jMinus1 = 1; j <= nCharString2; j++, jMinus1++) {
				V_MATRIX(nCharString1, j) =
					MAX(F_MATRIX(nCharString1, jMinus1),
					MAX(H_MATRIX(nCharString1, jMinus1),
						V_MATRIX(nCharString1, jMinus1)));
			}
		}

		/* Step 3e:  Get the optimal score for non-local alignments */
		align1InfoPtr->startMatch = 1;
		align2InfoPtr->startMatch = 1;
		maxScore =
			MAX(F_MATRIX(nCharString1, nCharString2),
			MAX(H_MATRIX(nCharString1, nCharString2),
			    V_MATRIX(nCharString1, nCharString2)));
	}

	if (scoreOnly == 0) {
		if (typeCode != LOCAL_ALIGNMENT) {
			/* Step 3f:  Generate profile scores for non-local alignments */
			if (endGap1) {
				align1InfoPtr->profile[0] =
					MAX(F_MATRIX(nCharString1, nCharString2),
						V_MATRIX(nCharString1, nCharString2));
				for (i = 1, iElt = nCharString1Minus1; i < nCharString1; i++, iElt--) {
					float fvScore = MAX(F_MATRIX(i, nCharString2), V_MATRIX(i, nCharString2));
					align1InfoPtr->profile[iElt] =
						SAFE_SUM(fvScore, gapOpening + iElt * gapExtension);
					if (fvScore >= H_MATRIX(i, nCharString2)) {
						align1InfoPtr->profile[iElt] =
							MAX(align1InfoPtr->profile[iElt],
								SAFE_SUM(H_MATRIX(i, nCharString2), iElt * gapExtension));
					}
				}
			} else {
				for (i = 1, iElt = nCharString1Minus1; i <= nCharString1; i++, iElt--) {
					align1InfoPtr->profile[iElt] =
						MAX(F_MATRIX(i, nCharString2), V_MATRIX(i, nCharString2));
				}
			}
			if (MAX(F_MATRIX(nCharString1, nCharString2),
					V_MATRIX(nCharString1, nCharString2)) >=
						H_MATRIX(nCharString1, nCharString2)) {
				align1InfoPtr->profile[0] = maxScore;
			}

			if (endGap2) {
				align2InfoPtr->profile[0] =
					MAX(F_MATRIX(nCharString1, nCharString2),
						H_MATRIX(nCharString1, nCharString2));
				for (j = 1, jElt = nCharString2Minus1; j < nCharString2; j++, jElt--) {
					float fhScore = MAX(F_MATRIX(nCharString1, j), H_MATRIX(nCharString1, j));
					align2InfoPtr->profile[jElt] =
						SAFE_SUM(fhScore, gapOpening + jElt * gapExtension);
					if (fhScore >= V_MATRIX(nCharString1, j)) {
						align2InfoPtr->profile[jElt] =
							MAX(align2InfoPtr->profile[jElt],
								SAFE_SUM(V_MATRIX(nCharString1, j),
										 jElt * gapExtension));
					}
				}
			} else {
				for (j = 1, jElt = nCharString2Minus1; j <= nCharString2; j++, jElt--) {
					align2InfoPtr->profile[jElt] =
						MAX(F_MATRIX(nCharString1, j), H_MATRIX(nCharString1, j));
				}
			}
			if (MAX(F_MATRIX(nCharString1, nCharString2),
					H_MATRIX(nCharString1, nCharString2)) >=
						V_MATRIX(nCharString1, nCharString2)) {
				align2InfoPtr->profile[0] = maxScore;
			}

		}

		/* Step 4:  Traceback through the score matrix */
		char previousAction = '?';
		i = nCharString1Plus1 - align1InfoPtr->startMatch;
		j = nCharString2Plus1 - align2InfoPtr->startMatch;
		while ((i > 0 && j > 0) && !(typeCode == LOCAL_ALIGNMENT && F_MATRIX(i, j) == 0)) {
			float scoreSubstitution = F_MATRIX(i, j);
			float scoreInsertion    = V_MATRIX(i, j);
			float scoreDeletion     = H_MATRIX(i, j);

			char action;
			if (scoreSubstitution >= MAX(scoreInsertion, scoreDeletion))
				action = SUBSTITUTION;
			else if (scoreInsertion >= scoreDeletion)
				action = INSERTION;
			else
				action = DELETION;
			switch (action) {
		    	case DELETION:
		    		if (j == nCharString2) {
		    			align1InfoPtr->startMatch++;
		    		} else {
		    			align1InfoPtr->widthMatch++;
			    		if (previousAction != DELETION) {
							align2InfoPtr->startInserts--;
							align2InfoPtr->widthInserts--;
							align2InfoPtr->lengthInserts++;
							*align2InfoPtr->startInserts = nCharString2Plus1 - j;
			    		}
		    		}
		    		*align2InfoPtr->widthInserts += 1;
		    		i--;
		    		break;
		    	case INSERTION:
		    		if (i == nCharString1) {
		    			align2InfoPtr->startMatch++;
		    		} else {
		    			align2InfoPtr->widthMatch++;
			    		if (previousAction != INSERTION) {
							align1InfoPtr->startInserts--;
							align1InfoPtr->widthInserts--;
							align1InfoPtr->lengthInserts++;
							*align1InfoPtr->startInserts = nCharString1Plus1 - i;
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

		align1InfoPtr->profile[align1InfoPtr->startMatch - 1] = maxScore;
		align2InfoPtr->profile[align2InfoPtr->startMatch - 1] = maxScore;

		int offset1 = align1InfoPtr->startMatch - 1;
		if (offset1 > 0 && align1InfoPtr->lengthInserts > 0) {
			for (i = 0; i < align1InfoPtr->lengthInserts; i++)
				align1InfoPtr->startInserts[i] -= offset1;
		}
		int offset2 = align2InfoPtr->startMatch - 1;
		if (offset2 > 0 && align2InfoPtr->lengthInserts > 0) {
			for (j = 0; j < align2InfoPtr->lengthInserts; j++)
				align2InfoPtr->startInserts[j] -= offset2;
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
 *                            1 = 'global',   2 = 'local',    3 = 'overlap',
 *                            4 = 'overlap1', 5 = 'overlap2')
 * 'scoreOnly':              denotes whether or not to only return the scores
 *                           of the optimal pairwise alignment
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
	/* Get the strings for alignment */
	RoSeq stringElements1 = _get_XString_asRoSeq(string1);
	RoSeq stringElements2 = _get_XString_asRoSeq(string2);
	RoSeq qualityElements1 = _get_XString_asRoSeq(quality1);
	RoSeq qualityElements2 = _get_XString_asRoSeq(quality2);
	int nCharString1 = stringElements1.nelt;
	int nCharString2 = stringElements2.nelt;

	/* Allocate memory for temporary buffers */
	struct AlignInfo align1Info, align2Info;
	int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;
	int* startInserts1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* startInserts2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* widthInserts1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	int* widthInserts2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
	align1Info.profile = (float *) R_alloc((long) nCharString1, sizeof(float));
	align2Info.profile = (float *) R_alloc((long) nCharString2, sizeof(float));

	/* Initialize inserts values */
	for(int i = 0; i < alignmentBufferSize; i++) {
		startInserts1Buffer[i] = 0;
		startInserts2Buffer[i] = 0;
		widthInserts1Buffer[i] = 0;
		widthInserts2Buffer[i] = 0;
	}
	align1Info.lengthInserts = 0;
	align2Info.lengthInserts = 0;
	align1Info.startInserts = startInserts1Buffer + alignmentBufferSize;
	align2Info.startInserts = startInserts2Buffer + alignmentBufferSize;
	align1Info.widthInserts = widthInserts1Buffer + alignmentBufferSize;
	align2Info.widthInserts = widthInserts2Buffer + alignmentBufferSize;

	/* Perform pairwise alignment */
	float score = pairwiseAlignment(
			stringElements1,
			stringElements2,
			qualityElements1,
			qualityElements2,
			INTEGER(typeCode)[0],
			LOGICAL(scoreOnly)[0],
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

	/* Prepare the output */
	SEXP output, outputNames, outputElement1, outputElement2;
	PROTECT(output = NEW_LIST(11));
	/* Set the names */
	PROTECT(outputNames = NEW_CHARACTER(11));
	SET_STRING_ELT(outputNames, 0, mkChar("startMatch1"));
	SET_STRING_ELT(outputNames, 1, mkChar("widthMatch1"));
	SET_STRING_ELT(outputNames, 2, mkChar("startInserts1"));
	SET_STRING_ELT(outputNames, 3, mkChar("widthInserts1"));
	SET_STRING_ELT(outputNames, 4, mkChar("profile1"));
	SET_STRING_ELT(outputNames, 5, mkChar("startMatch2"));
	SET_STRING_ELT(outputNames, 6, mkChar("widthMatch2"));
	SET_STRING_ELT(outputNames, 7, mkChar("startInserts2"));
	SET_STRING_ELT(outputNames, 8, mkChar("widthInserts2"));
	SET_STRING_ELT(outputNames, 9, mkChar("profile2"));
	SET_STRING_ELT(outputNames, 10, mkChar("score"));
	SET_NAMES(output, outputNames);
	UNPROTECT(1);
	/* Set the "startMatch1" element */
	PROTECT(outputElement1 = NEW_INTEGER(1));
	INTEGER(outputElement1)[0] = align1Info.startMatch;
	SET_ELEMENT(output, 0, outputElement1);
	UNPROTECT(1);
	/* Set the "widthMatch1" element */
	PROTECT(outputElement1 = NEW_INTEGER(1));
	INTEGER(outputElement1)[0] = align1Info.widthMatch;
	SET_ELEMENT(output, 1, outputElement1);
	UNPROTECT(1);
	/* Set the "startInserts1" and "widthInserts1" elements */
	PROTECT(outputElement1 = NEW_INTEGER(align1Info.lengthInserts));
	PROTECT(outputElement2 = NEW_INTEGER(align1Info.lengthInserts));
	for(int i = 0; i < align1Info.lengthInserts; i++) {
		INTEGER(outputElement1)[i] = align1Info.startInserts[align1Info.lengthInserts - 1 - i];
		INTEGER(outputElement2)[i] = align1Info.widthInserts[align1Info.lengthInserts - 1 - i];
	}
	SET_ELEMENT(output, 2, outputElement1);
	SET_ELEMENT(output, 3, outputElement2);
	UNPROTECT(2);
	/* Set the "profile1" element */
	PROTECT(outputElement1 = NEW_NUMERIC(nCharString1));
	for(int i = 0; i < nCharString1; i++) {
		REAL(outputElement1)[i] = align1Info.profile[i];
	}
	SET_ELEMENT(output, 4, outputElement1);
	UNPROTECT(1);
	/* Set the "startMatch2" element */
	PROTECT(outputElement1 = NEW_INTEGER(1));
	INTEGER(outputElement1)[0] = align2Info.startMatch;
	SET_ELEMENT(output, 5, outputElement1);
	UNPROTECT(1);
	/* Set the "widthMatch2" element */
	PROTECT(outputElement1 = NEW_INTEGER(1));
	INTEGER(outputElement1)[0] = align2Info.widthMatch;
	SET_ELEMENT(output, 6, outputElement1);
	UNPROTECT(1);
	/* Set the "startInserts2" and "widthInserts2" elements */
	PROTECT(outputElement1 = NEW_INTEGER(align2Info.lengthInserts));
	PROTECT(outputElement2 = NEW_INTEGER(align2Info.lengthInserts));
	for(int j = 0; j < align2Info.lengthInserts; j++) {
		INTEGER(outputElement1)[j] = align2Info.startInserts[align2Info.lengthInserts - 1 - j];
		INTEGER(outputElement2)[j] = align2Info.widthInserts[align2Info.lengthInserts - 1 - j];
	}
	SET_ELEMENT(output, 7, outputElement1);
	SET_ELEMENT(output, 8, outputElement2);
	UNPROTECT(2);
	/* Set the "profile2" element */
	PROTECT(outputElement1 = NEW_NUMERIC(nCharString2));
	for(int j = 0; j < nCharString2; j++) {
		REAL(outputElement1)[j] = align2Info.profile[j];
	}
	SET_ELEMENT(output, 9, outputElement1);
	UNPROTECT(1);
	/* Set the "score" element */
	PROTECT(outputElement1 = NEW_NUMERIC(1));
	REAL(outputElement1)[0] = score;
	SET_ELEMENT(output, 10, outputElement1);
	UNPROTECT(1);
	/* Output is ready */
	UNPROTECT(1);

	return output;
}
