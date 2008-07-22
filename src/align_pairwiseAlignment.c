#include "Biostrings.h"

#include <float.h>
#include <stdlib.h>

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

#define NEGATIVE_INFINITY R_NegInf

#define          GLOBAL_ALIGNMENT 1
#define           LOCAL_ALIGNMENT 2
#define         OVERLAP_ALIGNMENT 3
#define PATTERN_OVERLAP_ALIGNMENT 4
#define SUBJECT_OVERLAP_ALIGNMENT 5

#define SUBSTITUTION 'S'
#define DELETION     'D'
#define INSERTION    'I'
#define TERMINATION  'T'

#define CURR_MATRIX(i, j) (currMatrix[nCharString2Plus1 * i + j])
#define PREV_MATRIX(i, j) (prevMatrix[nCharString2Plus1 * i + j])
#define S_TRACE_MATRIX(i, j) (sTraceMatrix[nCharString2 * i + j])
#define D_TRACE_MATRIX(i, j) (dTraceMatrix[nCharString2 * i + j])
#define I_TRACE_MATRIX(i, j) (iTraceMatrix[nCharString2 * i + j])

#define SET_LOOKUP_VALUE(lookupTable, length, key) \
{ \
	unsigned char lookupKey = (unsigned char) (key); \
	if (lookupKey >= (length) || (lookupValue = (lookupTable)[lookupKey]) == NA_INTEGER) { \
		error("key %d not in lookup table", (int) lookupKey); \
	} \
}

#define CHECK_MALLOC(ptr, name) \
{ \
	if (ptr == NULL) { \
		error("pairwiseAlignment():  could not allocate memory for %s", name); \
	} \
}


/* Structure to hold alignment information */
struct AlignInfo {
	RoSeq string;
	RoSeq quality;
	int endGap;
	int startRange;
	int widthRange;
	int* mismatch;
	int lengthMismatch;
	int* startIndel;
	int* widthIndel;
	int lengthIndel;
};
void function(struct AlignInfo *);


/* Returns the score of the optimal pairwise alignment */
static double pairwiseAlignment(
		struct AlignInfo *align1InfoPtr,
		struct AlignInfo *align2InfoPtr,
		int localAlignment,
		int scoreOnly,
		float gapOpening,
		float gapExtension,
		int useQuality,
		const int *qualityLookupTable,
		int qualityLookupTableLength,
		const double *qualityMatchMatrix,
		const double *qualityMismatchMatrix,
		const int *qualityMatrixDim,
		const int *constantLookupTable,
		int constantLookupTableLength,
		const double *constantMatrix,
		const int *constantMatrixDim)
{
	int i, j, iMinus1, jMinus1;

	/* Step 0:  Make sure smaller strings are second for memory efficiencies */
	int swappedOrder = align1InfoPtr->string.nelt < align2InfoPtr->string.nelt;
	if (swappedOrder) {
		struct AlignInfo *tempInfoPtr = align1InfoPtr;
		align1InfoPtr = align2InfoPtr;
		align2InfoPtr = tempInfoPtr;
	}

	/* Step 1:  Get information on input XString objects */
	int nCharString1 = align1InfoPtr->string.nelt;
	int nCharString2 = align2InfoPtr->string.nelt;
	int nCharString2Plus1 = nCharString2 + 1;
	int nCharString1Minus1 = nCharString1 - 1;
	int nCharString2Minus1 = nCharString2 - 1;

	if (nCharString1 < 1 || nCharString2 < 1)
		return NA_REAL;

	/* Step 2:  Create objects for scores values */
	/* Rows of currMatrix and prevMatrix = (0) substitution, (1) deletion, and (2) insertion */
	float *currMatrix, *prevMatrix, *curr, *currMinus1, *prev, *prevMinus1;
	if (scoreOnly && gapOpening == 0.0) {
		currMatrix = (float *) malloc(nCharString2Plus1 * sizeof(float));
		CHECK_MALLOC(currMatrix, "currMatrix");
		prevMatrix = (float *) malloc(nCharString2Plus1 * sizeof(float));
		CHECK_MALLOC(prevMatrix, "prevMatrix");

		if (align2InfoPtr->endGap) {
			for (j = 0, curr = currMatrix; j <= nCharString2; j++, curr++)
				*curr = j * gapExtension;
		} else {
			for (j = 0, curr = currMatrix; j <= nCharString2; j++, curr++)
				*curr = 0.0;
		}
	} else {
		currMatrix = (float *) malloc((3 * nCharString2Plus1) * sizeof(float));
		CHECK_MALLOC(currMatrix, "currMatrix");
		prevMatrix = (float *) malloc((3 * nCharString2Plus1) * sizeof(float));
		CHECK_MALLOC(prevMatrix, "prevMatrix");

		CURR_MATRIX(0, 0) = 0.0;
		CURR_MATRIX(1, 0) = (align1InfoPtr->endGap ? gapOpening : 0.0);
		for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
			CURR_MATRIX(0, j) = NEGATIVE_INFINITY;
			CURR_MATRIX(1, j) = NEGATIVE_INFINITY;
		}
		if (align2InfoPtr->endGap) {
			for (j = 0; j <= nCharString2; j++)
				CURR_MATRIX(2, j) = gapOpening + j * gapExtension;
		} else {
			for (j = 0; j <= nCharString2; j++)
				CURR_MATRIX(2, j) = 0.0;
		}
	}

	/* Step 3:  Perform main alignment operations */
	RoSeq sequence1, sequence2;
	int scalar1, scalar2;
	const int *lookupTable;
	int lookupTableLength;
	const double *matchMatrix, *mismatchMatrix;
	const int *matrixDim;
	if (useQuality) {
		sequence1 = align1InfoPtr->quality;
		sequence2 = align2InfoPtr->quality;
		scalar1 = (align1InfoPtr->quality.nelt == 1);
		scalar2 = (align2InfoPtr->quality.nelt == 1);
		lookupTable = qualityLookupTable;
		lookupTableLength = qualityLookupTableLength;
		matchMatrix = qualityMatchMatrix;
		mismatchMatrix = qualityMismatchMatrix;
		matrixDim = qualityMatrixDim;
	} else {
		sequence1 = align1InfoPtr->string;
		sequence2 = align2InfoPtr->string;
		scalar1 = (nCharString1 == 1);
		scalar2 = (nCharString2 == 1);
		lookupTable = constantLookupTable;
		lookupTableLength = constantLookupTableLength;
		matchMatrix = constantMatrix;
		mismatchMatrix = constantMatrix;
		matrixDim = constantMatrixDim;
	}
	int lookupValue = 0, elements[2], iElt, jElt;
	int notSwappedOrder = 1 - swappedOrder;
	int noEndGap1 = !align1InfoPtr->endGap;
	int noEndGap2 = !align2InfoPtr->endGap;
	float gapOpeningPlusExtension = gapOpening + gapExtension;
	float endGapAddend = (align1InfoPtr->endGap ? gapExtension : 0.0);
	float *tempMatrix, substitutionValue;
	double maxScore = NEGATIVE_INFINITY;
	if (scoreOnly) {
		/* Simplified calculations when only need the alignment score */
		if (gapOpening == 0.0) {
			for (i = 1, iElt = nCharString1Minus1; i <= nCharString1; i++, iElt--) {
				tempMatrix = prevMatrix;
				prevMatrix = currMatrix;
				currMatrix = tempMatrix;

				CURR_MATRIX(0, 0) = PREV_MATRIX(0, 0) + endGapAddend;

				SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iElt]);
				elements[swappedOrder] = lookupValue;
				if (localAlignment) {
					for (j = 1, jElt = nCharString2Minus1,
							curr = (currMatrix+1), currMinus1 = currMatrix,
							prev = (prevMatrix+1), prevMinus1 = prevMatrix;
							j <= nCharString2; j++, jElt--, curr++, currMinus1++, prev++, prevMinus1++) {
						SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
						elements[notSwappedOrder] = lookupValue;
						if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
							substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
						else
							substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

						*curr =
							MAX(0.0,
							MAX(*prevMinus1 + substitutionValue,
							MAX(*prev, *currMinus1) + gapExtension));

						maxScore = MAX(*curr, maxScore);
					}
				} else {
					for (j = 1, jElt = nCharString2Minus1,
							curr = (currMatrix+1), currMinus1 = currMatrix,
							prev = (prevMatrix+1), prevMinus1 = prevMatrix;
							j <= nCharString2; j++, jElt--, curr++, currMinus1++, prev++, prevMinus1++) {
						SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
						elements[notSwappedOrder] = lookupValue;
						if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
							substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
						else
							substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

						*curr =
							MAX(*prevMinus1 + substitutionValue,
							MAX(*prev, *currMinus1) + gapExtension);
					}
					if (noEndGap1) {
						currMatrix[nCharString2] =
							MAX(currMatrix[nCharString2],
							MAX(prevMatrix[nCharString2], currMatrix[nCharString2Minus1]));
					}
					if (noEndGap2 && i == nCharString1) {
						for (j = 1, curr = (currMatrix+1), currMinus1 = currMatrix, prev = (prevMatrix+1);
						     j <= nCharString2; j++, curr++, currMinus1++, prev++) {
							*curr = MAX(*curr, MAX(*prev, *currMinus1));
						}
					}
				}
			}

			if (!localAlignment) {
				maxScore = currMatrix[nCharString2];
			}
		} else {
			for (i = 1, iElt = nCharString1Minus1; i <= nCharString1; i++, iElt--) {
				tempMatrix = prevMatrix;
				prevMatrix = currMatrix;
				currMatrix = tempMatrix;

				CURR_MATRIX(0, 0) = NEGATIVE_INFINITY;
				CURR_MATRIX(1, 0) = PREV_MATRIX(1, 0) + endGapAddend;
				CURR_MATRIX(2, 0) = NEGATIVE_INFINITY;

				SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iElt]);
				elements[swappedOrder] = lookupValue;
				if (localAlignment) {
					for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
						SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
						elements[notSwappedOrder] = lookupValue;
						if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
							substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
						else
							substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

						CURR_MATRIX(0, j) =
							MAX(0.0,
								MAX(PREV_MATRIX(0, jMinus1), MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) + substitutionValue);
						CURR_MATRIX(1, j) =
							MAX(MAX(PREV_MATRIX(0, j), PREV_MATRIX(2, j)) + gapOpeningPlusExtension,
							    PREV_MATRIX(1, j) + gapExtension);
						CURR_MATRIX(2, j) =
							MAX(MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)) + gapOpeningPlusExtension,
							    CURR_MATRIX(2, jMinus1) + gapExtension);

						maxScore = MAX(CURR_MATRIX(0, j), maxScore);
					}
				} else {
					for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
						SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
						elements[notSwappedOrder] = lookupValue;
						if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
							substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
						else
							substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

						CURR_MATRIX(0, j) =
							MAX(PREV_MATRIX(0, jMinus1), MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) + substitutionValue;
						CURR_MATRIX(1, j) =
							MAX(MAX(PREV_MATRIX(0, j), PREV_MATRIX(2, j)) + gapOpeningPlusExtension,
							    PREV_MATRIX(1, j) + gapExtension);
						CURR_MATRIX(2, j) =
							MAX(MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)) + gapOpeningPlusExtension,
							    CURR_MATRIX(2, jMinus1) + gapExtension);
					}
					if (noEndGap1) {
						CURR_MATRIX(1, nCharString2) =
							MAX(PREV_MATRIX(0, nCharString2), MAX(PREV_MATRIX(1, nCharString2), PREV_MATRIX(2, nCharString2)));
					}
					if (noEndGap2 && i == nCharString1) {
						for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
							CURR_MATRIX(2, j) =
								MAX(MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)), CURR_MATRIX(2, jMinus1));
						}
					}
				}
			}

			if (!localAlignment) {
				maxScore =
					MAX(CURR_MATRIX(0, nCharString2),
					MAX(CURR_MATRIX(1, nCharString2),
					    CURR_MATRIX(2, nCharString2)));
			}
		}
	} else {
		/* Step 3a:  Create objects for traceback values */
		char *sTraceMatrix = (char *) malloc((nCharString1 * nCharString2) * sizeof(char));
		CHECK_MALLOC(sTraceMatrix, "sTraceMatrix");
		char *iTraceMatrix = (char *) malloc((nCharString1 * nCharString2) * sizeof(char));
		CHECK_MALLOC(iTraceMatrix, "iTraceMatrix");
		char *dTraceMatrix = (char *) malloc((nCharString1 * nCharString2) * sizeof(char));
		CHECK_MALLOC(dTraceMatrix, "dTraceMatrix");

		/* Step 3b:  Prepare the alignment info object for alignment */
		int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;
		align1InfoPtr->mismatch   = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align1InfoPtr->mismatch, "align1InfoPtr->mismatch");
		align2InfoPtr->mismatch   = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align2InfoPtr->mismatch, "align2InfoPtr->mismatch");
		align1InfoPtr->startIndel = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align1InfoPtr->startIndel, "align1InfoPtr->startIndel");
		align2InfoPtr->startIndel = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align2InfoPtr->startIndel, "align2InfoPtr->startIndel");
		align1InfoPtr->widthIndel = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align1InfoPtr->widthIndel, "align1InfoPtr->widthIndel");
		align2InfoPtr->widthIndel = (int *) malloc(alignmentBufferSize * sizeof(int));
		CHECK_MALLOC(align2InfoPtr->widthIndel, "align2InfoPtr->widthIndel");

		align1InfoPtr->lengthMismatch = 0;
		align2InfoPtr->lengthMismatch = 0;
		align1InfoPtr->lengthIndel = 0;
		align2InfoPtr->lengthIndel = 0;

		memset(align1InfoPtr->mismatch,   0, alignmentBufferSize * sizeof(int));
		memset(align2InfoPtr->mismatch,   0, alignmentBufferSize * sizeof(int));
		memset(align1InfoPtr->startIndel, 0, alignmentBufferSize * sizeof(int));
		memset(align2InfoPtr->startIndel, 0, alignmentBufferSize * sizeof(int));
		memset(align1InfoPtr->widthIndel, 0, alignmentBufferSize * sizeof(int));
		memset(align2InfoPtr->widthIndel, 0, alignmentBufferSize * sizeof(int));

		align1InfoPtr->startRange = -1;
		align2InfoPtr->startRange = -1;
		align1InfoPtr->widthRange = 0;
		align2InfoPtr->widthRange = 0;

		for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
			tempMatrix = prevMatrix;
			prevMatrix = currMatrix;
			currMatrix = tempMatrix;

			CURR_MATRIX(0, 0) = NEGATIVE_INFINITY;
			CURR_MATRIX(1, 0) = PREV_MATRIX(1, 0) + endGapAddend;
			CURR_MATRIX(2, 0) = NEGATIVE_INFINITY;

			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iElt]);
			elements[swappedOrder] = lookupValue;
			if (localAlignment) {
				for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
					SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
					elements[notSwappedOrder] = lookupValue;
					if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
						substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
					else
						substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

					/* Step 3c:  Generate (0) substitution, (1) deletion, and (2) insertion scores 
					 *           and traceback values 
					 */
					if (PREV_MATRIX(0, jMinus1) >= MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(0, j) = PREV_MATRIX(0, jMinus1) + substitutionValue;
					} else if (PREV_MATRIX(2, jMinus1) >= PREV_MATRIX(1, jMinus1)) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(0, j) = PREV_MATRIX(2, jMinus1) + substitutionValue;
					} else {
						S_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(0, j) = PREV_MATRIX(1, jMinus1) + substitutionValue;
					}
					if (PREV_MATRIX(1, j) >= (MAX(PREV_MATRIX(0, j), PREV_MATRIX(2, j)) + gapOpening)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(1, j) = PREV_MATRIX(1, j) + gapExtension;
					} else if (PREV_MATRIX(0, j) >= PREV_MATRIX(2, j)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(1, j) = PREV_MATRIX(0, j) + gapOpeningPlusExtension;
					} else {
						D_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(1, j) = PREV_MATRIX(2, j) + gapOpeningPlusExtension;
					}
					if (CURR_MATRIX(2, jMinus1) >= (MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)) + gapOpening)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(2, jMinus1) + gapExtension;
					} else if (CURR_MATRIX(0, jMinus1) >= CURR_MATRIX(1, jMinus1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(0, jMinus1) + gapOpeningPlusExtension;
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(2, j) = CURR_MATRIX(1, jMinus1) + gapOpeningPlusExtension;
					}

					CURR_MATRIX(0, j) = MAX(0.0, CURR_MATRIX(0, j));
					if (CURR_MATRIX(0, j) == 0.0) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;
						D_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;
						I_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;
					}

					/* Step 3d:  Get the optimal score for local alignments */
					if (CURR_MATRIX(0, j) > maxScore) {
						align1InfoPtr->startRange = iElt + 1;
						align2InfoPtr->startRange = jElt + 1;
						maxScore = CURR_MATRIX(0, j);
					}
				}
			} else {
				for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
					SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
					elements[notSwappedOrder] = lookupValue;
					if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
						substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
					else
						substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

					/* Step 3c:  Generate (0) substitution, (1) deletion, and (2) insertion scores 
					 *           and traceback values 
					 */
					if (PREV_MATRIX(0, jMinus1) >= MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(0, j) = PREV_MATRIX(0, jMinus1) + substitutionValue;
					} else if (PREV_MATRIX(2, jMinus1) >= PREV_MATRIX(1, jMinus1)) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(0, j) = PREV_MATRIX(2, jMinus1) + substitutionValue;
					} else {
						S_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(0, j) = PREV_MATRIX(1, jMinus1) + substitutionValue;
					}
					if (PREV_MATRIX(1, j) >= (MAX(PREV_MATRIX(0, j), PREV_MATRIX(2, j)) + gapOpening)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(1, j) = PREV_MATRIX(1, j) + gapExtension;
					} else if (PREV_MATRIX(0, j) >= PREV_MATRIX(2, j)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(1, j) = PREV_MATRIX(0, j) + gapOpeningPlusExtension;
					} else {
						D_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(1, j) = PREV_MATRIX(2, j) + gapOpeningPlusExtension;
					}
					if (CURR_MATRIX(2, jMinus1) >= (MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)) + gapOpening)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(2, jMinus1) + gapExtension;
					} else if (CURR_MATRIX(0, jMinus1) >= CURR_MATRIX(1, jMinus1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(0, jMinus1) + gapOpeningPlusExtension;
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(2, j) = CURR_MATRIX(1, jMinus1) + gapOpeningPlusExtension;
					}
				}
			}

			if (noEndGap1) {
				if (PREV_MATRIX(1, nCharString2) >= MAX(PREV_MATRIX(0, nCharString2), PREV_MATRIX(2, nCharString2))) {
					D_TRACE_MATRIX(iMinus1, nCharString2Minus1) = DELETION;
					CURR_MATRIX(1, nCharString2) = PREV_MATRIX(1, nCharString2);
				} else if (PREV_MATRIX(0, nCharString2) >= PREV_MATRIX(2, nCharString2)) {
					D_TRACE_MATRIX(iMinus1, nCharString2Minus1) = SUBSTITUTION;
					CURR_MATRIX(1, nCharString2) = PREV_MATRIX(0, nCharString2);
				} else {
					D_TRACE_MATRIX(iMinus1, nCharString2Minus1) = INSERTION;
					CURR_MATRIX(1, nCharString2) = PREV_MATRIX(2, nCharString2);
				}
			}
			if (noEndGap2 && i == nCharString1) {
				for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++) {
					if (CURR_MATRIX(2, jMinus1) >= MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1))) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(2, jMinus1);
					} else if (CURR_MATRIX(0, jMinus1) >= CURR_MATRIX(1, jMinus1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(2, j) = CURR_MATRIX(0, jMinus1);
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(2, j) = CURR_MATRIX(1, jMinus1);
					}
				}
			}
		}

		char currTraceMatrix = '?';
		if (localAlignment) {
			if (maxScore == 0.0)
				currTraceMatrix = TERMINATION;
			else
				currTraceMatrix = SUBSTITUTION;
		} else {
			/* Step 3g:  Get the optimal score for non-local alignments */
			align1InfoPtr->startRange = 1;
			align2InfoPtr->startRange = 1;
			if (CURR_MATRIX(0, nCharString2) >=
					MAX(CURR_MATRIX(1, nCharString2), CURR_MATRIX(2, nCharString2))) {
				currTraceMatrix = SUBSTITUTION;
				maxScore = CURR_MATRIX(0, nCharString2);
			} else if (CURR_MATRIX(2, nCharString2) >= CURR_MATRIX(1, nCharString2)) {
				currTraceMatrix = INSERTION;
				maxScore = CURR_MATRIX(2, nCharString2);
			} else {
				currTraceMatrix = DELETION;
				maxScore = CURR_MATRIX(1, nCharString2);
			}
		}

		/* Step 4:  Traceback through the score matrices */
		i = nCharString1 - align1InfoPtr->startRange;
		j = nCharString2 - align2InfoPtr->startRange;
		char prevTraceMatrix = '?';
		while (currTraceMatrix != TERMINATION && i >= 0 && j >= 0) {
			switch (currTraceMatrix) {
	    		case DELETION:
	    			if (D_TRACE_MATRIX(i, j) != TERMINATION) {
		    			if (j == nCharString2Minus1) {
		    				align1InfoPtr->startRange++;
		    			} else {
		    				align1InfoPtr->widthRange++;
		    				if (prevTraceMatrix != DELETION) {
		    					align2InfoPtr->startIndel[align2InfoPtr->lengthIndel] = nCharString2 - j;
		    					align2InfoPtr->lengthIndel++;
		    				}
		    				align2InfoPtr->widthIndel[align2InfoPtr->lengthIndel - 1] += 1;
		    			}
	    			}
	    			prevTraceMatrix = currTraceMatrix;
					currTraceMatrix = D_TRACE_MATRIX(i, j);
	    			i--;
	    			break;
	    		case INSERTION:
	    			if (I_TRACE_MATRIX(i, j) != TERMINATION) {
		    			if (i == nCharString1Minus1) {
		    				align2InfoPtr->startRange++;
		    			} else {
		    				align2InfoPtr->widthRange++;
		    				if (prevTraceMatrix != INSERTION) {
		    					align1InfoPtr->startIndel[align1InfoPtr->lengthIndel] = nCharString1 - i;
		    					align1InfoPtr->lengthIndel++;
		    				}
		    				align1InfoPtr->widthIndel[align1InfoPtr->lengthIndel - 1] += 1;
		    			}
	    			}
	    			prevTraceMatrix = currTraceMatrix;
					currTraceMatrix = I_TRACE_MATRIX(i, j);
	    			j--;
	    			break;
	    		case SUBSTITUTION:
	    			if (S_TRACE_MATRIX(i, j) != TERMINATION) {
		    			align1InfoPtr->widthRange++;
		    			align2InfoPtr->widthRange++;
	    			}
	    			prevTraceMatrix = currTraceMatrix;
					currTraceMatrix = S_TRACE_MATRIX(i, j);
					if (align1InfoPtr->string.elts[nCharString1Minus1 - i] !=
						align2InfoPtr->string.elts[nCharString2Minus1 - j]) {
						align1InfoPtr->mismatch[align1InfoPtr->lengthMismatch] = nCharString1 - i;
						align2InfoPtr->mismatch[align2InfoPtr->lengthMismatch] = nCharString2 - j;
						align1InfoPtr->lengthMismatch++;
						align2InfoPtr->lengthMismatch++;
					}
					i--;
	    			j--;
	    			break;
	    		default:
	    			error("unknown traceback code %d", currTraceMatrix);
	    			break;
			}
		}

		int offset1 = align1InfoPtr->startRange - 1;
		if (offset1 > 0 && align1InfoPtr->lengthIndel > 0) {
			for (i = 0; i < align1InfoPtr->lengthIndel; i++)
				align1InfoPtr->startIndel[i] -= offset1;
		}
		int offset2 = align2InfoPtr->startRange - 1;
		if (offset2 > 0 && align2InfoPtr->lengthIndel > 0) {
			for (j = 0; j < align2InfoPtr->lengthIndel; j++)
				align2InfoPtr->startIndel[j] -= offset2;
		}

		free(sTraceMatrix);
		free(iTraceMatrix);
		free(dTraceMatrix);
	}

	free(currMatrix);
	free(prevMatrix);

	return (double) maxScore;
}

/*
 * INPUTS
 * 'pattern':                XString object for patterns
 * 'subject':                XString object for subject
 * 'patternQuality':         BString object for quality scores for pattern
 * 'subjectQuality':         BString object for quality scores for subject
 * 'type':                   type of pairwise alignment
 *                           (character vector of length 1;
 *                            'global', 'local', 'overlap',
 *                            'patternOverlap', 'subjectOverlap')
 * 'typeCode':               type of pairwise alignment
 *                           (integer vector of length 1;
 *                            1 = 'global', 2 = 'local', 3 = 'overlap',
 *                            4 = 'patternOverlap', 5 = 'subjectOverlap')
 * 'qualityType':            type of quality scores
 *                           (character vector of length 1)
 * 'scoreOnly':              denotes whether or not to only return the scores
 *                           of the optimal pairwise alignment
 *                           (logical vector of length 1)
 * 'gapOpening':             gap opening cost or penalty
 *                           (double vector of length 1)
 * 'gapExtension':           gap extension cost or penalty
 *                           (double vector of length 1)
 * 'useQuality':             denotes whether or not to use quality measures
 *                           in the optimal pairwise alignment
 *                           (logical vector of length 1)
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

SEXP XStringSet_align_pairwiseAlignment(
		SEXP pattern,
		SEXP subject,
		SEXP patternQuality,
		SEXP subjectQuality,
		SEXP type,
		SEXP typeCode,
		SEXP qualityType,
		SEXP scoreOnly,
		SEXP gapOpening,
		SEXP gapExtension,
		SEXP useQuality,
		SEXP qualityLookupTable,
		SEXP qualityMatchMatrix,
		SEXP qualityMismatchMatrix,
		SEXP qualityMatrixDim,
		SEXP constantLookupTable,
		SEXP constantMatrix,
		SEXP constantMatrixDim)
{
	int scoreOnlyValue = LOGICAL(scoreOnly)[0];
	int useQualityValue = LOGICAL(useQuality)[0];
	float gapOpeningValue = REAL(gapOpening)[0];
	float gapExtensionValue = REAL(gapExtension)[0];
	if (gapOpeningValue == NEGATIVE_INFINITY || gapExtensionValue == NEGATIVE_INFINITY) {
		gapOpeningValue = 0.0;
		gapExtensionValue = NEGATIVE_INFINITY;
	}
	int localAlignment = (INTEGER(typeCode)[0] == LOCAL_ALIGNMENT);

	/* Create the alignment info objects */
	struct AlignInfo align1Info, align2Info;
	align2Info.string = _get_XString_asRoSeq(subject);
	align2Info.quality = _get_XString_asRoSeq(subjectQuality);
	align1Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == SUBJECT_OVERLAP_ALIGNMENT);
	align2Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == PATTERN_OVERLAP_ALIGNMENT);

	int numberOfStrings = _get_XStringSet_length(pattern);
	CachedXStringSet cachedPattern = _new_CachedXStringSet(pattern);
	CachedXStringSet cachedPatternQuality = _new_CachedXStringSet(patternQuality);

	SEXP output;

	int i, qualityElement = 0;
	int qualityIncrement = ((_get_XStringSet_length(patternQuality) < numberOfStrings) ? 0 : 1);
	double *score;
	if (scoreOnlyValue) {
		PROTECT(output = NEW_NUMERIC(numberOfStrings));
		for (i = 0, score = REAL(output); i < numberOfStrings; i++, score++) {
			align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
			align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, qualityElement);
			*score = pairwiseAlignment(
					&align1Info,
					&align2Info,
					localAlignment,
					scoreOnlyValue,
					gapOpeningValue,
					gapExtensionValue,
					useQualityValue,
					INTEGER(qualityLookupTable),
					LENGTH(qualityLookupTable),
					REAL(qualityMatchMatrix),
					REAL(qualityMismatchMatrix),
					INTEGER(qualityMatrixDim),
					INTEGER(constantLookupTable),
					LENGTH(constantLookupTable),
					REAL(constantMatrix),
					INTEGER(constantMatrixDim));
			qualityElement += qualityIncrement;
		}
		UNPROTECT(1);
	} else {
		SEXP alignedPattern;
		SEXP alignedPatternRange, alignedPatternRangeStart, alignedPatternRangeWidth;
		SEXP alignedPatternMismatch, alignedPatternMismatchElt;
		SEXP alignedPatternIndel, alignedPatternIndelRange;
		SEXP alignedPatternIndelRangeStart, alignedPatternIndelRangeWidth;
		SEXP alignedSubject;
		SEXP alignedSubjectRange, alignedSubjectRangeStart, alignedSubjectRangeWidth;
		SEXP alignedSubjectMismatch, alignedSubjectMismatchElt;
		SEXP alignedSubjectIndel, alignedSubjectIndelRange;
		SEXP alignedSubjectIndelRangeStart, alignedSubjectIndelRangeWidth;
		SEXP alignedScore;

		PROTECT(output = NEW_OBJECT(MAKE_CLASS("PairwiseAlignment")));

		/* Set the "pattern" slot */
		if (useQualityValue) {
			PROTECT(alignedPattern = NEW_OBJECT(MAKE_CLASS("QualityAlignedXStringSet")));
			SET_SLOT(alignedPattern, mkChar("quality"), patternQuality);
			SET_SLOT(alignedPattern, mkChar("qualityType"), qualityType);
		} else {
			PROTECT(alignedPattern = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		}
		SET_SLOT(alignedPattern, mkChar("unaligned"), pattern);
		/* Set the "range" sub-slot */
		PROTECT(alignedPatternRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
		PROTECT(alignedPatternRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternRangeWidth = NEW_INTEGER(numberOfStrings));
		SET_SLOT(alignedPatternRange, mkChar("start"), alignedPatternRangeStart);
		SET_SLOT(alignedPatternRange, mkChar("width"), alignedPatternRangeWidth);
		SET_SLOT(alignedPattern, mkChar("range"), alignedPatternRange);
		/* Set the "mismatch" sub-slot */
		PROTECT(alignedPatternMismatch = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedPattern, mkChar("mismatch"), alignedPatternMismatch);
		/* Set the "indel" sub-slot */
		PROTECT(alignedPatternIndel = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedPattern, mkChar("indel"), alignedPatternIndel);
		SET_SLOT(output, mkChar("pattern"), alignedPattern);

		/* Set the "subject" slot */
		if (useQualityValue) {
			PROTECT(alignedSubject = NEW_OBJECT(MAKE_CLASS("QualityAlignedXStringSet")));
			SET_SLOT(alignedSubject, mkChar("quality"), subjectQuality);
			SET_SLOT(alignedSubject, mkChar("qualityType"), qualityType);
		} else {
			PROTECT(alignedSubject = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		}
		SET_SLOT(alignedSubject, mkChar("unaligned"), subject);
		/* Set the "range" sub-slot */
		PROTECT(alignedSubjectRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
		PROTECT(alignedSubjectRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectRangeWidth = NEW_INTEGER(numberOfStrings));
		SET_SLOT(alignedSubjectRange, mkChar("start"), alignedSubjectRangeStart);
		SET_SLOT(alignedSubjectRange, mkChar("width"), alignedSubjectRangeWidth);
		SET_SLOT(alignedSubject, mkChar("range"), alignedSubjectRange);
		/* Set the "mismatch" sub-slot */
		PROTECT(alignedSubjectMismatch = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedSubject, mkChar("mismatch"), alignedSubjectMismatch);
		/* Set the "indel" sub-slot */
		PROTECT(alignedSubjectIndel = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedSubject, mkChar("indel"), alignedSubjectIndel);
		SET_SLOT(output, mkChar("subject"), alignedSubject);

		/* Set the "score" slot */
		PROTECT(alignedScore = NEW_NUMERIC(numberOfStrings));
		SET_SLOT(output, mkChar("score"), alignedScore);

		/* Set the "type" slot */
		SET_SLOT(output, mkChar("type"), type);
		/* Set the "constantMatrix" slot */
		SET_SLOT(output, mkChar("constantMatrix"), constantMatrix);
		/* Set the "gapOpening" slot */
		SET_SLOT(output, mkChar("gapOpening"), gapOpening);
		/* Set the "gapExtension" slot */
		SET_SLOT(output, mkChar("gapExtension"), gapExtension);

		int *align1RangeStart, *align1RangeWidth, *align2RangeStart, *align2RangeWidth;
		for (i = 0, score = REAL(alignedScore),
				align1RangeStart = INTEGER(alignedPatternRangeStart),
				align1RangeWidth = INTEGER(alignedPatternRangeWidth),
				align2RangeStart = INTEGER(alignedSubjectRangeStart),
				align2RangeWidth = INTEGER(alignedSubjectRangeWidth); i < numberOfStrings;
				i++, score++, align1RangeStart++, align1RangeWidth++, align2RangeStart++, align2RangeWidth++) {
			align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
			align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, qualityElement);
			*score = pairwiseAlignment(
					&align1Info,
					&align2Info,
					localAlignment,
					scoreOnlyValue,
					gapOpeningValue,
					gapExtensionValue,
					useQualityValue,
					INTEGER(qualityLookupTable),
					LENGTH(qualityLookupTable),
					REAL(qualityMatchMatrix),
					REAL(qualityMismatchMatrix),
					INTEGER(qualityMatrixDim),
					INTEGER(constantLookupTable),
					LENGTH(constantLookupTable),
					REAL(constantMatrix),
					INTEGER(constantMatrixDim));

			PROTECT(alignedPatternMismatchElt = NEW_INTEGER(align1Info.lengthMismatch));
			PROTECT(alignedSubjectMismatchElt = NEW_INTEGER(align2Info.lengthMismatch));
			memcpy(INTEGER(alignedPatternMismatchElt), align1Info.mismatch,
			       align1Info.lengthMismatch * sizeof(int));
			memcpy(INTEGER(alignedSubjectMismatchElt), align2Info.mismatch,
			       align2Info.lengthMismatch * sizeof(int));
		    SET_VECTOR_ELT(alignedPatternMismatch, i, alignedPatternMismatchElt);
		    SET_VECTOR_ELT(alignedSubjectMismatch, i, alignedSubjectMismatchElt);
			UNPROTECT(2);

			*align1RangeStart = align1Info.startRange;
			*align1RangeWidth = align1Info.widthRange;
			PROTECT(alignedPatternIndelRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
			PROTECT(alignedPatternIndelRangeStart = NEW_INTEGER(align1Info.lengthIndel));
			PROTECT(alignedPatternIndelRangeWidth = NEW_INTEGER(align1Info.lengthIndel));
			memcpy(INTEGER(alignedPatternIndelRangeStart), align1Info.startIndel,
			       align1Info.lengthIndel * sizeof(int));
			memcpy(INTEGER(alignedPatternIndelRangeWidth), align1Info.widthIndel,
			       align1Info.lengthIndel * sizeof(int));
			SET_SLOT(alignedPatternIndelRange, mkChar("start"), alignedPatternIndelRangeStart);
			SET_SLOT(alignedPatternIndelRange, mkChar("width"), alignedPatternIndelRangeWidth);
		    SET_VECTOR_ELT(alignedPatternIndel, i, alignedPatternIndelRange);
		    UNPROTECT(3);

			*align2RangeStart = align2Info.startRange;
			*align2RangeWidth = align2Info.widthRange;
			PROTECT(alignedSubjectIndelRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
			PROTECT(alignedSubjectIndelRangeStart = NEW_INTEGER(align2Info.lengthIndel));
			PROTECT(alignedSubjectIndelRangeWidth = NEW_INTEGER(align2Info.lengthIndel));
			memcpy(INTEGER(alignedSubjectIndelRangeStart), align2Info.startIndel,
			       align2Info.lengthIndel * sizeof(int));
			memcpy(INTEGER(alignedSubjectIndelRangeWidth), align2Info.widthIndel,
			       align2Info.lengthIndel * sizeof(int));
			SET_SLOT(alignedSubjectIndelRange, mkChar("start"), alignedSubjectIndelRangeStart);
			SET_SLOT(alignedSubjectIndelRange, mkChar("width"), alignedSubjectIndelRangeWidth);
		    SET_VECTOR_ELT(alignedSubjectIndel, i, alignedSubjectIndelRange);
		    UNPROTECT(3);

			qualityElement += qualityIncrement;

			/* Free the memory allocated in pairwiseAlignment */
			free(align1Info.mismatch);
			free(align2Info.mismatch);
			free(align1Info.startIndel);
			free(align2Info.startIndel);
			free(align1Info.widthIndel);
			free(align2Info.widthIndel);
		}

		/* Output is ready */
		UNPROTECT(14);
	}

	return output;
}
