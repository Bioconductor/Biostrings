#include "Biostrings.h"
#include "IRanges_interface.h"
#include <R_ext/Utils.h>        /* R_CheckUserInterrupt */

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
void function1(struct AlignInfo *);


/* Structure to hold alignment buffers */
struct AlignBuffer {
	float *currMatrix;
	float *prevMatrix;
	char *sTraceMatrix;
	char *iTraceMatrix;
	char *dTraceMatrix;
};
void function2(struct AlignBuffer *);


/* Returns the score of the optimal pairwise alignment */
static double pairwiseAlignment(
		struct AlignInfo *align1InfoPtr,
		struct AlignInfo *align2InfoPtr,
		const int localAlignment,
		const int scoreOnly,
		const float gapOpening,
		const float gapExtension,
		const int useQuality,
		const int *qualityLookupTable,
		const int qualityLookupTableLength,
		const double *qualityMatchMatrix,
		const double *qualityMismatchMatrix,
		const int *qualityMatrixDim,
		const int *constantLookupTable,
		const int constantLookupTableLength,
		const double *constantMatrix,
		const int *constantMatrixDim,
		struct AlignBuffer *alignBufferPtr)
{
	int i, j, iMinus1, jMinus1;

	/* Step 0:  Make sure smaller strings are second for memory efficiencies */
	const int swappedOrder = align1InfoPtr->string.nelt < align2InfoPtr->string.nelt;
	if (swappedOrder) {
		struct AlignInfo *tempInfoPtr = align1InfoPtr;
		align1InfoPtr = align2InfoPtr;
		align2InfoPtr = tempInfoPtr;
	}

	/* Step 1:  Get information on input XString objects */
	const int nCharString1 = align1InfoPtr->string.nelt;
	const int nCharString2 = align2InfoPtr->string.nelt;
	const int nCharString2Plus1 = nCharString2 + 1;
	const int nCharString1Minus1 = nCharString1 - 1;
	const int nCharString2Minus1 = nCharString2 - 1;

	if (nCharString1 < 1 || nCharString2 < 1) {
		double zeroCharScore;
		if (nCharString2 >= 1 && align2InfoPtr->endGap)
			zeroCharScore = gapOpening + nCharString2 * gapExtension;
		else if (nCharString1 >= 1 && align1InfoPtr->endGap)
			zeroCharScore = gapOpening + nCharString1 * gapExtension;
		else
			zeroCharScore = 0.0;
		return zeroCharScore;
	}

	/* Step 2:  Create objects for scores values */
	/* Rows of currMatrix and prevMatrix = (0) substitution, (1) deletion, and (2) insertion */
	float *currMatrix = alignBufferPtr->currMatrix;
	float *prevMatrix = alignBufferPtr->prevMatrix;
	float *curr, *currMinus1, *prev, *prevMinus1;
	if (scoreOnly && gapOpening == 0.0) {
		if (align2InfoPtr->endGap) {
			for (j = 0, curr = currMatrix; j <= nCharString2; j++, curr++)
				*curr = j * gapExtension;
		} else {
			for (j = 0, curr = currMatrix; j <= nCharString2; j++, curr++)
				*curr = 0.0;
		}
	} else {
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
	const int notSwappedOrder = 1 - swappedOrder;
	const int noEndGap1 = !align1InfoPtr->endGap;
	const int noEndGap2 = !align2InfoPtr->endGap;
	const float gapOpeningPlusExtension = gapOpening + gapExtension;
	const float endGapAddend = (align1InfoPtr->endGap ? gapExtension : 0.0);
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
								MAX(PREV_MATRIX(0, jMinus1),
								MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) + substitutionValue);
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
							MAX(PREV_MATRIX(0, jMinus1),
							MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))) + substitutionValue;
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
		char *sTraceMatrix = alignBufferPtr->sTraceMatrix;
		char *iTraceMatrix = alignBufferPtr->iTraceMatrix;
		char *dTraceMatrix = alignBufferPtr->dTraceMatrix;

		/* Step 3b:  Prepare the alignment info object for alignment */
		const int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;

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

		const int offset1 = align1InfoPtr->startRange - 1;
		if (offset1 > 0 && align1InfoPtr->lengthIndel > 0) {
			for (i = 0; i < align1InfoPtr->lengthIndel; i++)
				align1InfoPtr->startIndel[i] -= offset1;
		}
		const int offset2 = align2InfoPtr->startRange - 1;
		if (offset2 > 0 && align2InfoPtr->lengthIndel > 0) {
			for (j = 0; j < align2InfoPtr->lengthIndel; j++)
				align2InfoPtr->startIndel[j] -= offset2;
		}
	}

	return (double) maxScore;
}

/*
 * INPUTS
 * 'pattern':                XStringSet object for patterns
 * 'subject':                XString object for subject
 * 'patternQuality':         BStringSet object for quality scores for pattern
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
 * If scoreOnly = TRUE, returns either a vector of scores
 * If scoreOnly = FALSE, returns an S4 PairwiseAlignment object.
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
	const int scoreOnlyValue = LOGICAL(scoreOnly)[0];
	const int useQualityValue = LOGICAL(useQuality)[0];
	const int localAlignment = (INTEGER(typeCode)[0] == LOCAL_ALIGNMENT);
	float gapOpeningValue = REAL(gapOpening)[0];
	float gapExtensionValue = REAL(gapExtension)[0];
	if (gapOpeningValue == NEGATIVE_INFINITY || gapExtensionValue == NEGATIVE_INFINITY) {
		gapOpeningValue = 0.0;
		gapExtensionValue = NEGATIVE_INFINITY;
	}

	/* Create the alignment info objects */
	struct AlignInfo align1Info, align2Info;
	align2Info.string = _get_XString_asRoSeq(subject);
	align2Info.quality = _get_XString_asRoSeq(subjectQuality);
	align1Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == SUBJECT_OVERLAP_ALIGNMENT);
	align2Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == PATTERN_OVERLAP_ALIGNMENT);

	const int numberOfStrings = _get_XStringSet_length(pattern);
	CachedXStringSet cachedPattern = _new_CachedXStringSet(pattern);
	CachedXStringSet cachedPatternQuality = _new_CachedXStringSet(patternQuality);

	SEXP output;

	int i, qualityElement = 0;
	const int qualityIncrement = ((_get_XStringSet_length(patternQuality) < numberOfStrings) ? 0 : 1);

	/* Create the alignment buffer object */
	struct AlignBuffer alignBuffer;
	int nCharString1 = 0;
	int nCharString2 = _get_XString_asRoSeq(subject).nelt;
	for (i = 0; i < numberOfStrings; i++) {
		nCharString1 = MAX(nCharString1, _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i).nelt);
	}
	const int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;
	if (scoreOnlyValue && gapOpeningValue == 0.0) {
		alignBuffer.currMatrix = (float *) R_alloc((long) alignmentBufferSize, sizeof(float));
		alignBuffer.prevMatrix = (float *) R_alloc((long) alignmentBufferSize, sizeof(float));
	} else {
		alignBuffer.currMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
		alignBuffer.prevMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
	}
	if (!scoreOnlyValue) {
		align1Info.mismatch   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.mismatch   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align1Info.startIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.startIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align1Info.widthIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.widthIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		alignBuffer.sTraceMatrix = (char *) R_alloc((long) (nCharString1 * nCharString2), sizeof(char));
		alignBuffer.iTraceMatrix = (char *) R_alloc((long) (nCharString1 * nCharString2), sizeof(char));
		alignBuffer.dTraceMatrix = (char *) R_alloc((long) (nCharString1 * nCharString2), sizeof(char));
	}

	double *score;
	if (scoreOnlyValue) {
		PROTECT(output = NEW_NUMERIC(numberOfStrings));
		for (i = 0, score = REAL(output); i < numberOfStrings; i++, score++) {
	        R_CheckUserInterrupt();
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
					INTEGER(constantMatrixDim),
					&alignBuffer);
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

		PROTECT(alignedPatternRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternRangeWidth = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternMismatch = NEW_LIST(numberOfStrings));
		PROTECT(alignedPatternIndel = NEW_LIST(numberOfStrings));

		PROTECT(alignedSubjectRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectRangeWidth = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectMismatch = NEW_LIST(numberOfStrings));
		PROTECT(alignedSubjectIndel = NEW_LIST(numberOfStrings));

		PROTECT(alignedScore = NEW_NUMERIC(numberOfStrings));

		int *align1RangeStart, *align1RangeWidth, *align2RangeStart, *align2RangeWidth;
		for (i = 0, score = REAL(alignedScore),
				align1RangeStart = INTEGER(alignedPatternRangeStart),
				align1RangeWidth = INTEGER(alignedPatternRangeWidth),
				align2RangeStart = INTEGER(alignedSubjectRangeStart),
				align2RangeWidth = INTEGER(alignedSubjectRangeWidth); i < numberOfStrings;
				i++, score++, align1RangeStart++, align1RangeWidth++, align2RangeStart++, align2RangeWidth++) {
	        R_CheckUserInterrupt();
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
					INTEGER(constantMatrixDim),
					&alignBuffer);

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
			PROTECT(alignedPatternIndelRangeStart = NEW_INTEGER(align1Info.lengthIndel));
			PROTECT(alignedPatternIndelRangeWidth = NEW_INTEGER(align1Info.lengthIndel));
			memcpy(INTEGER(alignedPatternIndelRangeStart), align1Info.startIndel,
			       align1Info.lengthIndel * sizeof(int));
			memcpy(INTEGER(alignedPatternIndelRangeWidth), align1Info.widthIndel,
			       align1Info.lengthIndel * sizeof(int));
			PROTECT(alignedPatternIndelRange = 
				new_IRanges("IRanges", alignedPatternIndelRangeStart, alignedPatternIndelRangeWidth, R_NilValue));
		    SET_VECTOR_ELT(alignedPatternIndel, i, alignedPatternIndelRange);
		    UNPROTECT(3);

			*align2RangeStart = align2Info.startRange;
			*align2RangeWidth = align2Info.widthRange;
			PROTECT(alignedSubjectIndelRangeStart = NEW_INTEGER(align2Info.lengthIndel));
			PROTECT(alignedSubjectIndelRangeWidth = NEW_INTEGER(align2Info.lengthIndel));
			memcpy(INTEGER(alignedSubjectIndelRangeStart), align2Info.startIndel,
			       align2Info.lengthIndel * sizeof(int));
			memcpy(INTEGER(alignedSubjectIndelRangeWidth), align2Info.widthIndel,
			       align2Info.lengthIndel * sizeof(int));
			PROTECT(alignedSubjectIndelRange = 
				new_IRanges("IRanges", alignedSubjectIndelRangeStart, alignedSubjectIndelRangeWidth, R_NilValue));
		    SET_VECTOR_ELT(alignedSubjectIndel, i, alignedSubjectIndelRange);
		    UNPROTECT(3);

			qualityElement += qualityIncrement;
		}

		/* Create the output object */
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
		PROTECT(alignedPatternRange =
			new_IRanges("IRanges", alignedPatternRangeStart, alignedPatternRangeWidth, R_NilValue));
		SET_SLOT(alignedPattern, mkChar("range"), alignedPatternRange);
		/* Set the "mismatch" sub-slot */
		SET_SLOT(alignedPattern, mkChar("mismatch"), alignedPatternMismatch);
		/* Set the "indel" sub-slot */
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
		PROTECT(alignedSubjectRange =
			new_IRanges("IRanges", alignedSubjectRangeStart, alignedSubjectRangeWidth, R_NilValue));
		SET_SLOT(alignedSubject, mkChar("range"), alignedSubjectRange);
		/* Set the "mismatch" sub-slot */
		SET_SLOT(alignedSubject, mkChar("mismatch"), alignedSubjectMismatch);
		/* Set the "indel" sub-slot */
		SET_SLOT(alignedSubject, mkChar("indel"), alignedSubjectIndel);
		SET_SLOT(output, mkChar("subject"), alignedSubject);

		/* Set the "score" slot */
		SET_SLOT(output, mkChar("score"), alignedScore);

		/* Set the "type" slot */
		SET_SLOT(output, mkChar("type"), type);
		/* Set the "constantMatrix" slot */
		SET_SLOT(output, mkChar("constantMatrix"), constantMatrix);
		/* Set the "gapOpening" slot */
		SET_SLOT(output, mkChar("gapOpening"), gapOpening);
		/* Set the "gapExtension" slot */
		SET_SLOT(output, mkChar("gapExtension"), gapExtension);

		/* Output is ready */
		UNPROTECT(14);
	}

	return output;
}



/*
 * INPUTS
 * 'string':                 XStringSet object for strings
 * 'stringQuality':          BStringSet object for quality scores for strings
 * 'type':                   type of pairwise alignment
 *                           (character vector of length 1;
 *                            'global', 'local', 'overlap')
 * 'typeCode':               type of pairwise alignment
 *                           (integer vector of length 1;
 *                            1 = 'global', 2 = 'local', 3 = 'overlap')
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
 * Return a numeric vector containing the lower triangle of the score matrix.
 */

SEXP XStringSet_align_distance(
		SEXP string,
		SEXP stringQuality,
		SEXP type,
		SEXP typeCode,
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
	int scoreOnlyValue = 1;
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
	align1Info.endGap = (INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT);
	align2Info.endGap = (INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT);

	int numberOfStrings = _get_XStringSet_length(string);
	CachedXStringSet cachedPattern = _new_CachedXStringSet(string);
	CachedXStringSet cachedPatternQuality = _new_CachedXStringSet(stringQuality);

	SEXP output;

	int i, j, iQualityElement = 0, jQualityElement = 0;
	int qualityIncrement = ((_get_XStringSet_length(stringQuality) < numberOfStrings) ? 0 : 1);

	/* Create the alignment buffer object */
	struct AlignBuffer alignBuffer;
	int nCharString = 0;
	for (i = 0; i < numberOfStrings; i++) {
		nCharString = MAX(nCharString, _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i).nelt);
	}
	int alignmentBufferSize = nCharString + 1;
	if (gapOpeningValue == 0.0) {
		alignBuffer.currMatrix = (float *) R_alloc((long) alignmentBufferSize, sizeof(float));
		alignBuffer.prevMatrix = (float *) R_alloc((long) alignmentBufferSize, sizeof(float));
	} else {
		alignBuffer.currMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
		alignBuffer.prevMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
	}

	double *score;
	PROTECT(output = NEW_NUMERIC((numberOfStrings * (numberOfStrings - 1)) / 2));
	score = REAL(output);
	for (i = 0; i < numberOfStrings; i++) {
        R_CheckUserInterrupt();
		align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
		align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, iQualityElement);
		jQualityElement = iQualityElement + qualityIncrement;
		for (j = i + 1; j < numberOfStrings; j++) {
			align2Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, j);
			align2Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, jQualityElement);
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
					INTEGER(constantMatrixDim),
					&alignBuffer);
			jQualityElement += qualityIncrement;
			score++;
		}
		iQualityElement += qualityIncrement;
	}
	UNPROTECT(1);

	return output;
}

