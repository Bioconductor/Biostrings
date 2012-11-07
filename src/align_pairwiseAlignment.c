#include "Biostrings.h"
#include "IRanges_interface.h"
#include <R_ext/Utils.h>        /* R_CheckUserInterrupt */

#include <float.h>
#include <stdlib.h>

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

#define NEGATIVE_INFINITY R_NegInf
#define MAX_BUF_SIZE      1048576

#define       GLOBAL_ALIGNMENT 1
#define        LOCAL_ALIGNMENT 2
#define      OVERLAP_ALIGNMENT 3
#define GLOBAL_LOCAL_ALIGNMENT 4
#define LOCAL_GLOBAL_ALIGNMENT 5

#define SUBSTITUTION 'S'
#define DELETION     'D'
#define INSERTION    'I'
#define TERMINATION  'T'

#define CURR_MATRIX(i, j) (currMatrix[i + nCharString1Plus1 * j])
#define PREV_MATRIX(i, j) (prevMatrix[i + nCharString1Plus1 * j])
#define S_TRACE_MATRIX(i, j) (sTraceMatrix[i + nCharString1 * j])
#define D_TRACE_MATRIX(i, j) (dTraceMatrix[i + nCharString1 * j])
#define I_TRACE_MATRIX(i, j) (iTraceMatrix[i + nCharString1 * j])
#define FUZZY_MATRIX(i, j) (fuzzyMatrix[i + fuzzyMatrixDim[0] * j])
#define SUBSTITUTION_ARRAY(i, j, k) (substitutionArray[i + substitutionArrayDim[0] * (j + substitutionArrayDim[1] * k)])

#define SET_LOOKUP_VALUE(lookupTable, length, key) \
{ \
	unsigned char lookupKey = (unsigned char) (key); \
	if (lookupKey >= (length) || (lookupValue = (lookupTable)[lookupKey]) == NA_INTEGER) { \
		error("key %d not in lookup table", (int) lookupKey); \
	} \
}


/* Structure to hold alignment information */
struct AlignInfo {
	cachedCharSeq string;
	cachedCharSeq quality;
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


/* Structure to hold mismatch buffers */
struct MismatchBuffer {
	int *pattern;
	int *subject;
	int usedSpace;
	int totalSpace;
};
void function3(struct MismatchBuffer *);


/* Structure to hold indel buffers */
struct IndelBuffer {
	int *start;
	int *width;
	int usedSpace;
	int totalSpace;
};
void function4(struct IndelBuffer *);


/* Returns the score of the optimal pairwise alignment */
static double pairwiseAlignment(
		struct AlignInfo *align1InfoPtr,
		struct AlignInfo *align2InfoPtr,
		const int localAlignment,
		const int scoreOnly,
		const float gapOpening,
		const float gapExtension,
		const int useQuality,
		const double *substitutionArray,
		const int *substitutionArrayDim,
		const int *substitutionLookupTable,
		const int substitutionLookupTableLength,
		const int *fuzzyMatrix,
		const int *fuzzyMatrixDim,
		const int *fuzzyLookupTable,
		const int fuzzyLookupTableLength,
		struct AlignBuffer *alignBufferPtr)
{
	int i, j, iMinus1, jMinus1;

	/* Step 1:  Get information on input XString objects */
	const int nCharString1 = align1InfoPtr->string.length;
	const int nCharString2 = align2InfoPtr->string.length;
	const int nCharString1Plus1 = nCharString1 + 1;
	const int nCharString1Minus1 = nCharString1 - 1;
	const int nCharString2Minus1 = nCharString2 - 1;

	if (nCharString1 < 1 || nCharString2 < 1) {
		double zeroCharScore;
		if (nCharString1 >= 1 && align1InfoPtr->endGap)
			zeroCharScore = gapOpening + nCharString1 * gapExtension;
		else if (nCharString2 >= 1 && align2InfoPtr->endGap)
			zeroCharScore = gapOpening + nCharString2 * gapExtension;
		else
			zeroCharScore = 0.0;
		return zeroCharScore;
	}

	/* Step 2:  Create objects for scores values */
	/* Rows of currMatrix and prevMatrix = (0) substitution, (1) deletion, and (2) insertion */
	float *currMatrix = alignBufferPtr->currMatrix;
	float *prevMatrix = alignBufferPtr->prevMatrix;
	CURR_MATRIX(0, 0) = 0.0;
	CURR_MATRIX(0, 1) = (align2InfoPtr->endGap ? gapOpening : 0.0);
	for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
		CURR_MATRIX(i, 0) = NEGATIVE_INFINITY;
		CURR_MATRIX(i, 1) = NEGATIVE_INFINITY;
	}
	if (align1InfoPtr->endGap) {
		for (i = 0; i <= nCharString1; i++)
			CURR_MATRIX(i, 2) = gapOpening + i * gapExtension;
	} else {
		for (i = 0; i <= nCharString1; i++)
			CURR_MATRIX(i, 2) = 0.0;
	}

	/* Step 3:  Perform main alignment operations */
	cachedCharSeq sequence1, sequence2;
	int scalar1, scalar2;
	if (useQuality) {
		sequence1 = align1InfoPtr->quality;
		sequence2 = align2InfoPtr->quality;
		scalar1 = (align1InfoPtr->quality.length == 1);
		scalar2 = (align2InfoPtr->quality.length == 1);
	} else {
		sequence1 = align1InfoPtr->string;
		sequence2 = align2InfoPtr->string;
		scalar1 = (nCharString1 == 1);
		scalar2 = (nCharString2 == 1);
	}
	int lookupValue = 0, element1, element2, stringElt1, stringElt2, fuzzy, iElt, jElt;
	const int noEndGap1 = !align1InfoPtr->endGap;
	const int noEndGap2 = !align2InfoPtr->endGap;
	const float gapOpeningPlusExtension = gapOpening + gapExtension;
	const float endGapAddend = (align2InfoPtr->endGap ? gapExtension : 0.0);
	float *tempMatrix, substitutionValue;
	double maxScore = NEGATIVE_INFINITY;
	if (scoreOnly) {
		/* Simplified calculations when only need the alignment score */
		for (j = 1, jElt = nCharString2Minus1; j <= nCharString2; j++, jElt--) {
			tempMatrix = prevMatrix;
			prevMatrix = currMatrix;
			currMatrix = tempMatrix;

			CURR_MATRIX(0, 0) = NEGATIVE_INFINITY;
			CURR_MATRIX(0, 1) = PREV_MATRIX(0, 1) + endGapAddend;
			CURR_MATRIX(0, 2) = NEGATIVE_INFINITY;

			SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align2InfoPtr->string.seq[jElt]);
			stringElt2 = lookupValue;
			SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence2.seq[scalar2 ? 0 : jElt]);
			element2 = lookupValue;
			if (localAlignment) {
				for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
					SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align1InfoPtr->string.seq[iElt]);
					stringElt1 = lookupValue;
					SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence1.seq[scalar1 ? 0 : iElt]);
					element1 = lookupValue;
					fuzzy = FUZZY_MATRIX(stringElt1, stringElt2);
					substitutionValue = (float) SUBSTITUTION_ARRAY(element1, element2, fuzzy);

					CURR_MATRIX(i, 0) =
						MAX(0.0,
							MAX(PREV_MATRIX(iMinus1, 0),
							MAX(PREV_MATRIX(iMinus1, 1), PREV_MATRIX(iMinus1, 2))) + substitutionValue);
					CURR_MATRIX(i, 1) =
						MAX(MAX(PREV_MATRIX(i, 0), PREV_MATRIX(i, 2)) + gapOpeningPlusExtension,
						    PREV_MATRIX(i, 1) + gapExtension);
					CURR_MATRIX(i, 2) =
						MAX(MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1)) + gapOpeningPlusExtension,
						    CURR_MATRIX(iMinus1, 2) + gapExtension);

					maxScore = MAX(CURR_MATRIX(i, 0), maxScore);
				}
			} else {
				for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
					SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align1InfoPtr->string.seq[iElt]);
					stringElt1 = lookupValue;
					SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence1.seq[scalar1 ? 0 : iElt]);
					element1 = lookupValue;
					fuzzy = FUZZY_MATRIX(stringElt1, stringElt2);
					substitutionValue = (float) SUBSTITUTION_ARRAY(element1, element2, fuzzy);

					CURR_MATRIX(i, 0) =
						MAX(PREV_MATRIX(iMinus1, 0),
						MAX(PREV_MATRIX(iMinus1, 1), PREV_MATRIX(iMinus1, 2))) + substitutionValue;
					CURR_MATRIX(i, 1) =
						MAX(MAX(PREV_MATRIX(i, 0), PREV_MATRIX(i, 2)) + gapOpeningPlusExtension,
						    PREV_MATRIX(i, 1) + gapExtension);
					CURR_MATRIX(i, 2) =
						MAX(MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1)) + gapOpeningPlusExtension,
						    CURR_MATRIX(iMinus1, 2) + gapExtension);
				}
				if (noEndGap2) {
					CURR_MATRIX(nCharString1, 1) =
						MAX(PREV_MATRIX(nCharString1, 0), MAX(PREV_MATRIX(nCharString1, 1), PREV_MATRIX(nCharString1, 2)));
				}
				if (noEndGap1 && j == nCharString2) {
					for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
						CURR_MATRIX(i, 2) =
							MAX(MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1)), CURR_MATRIX(iMinus1, 2));
					}
				}
			}
		}

		if (!localAlignment) {
			maxScore =
				MAX(CURR_MATRIX(nCharString1, 0),
				MAX(CURR_MATRIX(nCharString1, 1),
				    CURR_MATRIX(nCharString1, 2)));
		}
	} else {
		/* Step 3a:  Create objects for traceback values */
		char *sTraceMatrix = alignBufferPtr->sTraceMatrix;
		char *iTraceMatrix = alignBufferPtr->iTraceMatrix;
		char *dTraceMatrix = alignBufferPtr->dTraceMatrix;

		/* Step 3b:  Prepare the alignment info object for alignment */
		const int alignmentBufferSize = nCharString1Plus1;

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

		for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
			tempMatrix = prevMatrix;
			prevMatrix = currMatrix;
			currMatrix = tempMatrix;

			CURR_MATRIX(0, 0) = NEGATIVE_INFINITY;
			CURR_MATRIX(0, 1) = PREV_MATRIX(0, 1) + endGapAddend;
			CURR_MATRIX(0, 2) = NEGATIVE_INFINITY;

			SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align2InfoPtr->string.seq[jElt]);
			stringElt2 = lookupValue;
			SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence2.seq[scalar2 ? 0 : jElt]);
			element2 = lookupValue;
			if (localAlignment) {
				for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
					SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align1InfoPtr->string.seq[iElt]);
					stringElt1 = lookupValue;
					SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence1.seq[scalar1 ? 0 : iElt]);
					element1 = lookupValue;
					fuzzy = FUZZY_MATRIX(stringElt1, stringElt2);
					substitutionValue = (float) SUBSTITUTION_ARRAY(element1, element2, fuzzy);

					/* Step 3c:  Generate (0) substitution, (1) deletion, and (2) insertion scores
					 *           and traceback values
					 */
					if (PREV_MATRIX(iMinus1, 0) >= MAX(PREV_MATRIX(iMinus1, 1), PREV_MATRIX(iMinus1, 2))) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 0) + substitutionValue;
					} else if (PREV_MATRIX(iMinus1, 1) >= PREV_MATRIX(iMinus1, 2)) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 1) + substitutionValue;
					} else {
						S_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 2) + substitutionValue;
					}
					if (PREV_MATRIX(i, 1) > (MAX(PREV_MATRIX(i, 0), PREV_MATRIX(i, 2)) + gapOpening)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 1) + gapExtension;
					} else if (PREV_MATRIX(i, 0) >= PREV_MATRIX(i, 2)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 0) + gapOpeningPlusExtension;
					} else {
						D_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 2) + gapOpeningPlusExtension;
					}
					if (CURR_MATRIX(iMinus1, 2) > (MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1)) + gapOpening)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 2) + gapExtension;
					} else if (CURR_MATRIX(iMinus1, 0) >= CURR_MATRIX(iMinus1, 1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 0) + gapOpeningPlusExtension;
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 1) + gapOpeningPlusExtension;
					}

					CURR_MATRIX(i, 0) = MAX(0.0, CURR_MATRIX(i, 0));
					if (CURR_MATRIX(i, 0) == 0.0)
						S_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;
					CURR_MATRIX(i, 1) = MAX(0.0, CURR_MATRIX(i, 1));
					if (CURR_MATRIX(i, 1) == 0.0)
						D_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;
					CURR_MATRIX(i, 2) = MAX(0.0, CURR_MATRIX(i, 2));
					if (CURR_MATRIX(i, 2) == 0.0)
						I_TRACE_MATRIX(iMinus1, jMinus1) = TERMINATION;

					/* Step 3d:  Get the optimal score for local alignments */
					if (CURR_MATRIX(i, 0) >= maxScore) {
						align1InfoPtr->startRange = iElt + 1;
						align2InfoPtr->startRange = jElt + 1;
						maxScore = CURR_MATRIX(i, 0);
					}
				}
			} else {
				for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
					SET_LOOKUP_VALUE(fuzzyLookupTable, fuzzyLookupTableLength, align1InfoPtr->string.seq[iElt]);
					stringElt1 = lookupValue;
					SET_LOOKUP_VALUE(substitutionLookupTable, substitutionLookupTableLength, sequence1.seq[scalar1 ? 0 : iElt]);
					element1 = lookupValue;
					fuzzy = FUZZY_MATRIX(stringElt1, stringElt2);
					substitutionValue = (float) SUBSTITUTION_ARRAY(element1, element2, fuzzy);

					/* Step 3c:  Generate (0) substitution, (1) deletion, and (2) insertion scores
					 *           and traceback values
					 */
					if (PREV_MATRIX(iMinus1, 0) >= MAX(PREV_MATRIX(iMinus1, 1), PREV_MATRIX(iMinus1, 2))) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 0) + substitutionValue;
					} else if (PREV_MATRIX(iMinus1, 1) >= PREV_MATRIX(iMinus1, 2)) {
						S_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 1) + substitutionValue;
					} else {
						S_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 0) = PREV_MATRIX(iMinus1, 2) + substitutionValue;
					}
					if (PREV_MATRIX(i, 1) > (MAX(PREV_MATRIX(i, 0), PREV_MATRIX(i, 2)) + gapOpening)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 1) + gapExtension;
					} else if (PREV_MATRIX(i, 0) >= PREV_MATRIX(i, 2)) {
						D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 0) + gapOpeningPlusExtension;
					} else {
						D_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 1) = PREV_MATRIX(i, 2) + gapOpeningPlusExtension;
					}
					if (CURR_MATRIX(iMinus1, 2) > (MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1)) + gapOpening)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 2) + gapExtension;
					} else if (CURR_MATRIX(iMinus1, 0) >= CURR_MATRIX(iMinus1, 1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 0) + gapOpeningPlusExtension;
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 1) + gapOpeningPlusExtension;
					}
				}
			}

			if (noEndGap2) {
				if (PREV_MATRIX(nCharString1, 1) >= MAX(PREV_MATRIX(nCharString1, 0), PREV_MATRIX(nCharString1, 2))) {
					D_TRACE_MATRIX(nCharString1Minus1, jMinus1) = DELETION;
					CURR_MATRIX(nCharString1, 1) = PREV_MATRIX(nCharString1, 1);
				} else if (PREV_MATRIX(nCharString1, 0) >= PREV_MATRIX(nCharString1, 2)) {
					D_TRACE_MATRIX(nCharString1Minus1, jMinus1) = SUBSTITUTION;
					CURR_MATRIX(nCharString1, 1) = PREV_MATRIX(nCharString1, 0);
				} else {
					D_TRACE_MATRIX(nCharString1Minus1, jMinus1) = INSERTION;
					CURR_MATRIX(nCharString1, 1) = PREV_MATRIX(nCharString1, 2);
				}
			}
			if (noEndGap1 && j == nCharString2) {
				for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++) {
					if (CURR_MATRIX(iMinus1, 2) >= MAX(CURR_MATRIX(iMinus1, 0), CURR_MATRIX(iMinus1, 1))) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 2);
					} else if (CURR_MATRIX(iMinus1, 0) >= CURR_MATRIX(iMinus1, 1)) {
						I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 0);
					} else {
						I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
						CURR_MATRIX(i, 2) = CURR_MATRIX(iMinus1, 1);
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
			if (CURR_MATRIX(nCharString1, 0) >=
					MAX(CURR_MATRIX(nCharString1, 1), CURR_MATRIX(nCharString1, 2))) {
				currTraceMatrix = SUBSTITUTION;
				maxScore = CURR_MATRIX(nCharString1, 0);
			} else if (CURR_MATRIX(nCharString1, 1) >= CURR_MATRIX(nCharString1, 2)) {
				currTraceMatrix = DELETION;
				maxScore = CURR_MATRIX(nCharString1, 1);
			} else {
				currTraceMatrix = INSERTION;
				maxScore = CURR_MATRIX(nCharString1, 2);
			}
		}

		/* Step 4:  Traceback through the score matrices */
		i = nCharString1 - align1InfoPtr->startRange;
		j = nCharString2 - align2InfoPtr->startRange;
		char prevTraceMatrix = '?';
		while (currTraceMatrix != TERMINATION && i >= 0 && j >= 0) {
			switch (currTraceMatrix) {
	    		case INSERTION:
	    			if (I_TRACE_MATRIX(i, j) != TERMINATION) {
		    			if (j == nCharString2Minus1) {
		    				align1InfoPtr->startRange++;
		    			} else {
		    				align1InfoPtr->widthRange++;
		    				if (prevTraceMatrix != INSERTION) {
		    					align2InfoPtr->startIndel[align2InfoPtr->lengthIndel] = nCharString2 - j;
		    					align2InfoPtr->lengthIndel++;
		    				}
		    				align2InfoPtr->widthIndel[align2InfoPtr->lengthIndel - 1] += 1;
		    			}
	    			}
	    			prevTraceMatrix = currTraceMatrix;
					currTraceMatrix = I_TRACE_MATRIX(i, j);
	    			i--;
	    			break;
	    		case DELETION:
	    			if (D_TRACE_MATRIX(i, j) != TERMINATION) {
		    			if (i == nCharString1Minus1) {
		    				align2InfoPtr->startRange++;
		    			} else {
		    				align2InfoPtr->widthRange++;
		    				if (prevTraceMatrix != DELETION) {
		    					align1InfoPtr->startIndel[align1InfoPtr->lengthIndel] = nCharString1 - i;
		    					align1InfoPtr->lengthIndel++;
		    				}
		    				align1InfoPtr->widthIndel[align1InfoPtr->lengthIndel - 1] += 1;
		    			}
	    			}
	    			prevTraceMatrix = currTraceMatrix;
					currTraceMatrix = D_TRACE_MATRIX(i, j);
	    			j--;
	    			break;
	    		case SUBSTITUTION:
				prevTraceMatrix = currTraceMatrix;
				currTraceMatrix = S_TRACE_MATRIX(i, j);
	    			if (currTraceMatrix != TERMINATION) {
		    			align1InfoPtr->widthRange++;
		    			align2InfoPtr->widthRange++;
	    			}
					if (currTraceMatrix != TERMINATION && align1InfoPtr->string.seq[nCharString1Minus1 - i] !=
						align2InfoPtr->string.seq[nCharString2Minus1 - j]) {
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
 * 'pattern':                XStringSet or QualityScaledXStringSet object for patterns
 * 'subject':                XStringSet or QualityScaledXStringSet object for subject
 * 'type':                   type of pairwise alignment
 *                           (character vector of length 1;
 *                            'global', 'local', 'overlap', 'global-local',
 *                            'local-global')
 * 'typeCode':               type of pairwise alignment
 *                           (integer vector of length 1;
 *                            1 = 'global', 2 = 'local', 3 = 'overlap',
 *                            4 = 'global-local', 5 = 'local-global')
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
 * 'substitutionArray':        a three-dimensional double array where the first two
 *                             dimensions are for substitutions and the
 *                             third is for fuzziness of matches
 * 'substitutionArrayDim':     dimension of 'substitutionArray'
 *                             (integer vector of length 3)
 * 'substitutionLookupTable':  lookup table for translating XString bytes to
 *                             substitution indices
 *                             (integer vector)
 * 'fuzzyMatrix':              fuzzy matrix for matches
 *                             (double matrix)
 * 'fuzzyMatrixDim':           dimension of 'fuzzyMatrix'
 *                             (integer vector of length 2)
 * 'fuzzyLookupTable':         lookup table for translating XString bytes to
 *                             fuzzy indices
 *                             (integer vector)
 *
 * OUTPUT
 * If scoreOnly = TRUE, returns either a vector of scores
 * If scoreOnly = FALSE, returns an S4 PairwiseAlignments or
 * PairwiseAlignmentsSingleSubject object.
 */

SEXP XStringSet_align_pairwiseAlignment(
		SEXP pattern,
		SEXP subject,
		SEXP type,
		SEXP typeCode,
		SEXP scoreOnly,
		SEXP gapOpening,
		SEXP gapExtension,
		SEXP useQuality,
		SEXP substitutionArray,
		SEXP substitutionArrayDim,
		SEXP substitutionLookupTable,
		SEXP fuzzyMatrix,
		SEXP fuzzyMatrixDim,
		SEXP fuzzyLookupTable)
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

	cachedXStringSet cached_pattern = _cache_XStringSet(pattern);
	cachedXStringSet cached_subject = _cache_XStringSet(subject);
	const int numberOfStrings = _get_cachedXStringSet_length(&cached_pattern);
	const int multipleSubjects = _get_cachedXStringSet_length(&cached_subject) > 1;
	int lengthOfPatternQualitySet = 0;
	int lengthOfSubjectQualitySet = 0;

	SEXP patternQuality, subjectQuality;
	cachedXStringSet cached_patternQuality, cached_subjectQuality;
	if (useQualityValue) {
		patternQuality = GET_SLOT(pattern, install("quality"));
		subjectQuality = GET_SLOT(subject, install("quality"));
		cached_patternQuality = _cache_XStringSet(patternQuality);
		lengthOfPatternQualitySet = _get_XStringSet_length(patternQuality);
		cached_subjectQuality = _cache_XStringSet(subjectQuality);
		lengthOfSubjectQualitySet = _get_cachedXStringSet_length(&cached_subjectQuality);
	} else {
		patternQuality = R_NilValue;
		subjectQuality = R_NilValue;
	}

	/* Create the alignment info objects */
	struct AlignInfo align1Info, align2Info;
	align2Info.string = _get_cachedXStringSet_elt(&cached_subject, 0);
	if (useQualityValue)
		align2Info.quality = _get_cachedXStringSet_elt(&cached_subjectQuality, 0);
	align1Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == GLOBAL_LOCAL_ALIGNMENT);
	align2Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == LOCAL_GLOBAL_ALIGNMENT);

	SEXP output;

	int i, quality1Element = 0, quality2Element = 0;
	const int quality1Increment = ((lengthOfPatternQualitySet < numberOfStrings) ? 0 : 1);
	const int quality2Increment = ((lengthOfSubjectQualitySet < numberOfStrings) ? 0 : 1);

	/* Create the alignment buffer object */
	struct AlignBuffer alignBuffer;
	int nCharString1 = 0, nCharString2 = 0, nCharProduct = 0;
	if (multipleSubjects) {
		for (i = 0; i < numberOfStrings; i++) {
			int nchar1 = _get_cachedXStringSet_elt(&cached_pattern, i).length;
			int nchar2 = _get_cachedXStringSet_elt(&cached_subject, i).length;
			nCharString1 = MAX(nCharString1, nchar1);
			nCharString2 = MAX(nCharString2, nchar2);
			nCharProduct = MAX(nCharProduct, nchar1 * nchar2);
		}
	} else {
		for (i = 0; i < numberOfStrings; i++) {
			nCharString1 = MAX(nCharString1, _get_cachedXStringSet_elt(&cached_pattern, i).length);
		}
		nCharString2 = align2Info.string.length;
		nCharProduct = nCharString1 * nCharString2;
	}
	const int alignmentBufferSize = nCharString1 + 1;
	alignBuffer.currMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
	alignBuffer.prevMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));

	struct MismatchBuffer mismatchBuffer;
	struct IndelBuffer indel1Buffer;
	struct IndelBuffer indel2Buffer;
	int mismatchBufferSize = 0, indelBufferSize = 0;
	if (!scoreOnlyValue) {
		align1Info.mismatch   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.mismatch   = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align1Info.startIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.startIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align1Info.widthIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		align2Info.widthIndel = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		alignBuffer.sTraceMatrix = (char *) R_alloc((long) nCharProduct, sizeof(char));
		alignBuffer.iTraceMatrix = (char *) R_alloc((long) nCharProduct, sizeof(char));
		alignBuffer.dTraceMatrix = (char *) R_alloc((long) nCharProduct, sizeof(char));

		mismatchBufferSize = MIN(MAX_BUF_SIZE, alignmentBufferSize + numberOfStrings * (alignmentBufferSize/4));
		mismatchBuffer.pattern = (int *) R_alloc((long) mismatchBufferSize, sizeof(int));
		mismatchBuffer.subject = (int *) R_alloc((long) mismatchBufferSize, sizeof(int));
		mismatchBuffer.usedSpace = 0;
		mismatchBuffer.totalSpace = mismatchBufferSize;

		indelBufferSize = MIN(MAX_BUF_SIZE, alignmentBufferSize + numberOfStrings * (alignmentBufferSize/12));
		indel1Buffer.start = (int *) R_alloc((long) indelBufferSize, sizeof(int));
		indel1Buffer.width = (int *) R_alloc((long) indelBufferSize, sizeof(int));
		indel1Buffer.usedSpace = 0;
		indel1Buffer.totalSpace = indelBufferSize;
		indel2Buffer.start = (int *) R_alloc((long) indelBufferSize, sizeof(int));
		indel2Buffer.width = (int *) R_alloc((long) indelBufferSize, sizeof(int));
		indel2Buffer.usedSpace = 0;
		indel2Buffer.totalSpace = indelBufferSize;
	}

	double *score;
	if (scoreOnlyValue) {
		PROTECT(output = NEW_NUMERIC(numberOfStrings));
		for (i = 0, score = REAL(output); i < numberOfStrings; i++, score++) {
	        R_CheckUserInterrupt();
			align1Info.string = _get_cachedXStringSet_elt(&cached_pattern, i);
			if (useQualityValue) {
				align1Info.quality = _get_cachedXStringSet_elt(&cached_patternQuality, quality1Element);
				quality1Element += quality1Increment;
			}
			if (multipleSubjects) {
				align2Info.string = _get_cachedXStringSet_elt(&cached_subject, i);
				if (useQualityValue) {
					align2Info.quality = _get_cachedXStringSet_elt(&cached_subjectQuality, quality2Element);
					quality2Element += quality2Increment;
				}
			}
			*score = pairwiseAlignment(
					&align1Info,
					&align2Info,
					localAlignment,
					scoreOnlyValue,
					gapOpeningValue,
					gapExtensionValue,
					useQualityValue,
					REAL(substitutionArray),
					INTEGER(substitutionArrayDim),
					INTEGER(substitutionLookupTable),
					LENGTH(substitutionLookupTable),
					INTEGER(fuzzyMatrix),
					INTEGER(fuzzyMatrixDim),
					INTEGER(fuzzyLookupTable),
					LENGTH(fuzzyLookupTable),
					&alignBuffer);
		}
		UNPROTECT(1);
	} else {
		SEXP alignedPattern;
		SEXP alignedPatternRange, alignedPatternRangeStart, alignedPatternRangeWidth;
		SEXP alignedPatternMismatch;
		SEXP alignedPatternMismatchPartitioning, alignedPatternMismatchValues;
		SEXP alignedPatternMismatchEnds;
		SEXP alignedPatternIndel;
		SEXP alignedPatternIndelPartitioning, alignedPatternIndelRange;
		SEXP alignedPatternIndelRangeStart, alignedPatternIndelRangeWidth;
		SEXP alignedPatternIndelEnds;

		SEXP alignedSubject;
		SEXP alignedSubjectRange, alignedSubjectRangeStart, alignedSubjectRangeWidth;
		SEXP alignedSubjectMismatch;
		SEXP alignedSubjectMismatchPartitioning, alignedSubjectMismatchValues;
		SEXP alignedSubjectMismatchEnds;
		SEXP alignedSubjectIndel;
		SEXP alignedSubjectIndelPartitioning, alignedSubjectIndelRange;
		SEXP alignedSubjectIndelRangeStart, alignedSubjectIndelRangeWidth;
		SEXP alignedSubjectIndelEnds;

		SEXP alignedScore;

		PROTECT(alignedPatternRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternRangeWidth = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternMismatchEnds = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternIndelEnds = NEW_INTEGER(numberOfStrings));

		PROTECT(alignedSubjectRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectRangeWidth = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectMismatchEnds = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectIndelEnds = NEW_INTEGER(numberOfStrings));

		PROTECT(alignedScore = NEW_NUMERIC(numberOfStrings));

		int align1MismatchPrevEnd = 0, align1IndelPrevEnd = 0;
		int align2MismatchPrevEnd = 0, align2IndelPrevEnd = 0;
		int *tempIntPtr;
		int *align1RangeStart, *align1RangeWidth, *align1MismatchEnds, *align1IndelEnds;
		int *align2RangeStart, *align2RangeWidth, *align2MismatchEnds, *align2IndelEnds;
		for (i = 0, score = REAL(alignedScore),
				align1RangeStart = INTEGER(alignedPatternRangeStart),
				align1RangeWidth = INTEGER(alignedPatternRangeWidth),
				align1MismatchEnds = INTEGER(alignedPatternMismatchEnds),
				align1IndelEnds = INTEGER(alignedPatternIndelEnds),
				align2RangeStart = INTEGER(alignedSubjectRangeStart),
				align2RangeWidth = INTEGER(alignedSubjectRangeWidth),
				align2MismatchEnds = INTEGER(alignedSubjectMismatchEnds),
				align2IndelEnds = INTEGER(alignedSubjectIndelEnds);
		        i < numberOfStrings; i++, score++,
				align1RangeStart++, align1RangeWidth++, align1MismatchEnds++, align1IndelEnds++,
				align2RangeStart++, align2RangeWidth++, align2MismatchEnds++, align2IndelEnds++) {
	        R_CheckUserInterrupt();
			align1Info.string = _get_cachedXStringSet_elt(&cached_pattern, i);
			if (useQualityValue) {
				align1Info.quality = _get_cachedXStringSet_elt(&cached_patternQuality, quality1Element);
				quality1Element += quality1Increment;
			}
			if (multipleSubjects) {
				align2Info.string = _get_cachedXStringSet_elt(&cached_subject, i);
				if (useQualityValue) {
					align2Info.quality = _get_cachedXStringSet_elt(&cached_subjectQuality, quality2Element);
					quality2Element += quality2Increment;
				}
			}
			*score = pairwiseAlignment(
					&align1Info,
					&align2Info,
					localAlignment,
					scoreOnlyValue,
					gapOpeningValue,
					gapExtensionValue,
					useQualityValue,
					REAL(substitutionArray),
					INTEGER(substitutionArrayDim),
					INTEGER(substitutionLookupTable),
					LENGTH(substitutionLookupTable),
					INTEGER(fuzzyMatrix),
					INTEGER(fuzzyMatrixDim),
					INTEGER(fuzzyLookupTable),
					LENGTH(fuzzyLookupTable),
					&alignBuffer);

			*align1MismatchEnds = align1Info.lengthMismatch + align1MismatchPrevEnd;
			*align2MismatchEnds = align2Info.lengthMismatch + align2MismatchPrevEnd;
			if (align1Info.lengthMismatch > 0) {
				if ((mismatchBuffer.usedSpace + align1Info.lengthMismatch) > mismatchBuffer.totalSpace) {
					mismatchBuffer.totalSpace =
						mismatchBuffer.totalSpace +
							MIN(MAX_BUF_SIZE,
							    alignmentBufferSize + (numberOfStrings - (i+1)) * (alignmentBufferSize/4));
					tempIntPtr = (int *) R_alloc((long) mismatchBuffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, mismatchBuffer.pattern, mismatchBuffer.usedSpace * sizeof(int));
					mismatchBuffer.pattern = tempIntPtr;
					tempIntPtr = (int *) R_alloc((long) mismatchBuffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, mismatchBuffer.subject, mismatchBuffer.usedSpace * sizeof(int));
					mismatchBuffer.subject = tempIntPtr;
				}

				memcpy(&mismatchBuffer.pattern[mismatchBuffer.usedSpace], align1Info.mismatch,
					   align1Info.lengthMismatch * sizeof(int));

				memcpy(&mismatchBuffer.subject[mismatchBuffer.usedSpace], align2Info.mismatch,
					   align1Info.lengthMismatch * sizeof(int));
				mismatchBuffer.usedSpace = mismatchBuffer.usedSpace + align1Info.lengthMismatch;
			}

			*align1RangeStart = align1Info.startRange;
			*align1RangeWidth = align1Info.widthRange;
			*align1IndelEnds = align1Info.lengthIndel + align1IndelPrevEnd;
			if (align1Info.lengthIndel > 0) {
				if ((indel1Buffer.usedSpace + align1Info.lengthIndel) > indel1Buffer.totalSpace) {
					indel1Buffer.totalSpace =
						indel1Buffer.totalSpace +
							MIN(MAX_BUF_SIZE,
							    alignmentBufferSize + (numberOfStrings - (i+1)) * (alignmentBufferSize/12));
					tempIntPtr = (int *) R_alloc((long) indel1Buffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, indel1Buffer.start, indel1Buffer.usedSpace * sizeof(int));
					indel1Buffer.start = tempIntPtr;
					tempIntPtr = (int *) R_alloc((long) indel1Buffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, indel1Buffer.width, indel1Buffer.usedSpace * sizeof(int));
					indel1Buffer.width = tempIntPtr;
				}
				memcpy(&indel1Buffer.start[indel1Buffer.usedSpace], align1Info.startIndel,
					   align1Info.lengthIndel * sizeof(int));
				memcpy(&indel1Buffer.width[indel1Buffer.usedSpace], align1Info.widthIndel,
					   align1Info.lengthIndel * sizeof(int));
				indel1Buffer.usedSpace = indel1Buffer.usedSpace + align1Info.lengthIndel;
			}

			*align2RangeStart = align2Info.startRange;
			*align2RangeWidth = align2Info.widthRange;
			*align2IndelEnds = align2Info.lengthIndel + align2IndelPrevEnd;
			if (align2Info.lengthIndel > 0) {
				if ((indel2Buffer.usedSpace + align2Info.lengthIndel) > indel2Buffer.totalSpace) {
					indel2Buffer.totalSpace =
						indel2Buffer.totalSpace +
							MIN(MAX_BUF_SIZE,
							    alignmentBufferSize + (numberOfStrings - (i+1)) * (alignmentBufferSize/12));
					tempIntPtr = (int *) R_alloc((long) indel2Buffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, indel2Buffer.start, indel2Buffer.usedSpace * sizeof(int));
					indel2Buffer.start = tempIntPtr;
					tempIntPtr = (int *) R_alloc((long) indel2Buffer.totalSpace, sizeof(int));
					memcpy(tempIntPtr, indel2Buffer.width, indel2Buffer.usedSpace * sizeof(int));
					indel2Buffer.width = tempIntPtr;
				}
				memcpy(&indel2Buffer.start[indel2Buffer.usedSpace], align2Info.startIndel,
					   align2Info.lengthIndel * sizeof(int));
				memcpy(&indel2Buffer.width[indel2Buffer.usedSpace], align2Info.widthIndel,
					   align2Info.lengthIndel * sizeof(int));
				indel2Buffer.usedSpace = indel2Buffer.usedSpace + align2Info.lengthIndel;
			}

			align1MismatchPrevEnd = *align1MismatchEnds;
			align2MismatchPrevEnd = *align2MismatchEnds;
			align1IndelPrevEnd = *align1IndelEnds;
			align2IndelPrevEnd = *align2IndelEnds;
		}

		/* Create the output object */
		if (multipleSubjects) {
			PROTECT(output = NEW_OBJECT(MAKE_CLASS("PairwiseAlignments")));
		} else {
			PROTECT(output = NEW_OBJECT(MAKE_CLASS("PairwiseAlignmentsSingleSubject")));
		}

		/* Set the "pattern" slot */
		if (useQualityValue) {
			PROTECT(alignedPattern = NEW_OBJECT(MAKE_CLASS("QualityAlignedXStringSet")));
		} else {
			PROTECT(alignedPattern = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		}
		SET_SLOT(alignedPattern, mkChar("unaligned"), pattern);
		/* Set the "range" sub-slot */
		PROTECT(alignedPatternRange =
			new_IRanges("IRanges", alignedPatternRangeStart, alignedPatternRangeWidth, R_NilValue));
		SET_SLOT(alignedPattern, mkChar("range"), alignedPatternRange);
		/* Set the "mismatch" sub-slot */
		PROTECT(alignedPatternMismatch = NEW_OBJECT(MAKE_CLASS("CompressedIntegerList")));
		PROTECT(alignedPatternMismatchPartitioning = NEW_OBJECT(MAKE_CLASS("PartitioningByEnd")));
		PROTECT(alignedPatternMismatchValues = NEW_INTEGER(mismatchBuffer.usedSpace));
		memcpy(INTEGER(alignedPatternMismatchValues), mismatchBuffer.pattern,
			   mismatchBuffer.usedSpace * sizeof(int));
		SET_SLOT(alignedPatternMismatchPartitioning, mkChar("end"), alignedPatternMismatchEnds);
		SET_SLOT(alignedPatternMismatch, mkChar("partitioning"), alignedPatternMismatchPartitioning);
		SET_SLOT(alignedPatternMismatch, mkChar("unlistData"), alignedPatternMismatchValues);
		SET_SLOT(alignedPattern, mkChar("mismatch"), alignedPatternMismatch);
		/* Set the "indel" sub-slot */
		PROTECT(alignedPatternIndel = NEW_OBJECT(MAKE_CLASS("CompressedIRangesList")));
		PROTECT(alignedPatternIndelPartitioning = NEW_OBJECT(MAKE_CLASS("PartitioningByEnd")));
		PROTECT(alignedPatternIndelRangeStart = NEW_INTEGER(indel1Buffer.usedSpace));
		PROTECT(alignedPatternIndelRangeWidth = NEW_INTEGER(indel1Buffer.usedSpace));
		memcpy(INTEGER(alignedPatternIndelRangeStart), indel1Buffer.start,
			   indel1Buffer.usedSpace * sizeof(int));
		memcpy(INTEGER(alignedPatternIndelRangeWidth), indel1Buffer.width,
			   indel1Buffer.usedSpace * sizeof(int));
		PROTECT(alignedPatternIndelRange =
			new_IRanges("IRanges", alignedPatternIndelRangeStart, alignedPatternIndelRangeWidth, R_NilValue));
		SET_SLOT(alignedPatternIndelPartitioning, mkChar("end"), alignedPatternIndelEnds);
		SET_SLOT(alignedPatternIndel, mkChar("partitioning"), alignedPatternIndelPartitioning);
	    SET_SLOT(alignedPatternIndel, mkChar("unlistData"), alignedPatternIndelRange);
		SET_SLOT(alignedPattern, mkChar("indel"), alignedPatternIndel);
		SET_SLOT(output, mkChar("pattern"), alignedPattern);

		/* Set the "subject" slot */
		if (useQualityValue) {
			PROTECT(alignedSubject = NEW_OBJECT(MAKE_CLASS("QualityAlignedXStringSet")));
		} else {
			PROTECT(alignedSubject = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		}
		SET_SLOT(alignedSubject, mkChar("unaligned"), subject);
		/* Set the "range" sub-slot */
		PROTECT(alignedSubjectRange =
			new_IRanges("IRanges", alignedSubjectRangeStart, alignedSubjectRangeWidth, R_NilValue));
		SET_SLOT(alignedSubject, mkChar("range"), alignedSubjectRange);
		/* Set the "mismatch" sub-slot */
		PROTECT(alignedSubjectMismatch = NEW_OBJECT(MAKE_CLASS("CompressedIntegerList")));
		PROTECT(alignedSubjectMismatchPartitioning = NEW_OBJECT(MAKE_CLASS("PartitioningByEnd")));
		PROTECT(alignedSubjectMismatchValues = NEW_INTEGER(mismatchBuffer.usedSpace));
		memcpy(INTEGER(alignedSubjectMismatchValues), mismatchBuffer.subject,
			   mismatchBuffer.usedSpace * sizeof(int));
	    SET_SLOT(alignedSubjectMismatchPartitioning, mkChar("end"), alignedSubjectMismatchEnds);
		SET_SLOT(alignedSubjectMismatch, mkChar("partitioning"), alignedSubjectMismatchPartitioning);
		SET_SLOT(alignedSubjectMismatch, mkChar("unlistData"), alignedSubjectMismatchValues);
		SET_SLOT(alignedSubject, mkChar("mismatch"), alignedSubjectMismatch);
		/* Set the "indel" sub-slot */
		PROTECT(alignedSubjectIndel = NEW_OBJECT(MAKE_CLASS("CompressedIRangesList")));
		PROTECT(alignedSubjectIndelPartitioning = NEW_OBJECT(MAKE_CLASS("PartitioningByEnd")));
		PROTECT(alignedSubjectIndelRangeStart = NEW_INTEGER(indel2Buffer.usedSpace));
		PROTECT(alignedSubjectIndelRangeWidth = NEW_INTEGER(indel2Buffer.usedSpace));
		memcpy(INTEGER(alignedSubjectIndelRangeStart), indel2Buffer.start,
			   indel2Buffer.usedSpace * sizeof(int));
		memcpy(INTEGER(alignedSubjectIndelRangeWidth), indel2Buffer.width,
			   indel2Buffer.usedSpace * sizeof(int));
		PROTECT(alignedSubjectIndelRange =
			new_IRanges("IRanges", alignedSubjectIndelRangeStart, alignedSubjectIndelRangeWidth, R_NilValue));
	    SET_SLOT(alignedSubjectIndelPartitioning, mkChar("end"), alignedSubjectIndelEnds);
		SET_SLOT(alignedSubjectIndel, mkChar("partitioning"), alignedSubjectIndelPartitioning);
		SET_SLOT(alignedSubjectIndel, mkChar("unlistData"), alignedSubjectIndelRange);
		SET_SLOT(alignedSubject, mkChar("indel"), alignedSubjectIndel);
		SET_SLOT(output, mkChar("subject"), alignedSubject);

		/* Set the "score" slot */
		SET_SLOT(output, mkChar("score"), alignedScore);

		/* Set the "type" slot */
		SET_SLOT(output, mkChar("type"), type);

		/* Set the "gapOpening" slot */
		SET_SLOT(output, mkChar("gapOpening"), gapOpening);
		/* Set the "gapExtension" slot */
		SET_SLOT(output, mkChar("gapExtension"), gapExtension);

		/* Output is ready */
		UNPROTECT(30);
	}

	return output;
}



/*
 * INPUTS
 * 'string':                   XStringSet object for strings
 * 'type':                     type of pairwise alignment
 *                             (character vector of length 1;
 *                              'global', 'local', 'overlap')
 * 'typeCode':                 type of pairwise alignment
 *                             (integer vector of length 1;
 *                              1 = 'global', 2 = 'local', 3 = 'overlap')
 * 'gapOpening':               gap opening cost or penalty
 *                             (double vector of length 1)
 * 'gapExtension':             gap extension cost or penalty
 *                             (double vector of length 1)
 * 'useQuality':               denotes whether or not to use quality measures
 *                             in the optimal pairwise alignment
 *                             (logical vector of length 1)
 * 'substitutionArray':        a three-dimensional double array where the first two
 *                             dimensions are for substitutions and the
 *                             third is for fuzziness of matches
 * 'substitutionArrayDim':     dimension of 'substitutionArray'
 *                             (integer vector of length 3)
 * 'substitutionLookupTable':  lookup table for translating XString bytes to
 *                             substitution indices
 *                             (integer vector)
 * 'fuzzyMatrix':              fuzzy matrix for matches
 *                             (double matrix)
 * 'fuzzyMatrixDim':           dimension of 'fuzzyMatrix'
 *                             (integer vector of length 2)
 * 'fuzzyLookupTable':         lookup table for translating XString bytes to
 *                             fuzzy indices
 *                             (integer vector)
 *
 * OUTPUT
 * Return a numeric vector containing the lower triangle of the score matrix.
 */

SEXP XStringSet_align_distance(
		SEXP string,
		SEXP type,
		SEXP typeCode,
		SEXP gapOpening,
		SEXP gapExtension,
		SEXP useQuality,
		SEXP substitutionArray,
		SEXP substitutionArrayDim,
		SEXP substitutionLookupTable,
		SEXP fuzzyMatrix,
		SEXP fuzzyMatrixDim,
		SEXP fuzzyLookupTable)
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
	int lengthOfStringQualitySet = 0;

	SEXP stringQuality = R_NilValue;
	cachedXStringSet cached_string = _cache_XStringSet(string);
	cachedXStringSet cached_stringQuality;
	if (useQualityValue) {
		stringQuality = GET_SLOT(string, install("quality"));
		cached_stringQuality = _cache_XStringSet(stringQuality);
		lengthOfStringQualitySet = _get_XStringSet_length(stringQuality);
	}

	SEXP output;

	int i, j, iQualityElement = 0, jQualityElement = 0;
	int qualityIncrement = ((lengthOfStringQualitySet < numberOfStrings) ? 0 : 1);

	/* Create the alignment buffer object */
	struct AlignBuffer alignBuffer;
	int nCharString = 0;
	for (i = 0; i < numberOfStrings; i++) {
		nCharString = MAX(nCharString, _get_cachedXStringSet_elt(&cached_string, i).length);
	}
	int alignmentBufferSize = nCharString + 1;
	alignBuffer.currMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));
	alignBuffer.prevMatrix = (float *) R_alloc((long) 3 * alignmentBufferSize, sizeof(float));

	double *score;
	PROTECT(output = NEW_NUMERIC((numberOfStrings * (numberOfStrings - 1)) / 2));
	score = REAL(output);
	if (!useQualityValue) {
		for (i = 0; i < numberOfStrings; i++) {
	        R_CheckUserInterrupt();
			align1Info.string = _get_cachedXStringSet_elt(&cached_string, i);
			for (j = i + 1; j < numberOfStrings; j++) {
				align2Info.string = _get_cachedXStringSet_elt(&cached_string, j);
				*score = pairwiseAlignment(
						&align1Info,
						&align2Info,
						localAlignment,
						scoreOnlyValue,
						gapOpeningValue,
						gapExtensionValue,
						useQualityValue,
						REAL(substitutionArray),
						INTEGER(substitutionArrayDim),
						INTEGER(substitutionLookupTable),
						LENGTH(substitutionLookupTable),
						INTEGER(fuzzyMatrix),
						INTEGER(fuzzyMatrixDim),
						INTEGER(fuzzyLookupTable),
						LENGTH(fuzzyLookupTable),
						&alignBuffer);
				score++;
			}
		}
	} else {
		for (i = 0; i < numberOfStrings; i++) {
	        R_CheckUserInterrupt();
			align1Info.string = _get_cachedXStringSet_elt(&cached_string, i);
			align1Info.quality = _get_cachedXStringSet_elt(&cached_stringQuality, iQualityElement);
			jQualityElement = iQualityElement + qualityIncrement;
			iQualityElement += qualityIncrement;
			for (j = i + 1; j < numberOfStrings; j++) {
				align2Info.string = _get_cachedXStringSet_elt(&cached_string, j);
				align2Info.quality = _get_cachedXStringSet_elt(&cached_stringQuality, jQualityElement);
				jQualityElement += qualityIncrement;
				*score = pairwiseAlignment(
						&align1Info,
						&align2Info,
						localAlignment,
						scoreOnlyValue,
						gapOpeningValue,
						gapExtensionValue,
						useQualityValue,
						REAL(substitutionArray),
						INTEGER(substitutionArrayDim),
						INTEGER(substitutionLookupTable),
						LENGTH(substitutionLookupTable),
						INTEGER(fuzzyMatrix),
						INTEGER(fuzzyMatrixDim),
						INTEGER(fuzzyLookupTable),
						LENGTH(fuzzyLookupTable),
						&alignBuffer);
				score++;
			}
		}
	}
	UNPROTECT(1);

	return output;
}
