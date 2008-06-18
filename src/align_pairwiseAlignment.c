#include <float.h>
#include "Biostrings.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

#define NEGATIVE_INFINITY (- FLT_MAX)

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
	RoSeq string;
	RoSeq quality;
	int endGap;
	int startRange;
	int widthRange;
	int* startIndels;
	int* widthIndels;
	int lengthIndels;
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
	int nQuality1 = align1InfoPtr->quality.nelt;
	int nQuality2 = align2InfoPtr->quality.nelt;

	if (nCharString1 < 1 || nCharString2 < 1)
		return NA_REAL;

	/* Step 2:  Create objects for scores values */
	/* Rows of currMatrix and prevMatrix = (0) substitution, (1) deletion, and (2) insertion */
	float *currMatrix = (float *) R_alloc((long) (3 * nCharString2Plus1), sizeof(float));
	float *prevMatrix = (float *) R_alloc((long) (3 * nCharString2Plus1), sizeof(float));

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

	/* Step 3:  Perform main alignment operations */
	RoSeq sequence1, sequence2;
	int scalar1, scalar2;
	const int *lookupTable;
	int lookupTableLength;
	const double *matchMatrix, *mismatchMatrix;
	const int *matrixDim;
	if (nQuality1 == 0) {
		sequence1 = align1InfoPtr->string;
		sequence2 = align2InfoPtr->string;
		scalar1 = (nCharString1 == 1);
		scalar2 = (nCharString2 == 1);
		lookupTable = constantLookupTable;
		lookupTableLength = constantLookupTableLength;
		matchMatrix = constantMatrix;
		mismatchMatrix = constantMatrix;
		matrixDim = constantMatrixDim;
	} else {
		sequence1 = align1InfoPtr->quality;
		sequence2 = align2InfoPtr->quality;
		scalar1 = (nQuality1 == 1);
		scalar2 = (nQuality2 == 1);
		lookupTable = qualityLookupTable;
		lookupTableLength = qualityLookupTableLength;
		matchMatrix = qualityMatchMatrix;
		mismatchMatrix = qualityMismatchMatrix;
		matrixDim = qualityMatrixDim;
	}
	int lookupValue = 0, elements[2], iElt, jElt;
	int notSwappedOrder = 1 - swappedOrder;
	float gapOpeningPlusExtension = gapOpening + gapExtension;
	float endGapAddend = (align1InfoPtr->endGap ? gapExtension : 0.0);
	float *tempMatrix, substitutionValue;
	double maxScore = NEGATIVE_INFINITY;
	if (scoreOnly) {
		/* Simplified calculations when only need the alignment score */
		for (i = 1, iMinus1 = 0, iElt = nCharString1Minus1; i <= nCharString1; i++, iMinus1++, iElt--) {
			tempMatrix = prevMatrix;
			prevMatrix = currMatrix;
			currMatrix = tempMatrix;

			CURR_MATRIX(0, 0) = NEGATIVE_INFINITY;
			CURR_MATRIX(1, 0) = PREV_MATRIX(1, 0) + endGapAddend;
			CURR_MATRIX(2, 0) = NEGATIVE_INFINITY;

			SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence1.elts[scalar1 ? 0 : iElt]);
			elements[swappedOrder] = lookupValue;
			for (j = 1, jMinus1 = 0, jElt = nCharString2Minus1; j <= nCharString2; j++, jMinus1++, jElt--) {
				SET_LOOKUP_VALUE(lookupTable, lookupTableLength, sequence2.elts[scalar2 ? 0 : jElt]);
				elements[notSwappedOrder] = lookupValue;
				if (align1InfoPtr->string.elts[iElt] == align2InfoPtr->string.elts[jElt])
					substitutionValue = (float) matchMatrix[matrixDim[0] * elements[0] + elements[1]];
				else
					substitutionValue = (float) mismatchMatrix[matrixDim[0] * elements[0] + elements[1]];

				CURR_MATRIX(0, j) =
					SAFE_SUM(MAX(PREV_MATRIX(0, jMinus1), MAX(PREV_MATRIX(1, jMinus1), PREV_MATRIX(2, jMinus1))),
					         substitutionValue);
				CURR_MATRIX(1, j) = 
					MAX(SAFE_SUM(MAX(PREV_MATRIX(0, j), PREV_MATRIX(2, j)), gapOpeningPlusExtension),
					    SAFE_SUM(PREV_MATRIX(1, j), gapExtension));
				CURR_MATRIX(2, j) =
					MAX(SAFE_SUM(MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(1, jMinus1)), gapOpeningPlusExtension),
					    SAFE_SUM(CURR_MATRIX(2, jMinus1), gapExtension));

				if (localAlignment) {
					CURR_MATRIX(0, j) = MAX(0.0, CURR_MATRIX(0, j));
					maxScore = MAX(CURR_MATRIX(0, j), maxScore);
				} else {
					if (!align1InfoPtr->endGap && j == nCharString2)
						CURR_MATRIX(1, j) = MAX(PREV_MATRIX(0, j), PREV_MATRIX(1, j));

					if (!align2InfoPtr->endGap && i == nCharString1)
						CURR_MATRIX(2, j) = MAX(CURR_MATRIX(0, jMinus1), CURR_MATRIX(2, jMinus1));
				}
			}
		}

		if (!localAlignment) {
			maxScore =
				MAX(CURR_MATRIX(0, nCharString2),
				MAX(CURR_MATRIX(1, nCharString2),
				    CURR_MATRIX(2, nCharString2)));
		}
	} else {
		/* Step 3a:  Create objects for traceback values */
		char *sTraceMatrix = (char *) R_alloc((long) nCharString1 * nCharString2, sizeof(char));
		char *iTraceMatrix = (char *) R_alloc((long) nCharString1 * nCharString2, sizeof(char));
		char *dTraceMatrix = (char *) R_alloc((long) nCharString1 * nCharString2, sizeof(char));

		/* Step 3b:  Prepare the alignment info object for alignment */
		int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;
		int* startIndels1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* startIndels2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* widthIndels1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* widthIndels2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));

		for(int i = 0; i < alignmentBufferSize; i++) {
			startIndels1Buffer[i] = 0;
			startIndels2Buffer[i] = 0;
			widthIndels1Buffer[i] = 0;
			widthIndels2Buffer[i] = 0;
		}
		align1InfoPtr->lengthIndels = 0;
		align2InfoPtr->lengthIndels = 0;
		align1InfoPtr->startIndels = startIndels1Buffer + alignmentBufferSize;
		align2InfoPtr->startIndels = startIndels2Buffer + alignmentBufferSize;
		align1InfoPtr->widthIndels = widthIndels1Buffer + alignmentBufferSize;
		align2InfoPtr->widthIndels = widthIndels2Buffer + alignmentBufferSize;

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
					CURR_MATRIX(0, j) = SAFE_SUM(PREV_MATRIX(0, jMinus1), substitutionValue);
				} else if (PREV_MATRIX(2, jMinus1) >= PREV_MATRIX(1, jMinus1)) {
					S_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
					CURR_MATRIX(0, j) = SAFE_SUM(PREV_MATRIX(2, jMinus1), substitutionValue);
				} else {
					S_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
					CURR_MATRIX(0, j) = SAFE_SUM(PREV_MATRIX(1, jMinus1), substitutionValue);
				}
				if ((SAFE_SUM(PREV_MATRIX(0, j), gapOpening) >= PREV_MATRIX(1, j)) &&
					(PREV_MATRIX(0, j) >= PREV_MATRIX(2, j))) {
					D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
					CURR_MATRIX(1, j) = SAFE_SUM(PREV_MATRIX(0, j), gapOpeningPlusExtension);
				} else if (SAFE_SUM(PREV_MATRIX(2, j), gapOpening) > PREV_MATRIX(1, j)) {
					D_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
					CURR_MATRIX(1, j) = SAFE_SUM(PREV_MATRIX(2, j), gapOpeningPlusExtension);
				} else {
					D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
					CURR_MATRIX(1, j) = SAFE_SUM(PREV_MATRIX(1, j), gapExtension);
				}
				if ((SAFE_SUM(CURR_MATRIX(0, jMinus1), gapOpening) >= CURR_MATRIX(2, jMinus1)) &&
					(CURR_MATRIX(0, jMinus1) >= CURR_MATRIX(1, jMinus1))) {
					I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
					CURR_MATRIX(2, j) = SAFE_SUM(CURR_MATRIX(0, jMinus1), gapOpeningPlusExtension);
				} else if (SAFE_SUM(CURR_MATRIX(1, jMinus1), gapOpening) > CURR_MATRIX(2, jMinus1)) {
					I_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
					CURR_MATRIX(2, j) = SAFE_SUM(CURR_MATRIX(1, jMinus1), gapOpeningPlusExtension);
				} else {
					I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
					CURR_MATRIX(2, j) = SAFE_SUM(CURR_MATRIX(2, jMinus1), gapExtension);
				}

				if (localAlignment) {
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
				} else {

					if (!align1InfoPtr->endGap && j == nCharString2) {
						if (PREV_MATRIX(0, j) >= PREV_MATRIX(1, j)) {
							D_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
							CURR_MATRIX(1, j) = PREV_MATRIX(0, j);
						} else {
							D_TRACE_MATRIX(iMinus1, jMinus1) = DELETION;
							CURR_MATRIX(1, j) = PREV_MATRIX(1, j);
						}
					}

					if (!align2InfoPtr->endGap && i == nCharString1) {
						if (CURR_MATRIX(0, jMinus1) >= CURR_MATRIX(2, jMinus1)) {
							I_TRACE_MATRIX(iMinus1, jMinus1) = SUBSTITUTION;
							CURR_MATRIX(2, j) = CURR_MATRIX(0, jMinus1);
						} else {
							I_TRACE_MATRIX(iMinus1, jMinus1) = INSERTION;
							CURR_MATRIX(2, j) = CURR_MATRIX(2, jMinus1);
						}
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
		    					align2InfoPtr->startIndels--;
		    					align2InfoPtr->widthIndels--;
		    					align2InfoPtr->lengthIndels++;
		    					*align2InfoPtr->startIndels = nCharString2 - j;
		    				}
		    				*align2InfoPtr->widthIndels += 1;
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
		    					align1InfoPtr->startIndels--;
		    					align1InfoPtr->widthIndels--;
		    					align1InfoPtr->lengthIndels++;
		    					*align1InfoPtr->startIndels = nCharString1 - i;
		    				}
		    				*align1InfoPtr->widthIndels += 1;
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
	    			i--;
	    			j--;
	    			break;
	    		default:
	    			error("unknown traceback code %d", currTraceMatrix);
	    			break;
			}
		}

		int offset1 = align1InfoPtr->startRange - 1;
		if (offset1 > 0 && align1InfoPtr->lengthIndels > 0) {
			for (i = 0; i < align1InfoPtr->lengthIndels; i++)
				align1InfoPtr->startIndels[i] -= offset1;
		}
		int offset2 = align2InfoPtr->startRange - 1;
		if (offset2 > 0 && align2InfoPtr->lengthIndels > 0) {
			for (j = 0; j < align2InfoPtr->lengthIndels; j++)
				align2InfoPtr->startIndels[j] -= offset2;
		}
	}

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

SEXP XStringSet_align_pairwiseAlignment(
		SEXP pattern,
		SEXP subject,
		SEXP patternQuality,
		SEXP subjectQuality,
		SEXP type,
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
	int scoreOnlyValue = LOGICAL(scoreOnly)[0];

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

	int qualityIncrement = ((_get_XStringSet_length(patternQuality) < numberOfStrings) ? 0 : 1);
	int qualityElement = 0;
	double *score;
	if (scoreOnlyValue) {
		PROTECT(output = NEW_NUMERIC(numberOfStrings));
		score = REAL(output);
		for (int i = 0; i < numberOfStrings; i++) {
			align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
			align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, qualityElement);
			score[i] = pairwiseAlignment(
					&align1Info,
					&align2Info,
					(INTEGER(typeCode)[0] == LOCAL_ALIGNMENT),
					scoreOnlyValue,
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
					INTEGER(constantMatrixDim));
			qualityElement += qualityIncrement;
		}
		UNPROTECT(1);
	} else {
		SEXP alignedPattern;
		SEXP alignedPatternRange, alignedPatternRangeStart, alignedPatternRangeWidth;
		SEXP alignedPatternIndels, alignedPatternIndelsRange;
		SEXP alignedPatternIndelsRangeStart, alignedPatternIndelsRangeWidth;
		SEXP alignedSubject;
		SEXP alignedSubjectRange, alignedSubjectRangeStart, alignedSubjectRangeWidth;
		SEXP alignedSubjectIndels, alignedSubjectIndelsRange;
		SEXP alignedSubjectIndelsRangeStart, alignedSubjectIndelsRangeWidth;
		SEXP alignedScore;

		PROTECT(output = NEW_OBJECT(MAKE_CLASS("PairwiseAlignment")));

		/* Set the "pattern" slot */
		PROTECT(alignedPattern = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		SET_SLOT(alignedPattern, mkChar("unaligned"), pattern);
		SET_SLOT(alignedPattern, mkChar("quality"), patternQuality);
		/* Set the "range" sub-slot */
		PROTECT(alignedPatternRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
		PROTECT(alignedPatternRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedPatternRangeWidth = NEW_INTEGER(numberOfStrings));
		int *align1RangeStart = INTEGER(alignedPatternRangeStart);
		int *align1RangeWidth = INTEGER(alignedPatternRangeWidth);
		SET_SLOT(alignedPatternRange, mkChar("start"), alignedPatternRangeStart);
		SET_SLOT(alignedPatternRange, mkChar("width"), alignedPatternRangeWidth);
		SET_SLOT(alignedPattern, mkChar("range"), alignedPatternRange);
		/* Set the "indels" sub-slot */
		PROTECT(alignedPatternIndels = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedPattern, mkChar("indels"), alignedPatternIndels);
		SET_SLOT(output, mkChar("pattern"), alignedPattern);

		/* Set the "subject" slot */
		PROTECT(alignedSubject = NEW_OBJECT(MAKE_CLASS("AlignedXStringSet")));
		SET_SLOT(alignedSubject, mkChar("unaligned"), subject);
		SET_SLOT(alignedSubject, mkChar("quality"), subjectQuality);
		/* Set the "range" sub-slot */
		PROTECT(alignedSubjectRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
		PROTECT(alignedSubjectRangeStart = NEW_INTEGER(numberOfStrings));
		PROTECT(alignedSubjectRangeWidth = NEW_INTEGER(numberOfStrings));
		int *align2RangeStart = INTEGER(alignedSubjectRangeStart);
		int *align2RangeWidth = INTEGER(alignedSubjectRangeWidth);
		SET_SLOT(alignedSubjectRange, mkChar("start"), alignedSubjectRangeStart);
		SET_SLOT(alignedSubjectRange, mkChar("width"), alignedSubjectRangeWidth);
		SET_SLOT(alignedSubject, mkChar("range"), alignedSubjectRange);
		/* Set the "indels" sub-slot */
		PROTECT(alignedSubjectIndels = NEW_LIST(numberOfStrings));
		SET_SLOT(alignedSubject, mkChar("indels"), alignedSubjectIndels);
		SET_SLOT(output, mkChar("subject"), alignedSubject);

		/* Set the "score" slot */
		PROTECT(alignedScore = NEW_NUMERIC(numberOfStrings));
		score = REAL(alignedScore);
		SET_SLOT(output, mkChar("score"), alignedScore);

		/* Set the "type" slot */
		SET_SLOT(output, mkChar("type"), type);
		/* Set the "constantMatrix" slot */
		SET_SLOT(output, mkChar("constantMatrix"), constantMatrix);
		/* Set the "gapOpening" slot */
		SET_SLOT(output, mkChar("gapOpening"), gapOpening);
		/* Set the "gapExtension" slot */
		SET_SLOT(output, mkChar("gapExtension"), gapExtension);

		for (int i = 0; i < numberOfStrings; i++) {
			align1Info.lengthIndels = 0;
			align2Info.lengthIndels = 0;

			align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
			align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, qualityElement);
			score[i] = pairwiseAlignment(
					&align1Info,
					&align2Info,
					(INTEGER(typeCode)[0] == LOCAL_ALIGNMENT),
					scoreOnlyValue,
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
					INTEGER(constantMatrixDim));

			align1RangeStart[i] = align1Info.startRange;
			align1RangeWidth[i] = align1Info.widthRange;
			PROTECT(alignedPatternIndelsRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
			PROTECT(alignedPatternIndelsRangeStart = NEW_INTEGER(align1Info.lengthIndels));
			PROTECT(alignedPatternIndelsRangeWidth = NEW_INTEGER(align1Info.lengthIndels));
			for (int j = 0; j < align1Info.lengthIndels; j++) {
				int infoIndex = align1Info.lengthIndels - 1 - j;
				INTEGER(alignedPatternIndelsRangeStart)[j] = align1Info.startIndels[infoIndex];
				INTEGER(alignedPatternIndelsRangeWidth)[j] = align1Info.widthIndels[infoIndex];
			}
			SET_SLOT(alignedPatternIndelsRange, mkChar("start"), alignedPatternIndelsRangeStart);
			SET_SLOT(alignedPatternIndelsRange, mkChar("width"), alignedPatternIndelsRangeWidth);
		    SET_VECTOR_ELT(alignedPatternIndels, i, alignedPatternIndelsRange);
		    UNPROTECT(3);

			align2RangeStart[i] = align2Info.startRange;
			align2RangeWidth[i] = align2Info.widthRange;
			PROTECT(alignedSubjectIndelsRange = NEW_OBJECT(MAKE_CLASS("IRanges")));
			PROTECT(alignedSubjectIndelsRangeStart = NEW_INTEGER(align2Info.lengthIndels));
			PROTECT(alignedSubjectIndelsRangeWidth = NEW_INTEGER(align2Info.lengthIndels));
			for (int j = 0; j < align2Info.lengthIndels; j++) {
				int infoIndex = align2Info.lengthIndels - 1 - j;
				INTEGER(alignedSubjectIndelsRangeStart)[j] = align2Info.startIndels[infoIndex];
				INTEGER(alignedSubjectIndelsRangeWidth)[j] = align2Info.widthIndels[infoIndex];
			}
			SET_SLOT(alignedSubjectIndelsRange, mkChar("start"), alignedSubjectIndelsRangeStart);
			SET_SLOT(alignedSubjectIndelsRange, mkChar("width"), alignedSubjectIndelsRangeWidth);
		    SET_VECTOR_ELT(alignedSubjectIndels, i, alignedSubjectIndelsRange);
		    UNPROTECT(3);

			qualityElement += qualityIncrement;
		}

		/* Output is ready */
		UNPROTECT(12);
	}

	return output;
}
