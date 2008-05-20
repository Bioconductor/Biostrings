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
#define INSERTION    'I'
#define DELETION     'D'
#define TERMINATION  'T'

#define CURR_MATRIX(i, j) (currMatrix[nCharString2Plus1 * i + j])
#define PREV_MATRIX(i, j) (prevMatrix[nCharString2Plus1 * i + j])
#define TRACE_MATRIX(i, j) (traceMatrix[nCharString2Plus1 * i + j])

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
	int* startInserts;
	int* widthInserts;
	int lengthInserts;
	double* profile;
};
void function(struct AlignInfo *);


/* Returns the score of the optimal pairwise alignment */
static float pairwiseAlignment(
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
	int nCharString1Plus1 = nCharString1 + 1;
	int nCharString2Plus1 = nCharString2 + 1;
	int nCharString1Minus1 = nCharString1 - 1;
	int nCharString2Minus1 = nCharString2 - 1;
	int nQuality1 = align1InfoPtr->quality.nelt;
	int nQuality2 = align2InfoPtr->quality.nelt;

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
						CURR_MATRIX(1, j) = 
							MAX(PREV_MATRIX(0, j),
							MAX(PREV_MATRIX(1, j),
								PREV_MATRIX(2, j)));

					if (!align2InfoPtr->endGap && i == nCharString1)
						CURR_MATRIX(2, j) =
							MAX(CURR_MATRIX(0, jMinus1),
							MAX(CURR_MATRIX(1, jMinus1),
								CURR_MATRIX(2, jMinus1)));
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
		char *traceMatrix = (char *) R_alloc((long) (nCharString1Plus1 * nCharString2Plus1), sizeof(char));
		for (i = 0; i <= nCharString1; i++)
			TRACE_MATRIX(i, 0) = TERMINATION;
		for (j = 0; j <= nCharString2; j++)
			TRACE_MATRIX(0, j) = TERMINATION;

		/* Step 3b:  Prepare the alignment info object for alignment */
		int alignmentBufferSize = MIN(nCharString1, nCharString2) + 1;
		int* startInserts1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* startInserts2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* widthInserts1Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));
		int* widthInserts2Buffer = (int *) R_alloc((long) alignmentBufferSize, sizeof(int));

		for(int i = 0; i < alignmentBufferSize; i++) {
			startInserts1Buffer[i] = 0;
			startInserts2Buffer[i] = 0;
			widthInserts1Buffer[i] = 0;
			widthInserts2Buffer[i] = 0;
		}
		align1InfoPtr->lengthInserts = 0;
		align2InfoPtr->lengthInserts = 0;
		align1InfoPtr->startInserts = startInserts1Buffer + alignmentBufferSize;
		align2InfoPtr->startInserts = startInserts2Buffer + alignmentBufferSize;
		align1InfoPtr->widthInserts = widthInserts1Buffer + alignmentBufferSize;
		align2InfoPtr->widthInserts = widthInserts2Buffer + alignmentBufferSize;

		align1InfoPtr->startRange = -1;
		align2InfoPtr->startRange = -1;
		align1InfoPtr->widthRange = 0;
		align2InfoPtr->widthRange = 0;
		for (j = 1, jMinus1 = 0; j <= nCharString2; j++, jMinus1++)
			align2InfoPtr->profile[jMinus1] = NEGATIVE_INFINITY;
		for (i = 1, iMinus1 = 0; i <= nCharString1; i++, iMinus1++)
			align1InfoPtr->profile[iMinus1] = NEGATIVE_INFINITY;

		double scoreSubstitution = 0, scoreDeletion = 0, scoreInsertion = 0;
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

				/* Step 3c:  Generate (0) substitution, (1) deletion, and (2) insertion scores */
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

					/* Step 3d:  Get the optimal score for local alignments */
					if (CURR_MATRIX(0, j) > maxScore) {
						align1InfoPtr->startRange = iElt + 1;
						align2InfoPtr->startRange = jElt + 1;
						maxScore = CURR_MATRIX(0, j);
					}
				} else {
					if (!align1InfoPtr->endGap && j == nCharString2)
						CURR_MATRIX(1, j) = 
							MAX(PREV_MATRIX(0, j),
							MAX(PREV_MATRIX(1, j),
								PREV_MATRIX(2, j)));

					if (!align2InfoPtr->endGap && i == nCharString1)
						CURR_MATRIX(2, j) =
							MAX(CURR_MATRIX(0, jMinus1),
							MAX(CURR_MATRIX(1, jMinus1),
								CURR_MATRIX(2, jMinus1)));
				}

				/* Step 3e:  Generate the traceback values */
				scoreSubstitution = CURR_MATRIX(0, j);
				scoreDeletion     = CURR_MATRIX(1, j);
				scoreInsertion    = CURR_MATRIX(2, j);

				if (localAlignment && scoreSubstitution == 0.0)
					TRACE_MATRIX(i, j) = TERMINATION;
				else if (scoreSubstitution >= MAX(scoreInsertion, scoreDeletion))
					TRACE_MATRIX(i, j) = SUBSTITUTION;
				else if (scoreInsertion >= scoreDeletion)
					TRACE_MATRIX(i, j) = INSERTION;
				else
					TRACE_MATRIX(i, j) = DELETION;

				/* Step 3f:  Generate profile scores for local alignments */
				if (localAlignment) {
					align1InfoPtr->profile[iElt] =
						MAX(scoreSubstitution, align1InfoPtr->profile[iElt]);
					align2InfoPtr->profile[jElt] =
						MAX(scoreSubstitution, align2InfoPtr->profile[jElt]);
				}
			}

			if (!localAlignment) {
				double profile1Score = MAX(scoreSubstitution, scoreInsertion);
				if (!align1InfoPtr->endGap || i == nCharString1) {
						align1InfoPtr->profile[iElt] = profile1Score;
				} else {
					align1InfoPtr->profile[iElt] =
						SAFE_SUM(profile1Score, gapOpening + iElt * gapExtension);
					if (profile1Score >= scoreDeletion) {
						align1InfoPtr->profile[iElt] =
							MAX(align1InfoPtr->profile[iElt],
								SAFE_SUM(scoreDeletion, iElt * gapExtension));
					}
				}
			}
		}

		if (!localAlignment) {
			/* Step 3g:  Get the optimal score for non-local alignments */
			align1InfoPtr->startRange = 1;
			align2InfoPtr->startRange = 1;
			maxScore =
				MAX(CURR_MATRIX(0, nCharString2),
				MAX(CURR_MATRIX(1, nCharString2),
				    CURR_MATRIX(2, nCharString2)));

			if (align2InfoPtr->endGap) {
				align2InfoPtr->profile[0] =
					MAX(CURR_MATRIX(0, nCharString2), CURR_MATRIX(1, nCharString2));
				for (j = 1, jElt = nCharString2Minus1; j < nCharString2; j++, jElt--) {
					double profile2Score = MAX(CURR_MATRIX(0, j), CURR_MATRIX(1, j));
					align2InfoPtr->profile[jElt] =
						SAFE_SUM(profile2Score, gapOpening + jElt * gapExtension);
					if (profile2Score >= CURR_MATRIX(2, j)) {
						align2InfoPtr->profile[jElt] =
							MAX(align2InfoPtr->profile[jElt],
								SAFE_SUM(CURR_MATRIX(2, j), jElt * gapExtension));
					}
				}
			} else {
				for (j = 1, jElt = nCharString2Minus1; j <= nCharString2; j++, jElt--) {
					align2InfoPtr->profile[jElt] = MAX(CURR_MATRIX(0, j), CURR_MATRIX(1, j));
				}
			}

			if (MAX(CURR_MATRIX(0, nCharString2),
					CURR_MATRIX(2, nCharString2)) >=
						CURR_MATRIX(1, nCharString2)) {
				align1InfoPtr->profile[0] = maxScore;
			}

			if (MAX(CURR_MATRIX(0, nCharString2),
					CURR_MATRIX(1, nCharString2)) >=
						CURR_MATRIX(2, nCharString2)) {
				align2InfoPtr->profile[0] = maxScore;
			}
		}

		/* Step 4:  Traceback through the score matrix */
		char previousAction = '?';
		i = nCharString1Plus1 - align1InfoPtr->startRange;
		j = nCharString2Plus1 - align2InfoPtr->startRange;
		while (TRACE_MATRIX(i, j) != TERMINATION) {
			char action = TRACE_MATRIX(i, j);
			switch (action) {
		    	case DELETION:
		    		if (j == nCharString2) {
		    			align1InfoPtr->startRange++;
		    		} else {
		    			align1InfoPtr->widthRange++;
			    		if (previousAction != DELETION) {
							align2InfoPtr->startInserts--;
							align2InfoPtr->widthInserts--;
							align2InfoPtr->lengthInserts++;
							*align2InfoPtr->startInserts = nCharString2Plus1 - j;
			    		}
			    		*align2InfoPtr->widthInserts += 1;
		    		}
		    		i--;
		    		break;
		    	case INSERTION:
		    		if (i == nCharString1) {
		    			align2InfoPtr->startRange++;
		    		} else {
		    			align2InfoPtr->widthRange++;
			    		if (previousAction != INSERTION) {
							align1InfoPtr->startInserts--;
							align1InfoPtr->widthInserts--;
							align1InfoPtr->lengthInserts++;
							*align1InfoPtr->startInserts = nCharString1Plus1 - i;
			    		}
			    		*align1InfoPtr->widthInserts += 1;
		    		}
		    		j--;
		    		break;
		    	case SUBSTITUTION:
		    		align1InfoPtr->widthRange++;
		    		align2InfoPtr->widthRange++;
		    		i--;
		    		j--;
		    		break;
		    	default:
		    		error("unknown traceback code %d", action);
		    		break;
			}
			previousAction = action;
		}

		align1InfoPtr->profile[align1InfoPtr->startRange - 1] = maxScore;
		align2InfoPtr->profile[align2InfoPtr->startRange - 1] = maxScore;

		int offset1 = align1InfoPtr->startRange - 1;
		if (offset1 > 0 && align1InfoPtr->lengthInserts > 0) {
			for (i = 0; i < align1InfoPtr->lengthInserts; i++)
				align1InfoPtr->startInserts[i] -= offset1;
		}
		int offset2 = align2InfoPtr->startRange - 1;
		if (offset2 > 0 && align2InfoPtr->lengthInserts > 0) {
			for (j = 0; j < align2InfoPtr->lengthInserts; j++)
				align2InfoPtr->startInserts[j] -= offset2;
		}
	}

	return maxScore;
}

/*
 * INPUTS
 * 'pattern':                XString object for patterns
 * 'subject':                XString object for subject
 * 'patternQuality':         BString object for quality scores for pattern
 * 'subjectQuality':         BString object for quality scores for subject
 * 'patternProfile':         XNumeric object for profile of pattern
 * 'subjectProfile':         XNumeric object for profile of subject
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
SEXP XString_align_pairwiseAlignment(
		SEXP pattern,
		SEXP subject,
		SEXP patternQuality,
		SEXP subjectQuality,
		SEXP patternProfile,
		SEXP subjectProfile,
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
	align1Info.string = _get_XString_asRoSeq(pattern);
	align1Info.quality = _get_XString_asRoSeq(patternQuality);
	align1Info.profile = REAL(R_ExternalPtrTag(patternProfile));
	align1Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == SUBJECT_OVERLAP_ALIGNMENT);

	align2Info.string = _get_XString_asRoSeq(subject);
	align2Info.quality = _get_XString_asRoSeq(subjectQuality);
	align2Info.profile = REAL(R_ExternalPtrTag(subjectProfile));
	align2Info.endGap =
		(INTEGER(typeCode)[0] == GLOBAL_ALIGNMENT || INTEGER(typeCode)[0] == PATTERN_OVERLAP_ALIGNMENT);

	/* Perform pairwise alignment */
	double score = pairwiseAlignment(
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

	/* Prepare the output */
	SEXP output, outputNames, outputElement1, outputElement2;
	if (scoreOnlyValue) {
		PROTECT(output = NEW_NUMERIC(1));
		REAL(output)[0] = score;
		UNPROTECT(1);
	} else {
		PROTECT(output = NEW_LIST(9));
		/* Set the names */
		PROTECT(outputNames = NEW_CHARACTER(9));
		SET_STRING_ELT(outputNames, 0, mkChar("startPatternRange"));
		SET_STRING_ELT(outputNames, 1, mkChar("widthPatternRange"));
		SET_STRING_ELT(outputNames, 2, mkChar("startPatternInserts"));
		SET_STRING_ELT(outputNames, 3, mkChar("widthPatternInserts"));
		SET_STRING_ELT(outputNames, 4, mkChar("startSubjectRange"));
		SET_STRING_ELT(outputNames, 5, mkChar("widthSubjectRange"));
		SET_STRING_ELT(outputNames, 6, mkChar("startSubjectInserts"));
		SET_STRING_ELT(outputNames, 7, mkChar("widthSubjectInserts"));
		SET_STRING_ELT(outputNames, 8, mkChar("score"));
		SET_NAMES(output, outputNames);
		UNPROTECT(1);
		/* Set the "startPatternRange" element */
		PROTECT(outputElement1 = NEW_INTEGER(1));
		INTEGER(outputElement1)[0] = align1Info.startRange;
		SET_ELEMENT(output, 0, outputElement1);
		UNPROTECT(1);
		/* Set the "widthPatternRange" element */
		PROTECT(outputElement1 = NEW_INTEGER(1));
		INTEGER(outputElement1)[0] = align1Info.widthRange;
		SET_ELEMENT(output, 1, outputElement1);
		UNPROTECT(1);
		/* Set the "startPatternInserts" and "widthPatternInserts" elements */
		PROTECT(outputElement1 = NEW_INTEGER(align1Info.lengthInserts));
		PROTECT(outputElement2 = NEW_INTEGER(align1Info.lengthInserts));
		for(int i = 0; i < align1Info.lengthInserts; i++) {
			int infoIndex = align1Info.lengthInserts - 1 - i;
			INTEGER(outputElement1)[i] = align1Info.startInserts[infoIndex];
			INTEGER(outputElement2)[i] = align1Info.widthInserts[infoIndex];
		}
		SET_ELEMENT(output, 2, outputElement1);
		SET_ELEMENT(output, 3, outputElement2);
		UNPROTECT(2);
		/* Set the "startSubjectRange" element */
		PROTECT(outputElement1 = NEW_INTEGER(1));
		INTEGER(outputElement1)[0] = align2Info.startRange;
		SET_ELEMENT(output, 4, outputElement1);
		UNPROTECT(1);
		/* Set the "widthSubjectRange" element */
		PROTECT(outputElement1 = NEW_INTEGER(1));
		INTEGER(outputElement1)[0] = align2Info.widthRange;
		SET_ELEMENT(output, 5, outputElement1);
		UNPROTECT(1);
		/* Set the "startSubjectInserts" and "widthSubjectInserts" elements */
		PROTECT(outputElement1 = NEW_INTEGER(align2Info.lengthInserts));
		PROTECT(outputElement2 = NEW_INTEGER(align2Info.lengthInserts));
		for(int j = 0; j < align2Info.lengthInserts; j++) {
			int infoIndex = align2Info.lengthInserts - 1 - j;
			INTEGER(outputElement1)[j] = align2Info.startInserts[infoIndex];
			INTEGER(outputElement2)[j] = align2Info.widthInserts[infoIndex];
		}
		SET_ELEMENT(output, 6, outputElement1);
		SET_ELEMENT(output, 7, outputElement2);
		UNPROTECT(2);
		/* Set the "score" element */
		PROTECT(outputElement1 = NEW_NUMERIC(1));
		REAL(outputElement1)[0] = score;
		SET_ELEMENT(output, 8, outputElement1);
		UNPROTECT(1);
		/* Output is ready */
		UNPROTECT(1);
	}

	return output;
}



SEXP XStringSet_align_pairwiseAlignment(
		SEXP pattern,
		SEXP subject,
		SEXP patternQuality,
		SEXP subjectQuality,
		SEXP typeCode,
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
	PROTECT(output = NEW_NUMERIC(numberOfStrings));
	double *score = REAL(output);

	int qualityIncrement = ((_get_XStringSet_length(patternQuality) < numberOfStrings) ? 0 : 1);
	int qualityElement = 0;
	for(int i = 0; i < numberOfStrings; i++) {
		align1Info.string = _get_CachedXStringSet_elt_asRoSeq(&cachedPattern, i);
		align1Info.quality = _get_CachedXStringSet_elt_asRoSeq(&cachedPatternQuality, qualityElement);
		score[i] = pairwiseAlignment(
				&align1Info,
				&align2Info,
				(INTEGER(typeCode)[0] == LOCAL_ALIGNMENT),
				1,
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

	return output;
}
