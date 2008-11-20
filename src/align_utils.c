#include "Biostrings.h"
#include "IRanges_interface.h"


const char* get_qualityless_classname(SEXP object)
{
	const char *classname = get_classname(object);
	const char *outputClassname;
	if (strcmp(classname, "QualityScaledBStringSet") == 0) {
		outputClassname = "BStringSet";
	} else if (strcmp(classname, "QualityScaledDNAStringSet") == 0) {
		outputClassname = "DNAStringSet";
	} else if (strcmp(classname, "QualityScaledRNAStringSet") == 0) {
		outputClassname = "RNAStringSet";
	} else {
		outputClassname = classname;
	}
	return outputClassname;
}


/*
 * --- .Call ENTRY POINT ---
 */
SEXP PairwiseAlignment_nmatch(SEXP nchar, SEXP nmismatch, SEXP ninsertion,
                SEXP ndeletion)
{
	int ans_len, i, *ans_elt;
	const int *nchar_elt, *nmismatch_elt, *ninsertion_elt, *ndeletion_elt;
	SEXP ans;

	ans_len = LENGTH(nchar);
	PROTECT(ans = NEW_INTEGER(ans_len));
	for (i = 0, nchar_elt = INTEGER(nchar), nmismatch_elt = INTEGER(nmismatch),
	     ninsertion_elt = INTEGER(ninsertion), ndeletion_elt = INTEGER(ndeletion),
	     ans_elt = INTEGER(ans);
	     i < ans_len;
	     i++, nchar_elt++, nmismatch_elt++, ninsertion_elt++, ndeletion_elt++, ans_elt++)
	{
		*ans_elt = *nchar_elt - *nmismatch_elt - *ninsertion_elt - *ndeletion_elt;
	}
	UNPROTECT(1);
	return ans;
}

SEXP AlignedXStringSet_nchar(SEXP alignedXStringSet)
{
	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	SEXP indel = GET_SLOT(GET_SLOT(alignedXStringSet, install("indel")), install("elements"));
	int numberOfAlignments = LENGTH(indel);

	SEXP output;
	PROTECT(output = NEW_INTEGER(numberOfAlignments));
	int i, j, *outputPtr;
	const int *rangeWidth, *indelWidth;
	for (i = 0, rangeWidth = INTEGER(get_IRanges_width(range)), outputPtr = INTEGER(output);
			i < numberOfAlignments; i++, rangeWidth++, outputPtr++) {
		SEXP indelElement = VECTOR_ELT(indel, i);
		int numberOfIndels = LENGTH(get_IRanges_width(indelElement));
		*outputPtr = *rangeWidth;
		for (j = 0, indelWidth = INTEGER(get_IRanges_width(indelElement));
		     j < numberOfIndels; j++, indelWidth++) {
			*outputPtr += *indelWidth;
		}
	}
	UNPROTECT(1);

	return output;
}


SEXP AlignedXStringSet_align_aligned(SEXP alignedXStringSet, SEXP gapCode)
{
	int i, j;
	char gapCodeValue = (char) RAW(gapCode)[0];

	SEXP unaligned = GET_SLOT(alignedXStringSet, install("unaligned"));
	CachedXStringSet cachedAlignedXStringSet = _new_CachedXStringSet(unaligned);
	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	SEXP indel = GET_SLOT(GET_SLOT(alignedXStringSet, install("indel")), install("elements"));

	const char *stringSetClass = get_qualityless_classname(unaligned);
	const char *stringClass = get_classname(_get_XStringSet_super(unaligned));

	int numberOfStrings = _get_XStringSet_length(unaligned);
	int numberOfAlignments = LENGTH(indel);

	SEXP output;

	SEXP alignedString, alignedRanges, alignedStart, alignedWidth;
	PROTECT(alignedWidth = AlignedXStringSet_nchar(alignedXStringSet));
	PROTECT(alignedStart = NEW_INTEGER(LENGTH(alignedWidth)));
	int totalNChars = 0;
	const int *width_i, *start_i;
	int *start_iPlus1;
	for (i = 0, width_i = INTEGER(alignedWidth); i < LENGTH(alignedWidth); i++, width_i++) {
		totalNChars += *width_i;
	}
	if (totalNChars > 0) {
		INTEGER(alignedStart)[0] = 1;
		for (i = 0, start_i = INTEGER(alignedStart), width_i = INTEGER(alignedWidth),
				start_iPlus1 = INTEGER(alignedStart) + 1; i < LENGTH(alignedWidth) - 1;
				i++, start_i++, width_i++, start_iPlus1++) {
			*start_iPlus1 = *start_i + *width_i;
		}
	}
	SEXP alignedStringTag, alignedStringXData;
	PROTECT(alignedStringTag = NEW_RAW(totalNChars));
	PROTECT(alignedStringXData = new_SequencePtr("RawPtr", alignedStringTag));
	PROTECT(alignedString = new_XSequence(stringClass, alignedStringXData, 0, LENGTH(alignedStringTag)));
	PROTECT(alignedRanges = new_IRanges("IRanges", alignedStart, alignedWidth, R_NilValue));
	char *alignedStringPtr = (char *) RAW(alignedStringTag);
	PROTECT(output = _new_XStringSet(stringSetClass, alignedString, alignedRanges));

	int stringIncrement = (numberOfStrings == 1 ? 0 : 1);
	int index = 0, stringElement = 0;
	const int *rangeStart, *rangeWidth;
	for (i = 0, rangeStart = INTEGER(get_IRanges_start(range)), rangeWidth = INTEGER(get_IRanges_width(range));
	         i < numberOfAlignments; i++, rangeStart++, rangeWidth++) {
		RoSeq origString = _get_CachedXStringSet_elt_asRoSeq(&cachedAlignedXStringSet, stringElement);
		char *origStringPtr = (char *) (origString.elts + (*rangeStart - 1));
		SEXP indelElement = VECTOR_ELT(indel, i);
		int numberOfIndel = LENGTH(get_IRanges_start(indelElement));
		if (numberOfIndel == 0) {
			memcpy(&alignedStringPtr[index], origStringPtr, *rangeWidth * sizeof(char));
			index += *rangeWidth;
		} else {
			int prevStart = 0;
			const int *indelStart, *indelWidth;
			for (j = 0, indelStart = INTEGER(get_IRanges_start(indelElement)),
					indelWidth = INTEGER(get_IRanges_width(indelElement));
			        j < numberOfIndel; j++, indelStart++, indelWidth++) {
				int currStart = *indelStart - 1;
				int currWidth = *indelWidth;
				int copyElements = currStart - prevStart;
				if (copyElements > 0) {
					memcpy(&alignedStringPtr[index], origStringPtr, copyElements * sizeof(char));
					index += copyElements;
					origStringPtr += copyElements;
				}
				for (int k = 0; k < currWidth; k++) {
					alignedStringPtr[index] = gapCodeValue;
					index++;
				}
				prevStart = currStart;
			}
			int copyElements = *rangeWidth - prevStart;
			memcpy(&alignedStringPtr[index], origStringPtr, copyElements * sizeof(char));
			index += copyElements;
		}
		stringElement += stringIncrement;
	}
	UNPROTECT(7);

	return output;
}


SEXP PairwiseAlignment_align_aligned(SEXP alignment, SEXP gapCode)
{
	int i, j;
	char gapCodeValue = (char) RAW(gapCode)[0];

	SEXP pattern = GET_SLOT(alignment, install("pattern"));
	SEXP unalignedPattern = GET_SLOT(pattern, install("unaligned"));
	CachedXStringSet cachedUnalignedPattern = _new_CachedXStringSet(unalignedPattern);
	SEXP rangePattern = GET_SLOT(pattern, install("range"));
	SEXP namesPattern = get_IRanges_names(rangePattern);
	SEXP indelPattern = GET_SLOT(GET_SLOT(pattern, install("indel")), install("elements"));

	SEXP subject = GET_SLOT(alignment, install("subject"));
	SEXP rangeSubject = GET_SLOT(subject, install("range"));
	SEXP indelSubject = GET_SLOT(GET_SLOT(subject, install("indel")), install("elements"));

	const char *stringSetClass = get_qualityless_classname(unalignedPattern);
	const char *stringClass = get_classname(_get_XStringSet_super(unalignedPattern));

	int numberOfAlignments = LENGTH(indelPattern);
	int numberOfChars = INTEGER(_get_XStringSet_width(GET_SLOT(subject, install("unaligned"))))[0];

	SEXP output;

	SEXP mappedString, mappedRanges, mappedStart, mappedWidth;
	PROTECT(mappedWidth = NEW_INTEGER(numberOfAlignments));
	PROTECT(mappedStart = NEW_INTEGER(numberOfAlignments));
	int totalNChars = numberOfAlignments * numberOfChars;
	int *width_i, *start_i;
	if (totalNChars > 0) {
		for (i = 0, start_i = INTEGER(mappedStart), width_i = INTEGER(mappedWidth);
				i < numberOfAlignments; i++, start_i++, width_i++) {
			*start_i = i * numberOfChars + 1;
			*width_i = numberOfChars;
		}
	}
	SEXP mappedStringTag, mappedStringXData;
	PROTECT(mappedStringTag = NEW_RAW(totalNChars));
	PROTECT(mappedStringXData = new_SequencePtr("RawPtr", mappedStringTag));
	PROTECT(mappedString = new_XSequence(stringClass, mappedStringXData, 0, LENGTH(mappedStringTag)));
	PROTECT(mappedRanges = new_IRanges("IRanges", mappedStart, mappedWidth, namesPattern));
	char *mappedStringPtr = (char *) RAW(mappedStringTag);
	PROTECT(output = _new_XStringSet(stringSetClass, mappedString, mappedRanges));

	int index = 0;
	const int *rangeStartPattern, *rangeWidthPattern, *rangeStartSubject, *rangeWidthSubject;
	for (i = 0,
			rangeStartPattern = INTEGER(get_IRanges_start(rangePattern)),
			rangeWidthPattern = INTEGER(get_IRanges_width(rangePattern)),
			rangeStartSubject = INTEGER(get_IRanges_start(rangeSubject)),
			rangeWidthSubject = INTEGER(get_IRanges_width(rangeSubject));
	         i < numberOfAlignments;
	         i++, rangeStartPattern++, rangeWidthPattern++, rangeStartSubject++, rangeWidthSubject++) {
		RoSeq origString = _get_CachedXStringSet_elt_asRoSeq(&cachedUnalignedPattern, i);
		char *origStringPtr = (char *) (origString.elts + (*rangeStartPattern - 1));
		SEXP indelElementPattern = VECTOR_ELT(indelPattern, i);
		SEXP indelElementSubject = VECTOR_ELT(indelSubject, i);
		int numberOfIndelPattern = LENGTH(get_IRanges_start(indelElementPattern));
		int numberOfIndelSubject = LENGTH(get_IRanges_start(indelElementSubject));

		for (j = 0; j < *rangeStartSubject - 1; j++) {
			mappedStringPtr[index] = gapCodeValue;
			index++;
		}
		int jPattern = 1;
		const int *indelStartPattern, *indelWidthPattern, *indelStartSubject, *indelWidthSubject;
		if (numberOfIndelPattern > 0) {
			indelStartPattern = INTEGER(get_IRanges_start(indelElementPattern));
			indelWidthPattern = INTEGER(get_IRanges_width(indelElementPattern));
		}
		if (numberOfIndelSubject > 0) {
			indelStartSubject = INTEGER(get_IRanges_start(indelElementSubject));
			indelWidthSubject = INTEGER(get_IRanges_width(indelElementSubject));
		}
		for (j = 1; j <= *rangeWidthSubject; j++) {
			if ((numberOfIndelSubject == 0) || (j < *indelStartSubject)) {
				if ((numberOfIndelPattern == 0) || (jPattern < *indelStartPattern)) {
					mappedStringPtr[index] = *origStringPtr;
					index++;
					origStringPtr++;
					jPattern++;
				} else {
					for (int k = 0; k < *indelWidthPattern; k++) {
						mappedStringPtr[index] = gapCodeValue;
						index++;
					}
					j += *indelWidthPattern - 1;
					indelStartPattern++;
					indelWidthPattern++;
					numberOfIndelPattern--;
				}
			} else {
				origStringPtr += *indelWidthSubject;
				jPattern += *indelWidthSubject;
				j--;
				indelStartSubject++;
				indelWidthSubject++;
				numberOfIndelSubject--;
			}
		}
		for (j = *rangeStartSubject + (*rangeWidthSubject - 1); j < numberOfChars; j++) {
			mappedStringPtr[index] = gapCodeValue;
			index++;
		}
	}
	UNPROTECT(7);

	return(output);
}


SEXP align_compareStrings(SEXP patternStrings, SEXP subjectStrings, SEXP maxNChar,
                          SEXP insertionCode, SEXP deletionCode, SEXP mismatchCode)
{
	char insertionChar = CHAR(STRING_ELT(insertionCode, 0))[0];
	char deletionChar = CHAR(STRING_ELT(deletionCode, 0))[0];
	char mismatchChar = CHAR(STRING_ELT(mismatchCode, 0))[0];
	int numberOfStrings = LENGTH(patternStrings);
	char *outputPtr = (char *) R_alloc((long) (INTEGER(maxNChar)[0] + 1), sizeof(char));
	SEXP output;
	PROTECT(output = NEW_CHARACTER(numberOfStrings));
	int i, j;
	char *output_j;
	const char *subject_j;
	for (i = 0; i < numberOfStrings; i++) {
		const char *patternPtr = (char *) CHAR(STRING_ELT(patternStrings, i));
		const char *subjectPtr = (char *) CHAR(STRING_ELT(subjectStrings, i));
		int numberOfChars = strlen(patternPtr);
		memcpy(outputPtr, patternPtr, numberOfChars * sizeof(char));
		outputPtr[numberOfChars] = '\0';
		for (j = 0, output_j = outputPtr, subject_j = subjectPtr;
		     j < numberOfChars; j++, output_j++, subject_j++) {
			if (*output_j != deletionChar) {
				if (*subject_j == deletionChar) {
					*output_j = insertionChar;
				} else if (*subject_j != *output_j) {
					*output_j = mismatchChar;
				}
			}
		}
		SET_STRING_ELT(output, i, mkChar(outputPtr));
	}
	UNPROTECT(1);
	return(output);
}
