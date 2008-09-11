#include "Biostrings.h"
#include "IRanges_interface.h"

SEXP AlignedXStringSet_nchar(SEXP alignedXStringSet)
{
	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	SEXP indel = GET_SLOT(alignedXStringSet, install("indel"));
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
	SEXP indel = GET_SLOT(alignedXStringSet, install("indel"));

	const char *stringSetClass = get_class(unaligned);
	const char *stringClass = get_class(GET_SLOT(unaligned, install("super")));

	int numberOfStrings = _get_XStringSet_length(unaligned);
	int numberOfAlignments = LENGTH(indel);

	SEXP output;
	PROTECT(output = NEW_OBJECT(MAKE_CLASS(stringSetClass)));

	SEXP alignedString, alignedStart, alignedWidth;
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
	PROTECT(alignedString = _new_XString(stringClass, alignedStringXData, 0, LENGTH(alignedStringTag)));
	char *alignedStringPtr = (char *) RAW(alignedStringTag);
	SET_SLOT(output, mkChar("super"), alignedString);
	SET_SLOT(output, mkChar("start"), alignedStart);
	SET_SLOT(output, mkChar("width"), alignedWidth);

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
	UNPROTECT(6);

	return output;
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
