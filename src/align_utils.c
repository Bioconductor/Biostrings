#include "Biostrings.h"


SEXP AlignedXStringSet_nchar(SEXP alignedXStringSet)
{
	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	int *rangeWidth = INTEGER(_get_IRanges_width(range));
	SEXP indels = GET_SLOT(alignedXStringSet, install("indels"));
	int numberOfAlignments = LENGTH(indels);

	SEXP output;
	PROTECT(output = NEW_INTEGER(numberOfAlignments));
	int *outputPtr = INTEGER(output);
	for (int i = 0; i < numberOfAlignments; i++) {
		SEXP indelsElement = VECTOR_ELT(indels, i);
		int numberOfIndels = LENGTH(_get_IRanges_width(indelsElement));
		int *indelsWidth = INTEGER(_get_IRanges_width(indelsElement));
		outputPtr[i] = rangeWidth[i];
		for (int j = 0; j < numberOfIndels; j++) {
			outputPtr[i] += indelsWidth[j];
		}
	}
	UNPROTECT(1);

	return output;
}


SEXP AlignedXStringSet_align_aligned(SEXP alignedXStringSet, SEXP gapCode)
{
	char gapCodeValue = (char) RAW(gapCode)[0];

	SEXP unaligned = GET_SLOT(alignedXStringSet, install("unaligned"));
	CachedXStringSet cachedAlignedXStringSet = _new_CachedXStringSet(unaligned);
	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	int *rangeStart = INTEGER(_get_IRanges_start(range));
	int *rangeWidth = INTEGER(_get_IRanges_width(range));
	SEXP indels = GET_SLOT(alignedXStringSet, install("indels"));

	const char *stringSetClass = _get_class(unaligned);
	const char *stringClass = _get_class(GET_SLOT(unaligned, install("super")));

	int numberOfStrings = _get_XStringSet_length(unaligned);
	int numberOfAlignments = LENGTH(indels);

	SEXP output;
	PROTECT(output = NEW_OBJECT(MAKE_CLASS(stringSetClass)));

	SEXP alignedString, alignedStart, alignedWidth;
	PROTECT(alignedWidth = AlignedXStringSet_nchar(alignedXStringSet));
	PROTECT(alignedStart = NEW_INTEGER(LENGTH(alignedWidth)));
	int totalNChars = 0;
	for (int i = 0; i < LENGTH(alignedWidth); i++) {
		totalNChars += INTEGER(alignedWidth)[i];
	}
	if (totalNChars > 0) {
		INTEGER(alignedStart)[0] = 1;
		for (int i = 0; i < LENGTH(alignedWidth) - 1; i++) {
			INTEGER(alignedStart)[i+1] = INTEGER(alignedStart)[i] + INTEGER(alignedWidth)[i];
		}
	}
	SEXP alignedStringTag, alignedStringXData;
	PROTECT(alignedStringTag = NEW_RAW(totalNChars));
	PROTECT(alignedStringXData = _new_XRaw(alignedStringTag));
	PROTECT(alignedString = _new_XString(stringClass, alignedStringXData, 0, LENGTH(alignedStringTag)));
	char *alignedStringPtr = (char *) RAW(alignedStringTag);
	SET_SLOT(output, mkChar("super"), alignedString);
	SET_SLOT(output, mkChar("start"), alignedStart);
	SET_SLOT(output, mkChar("width"), alignedWidth);

	int stringIncrement = (numberOfStrings == 1 ? 0 : 1);
	int index = 0, stringElement = 0;
	for (int i = 0; i < numberOfAlignments; i++) {
		RoSeq origString = _get_CachedXStringSet_elt_asRoSeq(&cachedAlignedXStringSet, stringElement);
		char *origStringPtr = (char *) (origString.elts + (rangeStart[i] - 1));
		SEXP indelsElement = VECTOR_ELT(indels, i);
		int numberOfIndels = LENGTH(_get_IRanges_start(indelsElement));
		int *indelsStart = INTEGER(_get_IRanges_start(indelsElement));
		int *indelsWidth = INTEGER(_get_IRanges_width(indelsElement));
		if (numberOfIndels == 0) {
			memcpy(&alignedStringPtr[index], origStringPtr, rangeWidth[i] * sizeof(char));
			index += rangeWidth[i];
		} else {
			int prevStart = 0;
			for (int j = 0; j < numberOfIndels; j++) {
				int currStart = indelsStart[j] - 1;
				int currWidth = indelsWidth[j];
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
			int copyElements = rangeWidth[i] - prevStart;
			memcpy(&alignedStringPtr[index], origStringPtr, copyElements * sizeof(char));
			index += copyElements;
		}
		stringElement += stringIncrement;
	}
	UNPROTECT(6);

	return output;
}


void align_compareStrings(SEXP patternStrings, SEXP subjectStrings,
                          SEXP insertionCode, SEXP deletionCode, SEXP mismatchCode)
{
	char insertionChar = CHAR(STRING_ELT(insertionCode, 0))[0];
	char deletionChar = CHAR(STRING_ELT(deletionCode, 0))[0];
	char mismatchChar = CHAR(STRING_ELT(mismatchCode, 0))[0];
	int numberOfStrings = LENGTH(patternStrings);
	for (int i = 0; i < numberOfStrings; i++) {
		char *patternPtr = (char *) CHAR(STRING_ELT(patternStrings, i));
		char *subjectPtr = (char *) CHAR(STRING_ELT(subjectStrings, i));
		int numberOfChars = strlen(patternPtr);
		for (int j = 0; j < numberOfChars; j++) {
			if (patternPtr[j] != deletionChar) {
				if (subjectPtr[j] == deletionChar) {
					patternPtr[j] = insertionChar;
				} else if (subjectPtr[j] != patternPtr[j]) {
					patternPtr[j] = mismatchChar;
				}
			}
		}
		SET_STRING_ELT(patternStrings, i, mkChar(patternPtr));
	}
}
