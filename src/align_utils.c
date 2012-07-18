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
SEXP PairwiseAlignments_nmatch(SEXP nchar, SEXP nmismatch, SEXP ninsertion,
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
	int numberOfAlignments = get_IRanges_length(range);

	SEXP indel = GET_SLOT(alignedXStringSet, install("indel"));
	cachedCompressedIRangesList cached_indel = cache_CompressedIRangesList(indel);

	SEXP output;
	PROTECT(output = NEW_INTEGER(numberOfAlignments));
	int i, j, *outputPtr;
	const int *rangeWidth;
	for (i = 0, rangeWidth = INTEGER(get_IRanges_width(range)), outputPtr = INTEGER(output);
			i < numberOfAlignments; i++, rangeWidth++, outputPtr++) {
		cachedIRanges indelElement = get_cachedCompressedIRangesList_elt(&cached_indel, i);
		int numberOfIndels = get_cachedIRanges_length(&indelElement);
		*outputPtr = *rangeWidth;
		for (j = 0; j < numberOfIndels; j++)
			*outputPtr += get_cachedIRanges_elt_width(&indelElement, j);
	}
	UNPROTECT(1);

	return output;
}


SEXP AlignedXStringSet_align_aligned(SEXP alignedXStringSet, SEXP gapCode)
{
	int i, j;
	char gapCodeValue = (char) RAW(gapCode)[0];

	SEXP unaligned = GET_SLOT(alignedXStringSet, install("unaligned"));
	cachedXStringSet cached_unaligned = _cache_XStringSet(unaligned);

	SEXP range = GET_SLOT(alignedXStringSet, install("range"));
	int numberOfAlignments = get_IRanges_length(range);

	SEXP indel = GET_SLOT(alignedXStringSet, install("indel"));
	cachedCompressedIRangesList cached_indel = cache_CompressedIRangesList(indel);

	const char *stringSetClass = get_qualityless_classname(unaligned);
	const char *stringClass = _get_XStringSet_xsbaseclassname(unaligned);

	int numberOfStrings = _get_XStringSet_length(unaligned);

	SEXP output;

	SEXP alignedRanges, alignedStart, alignedWidth;
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
	SEXP alignedStringTag;
	PROTECT(alignedStringTag = NEW_RAW(totalNChars));
	PROTECT(alignedRanges = new_IRanges("IRanges", alignedStart, alignedWidth, R_NilValue));
	char *alignedStringPtr = (char *) RAW(alignedStringTag);
	PROTECT(output = new_XRawList_from_tag(stringSetClass, stringClass, alignedStringTag, alignedRanges));

	int stringIncrement = (numberOfStrings == 1 ? 0 : 1);
	int index = 0, stringElement = 0;
	const int *rangeStart, *rangeWidth;
	for (i = 0, rangeStart = INTEGER(get_IRanges_start(range)), rangeWidth = INTEGER(get_IRanges_width(range));
	         i < numberOfAlignments; i++, rangeStart++, rangeWidth++) {
		cachedCharSeq origString = _get_cachedXStringSet_elt(&cached_unaligned, stringElement);
		char *origStringPtr = (char *) (origString.seq + (*rangeStart - 1));
		cachedIRanges indelElement = get_cachedCompressedIRangesList_elt(&cached_indel, i);
		int numberOfIndel = get_cachedIRanges_length(&indelElement);
		if (numberOfIndel == 0) {
			memcpy(&alignedStringPtr[index], origStringPtr, *rangeWidth * sizeof(char));
			index += *rangeWidth;
		} else {
			int prevStart = 0;
			for (j = 0; j < numberOfIndel; j++) {
				int currStart = get_cachedIRanges_elt_start(&indelElement, j) - 1;
				int currWidth = get_cachedIRanges_elt_width(&indelElement, j);
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
	UNPROTECT(5);

	return output;
}


SEXP PairwiseAlignmentsSingleSubject_align_aligned(SEXP alignment, SEXP gapCode, SEXP endgapCode)
{
	int i, j;
	char gapCodeValue = (char) RAW(gapCode)[0];
	char endgapCodeValue = (char) RAW(endgapCode)[0];

	SEXP pattern = GET_SLOT(alignment, install("pattern"));
	SEXP unalignedPattern = GET_SLOT(pattern, install("unaligned"));
	cachedXStringSet cached_unalignedPattern = _cache_XStringSet(unalignedPattern);
	SEXP rangePattern = GET_SLOT(pattern, install("range"));
	SEXP namesPattern = get_IRanges_names(rangePattern);
	SEXP indelPattern = GET_SLOT(pattern, install("indel"));
	cachedCompressedIRangesList cached_indelPattern = cache_CompressedIRangesList(indelPattern);

	SEXP subject = GET_SLOT(alignment, install("subject"));
	SEXP rangeSubject = GET_SLOT(subject, install("range"));
	SEXP indelSubject = GET_SLOT(subject, install("indel"));
	cachedCompressedIRangesList cached_indelSubject = cache_CompressedIRangesList(indelSubject);

	const char *stringSetClass = get_qualityless_classname(unalignedPattern);
	const char *stringClass = _get_XStringSet_xsbaseclassname(unalignedPattern);

	int numberOfAlignments = get_IRanges_length(rangePattern);
	int numberOfChars = INTEGER(_get_XStringSet_width(GET_SLOT(subject, install("unaligned"))))[0];

	SEXP output;

	SEXP mappedRanges, mappedStart, mappedWidth;
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
	SEXP mappedStringTag;
	PROTECT(mappedStringTag = NEW_RAW(totalNChars));
	PROTECT(mappedRanges = new_IRanges("IRanges", mappedStart, mappedWidth, namesPattern));
	char *mappedStringPtr = (char *) RAW(mappedStringTag);
	PROTECT(output = new_XRawList_from_tag(stringSetClass, stringClass, mappedStringTag, mappedRanges));

	int index = 0;
	const int *rangeStartPattern, *rangeWidthPattern, *rangeStartSubject, *rangeWidthSubject;
	for (i = 0,
			rangeStartPattern = INTEGER(get_IRanges_start(rangePattern)),
			rangeWidthPattern = INTEGER(get_IRanges_width(rangePattern)),
			rangeStartSubject = INTEGER(get_IRanges_start(rangeSubject)),
			rangeWidthSubject = INTEGER(get_IRanges_width(rangeSubject));
	         i < numberOfAlignments;
	         i++, rangeStartPattern++, rangeWidthPattern++, rangeStartSubject++, rangeWidthSubject++) {
		cachedCharSeq origString = _get_cachedXStringSet_elt(&cached_unalignedPattern, i);
		char *origStringPtr = (char *) (origString.seq + (*rangeStartPattern - 1));
		cachedIRanges indelElementPattern = get_cachedCompressedIRangesList_elt(&cached_indelPattern, i);
		cachedIRanges indelElementSubject = get_cachedCompressedIRangesList_elt(&cached_indelSubject, i);
		int numberOfIndelPattern = get_cachedIRanges_length(&indelElementPattern);
		int numberOfIndelSubject = get_cachedIRanges_length(&indelElementSubject);

		for (j = 0; j < *rangeStartSubject - 1; j++) {
			mappedStringPtr[index] = endgapCodeValue;
			index++;
		}
		int jPattern = 1, jp = 0, js = 0;
		int indelStartPattern, indelWidthPattern, indelStartSubject, indelWidthSubject;
		if (numberOfIndelPattern > 0) {
			indelStartPattern = get_cachedIRanges_elt_start(&indelElementPattern, jp);
			indelWidthPattern = get_cachedIRanges_elt_width(&indelElementPattern, jp);
		}
		if (numberOfIndelSubject > 0) {
			indelStartSubject = get_cachedIRanges_elt_start(&indelElementSubject, js);
			indelWidthSubject = get_cachedIRanges_elt_width(&indelElementSubject, js);
		}
		for (j = 1; j <= *rangeWidthSubject; j++) {
			if ((numberOfIndelSubject == 0) || (j < indelStartSubject)) {
				if ((numberOfIndelPattern == 0) || (jPattern < indelStartPattern)) {
					mappedStringPtr[index] = *origStringPtr;
					index++;
					origStringPtr++;
					jPattern++;
				} else {
					for (int k = 0; k < indelWidthPattern; k++) {
						mappedStringPtr[index] = gapCodeValue;
						index++;
					}
					j += indelWidthPattern - 1;
					jp++;
					indelStartPattern = get_cachedIRanges_elt_start(&indelElementPattern, jp);
					indelWidthPattern = get_cachedIRanges_elt_width(&indelElementPattern, jp);
					numberOfIndelPattern--;
				}
			} else {
				origStringPtr += indelWidthSubject;
				jPattern += indelWidthSubject;
				j--;
				js++;
				indelStartSubject = get_cachedIRanges_elt_start(&indelElementSubject, js);
				indelWidthSubject = get_cachedIRanges_elt_width(&indelElementSubject, js);
				numberOfIndelSubject--;
			}
		}
		for (j = *rangeStartSubject + (*rangeWidthSubject - 1); j < numberOfChars; j++) {
			mappedStringPtr[index] = endgapCodeValue;
			index++;
		}
	}
	UNPROTECT(5);

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
