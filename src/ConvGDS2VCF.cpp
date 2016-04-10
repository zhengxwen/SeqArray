// ===========================================================
//
// ConvGDS2VCF.cpp: format conversion from GDS to VCF
//
// Copyright (C) 2013-2016    Xiuwen Zheng
//
// This file is part of SeqArray.
//
// SeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.

#include "Common.h"
#include "vectorization.h"
#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;

namespace SeqArray
{

// double quote the text if it is needed
static string QuoteText(const char *p)
{
	string rv;

	rv.clear();
	bool flag = false;
	for (; *p != 0; p ++)
	{
		switch (*p)
		{
			case ',': case ';':
				flag = true; rv.push_back(*p); break;
			case '\"':
				flag = true; rv.append("\\\""); break;
			case '\'':
				flag = true; rv.append("\\\'"); break;
			case ' ':
				flag = true; rv.push_back(' '); break;
			default:
				rv.push_back(*p);
		}
	}
	if (flag) // add double quote
	{
		rv.insert(0, "\"");
		rv.push_back('\"');
	}

	return rv;
}



// ========================================================================

/// used in SEQ_OutVCF4
static vector<int> VCF_INFO_Number;    ///<
static vector<int> VCF_FORMAT_Number;  ///<
static vector<SEXP> VCF_FORMAT_List;

static size_t VCF_NumAllele;  ///< the number of alleles
static size_t VCF_NumSample;  ///< the number of samples
static Rconnection VCF_File = NULL;  ///< R connection object


static const size_t LINE_BUFFER_SIZE = 4096;
static vector<char> LineBuffer;
static char *LineBegin = NULL;
static char *LinePtr   = NULL;
static char *LineEnd   = NULL;


inline static void LineBuf_Init()
{
	LineBuffer.resize(LINE_BUFFER_SIZE);
	LinePtr = LineBegin = &LineBuffer[0];
	LineEnd = LinePtr + LINE_BUFFER_SIZE;
}

inline static void LineBuf_Done()
{
	LineBuffer.clear();
	LineBegin = LinePtr = LineEnd = NULL;
	VCF_FORMAT_List.clear();
}

inline static void LineBuf_InitPtr()
{
	LinePtr = LineBegin = &LineBuffer[0];
}

inline static void LineBuf_NeedSize(size_t st)
{
	if (LinePtr + st > LineEnd)
	{
		size_t p = LinePtr - LineBegin;
		size_t n = p + st;
		n = (n / LINE_BUFFER_SIZE + 1) * LINE_BUFFER_SIZE;
		LineBuffer.resize(n);
		LineBegin = &LineBuffer[0];
		LinePtr = LineBegin + p;
		LineEnd = LineBegin + n;
	}
}

inline static char *fast_itoa(char *p, int32_t val)
{
	static int base[10] = {
		1000000000, 100000000, 10000000, 1000000, 100000,
		10000, 1000, 100, 10, 1
	};
	if (val < 0)
		{ *p++ = '-'; val = -val; }

	size_t n=9;
	for (int *b = base+8; n > 0; n--)
		if (val < (*b--)) break;

	for (int *b = base+n; n < 9; n++, b++)
	{
		*p++ = val / (*b) + '0';
		val %= (*b);
	}

	*p++ = val + '0';
	return p;
}

inline static void _Line_Append(int val)
{
	if (val != NA_INTEGER)
		LinePtr = fast_itoa(LinePtr, val);
	else
		*LinePtr++ = '.';
}

inline static void _Line_Append_Geno(int val)
{
	if (val >= 0)
	{
		if (val < 10)
			*LinePtr++ = val + '0';
		else
			LinePtr = fast_itoa(LinePtr, val);
	} else 
		*LinePtr++ = '.';
}

inline static void _Line_Append(double val)
{
	if (R_FINITE(val))
		LinePtr += sprintf(LinePtr, "%0.6g", val);
	else
		*LinePtr++ = '.';
}


inline static void LineBuf_Append(int val)
{
	LineBuf_NeedSize(32);
	_Line_Append(val);
}

inline static void LineBuf_Append(double val)
{
	LineBuf_NeedSize(32);
	_Line_Append(val);
}

inline static void LineBuf_Append(const char *txt)
{
	const size_t n = strlen(txt);
	LineBuf_NeedSize(n + 16);
	memcpy(LinePtr, txt, n);
	LinePtr += n;
}


inline static void put_text(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*VCF_File->vfprintf)(VCF_File, fmt, args);
	va_end(args);
}



// ========================================================================

/// return the number in the INFO field
inline static int INFO_GetNum(SEXP X, int n)
{
	if (n < 0) n = Rf_length(X);

	if (IS_INTEGER(X))
	{
		int *p = INTEGER(X) + n;
		for (; n > 0; n--)
			if (*(--p) != NA_INTEGER) break;
	} else if (IS_NUMERIC(X))
	{
		double *p = REAL(X) + n;
		for (; n > 0; n--)
			if (R_finite(*(--p))) break;
	} else if (IS_CHARACTER(X))
	{
		for (; n > 0; n--)
		{
			SEXP s = STRING_ELT(X, n-1);
			if ((s != NA_STRING) && (CHAR(s)[0] != 0)) break;
		}
	}

	return n;
}


/// write info values
inline static void INFO_Write(SEXP X, size_t n)
{
	size_t i = 0;
	if (IS_INTEGER(X))
	{
		LineBuf_NeedSize(12*n + 32);
		for (int *p = INTEGER(X); i < n; i++)
		{
			if (i > 0) *LinePtr++ = ',';
			_Line_Append(*p++);
		}
	} else if (IS_NUMERIC(X))
	{
		LineBuf_NeedSize(12*n + 32);
		for (double *p = REAL(X); i < n; i++)
		{
			if (i > 0) *LinePtr++ = ',';
			_Line_Append(*p++);
		}
	} else if (IS_CHARACTER(X))
	{
		for (; i < n; i++)
		{
			if (i > 0) *LinePtr++ = ',';
			LineBuf_Append(CHAR(STRING_ELT(X, i)));
		}
	}
}


/// write format values
inline static void FORMAT_Write(SEXP X, size_t n, size_t Start, size_t Step)
{
	if (IS_INTEGER(X) || IS_LOGICAL(X))
	{
		int *base = (IS_INTEGER(X) ? INTEGER(X) : LOGICAL(X)) + Start;
		int *p = base + (n - 1)*Step;
		for (; n > 0; n--, p-=Step)
			if (*p != NA_INTEGER) break;
		LineBuf_NeedSize(12*n + 32);
		p = base;
		for (size_t i=0; i < n; i++)
		{
			if (i > 0) *LinePtr++ = ',';
			_Line_Append(*p);
			p += Step;
		}
	} else if (IS_NUMERIC(X))
	{
		double *base = REAL(X) + Start;
		double *p = base + (n - 1)*Step;
		for (; n > 0; n--, p-=Step)
			if (R_finite(*p)) break;
		LineBuf_NeedSize(12*n + 32);
		p = base;
		for (size_t i=0; i < n; i++)
		{
			if (i > 0) *LinePtr++ = ',';
			_Line_Append(*p);
			p += Step;
		}
	} else if (IS_CHARACTER(X) || Rf_isFactor(X))
	{
		if (Rf_isFactor(X))
			X = Rf_asCharacterFactor(X);
		for (; n > 0; n--)
		{
			SEXP s = STRING_ELT(X, Start + (n - 1)*Step);
			if ((s != NA_STRING) && (CHAR(s)[0] != 0)) break;
		}
		for (size_t i=0; i < n; i++, Start += Step)
		{
			if (i > 0) *LinePtr++ = ',';
			SEXP s = STRING_ELT(X, Start);
			if (s != NA_STRING)
				LineBuf_Append(CHAR(s));
			else
				*LinePtr++ = '.';
		}
	}

	if (n <= 0) *LinePtr++ = '.';
}





// ========================================================================

/// export the first seven columns: chr, pos, id, allele (REF/ALT), qual, filter
inline static void ExportHead(SEXP X)
{
	// CHROM
	LineBuf_Append(CHAR(STRING_ELT(VECTOR_ELT(X, 0), 0)));
	*LinePtr++ = '\t';

	// POS
	LineBuf_Append(Rf_asInteger(VECTOR_ELT(X, 1)));
	*LinePtr++ = '\t';

	// ID
	LineBuf_Append(CHAR(STRING_ELT(VECTOR_ELT(X, 2), 0)));
	*LinePtr++ = '\t';

	// allele -- REF/ALT
	size_t n = LinePtr - LineBegin;
	LineBuf_Append(CHAR(STRING_ELT(VECTOR_ELT(X, 3), 0)));

	char *s = LineBegin + n;
	for (; s < LinePtr; s++)
	{
		if (*s == ',')
			{ *s = '\t'; break; }
	}
	if (s == LinePtr)
	{
		*LinePtr++ = '\t';
		*LinePtr++ = '.';
	}
	*LinePtr++ = '\t';

	// QUAL
	LineBuf_Append(Rf_asReal(VECTOR_ELT(X, 4)));
	*LinePtr++ = '\t';

	// FILTER
	SEXP tmp = VECTOR_ELT(X, 5);
	if (Rf_isFactor(tmp))
		tmp = Rf_asCharacterFactor(tmp);
	else
		tmp = AS_CHARACTER(tmp);
	LineBuf_Append(CHAR(STRING_ELT(tmp, 0)));
	*LinePtr++ = '\t';
}


/// export the INFO and FORMAT fields
inline static void ExportInfoFormat(SEXP X)
{
	// variable list
	SEXP VarNames = getAttrib(X, R_NamesSymbol);

	//====  INFO  ====//

	LineBuf_NeedSize(32);
	size_t cnt_info = VCF_INFO_Number.size();
	size_t n = 0;
	for (size_t i=0; i < cnt_info; i++)
	{
		// name, "info.*"
		const char *nm = CHAR(STRING_ELT(VarNames, i + 8)) + 5;
		// SEXP
		SEXP D = VECTOR_ELT(X, i + 8);

		if (IS_LOGICAL(D))  // FLAG type
		{
			if (Rf_asLogical(D) == TRUE)
			{
				if (n > 0) *LinePtr++ = ';';
				LineBuf_Append(nm);
				n ++;
			}
		} else {
			int m = INFO_GetNum(D, VCF_INFO_Number[i]);
			if (m > 0)
			{
				if (n > 0) *LinePtr++ = ';';
				LineBuf_Append(nm);
				*LinePtr++ = '=';
				INFO_Write(D, m);
				n ++;
			}
		}
	}

	if (n <= 0) *LinePtr++ = '.';
	*LinePtr++ = '\t';

	//====  FORMAT  ====//

	VCF_FORMAT_List.clear();
	LineBuf_NeedSize(32);
	LinePtr[0] = 'G', LinePtr[1] = 'T'; LinePtr += 2;

	size_t cnt_fmt = VCF_FORMAT_Number.size();
	for (size_t i=0; i < cnt_fmt; i++)
	{
		// name, "fmt.*"
		const char *nm = CHAR(STRING_ELT(VarNames, i + cnt_info + 8));
		SEXP D = VECTOR_ELT(X, i + cnt_info + 8);
		if (!isNull(D))
		{
			*LinePtr++ = ':';
			LineBuf_Append(nm + 4);
			VCF_FORMAT_List.push_back(D);
		}
	}
	*LinePtr++ = '\t';
}

}


extern "C"
{
using namespace SeqArray;

// ========================================================================
// Convert to VCF4: GDS -> VCF4
// ========================================================================

/// double quote text if needed
COREARRAY_DLL_EXPORT SEXP SEQ_Quote(SEXP text, SEXP dQuote)
{
	SEXP NewText, ans;
	PROTECT(NewText = AS_CHARACTER(text));
	PROTECT(ans = NEW_CHARACTER(Rf_length(NewText)));

	for (int i=0; i < Rf_length(NewText); i++)
	{
		string tmp = QuoteText(CHAR(STRING_ELT(NewText, i)));
		if (LOGICAL(dQuote)[0] == TRUE)
		{
			if ((tmp[0] != '\"') || (tmp[tmp.size()-1] != '\"'))
			{
				tmp.insert(0, "\"");
				tmp.push_back('\"');
			}
		}
		SET_STRING_ELT(ans, i, mkChar(tmp.c_str()));
	}

	UNPROTECT(2);
	return ans;
}




// ========================================================================


/// initialize
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Init(SEXP Sel, SEXP Info, SEXP Format,
	SEXP File)
{
	VCF_NumAllele = INTEGER(Sel)[0];
	VCF_NumSample = INTEGER(Sel)[1];
	VCF_File = R_GetConnection(File);

	int *pInfo = INTEGER(Info);
	VCF_INFO_Number.assign(pInfo, pInfo + Rf_length(Info));

	int *pFmt = INTEGER(Format);
	VCF_FORMAT_Number.assign(pFmt, pFmt + Rf_length(Format));

	VCF_FORMAT_List.reserve(256);
	LineBuf_Init();

	return R_NilValue;
}

/// finalize
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Done()
{
	LineBuf_Done();
	return R_NilValue;
}



/// convert to VCF4
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();

	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);

	// INFO, FORMAT
	ExportInfoFormat(X);

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);
	int *pSamp = INTEGER(geno);

	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	int *pAllele = INTEGER(phase);

	// for-loop of samples
	for (size_t i=0; i < VCF_NumSample; i ++)
	{
		// add '\t'
		if (i > 0) *LinePtr++ = '\t';

		// genotypes
		LineBuf_NeedSize(VCF_NumAllele << 4); // NumAllele*16
		if (VCF_NumAllele == 2)
		{
			_Line_Append_Geno(*pSamp++);
			*LinePtr++ = (*pAllele++) ? '|' : '/';
			_Line_Append_Geno(*pSamp++);
		} else {
			for (size_t j=0; j < VCF_NumAllele; j++)
			{
				if (j > 0)
					*LinePtr++ = (*pAllele++) ? '|' : '/';
				_Line_Append_Geno(*pSamp++);
			}
		}

		// annotation
		vector<SEXP>::iterator p;
		for (p=VCF_FORMAT_List.begin(); p != VCF_FORMAT_List.end(); p++)
		{
			*LinePtr++ = ':';
			size_t n = Rf_length(*p) / VCF_NumSample;
			FORMAT_Write(*p, n, i, VCF_NumSample);
		}
	}

	*LinePtr++ = '\n';

	// output
	if (VCF_File->text)
	{
		*LinePtr = 0;
		put_text("%s", LineBegin);
	} else {
		size_t size = LinePtr - LineBegin;
		size_t n = R_WriteConnection(VCF_File, LineBegin, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}



// --------------------------------------------------------------

#ifdef __SSE2__

static const __m128i char_unphased = _mm_set_epi8(
	'\t', '.', '/', '.',   '\t', '.', '/', '.',
	'\t', '.', '/', '.',   '\t', '.', '/', '.');
static const __m128i char_phased = _mm_set_epi8(
	'\t', '.', '|', '.',   '\t', '.', '|', '.',
	'\t', '.', '|', '.',   '\t', '.', '|', '.');
static const __m128i nines = _mm_set1_epi32(9);
static const __m128i char_zero = _mm_set1_epi16('0');
static const __m128i char_mask = _mm_set1_epi16(0xFF00);

#endif

#ifdef __AVX2__

static const __m256i char_unphased_256 = _mm256_set_epi8(
	'\t', '.', '/', '.',   '\t', '.', '/', '.',
	'\t', '.', '/', '.',   '\t', '.', '/', '.',
	'\t', '.', '/', '.',   '\t', '.', '/', '.',
	'\t', '.', '/', '.',   '\t', '.', '/', '.');
static const __m256i char_phased_256 = _mm256_set_epi8(
	'\t', '.', '|', '.',   '\t', '.', '|', '.',
	'\t', '.', '|', '.',   '\t', '.', '|', '.',
	'\t', '.', '|', '.',   '\t', '.', '|', '.',
	'\t', '.', '|', '.',   '\t', '.', '|', '.');
static const __m256i nines_256 = _mm256_set1_epi32(9);
static const __m256i char_zero_256 = _mm256_set1_epi16('0');
static const __m256i char_mask_256 = _mm256_set1_epi16(0xFF00);

#endif


/// convert to VCF4, diploid without FORMAT variables
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Di_WrtFmt(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();

	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);

	// INFO, FORMAT
	ExportInfoFormat(X);

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);
	int *pSamp = INTEGER(geno);

	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	int *pAllele = INTEGER(phase);

	// for-loop, genotypes
	size_t n = VCF_NumSample;
	size_t offset = 0;

#ifdef __SSE2__

	LineBuf_NeedSize(n*4 + 64);

	// 32-byte alignment
	offset = (size_t)LinePtr & 0x1F;
	if (offset > 0)
	{
		offset = 32 - offset;
		memmove(LineBegin+offset, LineBegin, LinePtr-LineBegin);
		LinePtr += offset;
	}

#ifdef __AVX2__
	for (; n >= 8; n -= 8)
	{
		__m256i v1 = MM_LOADU_256(pSamp);
		if (_mm256_movemask_epi8(_mm256_cmpgt_epi32(v1, nines_256))) goto tail;

		__m256i v2 = MM_LOADU_256(pSamp+8);
		if (_mm256_movemask_epi8(_mm256_cmpgt_epi32(v2, nines_256))) goto tail;
		pSamp += 16;

		__m256i phase = MM_LOADU_256(pAllele);
		pAllele += 8;
		__m256i m = _mm256_cmpeq_epi32(phase, _mm256_setzero_si256());
		__m256i bkg = MM_BLEND_256(char_unphased_256, char_phased_256, m);

		__m256i v = _mm256_add_epi16(_mm256_packs_epi32(v1, v2), char_zero_256);
		v = _mm256_permute4x64_epi64(v, _MM_SHUFFLE(3,1,2,0));
		m = _mm256_cmpgt_epi32(_mm256_setzero_si256(), v);
		m = _mm256_or_si256(m, char_mask_256);
		v = MM_BLEND_256(bkg, v, m);
		_mm256_store_si256((__m256i *)LinePtr, v);
		LinePtr += 32;
	}
#endif

	for (; n >= 4; n -= 4)
	{
		__m128i v1 = MM_LOADU_128(pSamp);
		if (_mm_movemask_epi8(_mm_cmpgt_epi32(v1, nines))) break;

		__m128i v2 = MM_LOADU_128(pSamp+4);
		if (_mm_movemask_epi8(_mm_cmpgt_epi32(v2, nines))) break;
		pSamp += 8;

		__m128i phase = MM_LOADU_128(pAllele);
		pAllele += 4;
		__m128i m = _mm_cmpeq_epi32(phase, _mm_setzero_si128());
		__m128i bkg = MM_BLEND_128(char_unphased, char_phased, m);

		__m128i v = _mm_add_epi16(_mm_packs_epi32(v1, v2), char_zero);
		m = _mm_cmplt_epi16(v, _mm_setzero_si128());
		m = _mm_or_si128(m, char_mask);
		v = MM_BLEND_128(bkg, v, m);
		_mm_store_si128((__m128i *)LinePtr, v);
		LinePtr += 16;
	}

#ifdef __AVX2__
tail:
#endif

#endif
	// tail
	for (; n > 0; n--)
	{
		LineBuf_NeedSize(32);
		_Line_Append_Geno(*pSamp++);
		*LinePtr++ = (*pAllele++) ? '|' : '/';
		_Line_Append_Geno(*pSamp++);
		*LinePtr++ = '\t';
	}
	LinePtr --;
	*LinePtr++ = '\n';

	// output
	if (VCF_File->text)
	{
		*LinePtr = 0;
		put_text("%s", LineBegin+offset);
	} else {
		size_t size = LinePtr - LineBegin - offset;
		size_t n = R_WriteConnection(VCF_File, LineBegin + offset, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}

} // extern "C"
