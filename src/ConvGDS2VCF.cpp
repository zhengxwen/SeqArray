// ===========================================================
//
// ConvGDS2VCF.cpp: Format conversion from GDS to VCF
//
// Copyright (C) 2013-2022    Xiuwen Zheng
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

#include <cstdio>
#include <cstring>
#include <vector>
#include "Index.h"
#include "vectorization.h"

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
static size_t VCF_ChrPrefix_NChar = 0;    ///< # of characters in chromosome prefix
static const char *VCF_ChrPrefix = NULL;  ///< pointer to chromosome prefix
static Rconnection VCF_File = NULL;  ///< R connection object

static const size_t LINE_BUFFER_SIZE = 4096;
static vector<char> LineBuffer;
static char *LineBegin = NULL;
static char *LineEnd   = NULL;
static char *pLine     = NULL;


inline static void LineBuf_Init()
{
	LineBuffer.resize(LINE_BUFFER_SIZE);
	pLine = LineBegin = &LineBuffer[0];
	LineEnd = pLine + LINE_BUFFER_SIZE;
}

inline static void LineBuf_Done()
{
	LineBuffer.clear();
	vector<char>().swap(LineBuffer);
	LineBegin = pLine = LineEnd = NULL;
	VCF_INFO_Number.clear();
	vector<int>().swap(VCF_INFO_Number);
	VCF_FORMAT_Number.clear();
	vector<int>().swap(VCF_FORMAT_Number);
	VCF_FORMAT_List.clear();
	vector<SEXP>().swap(VCF_FORMAT_List);
}

inline static void LineBuf_InitPtr()
{
	pLine = LineBegin = &LineBuffer[0];
}

inline static void LineBuf_NeedSize(size_t st)
{
	if (pLine + st > LineEnd)
	{
		size_t p = pLine - LineBegin;
		size_t n = p + st;
		n = (n / LINE_BUFFER_SIZE + 1) * LINE_BUFFER_SIZE;
		LineBuffer.resize(n);
		LineBegin = &LineBuffer[0];
		pLine = LineBegin + p;
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
		pLine = fast_itoa(pLine, val);
	else
		*pLine++ = '.';
}

inline static void _Line_Append_Geno(int val)
{
	if (val >= 0)
	{
		if (val < 10)
			*pLine++ = val + '0';
		else
			pLine = fast_itoa(pLine, val);
	} else 
		*pLine++ = '.';
}

inline static void _Line_Append_Geno_Raw(C_UInt8 val)
{
	if (val < 10)
		*pLine++ = val + '0';
	else if (val == NA_RAW)
		*pLine++ = '.';
	else
		pLine = fast_itoa(pLine, val);
}

inline static void _Line_Append(double val)
{
	if (R_FINITE(val))
		pLine += snprintf(pLine, 32, "%g", val);
	else
		*pLine++ = '.';
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
	memcpy(pLine, txt, n);
	pLine += n;
}

inline static void LineBuf_Append(const char *txt, const size_t n)
{
	LineBuf_NeedSize(n + 16);
	memcpy(pLine, txt, n);
	pLine += n;
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
	if (n < 0)
		n = !Rf_isNull(X) ? Rf_length(X) : 0;

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
	switch (TYPEOF(X))
	{
	case INTSXP:
		if (!Rf_isFactor(X))
		{
			LineBuf_NeedSize(12*n + 32);
			for (int *p = INTEGER(X); i < n; i++)
			{
				if (i > 0) *pLine++ = ',';
				_Line_Append(*p++);
			}
			break;
		} else {
			X = Rf_asCharacterFactor(X);
		}
	case STRSXP:
		for (; i < n; i++)
		{
			if (i > 0) *pLine++ = ',';
			SEXP s = STRING_ELT(X, i);
			LineBuf_Append(
				((s != NA_STRING) && (CHAR(s)[0] != 0)) ? CHAR(s) : ".");
		}
		break;
	case REALSXP:
		LineBuf_NeedSize(12*n + 32);
		for (double *p = REAL(X); i < n; i++)
		{
			if (i > 0) *pLine++ = ',';
			_Line_Append(*p++);
		}
		break;
	case LGLSXP:
		LineBuf_NeedSize(2*n + 32);
		for (int *p = LOGICAL(X); i < n; i++)
		{
			if (i > 0) *pLine++ = ',';
			if (*p == NA_INTEGER)
				*pLine++ = '.';
			else if (*p == 0)
				*pLine++ = '0';
			else
				*pLine++ = '1';
		}
		break;
	default:
		throw ErrSeqArray("INFO_Write: invalid data type.");
	}
}


/// write format values
inline static void FORMAT_Write(SEXP X, size_t n, size_t Start, size_t Step)
{
	switch (TYPEOF(X))
	{
	case INTSXP:
		if (!Rf_isFactor(X))
		{
			int *base = INTEGER(X) + Start, *p = base + (n - 1)*Step;
			for (; n > 0; n--, p-=Step)
				if (*p != NA_INTEGER) break;
			LineBuf_NeedSize(12*n + 32);
			p = base;
			for (size_t i=0; i < n; i++)
			{
				if (i > 0) *pLine++ = ',';
				_Line_Append(*p); p += Step;
			}
			break;
		} else {
			X = Rf_asCharacterFactor(X);
		}
	case STRSXP:
		{
			for (; n > 0; n--)
			{
				SEXP s = STRING_ELT(X, Start + (n - 1)*Step);
				if ((s != NA_STRING) && (CHAR(s)[0] != 0)) break;
			}
			for (size_t i=0; i < n; i++, Start += Step)
			{
				if (i > 0) *pLine++ = ',';
				SEXP s = STRING_ELT(X, Start);
				if ((s != NA_STRING) && (CHAR(s)[0] != 0))
					LineBuf_Append(CHAR(s));
				else
					*pLine++ = '.';
			}
			break;
		}
	case REALSXP:
		{
			double *base = REAL(X) + Start, *p = base + (n - 1)*Step;
			for (; n > 0; n--, p-=Step)
				if (R_finite(*p)) break;
			LineBuf_NeedSize(12*n + 32);
			p = base;
			for (size_t i=0; i < n; i++)
			{
				if (i > 0) *pLine++ = ',';
				_Line_Append(*p); p += Step;
			}
			break;
		}
	case LGLSXP:
		{
			int *base = LOGICAL(X) + Start, *p = base + (n - 1)*Step;
			for (; n > 0; n--, p-=Step)
				if (*p != NA_INTEGER) break;
			LineBuf_NeedSize(2*n + 32);
			p = base;
			for (size_t i=0; i < n; i++)
			{
				if (i > 0) *pLine++ = ',';
				if (*p == NA_INTEGER)
					*pLine++ = '.';
				else if (*p == 0)
					*pLine++ = '0';
				else
					*pLine++ = '1';
				p += Step;
			}
			break;
		}
	default:
		throw ErrSeqArray("FORMAT_Write: invalid data type.");
	}

	if (n <= 0) *pLine++ = '.';
}





// ========================================================================

/// export the first seven columns: chr, pos, id, allele (REF/ALT), qual, filter
inline static void ExportHead(SEXP X)
{
	// CHROM
	if (VCF_ChrPrefix_NChar > 0)
		LineBuf_Append(VCF_ChrPrefix, VCF_ChrPrefix_NChar);
	LineBuf_Append(CHAR(STRING_ELT(VECTOR_ELT(X, 0), 0)));
	*pLine++ = '\t';

	// POS
	LineBuf_Append(Rf_asInteger(VECTOR_ELT(X, 1)));
	*pLine++ = '\t';

	// ID
	char *s = (char*)CHAR(STRING_ELT(VECTOR_ELT(X, 2), 0));
	if (*s != 0)
		LineBuf_Append(s);
	else
		*pLine++ = '.';
	*pLine++ = '\t';

	// allele -- REF/ALT
	size_t n = pLine - LineBegin;
	LineBuf_Append(CHAR(STRING_ELT(VECTOR_ELT(X, 3), 0)));

	for (s = LineBegin+n; s < pLine; s++)
	{
		if (*s == ',')
			{ *s = '\t'; break; }
	}
	if (s == pLine)
	{
		*pLine++ = '\t';
		*pLine++ = '.';
	}
	*pLine++ = '\t';

	// QUAL
	LineBuf_Append(Rf_asReal(VECTOR_ELT(X, 4)));
	*pLine++ = '\t';

	// FILTER
	SEXP x = VECTOR_ELT(X, 5);
	bool is_na = false;
	switch (TYPEOF(x))
	{
	case LGLSXP:
		is_na = (LOGICAL_ELT(x, 0) == NA_LOGICAL); break;
	case INTSXP:
		is_na = (INTEGER_ELT(x, 0) == NA_INTEGER); break;
	case REALSXP:
		is_na = R_FINITE(REAL_ELT(x, 0)); break;
	case STRSXP:
		is_na = (STRING_ELT(x, 0) == NA_STRING); break;
	}
	if (!is_na)
	{
		if (Rf_isFactor(x))
			x = Rf_asCharacterFactor(x);
		else
			x = AS_CHARACTER(x);
		LineBuf_Append(CHAR(STRING_ELT(x, 0)));
	} else {
		*pLine++ = '.';
	}
	*pLine++ = '\t';
}


/// export the INFO and FORMAT fields
inline static void ExportInfoFormat(SEXP X, size_t info_st)
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
		const char *nm = CHAR(STRING_ELT(VarNames, i + info_st)) + 5;
		// SEXP
		SEXP D = VECTOR_ELT(X, i + info_st);

		if (IS_LOGICAL(D))  // FLAG type
		{
			if (Rf_asLogical(D) == TRUE)
			{
				if (n > 0) *pLine++ = ';';
				LineBuf_Append(nm);
				n ++;
			}
		} else {
			int m = INFO_GetNum(D, VCF_INFO_Number[i]);
			if (m > 0)
			{
				if (n > 0) *pLine++ = ';';
				LineBuf_Append(nm);
				*pLine++ = '=';
				INFO_Write(D, m);
				n ++;
			}
		}
	}

	if (n <= 0) *pLine++ = '.';
	*pLine++ = '\t';

	//====  FORMAT  ====//

	if (VCF_NumSample <= 0) return;  // no FORMAT column
	VCF_FORMAT_List.clear();
	size_t cnt_fmt = VCF_FORMAT_Number.size();

	LineBuf_NeedSize(32);
	if (info_st > 6)
	{
		pLine[0] = 'G'; pLine[1] = 'T';
		pLine += 2;
	} else if (cnt_fmt <= 0)
	{
		pLine[0] = '.';
		pLine ++;
	}

	for (size_t i=0; i < cnt_fmt; i++)
	{
		// name, "fmt.*"
		SEXP D = VECTOR_ELT(X, i + cnt_info + info_st);
		if (!isNull(D))
		{
			if (i > 0 || info_st > 6)
				*pLine++ = ':';
			const char *nm = CHAR(STRING_ELT(VarNames, i + cnt_info + info_st));
			LineBuf_Append(nm + 4);
			VCF_FORMAT_List.push_back(D);
		}
	}
	*pLine++ = '\t';
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
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Init(SEXP SelDim, SEXP ChrPrefix, SEXP Info,
	SEXP Format, SEXP File, SEXP Verbose)
{
	int num_allele = INTEGER(SelDim)[0];
	if (num_allele <= 0) num_allele = 2;  // diploid (notice NA < 0)
	VCF_NumAllele = num_allele;
	VCF_NumSample = INTEGER(SelDim)[1];

	SEXP chr = STRING_ELT(ChrPrefix, 0);
	if (chr != NA_STRING)
	{
		VCF_ChrPrefix = CHAR(chr);
		VCF_ChrPrefix_NChar = strlen(VCF_ChrPrefix);
	} else {
		VCF_ChrPrefix_NChar = 0;
		VCF_ChrPrefix = "";
	}

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



/// convert to VCF4 in general
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();
	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);
	// INFO, FORMAT
	ExportInfoFormat(X, 8);
	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	C_UInt8 *pAllele = (C_UInt8*)RAW(phase);

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);

	if (TYPEOF(geno) == RAWSXP)
	{
		C_UInt8 *pSamp = (C_UInt8*)RAW(geno);
		// for-loop of samples
		for (size_t i=0; i < VCF_NumSample; i++)
		{
			// add '\t'
			if (i > 0) *pLine++ = '\t';
			// genotypes
			LineBuf_NeedSize(VCF_NumAllele << 4); // NumAllele*16
			if (VCF_NumAllele == 2)
			{
				_Line_Append_Geno_Raw(*pSamp++);
				*pLine++ = (*pAllele++) ? '|' : '/';
				_Line_Append_Geno_Raw(*pSamp++);
			} else {
				for (size_t j=0; j < VCF_NumAllele; j++)
				{
					if (j > 0)
						*pLine++ = (*pAllele++) ? '|' : '/';
					_Line_Append_Geno_Raw(*pSamp++);
				}
			}
			// annotation
			vector<SEXP>::iterator p;
			for (p=VCF_FORMAT_List.begin(); p != VCF_FORMAT_List.end(); p++)
			{
				*pLine++ = ':';
				size_t n = Rf_length(*p) / VCF_NumSample;
				FORMAT_Write(*p, n, i, VCF_NumSample);
			}
		}
	} else {
		int *pSamp = INTEGER(geno);
		// for-loop of samples
		for (size_t i=0; i < VCF_NumSample; i++)
		{
			// add '\t'
			if (i > 0) *pLine++ = '\t';
			// genotypes
			LineBuf_NeedSize(VCF_NumAllele << 4); // NumAllele*16
			if (VCF_NumAllele == 2)
			{
				_Line_Append_Geno(*pSamp++);
				*pLine++ = (*pAllele++) ? '|' : '/';
				_Line_Append_Geno(*pSamp++);
			} else {
				for (size_t j=0; j < VCF_NumAllele; j++)
				{
					if (j > 0)
						*pLine++ = (*pAllele++) ? '|' : '/';
					_Line_Append_Geno(*pSamp++);
				}
			}
			// annotation
			vector<SEXP>::iterator p;
			for (p=VCF_FORMAT_List.begin(); p != VCF_FORMAT_List.end(); p++)
			{
				*pLine++ = ':';
				size_t n = Rf_length(*p) / VCF_NumSample;
				FORMAT_Write(*p, n, i, VCF_NumSample);
			}
		}
	}

	*pLine++ = '\n';

	// output
	if (VCF_File->text)
	{
		*pLine = 0;
		put_text("%s", LineBegin);
	} else {
		size_t size = pLine - LineBegin;
		size_t n = R_WriteConnection(VCF_File, LineBegin, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}



// --------------------------------------------------------------

/// convert to VCF4, diploid without FORMAT variables
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Di_WrtFmt(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();
	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);
	// INFO, FORMAT
	ExportInfoFormat(X, 8);
	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	C_UInt8 *pAllele = (C_UInt8*)RAW(phase);

	// for-loop, genotypes
	size_t n = VCF_NumSample;
	size_t offset = 0;

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);

	if (TYPEOF(geno) == RAWSXP)
	{
		C_UInt8 *pSamp = (C_UInt8*)RAW(geno);

	#ifdef COREARRAY_SIMD_SSE2

		// need buffer
		LineBuf_NeedSize(n*4 + 64);
		// 32-byte alignment
		offset = (size_t)pLine & 0x1F;
		if (offset > 0)
		{
			offset = 32 - offset;
			memmove(LineBegin+offset, LineBegin, pLine-LineBegin);
			pLine += offset;
		}

		static const __m128i char_unphased = _mm_set1_epi8('/');
		static const __m128i char_phased = _mm_set1_epi8('|');
		static const __m128i char_tab = _mm_set1_epi8('\t');

	#ifdef COREARRAY_SIMD_AVX2
		static const __m256i ten = _mm256_set1_epi8(10);
		static const __m256i na  = _mm256_set1_epi8(0xFF);
		static const __m256i char_zero = _mm256_set1_epi8('0');
		static const __m256i char_na = _mm256_set1_epi8('.');

		for (; n >= 16; n-=16)
		{
			__m256i v1 = MM_LOADU_256(pSamp);
			__m256i m1 = _mm256_cmpeq_epi8(v1, na);
			if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(ten, _mm256_min_epu8(
				_mm256_andnot_si256(m1, v1), ten)))) break;
			pSamp += 32;

			v1 = _mm256_add_epi8(v1, char_zero);
			v1 = MM_BLEND_256(char_na, v1, m1);

			__m128i v3 = MM_LOADU_128(pAllele);
			pAllele += 16;
			__m128i m = _mm_cmpeq_epi8(v3, _mm_setzero_si128());
			v3 = MM_BLEND_128(char_unphased, char_phased, m);
			__m256i phase = MM_SET_M128(_mm_unpackhi_epi8(v3, char_tab),
				_mm_unpacklo_epi8(v3, char_tab));

			__m256i w1 = _mm256_unpacklo_epi8(v1, phase);
			__m256i w2 = _mm256_unpackhi_epi8(v1, phase);
			_mm256_store_si256((__m256i *)pLine,
				_mm256_permute2x128_si256(w1, w2, 0x20));
			_mm256_store_si256((__m256i *)(pLine+32),
				_mm256_permute2x128_si256(w1, w2, 0x31));
			pLine += 64;
		}
	#else
		static const __m128i ten = _mm_set1_epi8(10);
		static const __m128i na  = _mm_set1_epi8(0xFF);
		static const __m128i char_zero = _mm_set1_epi8('0');
		static const __m128i char_na = _mm_set1_epi8('.');

		for (; n >= 16; n-=16)
		{
			__m128i v1 = MM_LOADU_128(pSamp);
			__m128i m1 = _mm_cmpeq_epi8(v1, na);
			if (_mm_movemask_epi8(_mm_cmpeq_epi8(ten, _mm_min_epu8(
				_mm_andnot_si128(m1, v1), ten)))) break;

			__m128i v2 = MM_LOADU_128((pSamp+16));
			__m128i m2 = _mm_cmpeq_epi8(v2, na);
			if (_mm_movemask_epi8(_mm_cmpeq_epi8(ten, _mm_min_epu8(
				_mm_andnot_si128(m2, v2), ten)))) break;
			pSamp += 32;

			v1 = _mm_add_epi8(v1, char_zero);
			v1 = MM_BLEND_128(char_na, v1, m1);
			v2 = _mm_add_epi8(v2, char_zero);
			v2 = MM_BLEND_128(char_na, v2, m2);

			__m128i v3 = MM_LOADU_128(pAllele);
			pAllele += 16;
			__m128i m = _mm_cmpeq_epi8(v3, _mm_setzero_si128());
			v3 = MM_BLEND_128(char_unphased, char_phased, m);

			__m128i p1 = _mm_unpacklo_epi8(v3, char_tab);
			__m128i p2 = _mm_unpackhi_epi8(v3, char_tab);

			_mm_store_si128((__m128i *)pLine, _mm_unpacklo_epi8(v1, p1));
			_mm_store_si128((__m128i *)(pLine+16), _mm_unpackhi_epi8(v1, p1));
			_mm_store_si128((__m128i *)(pLine+32), _mm_unpacklo_epi8(v2, p2));
			_mm_store_si128((__m128i *)(pLine+48), _mm_unpackhi_epi8(v2, p2));
			pLine += 64;
		}
	#endif

	#endif

		// tail
		for (; n > 0; n--)
		{
			LineBuf_NeedSize(32);
			_Line_Append_Geno_Raw(*pSamp++);
			*pLine++ = (*pAllele++) ? '|' : '/';
			_Line_Append_Geno_Raw(*pSamp++);
			*pLine++ = '\t';
		}
	} else {
		// integer vector for genotypes
		int *pSamp = INTEGER(geno);
		for (; n > 0; n--)
		{
			LineBuf_NeedSize(32);
			_Line_Append_Geno(*pSamp++);
			*pLine++ = (*pAllele++) ? '|' : '/';
			_Line_Append_Geno(*pSamp++);
			*pLine++ = '\t';
		}
	}

	pLine--; *pLine++ = '\n';

	// output
	if (VCF_File->text)
	{
		*pLine = 0;
		put_text("%s", LineBegin+offset);
	} else {
		size_t size = pLine - LineBegin - offset;
		size_t n = R_WriteConnection(VCF_File, LineBegin + offset, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}



// --------------------------------------------------------------

/// convert to haploid VCF4
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_Haploid(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();
	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);
	// INFO, FORMAT
	ExportInfoFormat(X, 7);

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);

	if (TYPEOF(geno) == RAWSXP)
	{
		C_UInt8 *pSamp = (C_UInt8*)RAW(geno);
		// for-loop of samples
		for (size_t i=0; i < VCF_NumSample; i++)
		{
			// add '\t'
			if (i > 0) *pLine++ = '\t';
			// genotypes
			LineBuf_NeedSize(VCF_NumAllele << 3); // NumAllele*8
			_Line_Append_Geno_Raw(*pSamp++);
			// annotation
			vector<SEXP>::iterator p;
			for (p=VCF_FORMAT_List.begin(); p != VCF_FORMAT_List.end(); p++)
			{
				*pLine++ = ':';
				size_t n = Rf_length(*p) / VCF_NumSample;
				FORMAT_Write(*p, n, i, VCF_NumSample);
			}
		}
	} else {
		int *pSamp = INTEGER(geno);
		// for-loop of samples
		for (size_t i=0; i < VCF_NumSample; i++)
		{
			// add '\t'
			if (i > 0) *pLine++ = '\t';
			// genotypes
			LineBuf_NeedSize(VCF_NumAllele << 3); // NumAllele*8
			_Line_Append_Geno(*pSamp++);
			// annotation
			vector<SEXP>::iterator p;
			for (p=VCF_FORMAT_List.begin(); p != VCF_FORMAT_List.end(); p++)
			{
				*pLine++ = ':';
				size_t n = Rf_length(*p) / VCF_NumSample;
				FORMAT_Write(*p, n, i, VCF_NumSample);
			}
		}
	}

	*pLine++ = '\n';

	// output
	if (VCF_File->text)
	{
		*pLine = 0;
		put_text("%s", LineBegin);
	} else {
		size_t size = pLine - LineBegin;
		size_t n = R_WriteConnection(VCF_File, LineBegin, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}



// --------------------------------------------------------------

/// convert to VCF4 without genotypes
COREARRAY_DLL_EXPORT SEXP SEQ_ToVCF_NoGeno(SEXP X)
{
	// initialize line pointer
	LineBuf_InitPtr();
	// CHROM, POS, ID, REF, ALT, QUAL, FILTER
	ExportHead(X);
	// INFO, FORMAT
	ExportInfoFormat(X, 6);

	// for-loop of samples
	LineBuf_NeedSize(VCF_NumSample + 16);
	for (size_t i=0; i < VCF_NumSample; i++)
	{
		// add '\t'
		if (i > 0) *pLine++ = '\t';
		// annotation
		vector<SEXP>::iterator p, st = VCF_FORMAT_List.begin();
		for (p=st; p != VCF_FORMAT_List.end(); p++)
		{
			if (p != st) *pLine++ = ':';
			size_t n = Rf_length(*p) / VCF_NumSample;
			FORMAT_Write(*p, n, i, VCF_NumSample);
		}
	}
	*pLine++ = '\n';

	// output
	if (VCF_File->text)
	{
		*pLine = 0;
		put_text("%s", LineBegin);
	} else {
		size_t size = pLine - LineBegin;
		size_t n = R_WriteConnection(VCF_File, LineBegin, size);
		if (size != n)
			throw ErrSeqArray("writing error.");
	}

	return R_NilValue;
}

} // extern "C"
