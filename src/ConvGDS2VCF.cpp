// ===========================================================
//
// ConvGDS2VCF.cpp: the C++ code for the conversion from GDS to VCF
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
#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;


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



extern "C"
{
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

/// used in SEQ_OutVCF4
static vector<int> _VCF4_INFO_Number;    ///<
static vector<int> _VCF4_FORMAT_Number;  ///<
static vector<SEXP> _VCF4_FORMAT_List;

static size_t _VCF4_NumAllele;
static size_t _VCF4_NumSample;


const size_t LINE_BUFFER_SIZE = 4096;
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
	_VCF4_FORMAT_List.clear();
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
	size_t cnt_info = _VCF4_INFO_Number.size();
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
			int m = INFO_GetNum(D, _VCF4_INFO_Number[i]);
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

	_VCF4_FORMAT_List.clear();
	LineBuf_NeedSize(32);
	LinePtr[0] = 'G', LinePtr[1] = 'T'; LinePtr += 2;

	size_t cnt_fmt = _VCF4_FORMAT_Number.size();
	for (size_t i=0; i < cnt_fmt; i++)
	{
		// name, "fmt.*"
		const char *nm = CHAR(STRING_ELT(VarNames, i + cnt_info + 8));
		SEXP D = VECTOR_ELT(X, i + cnt_info + 8);
		if (!isNull(D))
		{
			*LinePtr++ = ':';
			LineBuf_Append(nm + 4);
			_VCF4_FORMAT_List.push_back(D);
		}
	}
	*LinePtr++ = '\t';
}



// ========================================================================


/// initialize
COREARRAY_DLL_EXPORT SEXP SEQ_InitOutVCF4(SEXP Sel, SEXP Info, SEXP Format)
{
	_VCF4_NumAllele = INTEGER(Sel)[0];
	_VCF4_NumSample = INTEGER(Sel)[1];

	int *pInfo = INTEGER(Info);
	_VCF4_INFO_Number.assign(pInfo, pInfo + Rf_length(Info));

	int *pFmt = INTEGER(Format);
	_VCF4_FORMAT_Number.assign(pFmt, pFmt + Rf_length(Format));

	_VCF4_FORMAT_List.reserve(256);
	LineBuf_Init();

	return R_NilValue;
}

/// finalize
COREARRAY_DLL_EXPORT SEXP SEQ_DoneOutVCF4()
{
	LineBuf_Done();
	return R_NilValue;
}



/// convert to VCF4
COREARRAY_DLL_EXPORT SEXP SEQ_OutVCF4(SEXP X)
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
	for (size_t i=0; i < _VCF4_NumSample; i ++)
	{
		// add '\t'
		if (i > 0) *LinePtr++ = '\t';

		// genotypes
		LineBuf_NeedSize(_VCF4_NumAllele << 4); // NumAllele*16
		if (_VCF4_NumAllele == 2)
		{
			_Line_Append_Geno(*pSamp++);
			*LinePtr++ = (*pAllele++) ? '|' : '/';
			_Line_Append_Geno(*pSamp++);
		} else {
			for (size_t j=0; j < _VCF4_NumAllele; j++)
			{
				if (j > 0)
					*LinePtr++ = (*pAllele++) ? '|' : '/';
				_Line_Append(*pSamp++);
			}
		}

		// annotation
		vector<SEXP>::iterator p;
		for (p=_VCF4_FORMAT_List.begin(); p != _VCF4_FORMAT_List.end(); p++)
		{
			*LinePtr++ = ':';
			size_t n = Rf_length(*p) / _VCF4_NumSample;
			FORMAT_Write(*p, n, i, _VCF4_NumSample);
		}
	}

	// return
	SEXP ans = PROTECT(NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, Rf_mkCharLen(LineBegin, LinePtr-LineBegin));
	UNPROTECT(1);
	return ans;
}

} // extern "C"
