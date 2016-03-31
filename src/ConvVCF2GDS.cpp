// ===========================================================
//
// ConvVCF2GDS.cpp: format conversion from VCF to GDS
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
#include <vector>
#include <set>
#include <algorithm>

extern "C"
{
#define class xclass
#define private xprivate
#include <R_ext/Connections.h>
#undef class
#undef private
}


namespace SeqArray
{
using namespace std;

// ===========================================================
// define 
// ===========================================================

static Rconnection VCF_File = NULL;  ///< R connection object

static vector<char> VCF_Buffer;  ///< reading buffer
static char *VCF_Buffer_Ptr;     ///< the current pointer to reading buffer
static char *VCF_Buffer_EndPtr;  ///< the end pointer to reading buffer
static const size_t VCF_BUFFER_SIZE = 65536;  ///< reading buffer size
static const size_t VCF_BUFFER_SIZE_PLUS = 32;  ///< additional buffer is needed since *VCF_Buffer_EndPtr might be revised

/// read file buffer
inline static void Read_VCF_Buffer()
{
	VCF_Buffer_Ptr = &VCF_Buffer[0];
	size_t n = R_ReadConnection(VCF_File, VCF_Buffer_Ptr, VCF_BUFFER_SIZE);
	VCF_Buffer_EndPtr = VCF_Buffer_Ptr + n;
	if (n <= 0)
	{
		if (VCF_File->EOF_signalled)
			throw ErrSeqArray("read text error.");
		VCF_File->EOF_signalled = TRUE;
	}
}

/// test EOF
inline static bool VCF_EOF()
{
	if (VCF_File->EOF_signalled) return true;
	if (VCF_Buffer_Ptr >= VCF_Buffer_EndPtr)
		Read_VCF_Buffer();
	return (VCF_Buffer_Ptr >= VCF_Buffer_EndPtr);
}



// ===========================================================

static vector<char> Text_Buffer;  ///< text buffer
static char *Text_BeginPtr;  ///< the starting pointer for the text buffer
static char *Text_EndPtr;    ///< the end pointer for the text buffer
static size_t Text_Size;     ///< the buffer size in Text_Buffer

static C_Int64 VCF_LineNum;  ///< the current line number
static int VCF_ColumnNum;    ///< the current column number
static C_Int64 VCF_NextLineNum;  ///< the next line number
static int VCF_NextColumnNum;    ///< the next column number


/// get a string with a seperator '\t', which is saved in _Text_Buffer
static void GetText(int last_column)
{
	if (VCF_File->EOF_signalled)
		throw ErrSeqArray("it is the end of file.");

	VCF_ColumnNum = VCF_NextColumnNum;
	VCF_LineNum = VCF_NextLineNum;

	Text_EndPtr = Text_BeginPtr = &Text_Buffer[0];
	int ch = -1;

	while (true)
	{
		char *p = VCF_Buffer_Ptr;
		while (p < VCF_Buffer_EndPtr)
		{
			ch = *p;
			if (ch=='\t' || ch=='\n' || ch=='\r')
				break;
			p ++;
		}

		// copy to Text_Buffer
		size_t n = p - VCF_Buffer_Ptr;
		size_t m = Text_EndPtr - Text_BeginPtr;
		size_t nn = m + n;
		if (nn > Text_Size)
		{
			nn = ((nn / 1024) + 1) * 1024;
			Text_Buffer.resize(nn + VCF_BUFFER_SIZE_PLUS);
			Text_Size = nn;
			Text_BeginPtr = &Text_Buffer[0];
			Text_EndPtr = Text_BeginPtr + m;
		}
		memcpy(Text_EndPtr, VCF_Buffer_Ptr, n);
		VCF_Buffer_Ptr += n;
		Text_EndPtr += n;

		if (p < VCF_Buffer_EndPtr || VCF_File->EOF_signalled)
			break;
		else
			Read_VCF_Buffer();
	}

	if (ch == '\t')
	{
		if (last_column == TRUE)
			throw ErrSeqArray("more columns than what expected.");
		VCF_NextColumnNum ++;
		VCF_Buffer_Ptr ++;
	} else {
		if (last_column == FALSE)
			throw ErrSeqArray("fewer columns than what expected.");
		VCF_NextColumnNum = 1;
		VCF_NextLineNum ++;

		// skip '\n' and '\r'
		if (ch=='\n' || ch=='\r')
		{
			do {
				VCF_Buffer_Ptr ++;
				if (VCF_Buffer_Ptr >= VCF_Buffer_EndPtr)
				{
					if (VCF_File->EOF_signalled)
						break;
					Read_VCF_Buffer();
				}
				ch = *VCF_Buffer_Ptr;
			} while (ch=='\n' || ch=='\r');
		}
	}
}

/// skip the current line
static void SkipLine()
{
	VCF_ColumnNum = VCF_NextColumnNum;
	VCF_LineNum = VCF_NextLineNum;

	int ch = -1;
	while (true)
	{
		while (VCF_Buffer_Ptr < VCF_Buffer_EndPtr)
		{
			ch = *VCF_Buffer_Ptr;
			if (ch=='\n' || ch=='\r') break;
			VCF_Buffer_Ptr ++;
		}
		if (VCF_Buffer_Ptr < VCF_Buffer_EndPtr || VCF_File->EOF_signalled)
			break;
		else
			Read_VCF_Buffer();
	}

	// skip '\n' and '\r'
	if (ch=='\n' || ch=='\r')
	{
		do {
			VCF_Buffer_Ptr ++;
			if (VCF_Buffer_Ptr >= VCF_Buffer_EndPtr)
			{
				if (VCF_File->EOF_signalled)
					break;
				Read_VCF_Buffer();
			}
			ch = *VCF_Buffer_Ptr;
		} while (ch=='\n' || ch=='\r');
	}

	VCF_NextColumnNum = 1;
	VCF_NextLineNum ++;
}

/// skip white space
inline static void SkipWhiteSpace()
{
	// skip whitespace
	while ((VCF_Buffer_Ptr < VCF_Buffer_EndPtr) && (*VCF_Buffer_Ptr == ' '))
		VCF_Buffer_Ptr ++;
	while ((VCF_Buffer_Ptr < VCF_Buffer_EndPtr) && (*(VCF_Buffer_EndPtr-1) == ' '))
		VCF_Buffer_EndPtr --;
}

/// skip a dot
inline static void SkipDot()
{
	SkipWhiteSpace();
	if ((Text_EndPtr-Text_BeginPtr == 1) && (*Text_BeginPtr == '.'))
		Text_BeginPtr ++;
}



// ===========================================================
// VCF structure
// ===========================================================

const static int FIELD_TYPE_INT      = 1;
const static int FIELD_TYPE_FLOAT    = 2;
const static int FIELD_TYPE_FLAG     = 3;
const static int FIELD_TYPE_STRING   = 4;


/// the structure of INFO field
struct COREARRAY_DLL_LOCAL TVCF_Field_Info
{
	string name;         //< INFO ID
	int type;            //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;    //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;  //< according to 'Number=' field
	             // -1: variable-length (.), -2: # of alternate alleles (A)
	             // -3: # of possible genotypes (G), -4: # of alleles (R)
	bool used;   //< if TRUE, it has been parsed for the current line

	TVCF_Field_Info()
	{
		type = 0;
		data_obj = len_obj = NULL;
		number = 0;
		used = false;
	}

	template<typename TYPE> inline void Index(vector<TYPE> &array, int num_allele)
	{
		C_Int32 I32 = array.size();
		switch (number)
		{
		case -1:  // variable-length, .
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
			break;

		case -2:  // # of alternate alleles, A
			if (I32 != (num_allele-1))
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=A) should have %d value(s), but receives %d.",
					name.c_str(), num_allele-1, I32);
			}
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
			break;

		case -3:  // # of all possible genotypes, G
			if (I32 != (num_allele+1)*num_allele/2)
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=G) should have %d value(s), but receives %d.",
					name.c_str(), (num_allele+1)*num_allele/2, I32);
			}
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
			break;

		case -4:  // # of alleles, R
			if (I32 != num_allele)
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=R) should have %d value(s), but receives %d.",
					name.c_str(), num_allele, I32);
			}
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
			break;

		default:
			if (number >= 0)
			{
				if (number != (int)array.size())
				{
					throw ErrSeqArray(
						"INFO ID '%s' should have %d value(s), but receives %d.",
						name.c_str(), number, (int)array.size());
				}
			} else
				throw ErrSeqArray("Invalid value 'number' in TVCF_Field_Info.");
		}
	}

	template<typename TYPE> void Fill(vector<TYPE> &array, TYPE val)
	{
		if (number < 0)
		{
			C_Int32 I32 = 0;
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
		} else {
			array.clear();
			array.resize(number, val);
			GDS_Array_AppendData(data_obj, number, &(array[0]),
				TdTraits<TYPE>::SVType);
		}
	}
};


/// the structure of FORMAT field
struct COREARRAY_DLL_LOCAL TVCF_Field_Format
{
	string name;         //< FORMAT ID
	int type;            //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;    //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;  //< according to 'Number=' field
	             // -1: variable-length (.), -2: # of alternate alleles (A)
	             // -3: # of possible genotypes (G), -4: # of alleles (R)
	bool used;   //< if TRUE, it has been parsed for the current line

	/// data -- Int32
	vector< vector<C_Int32> > I32ss;
	/// data -- C_Float32
	vector< vector<C_Float64> > F64ss;
	/// data -- UTF8 string
	vector< vector<string> > UTF8ss;


	TVCF_Field_Format()
	{
		type = 0;
		data_obj = len_obj = NULL;
		number = 0;
		used = false;
	}

	// FORMAT field

	template<typename TYPE> void Check(vector<TYPE> &array, int num_allele,
		const TYPE &missing)
	{
		switch (number)
		{
			case -1:
				break;

			case -2:
				// # of alternate alleles, A
				if ((int)array.size() > (num_allele-1))
				{
					throw ErrSeqArray(
						"FORMAT ID '%s' (Number=A) should have %d value(s), but receives %d.",
						name.c_str(), num_allele-1, (int)array.size());
				} else
					array.resize(num_allele-1, missing);
				break;		

			case -3:
				// # of all possible genotypes, G
				if ((int)array.size() > (num_allele+1)*num_allele/2)
				{
					throw ErrSeqArray(
						"FORMAT ID '%s' (Number=G) should have %d value(s), but receives %d.",
						name.c_str(), (num_allele+1)*num_allele/2,
						(int)array.size());
				} else
					array.resize((num_allele+1)*num_allele/2, missing);
				break;		

			case -4:
				// # of alleles, R
				if ((int)array.size() > num_allele)
				{
					throw ErrSeqArray(
						"FORMAT ID '%s' (Number=R) should have %d value(s), but receives %d.",
						name.c_str(), num_allele, (int)array.size());
				} else
					array.resize(num_allele, missing);
				break;		

			default:
				if (number >= 0)
				{
					if ((int)array.size() > number)
					{
						throw ErrSeqArray(
							"FORMAT ID '%s' should have %d value(s), but receives %d.",
							name.c_str(), number, (int)array.size());
					} else {
						array.resize(number, missing);
					}
				} else
					throw ErrSeqArray("Invalid value 'number' in TVCF_Field_Format.");
		}
	}

	void WriteFixedLength()
	{
		if (number < 0)
		{
			throw ErrSeqArray(
				"Invalid call 'WriteFixedLength' in TVCF_Field_Format.");
		}
		switch (type)
		{
			case FIELD_TYPE_INT:
				for (vector< vector<C_Int32> >::iterator it = I32ss.begin();
					it != I32ss.end(); it ++)
				{
					GDS_Array_AppendData(data_obj, number, &((*it)[0]), svInt32);
				}
				break;

			case FIELD_TYPE_FLOAT:
				for (vector< vector<double> >::iterator it = F64ss.begin();
					it != F64ss.end(); it ++)
				{
					GDS_Array_AppendData(data_obj, number, &((*it)[0]), svFloat64);
				}
				break;

			case FIELD_TYPE_STRING:
				for (vector< vector<string> >::iterator it = UTF8ss.begin();
					it != UTF8ss.end(); it ++)
				{
					vector<string> &text = (*it);
					size_t n = text.size();
					for (size_t j=0; j < n; j ++)
						GDS_Array_AppendString(data_obj, text[j].c_str());
				}
				break;

			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
		}
	}

	int WriteVariableLength(int nTotalSample, vector<C_Int32> &I32s,
		vector<double> &F64s)
	{
		if (number >= 0)
		{
			throw ErrSeqArray(
				"Invalid call 'WriteVariableLength' in TVCF_Field_Format.");
		}
		int nMax = 0;
		switch (type)
		{
			case FIELD_TYPE_INT:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)I32ss[j].size())
						nMax = I32ss[j].size();
				}
				I32s.resize(nTotalSample);
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<C_Int32> &B = I32ss[j];
						I32s[j] = (i < (int)B.size()) ? B[i] : NA_INTEGER;
					}
					GDS_Array_AppendData(data_obj, nTotalSample,
						&(I32s[0]), svInt32);
				}
				break;

			case FIELD_TYPE_FLOAT:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)F64ss[j].size())
						nMax = F64ss[j].size();
				}
				F64s.resize(nTotalSample);
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<double> &L = F64ss[j];
						F64s[j] = (i < (int)L.size()) ? L[i] : R_NaN;
					}
					GDS_Array_AppendData(data_obj, nTotalSample,
						&(F64s[0]), svFloat64);
				}
				break;

			case FIELD_TYPE_STRING:
				for (int j=0; j < nTotalSample; j++)
				{
					if (nMax < (int)UTF8ss[j].size())
						nMax = UTF8ss[j].size();
				}
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<string> &B = UTF8ss[j];
						GDS_Array_AppendString(data_obj,
							(i < (int)B.size()) ? B[i].c_str() : "");
					}
				}
				break;

			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
		}
		return nMax;
	}
};



// ===========================================================
// Text conversion
// ===========================================================

inline static string SHORT(const char *p, const char *end)
{
	string s(p, end);
	if (s.size() > 32)
	{
		s.resize(32);
		s.append(" ...");
	}
	return s;
}



static const char *ERR_INT_CONV = "Invalid integer conversion '%s'";
static const char *ERR_INT_OUT_RANGE = "Integer conversion is out of range '%s'";
static const char *ERR_FLOAT_CONV = "Invalid float conversion '%s'";

/// get an integer from a string
inline static int getInt32(const char *p, const char *end, bool raise_error)
{
	while ((p < end) && (*p == ' '))
		p ++;
	const char *start = p;

	if ((p < end) && (*p == '.'))
	{
		p ++;
		while ((p < end) && (*p == ' ')) p ++;
		if ((p < end) && raise_error)
			throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
		return NA_INTEGER;
	}

	bool sign = ((p < end) && (*p == '-'));
	if (sign) p ++;

	C_Int64 val = 0;
	while (p < end)
	{
		char ch = *p ++;
		if ('0' <= ch && ch <= '9')
		{
			val = val*10 + (ch - '0');
			if (raise_error && (val > INT_MAX))
				throw ErrSeqArray(ERR_INT_OUT_RANGE, SHORT(start, end).c_str());
		} else {
			if (ch == ' ')
			{
				while ((p < end) && (*p == ' '))
					p ++;
				if (p >= end) break;
			}
			if (raise_error)
				throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
			return NA_INTEGER;
		}
	}

	return sign ? -val : val;
}

/// get multiple integers from a string
static void getInt32Array(const char *p, const char *end, vector<C_Int32> &I32s,
	bool raise_error)
{
	I32s.clear();

	while (p < end)
	{
		while ((p < end) && (*p == ' '))
			p ++;
		const char *start = p;

		if ((p < end) && (*p == '.'))
		{
			p ++;
			while ((p < end) && (*p == ' ')) p ++;
			if ((p < end) && (*p != ','))
			{
				if (raise_error)
					throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
				while ((p < end) && (*p != ',')) p ++;
			}
			I32s.push_back(NA_INTEGER);
			if ((p < end) && (*p == ',')) p ++;
			continue;
		}

		bool sign = ((p < end) && (*p == '-'));
		if (sign) p ++;

		C_Int64 val = 0;
		while (p < end)
		{
			char ch = *p;
			if ('0' <= ch && ch <= '9')
			{
				val = val*10 + (ch - '0');
				if (raise_error && (val > INT_MAX))
					throw ErrSeqArray(ERR_INT_OUT_RANGE, SHORT(start, end).c_str());
				p ++;
			} else {
				if (ch == ' ')
				{
					while ((p < end) && (*p == ' '))
						p ++;
					if (p >= end) break;
					ch = *p;
				}

				if (ch == ',')
				{
					p ++;
					break;
				} else {
					if (raise_error)
						throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
					else
						while ((p < end) && (*p != ',')) p ++;
				}

				val = NA_INTEGER; sign = false;
				if ((p < end) && (*p == ',')) p ++;
				break;
			}
		}

		I32s.push_back(sign ? -val : val);
	}
}


/// get a real number from a string
inline static double getFloat(char *p, char *end, bool raise_error)
{
	while ((p < end) && (*p == ' ')) p ++;
	while ((p < end) && (*(end-1) == ' ')) end --;
	const char *start = p;

	if (!((end-p == 1) && (*p == '.')))
	{
		char *endptr = (char*)p;
		*end = 0;  // no worry, see VCF_BUFFER_SIZE_PLUS
		double val = strtod(p, &endptr);

		if (endptr == p)
		{
			if (raise_error)
				throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
			val = R_NaN;
		} else {
			p = endptr;
			while ((p < end) && (*p == ' ')) p ++;
			if (p < end)
			{
				val = R_NaN;
				if (raise_error)
					throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
			}
		}

		return val;
	} else {
		return R_NaN;
	}
}

/// get multiple real numbers from a string
static void getFloatArray(char *p, char *end, vector<double> &F64s,
	bool raise_error)
{
	while ((p < end) && (*(end-1) == ' ')) end --;
	*end = 0;  // no worry, see VCF_BUFFER_SIZE_PLUS
	F64s.clear();

	while (p < end)
	{
		while ((p < end) && (*p == ' ')) p ++;
		const char *start = p;
		double val;

		bool is_dot = false;
		if ((p < end) && (*p == '.'))
		{
			const char *s = p + 1;
			while ((s < end) && (*s == ' ')) s ++;
			is_dot = (s >= end) || (*s == ',');
		}

		if (!is_dot)
		{
			char *endptr = (char*)p;
			val = strtod(p, &endptr);

			if (endptr == p)
			{
				if (raise_error)
					throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
				val = R_NaN;
				while ((p < end) && (*p != ',')) p ++;
			} else {
				p = endptr;
				while ((p < end) && (*p == ' ')) p ++;
				if ((p < end) && (*p != ','))
				{
					if (raise_error)
						throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
					val = R_NaN;
					while ((p < end) && (*p != ',')) p ++;
				}
			}
		} else {
			val = R_NaN;
		}

		F64s.push_back(val);
		if ((p < end) && (*p == ',')) p ++;
	}
}


/// get an integer from a string
static void getStringArray(char *p, char *end, vector<string> &UTF8s)
{
	UTF8s.clear();
	while (p < end)
	{
		while ((p < end) && (*p == ' ')) p ++;
		const char *s = p;

		while ((p < end) && (*p != ',')) p ++;
		const char *e = p;
		while ((s < e) && (*(e-1) == ' ')) e --;

		UTF8s.push_back(string(s, e));
		if ((p < end) && (*p == ',')) p ++;
	}
}


static const string BlackString;



}


extern "C"
{
using namespace SeqArray;

// ===========================================================
// Convert from VCF4: VCF4 -> GDS
// ===========================================================

/// return true, if matching
inline static bool StrCaseCmp(const char *prefix, const char *txt, size_t nmax)
{
	while (*prefix && *txt && nmax>0)
	{
		if (toupper(*prefix) != toupper(*txt))
			return false;
		prefix ++; txt ++; nmax --;
	}
	return (*prefix == 0);
}


/// VCF format --> SeqArray GDS format
COREARRAY_DLL_EXPORT SEXP SEQ_Parse_VCF(SEXP vcf_fn, SEXP header,
	SEXP gds_root, SEXP param, SEXP rho)
{
	const char *fn = CHAR(STRING_ELT(vcf_fn, 0));

	COREARRAY_TRY

		// the number of calling PROTECT
		int nProtected = 0;

		// =========================================================
		// initialize variables		

		// VCF file buffer
		VCF_Buffer.resize(VCF_BUFFER_SIZE + VCF_BUFFER_SIZE_PLUS);
		VCF_Buffer_EndPtr = VCF_Buffer_Ptr = &VCF_Buffer[0];

		// Text buffer
		Text_Buffer.resize(1024);
		Text_Size = 1024;
		Text_BeginPtr = Text_EndPtr = &Text_Buffer[0];

		// line and column numbers
		VCF_LineNum = VCF_ColumnNum = 0;
		VCF_NextLineNum = VCF_NextColumnNum = 1;

		// the total number of samples
		size_t nTotalSamp = Rf_asInteger(RGetListElement(param, "sample.num"));
		// the variable name for genotypic data
		string geno_id = CHAR(STRING_ELT(RGetListElement(param, "genotype.var.name"), 0));
		// raise an error
		bool RaiseError = (Rf_asLogical(RGetListElement(param, "raise.error")) == TRUE);
		// variant start
		int variant_start = Rf_asInteger(RGetListElement(param, "start"));
		// variant count
		int variant_count = Rf_asInteger(RGetListElement(param, "count"));
		// input file
		VCF_File = R_GetConnection(RGetListElement(param, "infile"));
		VCF_File->EOF_signalled = FALSE;
		// chromosome prefix
		SEXP ChrPrefix = RGetListElement(param, "chr.prefix");
		// verbose
		// bool Verbose = (LOGICAL(RGetListElement(param, "verbose"))[0] == TRUE);

		// the number of ploidy
		int num_ploidy = Rf_asInteger(RGetListElement(header, "ploidy"));
		if (num_ploidy <= 0)
			throw ErrSeqArray("Invalid header$ploidy: %d.", num_ploidy);

		// filter level list
		vector<string> filter_list;
		{
			SEXP level = RGetListElement(param, "filter.levels");
			const int n = RLength(level);
			for (int i=0; i < n; i++)
				filter_list.push_back(CHAR(STRING_ELT(level, i)));
		}

		// GDS nodes
		PdAbstractArray Root = GDS_R_SEXP2Obj(gds_root, FALSE);

		PdAbstractArray varIdx = GDS_Node_Path(Root, "variant.id", TRUE);
		PdAbstractArray varChr = GDS_Node_Path(Root, "chromosome", TRUE);
		PdAbstractArray varPos = GDS_Node_Path(Root, "position", TRUE);
		PdAbstractArray varRSID = GDS_Node_Path(Root, "annotation/id", TRUE);
		PdAbstractArray varAllele = GDS_Node_Path(Root, "allele", TRUE);

		PdAbstractArray varQual = GDS_Node_Path(Root, "annotation/qual", TRUE);
		PdAbstractArray varFilter = GDS_Node_Path(Root, "annotation/filter", TRUE);

		PdAbstractArray varGeno = GDS_Node_Path(Root, "genotype/data", TRUE);
		PdAbstractArray varGenoLen = GDS_Node_Path(Root, "genotype/@data", TRUE);
		PdAbstractArray varGenoExtraIdx = GDS_Node_Path(Root, "genotype/extra.index", TRUE);
		PdAbstractArray varGenoExtra = GDS_Node_Path(Root, "genotype/extra", TRUE);

		const int GenoNumBits= GDS_Array_GetBitOf(varGeno);
		if (GenoNumBits != 2)
			throw ErrSeqArray("Invalid data type in genotype/data, it should be bit2.");
		const int GenoBitMask = ~((-1) << GenoNumBits);

		PdAbstractArray varPhase, varPhaseExtraIdx, varPhaseExtra;
		if (num_ploidy > 1)
		{
			varPhase = GDS_Node_Path(Root, "phase/data", TRUE);
			varPhaseExtraIdx = GDS_Node_Path(Root, "phase/extra.index", TRUE);
			varPhaseExtra = GDS_Node_Path(Root, "phase/extra", TRUE);
		} else {
			varPhase = varPhaseExtraIdx = varPhaseExtra = NULL;
		}

		// INFO
		vector<TVCF_Field_Info> info_list;
		set<string> info_missing;
		{
			SEXP info = RGetListElement(header, "info");
			SEXP info_ID = RGetListElement(info, "ID");
			SEXP info_inttype = RGetListElement(info, "int_type");
			SEXP info_intnum = RGetListElement(info, "int_num");
			SEXP info_flag = RGetListElement(info, "import.flag");
			TVCF_Field_Info val;

			for (size_t i=0; i < RLength(info_ID); i++)
			{
				val.name = CHAR(STRING_ELT(info_ID, i));
				val.type = INTEGER(info_inttype)[i];
				val.import_flag = (LOGICAL(info_flag)[i] == TRUE);
				val.number = INTEGER(info_intnum)[i];
				val.data_obj = GDS_Node_Path(Root,
					(string("annotation/info/") + val.name).c_str(), FALSE);
				val.len_obj = GDS_Node_Path(Root,
					(string("annotation/info/@") + val.name).c_str(), FALSE);

				info_list.push_back(val);
			}
		}

		// FORMAT
		vector<TVCF_Field_Format> format_list;
		set<string> format_missing;
		{
			SEXP fmt = RGetListElement(header, "format");
			SEXP fmt_ID = RGetListElement(fmt, "ID");
			SEXP fmt_inttype = RGetListElement(fmt, "int_type");
			SEXP fmt_intnum = RGetListElement(fmt, "int_num");
			SEXP fmt_flag = RGetListElement(fmt, "import.flag");
			TVCF_Field_Format val;

			for (size_t i=0; i < RLength(fmt_ID); i++)
			{
				val.name = CHAR(STRING_ELT(fmt_ID, i));
				val.type = INTEGER(fmt_inttype)[i];
				val.import_flag = (LOGICAL(fmt_flag)[i] == TRUE);
				val.number = INTEGER(fmt_intnum)[i];
				val.data_obj = GDS_Node_Path(Root,
					(string("annotation/format/") + val.name + "/data").c_str(), FALSE);
				val.len_obj = GDS_Node_Path(Root,
					(string("annotation/format/") + val.name + "/@data").c_str(), FALSE);
				format_list.push_back(val);

				switch (val.type)
				{
					case FIELD_TYPE_INT:
						format_list.back().I32ss.resize(nTotalSamp);
						break;
					case FIELD_TYPE_FLOAT:
						format_list.back().F64ss.resize(nTotalSamp);
						break;
					case FIELD_TYPE_STRING:
						format_list.back().UTF8ss.resize(nTotalSamp);
						break;
					default:
						throw ErrSeqArray("Invalid FORMAT Type.");
				}
			}
		}

		// variant id (integer)
		C_Int32 variant_index = GDS_Array_GetTotalCount(varIdx);

		// the string buffer
		string cell, value;
		cell.reserve(1024);
		value.reserve(4096);

		// the numeric buffer
		C_Int32 I32;
		vector<C_Int32> I32s;
		I32s.reserve(nTotalSamp);

		C_Float64 F64;
		vector<C_Float64> F64s;
		F64s.reserve(nTotalSamp);

		// the string buffer
		vector<string> StrList;
		StrList.reserve(nTotalSamp);

		// genotypes
		vector< vector<C_Int16> > Geno;
		Geno.resize(nTotalSamp);
		vector<C_Int8> I8s;
		I8s.reserve(nTotalSamp * num_ploidy);

		vector< TVCF_Field_Info >::iterator pI;
		vector< TVCF_Field_Format* >::iterator pF;
		vector< TVCF_Field_Format* > fmt_ptr;
		fmt_ptr.reserve(format_list.size());

		// chromosome prefix
		vector<const char *> ChrPref;
		for (size_t i=0; i < RLength(ChrPrefix); i++)
			ChrPref.push_back(CHAR(STRING_ELT(ChrPrefix, i)));


		// =========================================================
		// skip the header

		while (!VCF_File->EOF_signalled)
		{
			GetText(NA_INTEGER);
			if (strncmp(Text_BeginPtr, "#CHROM", 6) == 0)
			{
				SkipLine();
				break;
			}
		}

		// =========================================================
		// parse the context

		while (!VCF_EOF())
		{
			// -----------------------------------------------------
			// column 1: CHROM
			GetText(FALSE);

			// -----------------------------------------------------
			// variant id
			variant_index ++;
			if (variant_index < variant_start)
			{
				SkipLine();
				continue;
			} else if (variant_count >= 0)
			{
				if (variant_index >= variant_start+variant_count)
					break;
			}
			GDS_Array_AppendData(varIdx, 1, &variant_index, svInt32);


			// column 1: CHROM
			for (vector<const char *>::iterator p=ChrPref.begin(); p != ChrPref.end(); p++)
			{
				if (StrCaseCmp(*p, Text_BeginPtr, Text_EndPtr-Text_BeginPtr))
				{
					Text_BeginPtr += strlen(*p);
					break;
				}
			}
			GDS_Array_AppendStrLen(varChr, Text_BeginPtr, Text_EndPtr-Text_BeginPtr);


			// -----------------------------------------------------
			// column 2: POS
			GetText(FALSE);
			I32 = getInt32(Text_BeginPtr, Text_EndPtr, RaiseError);
			GDS_Array_AppendData(varPos, 1, &I32, svInt32);


			// -----------------------------------------------------
			// column 3: ID
			GetText(FALSE);
			SkipDot();
			GDS_Array_AppendStrLen(varRSID, Text_BeginPtr, Text_EndPtr-Text_BeginPtr);


			// -----------------------------------------------------
			// column 4 & 5: REF + ALT 
			GetText(FALSE);  // REF
			SkipWhiteSpace();
			cell.assign(Text_BeginPtr, Text_EndPtr);

			GetText(FALSE);  // ALT
			SkipDot();
			if (Text_EndPtr > Text_BeginPtr)
			{
				cell.push_back(',');
				cell.append(Text_BeginPtr, Text_EndPtr);
			}

			GDS_Array_AppendData(varAllele, 1, &value, svStrUTF8);

			// determine how many alleles
			int num_allele = 0;
			for (const char *p = cell.c_str(); *p; )
			{
				num_allele ++;
				while (*p && (*p != ',')) p ++;
				if (*p == ',') p ++;
			}


			// -----------------------------------------------------
			// column 6: QUAL
			GetText(FALSE);
			F64 = getFloat(Text_BeginPtr, Text_EndPtr, RaiseError);
			GDS_Array_AppendData(varQual, 1, &F64, svFloat64);


			// -----------------------------------------------------
			// column 7: FILTER
			GetText(FALSE);
			SkipDot();
			if (Text_BeginPtr < Text_EndPtr)
			{
				cell.assign(Text_BeginPtr, Text_EndPtr);
				vector<string>::iterator p =
					std::find(filter_list.begin(), filter_list.end(), cell);
				if (p == filter_list.end())
				{
					filter_list.push_back(cell);
					I32 = filter_list.size();
				} else
					I32 = p - filter_list.begin() + 1;
			} else
				I32 = NA_INTEGER;
			GDS_Array_AppendData(varFilter, 1, &I32, svInt32);


			// -----------------------------------------------------
			// column 8: INFO

			// initialize
			for (pI = info_list.begin(); pI != info_list.end(); pI++)
				pI->used = false;
			GetText(FALSE);
			SkipDot();

			// parse
			while (Text_BeginPtr < Text_EndPtr)
			{
				// format: name=val | name
				char *s, *p;
				s = p = Text_BeginPtr;
				while ((p < Text_EndPtr) && (*p != ';') && (*p != '='))
					p ++;
				Text_BeginPtr = p;

				// variable name
				while ((s < p) && (*(p-1) == ' ')) p --;
				cell.assign(s, p);

				// variable value
				char *ValBegin, *ValEnd;
				ValBegin = ValEnd = p = Text_BeginPtr;
				if (p < Text_EndPtr)
				{
					if (*p == '=')
					{
						p ++;
						while ((p < Text_EndPtr) && (*p == ' ')) p ++;
						ValBegin = p;
						while ((p < Text_EndPtr) && (*p != ';')) p ++;
						Text_BeginPtr = p;
						if (p < Text_EndPtr) Text_BeginPtr ++;
						while ((ValBegin < p) && (*(p-1) == ' ')) p --;
						ValEnd = p;
					} else if (*p == ';')
						Text_BeginPtr = p + 1;
					else
						Text_BeginPtr = p;
				}

				for (pI=info_list.begin(); pI != info_list.end(); pI++)
				{
					if (pI->name == cell)
						break;
				}

				if (pI != info_list.end())
				{
					// it is in the list of INFO variables
					if (pI->used)
					{
						Rf_warning("LINE: %d, ignore duplicated INFO ID (%s).",
							VCF_LineNum, cell.c_str());
						continue;
					}

					if (pI->import_flag)
					{
						switch (pI->type)
						{
						case FIELD_TYPE_INT:
							getInt32Array(ValBegin, ValEnd, I32s, RaiseError);
							pI->Index(I32s, num_allele);
							GDS_Array_AppendData(pI->data_obj, I32s.size(),
								&(I32s[0]), svInt32);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(ValBegin, ValEnd, F64s, RaiseError);
							pI->Index(F64s, num_allele);
							GDS_Array_AppendData(pI->data_obj, F64s.size(),
								&(F64s[0]), svFloat64);
							break;

						case FIELD_TYPE_FLAG:
							if (ValBegin < ValEnd)
							{
								throw ErrSeqArray(
									"INFO ID '%s' should be a flag without values.",
									cell.c_str());
							}
							I32 = 1;
							GDS_Array_AppendData(pI->data_obj, 1, &I32, svInt32);
							break;

						case FIELD_TYPE_STRING:
							getStringArray(ValBegin, ValEnd, StrList);
							pI->Index(StrList, num_allele);
							GDS_Array_AppendData(pI->data_obj, StrList.size(),
								&(StrList[0]), svStrUTF8);
							break;

						default:
							throw ErrSeqArray("Invalid INFO Type.");
						}
					}

					pI->used = true;
				} else {
					set<string>::iterator it = info_missing.find(cell);
					if (it == info_missing.end())
					{
						info_missing.insert(cell);
						Rf_warning("Unknown INFO ID '%s' is ignored.", cell.c_str());
					}
				}
			}

			// for which does not exist
			for (pI = info_list.begin(); pI != info_list.end(); pI++)
			{
				if (!pI->used && pI->import_flag)
				{
					switch (pI->type)
					{
					case FIELD_TYPE_INT:
						pI->Fill(I32s, NA_INTEGER);
						break;
					case FIELD_TYPE_FLOAT:
						pI->Fill(F64s, R_NaN);
						break;
					case FIELD_TYPE_FLAG:
						I32 = 0;
						GDS_Array_AppendData(pI->data_obj, 1, &I32, svInt32);
						break;
					case FIELD_TYPE_STRING:
						pI->Fill(StrList, string());
						break;
					default:
						throw ErrSeqArray("Invalid INFO Type.");
					}
					pI->used = true;
				}
			}


SkipLine();


			// -----------------------------------------------------
			// column 9: FORMAT

/*			// initialize
			for (pF = fmt_ptr.begin(); pF != fmt_ptr.end(); pF++)
				(*pF)->used = false;
			GetText(false);

			// parse
			bool first_id_flag = true;
			fmt_ptr.clear();
			pCh = cell.c_str();
			while (*pCh)
			{
				name.clear();
				while ((*pCh != 0) && (*pCh != ':'))
					{ name.push_back(*pCh); pCh ++; }
				if (*pCh == ':') pCh ++;
				_Trim_(name);

				if (first_id_flag)
				{
					// genotype ID
					if (name != geno_id)
					{
						throw ErrSeqArray("The first FORMAT ID should be '%s'.",
							geno_id.c_str());
					}
					first_id_flag = false;

				} else {
					// find ID
					vector<TVCF_Field_Format>::iterator it;
					for (it = format_list.begin(); it != format_list.end(); it++)
					{
						if (it->name == name)
							{ it->used = true; break; }
					}
					if (it == format_list.end())
					{
						set<string>::iterator it = format_missing.find(name);
						if (it == format_missing.end())
						{
							format_missing.insert(name);
							warning("Unknown FORMAT ID '%s' is ignored.",
								name.c_str());
						}
					} else {
						// push
						fmt_ptr.push_back(&(*it));
					}
				}
			}


			// -----------------------------------------------------
			// columns for samples

			// for-loop
			for (int samp_idx=0; samp_idx < nTotalSamp; samp_idx ++)
			{
				const char *p;

				// read
				GetText(samp_idx >= (nTotalSamp-1));

				// -------------------------------------------------
				// the first field -- genotypes
				pCh = p = cell.c_str();
				while ((*p != 0) && (*p != ':')) p ++;
				value.assign(pCh, p);
				pCh = (*p == ':') ? (p + 1) : p;

				I32s.clear();
				vector<C_Int16> &pAllele = Geno[samp_idx];
				pAllele.clear();

				p = SKIP(value.c_str());
				while (*p)
				{
					char *endptr = (char*)p;
					C_Int32 val = strtol(p, &endptr, 10);

					if (endptr == p)
					{
						if ((*p != '.') && RaiseError)
						{
							throw ErrSeqArray(
								"Invalid integer conversion \"%s\".",
								SHORT_TEXT(p).c_str());
						}
						val = -1;
					} else {
						if (val < 0)
						{
							val = -1;
							if (RaiseError)
							{
								throw ErrSeqArray(
									"Genotype code should be non-negative \"%s\".",
									SHORT_TEXT(p).c_str());
							}
						} else if (val >= num_allele)
						{
							val = -1;
							if (RaiseError)
							{
								throw ErrSeqArray(
									"Genotype code is out of range \"%s\".",
									SHORT_TEXT(p).c_str());
							}
						}
						p = endptr;
					}

					pAllele.push_back(val);
					while ((*p != 0) && (*p != '|') && (*p != '/'))
						p ++;
					if (*p == '|')
					{
						I32s.push_back(1); p ++;
					} else if (*p == '/')
					{
						I32s.push_back(0); p ++;
					}
				}

				// check pAllele
				if ((int)pAllele.size() < num_ploidy)
					pAllele.resize(num_ploidy, -1);

				if (varPhase)
				{
					// write phasing information
					if ((int)I32s.size() < (num_ploidy-1))
						I32s.resize(num_ploidy-1, 0);
					GDS_Array_AppendData(varPhase, num_ploidy-1, &(I32s[0]), svInt32);
					if ((int)I32s.size() > (num_ploidy-1))
					{
						// E.g., triploid call: 0/0/1
						int Len = num_ploidy - (int)I32s.size() - 1;
						GDS_Array_AppendData(varPhaseExtra, Len, &(I32s[num_ploidy-1]), svInt32);
						I32 = samp_idx + 1;
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
						I32 = variant_index;
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
						I32 = Len;
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
					}
				}

				// -------------------------------------------------
				// the other field -- format id
				for (size_t i=0; i < fmt_ptr.size(); i++)
				{
					TVCF_Field_Format *pFmt = fmt_ptr[i];
					p = pCh;
					while ((*p != 0) && (*p != ':')) p ++;

					if ((pFmt!=NULL) && pFmt->import_flag)
					{
						// get the field context
						value.assign(pCh, p);

						// parse the field
						switch (pFmt->type)
						{
						case FIELD_TYPE_INT:
							getInt32Array(value, pFmt->I32ss[samp_idx], RaiseError);
							pFmt->Check(pFmt->I32ss[samp_idx], num_allele, NA_INTEGER);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(value, pFmt->F64ss[samp_idx], RaiseError);
							pFmt->Check(pFmt->F64ss[samp_idx], num_allele, R_NaN);
							break;

						case FIELD_TYPE_STRING:
							getStringArray(value, pFmt->UTF8ss[samp_idx]);
							pFmt->Check(pFmt->UTF8ss[samp_idx], num_allele, BlackString);
							break;

						default:
							throw ErrSeqArray("Invalid FORMAT Type.");
						}
					}

					pCh = (*p == ':') ? (p + 1) : p;
				}
			}

			// -------------------------------------------------
			// write genotypes

			// determine how many bits (GenoNumBits = 2)
			int num_bits = GenoNumBits;
			// plus ONE for missing value
			while ((num_allele + 1) > (1 << num_bits))
				num_bits += GenoNumBits;
			I32 = num_bits / GenoNumBits;
			GDS_Array_AppendData(varGenoLen, 1, &I32, svInt32);

			// write to the variable "genotype"
			for (int bits=0; bits < num_bits; bits += GenoNumBits)
			{
				I8s.clear();
				for (int i=0; i < nTotalSamp; i++)
				{
					vector<C_Int16> &pAllele = Geno[i];
					for (int j=0; j < num_ploidy; j++)
						I8s.push_back((pAllele[j] >> bits) & GenoBitMask);
				}
				GDS_Array_AppendData(varGeno, nTotalSamp*num_ploidy,
					&(I8s[0]), svInt8);
			}

			// write to "genotype/extra"
			for (int i=0; i < nTotalSamp; i++)
			{
				vector<C_Int16> &pGeno = Geno[i];
				if ((int)pGeno.size() > num_ploidy)
				{
					// E.g., triploid call: 0/0/1
					int Len = num_ploidy - (int)pGeno.size() - 1;
					GDS_Array_AppendData(varGenoExtra, Len, &(pGeno[num_ploidy-1]), svInt32);
					I32 = i + 1;
					GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
					I32 = variant_index;
					GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
					I32 = Len;
					GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
				}
			}


			// -------------------------------------------------
			// for-loop all format IDs: write
			for (vector<TVCF_Field_Format>::iterator it = format_list.begin();
				it != format_list.end(); it++)
			{
				if (it->import_flag)
				{
					if (it->used)
					{
						if (it->number > 0)
						{
							// fixed-length array
							it->WriteFixedLength();
							I32 = 1;
						} else if (it->number < 0)
						{
							// variable-length array
							I32 = it->WriteVariableLength(nTotalSamp, I32s, F64s);
						} else
							throw ErrSeqArray("Invalid FORMAT Number.");
						GDS_Array_AppendData(it->len_obj, 1, &I32, svInt32);
					} else {
						I32 = 0;
						GDS_Array_AppendData(it->len_obj, 1, &I32, svInt32);
					}
				}
			}
*/

		}

		// set returned value: levels(filter)
		PROTECT(rv_ans = NEW_CHARACTER(filter_list.size()));
		for (int i=0; i < (int)filter_list.size(); i++)
			SET_STRING_ELT(rv_ans, i, mkChar(filter_list[i].c_str()));
		nProtected ++;

		UNPROTECT(nProtected);

		VCF_File = NULL;
		VCF_Buffer.clear();
		Text_Buffer.clear();

	CORE_CATCH({
		char buf[4096];
		if (VCF_ColumnNum > 0)
		{
			snprintf(buf, sizeof(buf),
				"%s\nFILE: %s\nLINE: %lld, COLUMN: %d, %s\n",
				GDS_GetError(), fn, VCF_LineNum, VCF_ColumnNum,
				string(Text_BeginPtr, Text_EndPtr).c_str());
		} else {
			snprintf(buf, sizeof(buf), "%s\nFILE: %s\nLINE: %lld\n",
				GDS_GetError(), fn, VCF_LineNum);
		}
		GDS_SetError(buf);
		has_error = true;
	});
	if (has_error) error(GDS_GetError());

	// output
	return(rv_ans);
}

} // extern "C"
