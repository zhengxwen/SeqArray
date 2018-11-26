// ===========================================================
//
// ConvVCF2GDS.cpp: format conversion from VCF to GDS
//
// Copyright (C) 2013-2018    Xiuwen Zheng
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

#include "Index.h"
#include "vectorization.h"
#include <vector>
#include <set>
#include <algorithm>


// 1: load format field except GT
// 2: save format field to GDS file except GT (~70%)
// 3: load and save info field (< 3%)
#define GDS_TIMING	0

#if (GDS_TIMING > 0)
#   include <time.h>
#endif


namespace SeqArray
{
using namespace std;


#if (GDS_TIMING > 0)

static clock_t _timing_ = 0;
static clock_t _timing_last_point;

static inline void start_timing()
{
	_timing_last_point = clock();
}
static inline void end_timing()
{
	clock_t t = clock();
	_timing_ += t - _timing_last_point;
	_timing_last_point = t;
}

#endif


// ===========================================================
// define 
// ===========================================================

#if R_CONNECTIONS_VERSION != 1
#error "No support of R_CONNECTIONS_VERSION"
#endif

static Rconnection VCF_File = NULL;  ///< R connection object

static vector<char> VCF_Buffer;  ///< reading buffer
static char *VCF_Buffer_Ptr;     ///< the current pointer to reading buffer
static char *VCF_Buffer_EndPtr;  ///< the end pointer to reading buffer
static const size_t VCF_BUFFER_SIZE = 65536;  ///< reading buffer size
static const size_t VCF_BUFFER_SIZE_PLUS = 32;  ///< additional buffer is needed since *VCF_Buffer_EndPtr might be revised

/// initialize
inline static void Init_VCF_Buffer(SEXP File)
{
	VCF_File = R_GetConnection(File);
	VCF_File->EOF_signalled = FALSE;
	VCF_Buffer.resize(VCF_BUFFER_SIZE + VCF_BUFFER_SIZE_PLUS);
	VCF_Buffer_EndPtr = VCF_Buffer_Ptr = &VCF_Buffer[0];
}

/// finalize
inline static void Done_VCF_Buffer()
{
	VCF_File = NULL;
	VCF_Buffer.clear();
	vector<char>().swap(VCF_Buffer);
	VCF_Buffer_Ptr = VCF_Buffer_EndPtr = NULL;
}


/// read file buffer
inline static void Read_VCF_Buffer()
{
	VCF_Buffer_Ptr = &VCF_Buffer[0];
	size_t n = 0;
	size_t unread_len = VCF_File->buff_stored_len - VCF_File->buff_pos;
	if (unread_len > 0)
	{
		if (unread_len > VCF_BUFFER_SIZE) unread_len = VCF_BUFFER_SIZE;
		memcpy(VCF_Buffer_Ptr, VCF_File->buff + VCF_File->buff_pos, unread_len);
		VCF_Buffer_Ptr += unread_len;
		VCF_File->buff_pos += unread_len;
		n += unread_len;
	}
	if (n < VCF_BUFFER_SIZE)
	{
		size_t m = R_ReadConnection(VCF_File, VCF_Buffer_Ptr, VCF_BUFFER_SIZE-n);
		n += m;
	}
	VCF_Buffer_Ptr = &VCF_Buffer[0];
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
static char *Text_pBegin;  ///< the starting pointer for the text buffer
static char *Text_pEnd;    ///< the end pointer for the text buffer
static size_t Text_Size;     ///< the buffer size in Text_Buffer
static char *save_pBegin, *save_pEnd;  ///< save Text_pBegin, Text_pEnd

static C_Int64 VCF_LineNum;  ///< the current line number
static int VCF_ColumnNum;    ///< the current column number
static C_Int64 VCF_NextLineNum;  ///< the next line number
static int VCF_NextColumnNum;    ///< the next column number

static bool VCF_RaiseError;  ///< raise conversion error if true
static size_t SampleNum;  ///< the total number of samples

/// initialize text buffer
inline static void InitText()
{
	Text_Buffer.resize(1024);
	Text_Size = 1024;
	Text_pBegin = Text_pEnd = &Text_Buffer[0];
	save_pBegin = save_pEnd = Text_pBegin;
	VCF_LineNum = VCF_ColumnNum = 0;
	VCF_NextLineNum = VCF_NextColumnNum = 1;
}

/// finalize text buffer
inline static void DoneText()
{
	Text_Buffer.clear();
	vector<char>().swap(Text_Buffer);
	Text_pBegin = Text_pEnd;
	save_pBegin = save_pEnd = Text_pBegin;
}

/// get a string with a seperator '\t', which is saved in _Text_Buffer
inline static void GetText(int last_column)
{
	if (VCF_File->EOF_signalled)
		throw ErrSeqArray("it is the end of file.");

	VCF_ColumnNum = VCF_NextColumnNum;
	VCF_LineNum = VCF_NextLineNum;

	Text_pEnd = Text_pBegin = &Text_Buffer[0];
	int ch = -1;
	bool flag = true;

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

		if (flag && (p < VCF_Buffer_EndPtr))
		{
			Text_pBegin = VCF_Buffer_Ptr;
			Text_pEnd = p;
			VCF_Buffer_Ptr = p;
			break;
		} else
			flag = false;

		// copy to Text_Buffer
		size_t n = p - VCF_Buffer_Ptr;
		size_t m = Text_pEnd - Text_pBegin;
		size_t nn = m + n;
		if (nn > Text_Size)
		{
			nn = ((nn / 1024) + 1) * 1024;
			Text_Buffer.resize(nn + VCF_BUFFER_SIZE_PLUS);
			Text_Size = nn;
			Text_pBegin = &Text_Buffer[0];
			Text_pEnd = Text_pBegin + m;
		}
		memcpy(Text_pEnd, VCF_Buffer_Ptr, n);
		VCF_Buffer_Ptr += n;
		Text_pEnd += n;

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
					if (flag)
					{   // copy to Text_Buffer
						size_t n = Text_pEnd - Text_pBegin;
						if (n > Text_Size)
						{
							size_t nn = ((n / 1024) + 1) * 1024;
							Text_Buffer.resize(nn + VCF_BUFFER_SIZE_PLUS);
							Text_Size = nn;
						}
						memcpy(&Text_Buffer[0], Text_pBegin, n);
						Text_pBegin = &Text_Buffer[0];
						Text_pEnd = Text_pBegin + n;
						flag = false;
					}
					Read_VCF_Buffer();
				}
				ch = *VCF_Buffer_Ptr;
			} while (ch=='\n' || ch=='\r');
		}
	}

	save_pBegin = Text_pBegin;
	save_pEnd = Text_pEnd;
}

/// skip the current line
inline static void SkipLine()
{
	VCF_ColumnNum = VCF_NextColumnNum;
	VCF_LineNum = VCF_NextLineNum;

	int ch = -1;
	// search for '\n' or '\r'
	while (true)
	{
		VCF_Buffer_Ptr = (char*)vec_char_find_CRLF(VCF_Buffer_Ptr,
			VCF_Buffer_EndPtr - VCF_Buffer_Ptr);

		if (VCF_Buffer_Ptr < VCF_Buffer_EndPtr)
		{
			ch = *VCF_Buffer_Ptr;
			break;
		} else if (!VCF_File->EOF_signalled)
			Read_VCF_Buffer();
		else
			break;
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
	save_pBegin = save_pEnd = Text_pBegin;
}

/// skip white space
inline static void SkipWhiteSpace()
{
	// skip whitespace
	while ((Text_pBegin < Text_pEnd) && (*Text_pBegin == ' '))
		Text_pBegin ++;
	while ((Text_pBegin < Text_pEnd) && (*(Text_pEnd-1) == ' '))
		Text_pEnd --;
}

/// skip a dot
inline static void SkipTextWithDot()
{
	SkipWhiteSpace();
	if ((Text_pEnd-Text_pBegin == 1) && (*Text_pBegin == '.'))
		Text_pBegin ++;
}



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
static const char *ERR_GENO_CONV = "Invalid integer conversion for genotypes '%s'";
static const char *ERR_GENO_OUT_RANGE = "Genotype is out of range '%s'";
static const char *ERR_FLOAT_CONV = "Invalid float conversion '%s'";

/// get an integer from a string
inline static int getInt32(const char *p, const char *end)
{
	while ((p < end) && (*p == ' '))
		p ++;
	const char *start = p;

	if ((p < end) && (*p == '.'))
	{
		p ++;
		while ((p < end) && (*p == ' ')) p ++;
		if ((p < end) && VCF_RaiseError)
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
			if (VCF_RaiseError && (val > INT_MAX))
				throw ErrSeqArray(ERR_INT_OUT_RANGE, SHORT(start, end).c_str());
		} else {
			if (ch == ' ')
			{
				while ((p < end) && (*p == ' '))
					p ++;
				if (p >= end) break;
			}
			if (VCF_RaiseError)
				throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
			return NA_INTEGER;
		}
	}

	return sign ? -val : val;
}

/// get multiple integers from a string
inline static void getInt32Array(const char *p, const char *end,
	vector<C_Int32> &I32s)
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
				if (VCF_RaiseError)
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
				if (VCF_RaiseError && (val > INT_MAX))
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
					if (VCF_RaiseError)
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

/// get a genotype from a string
inline static C_Int16 getGeno(const char *p, const char *end, int num_allele)
{
	const char *start = p;
	if ((p < end) && (*p == '.'))
	{
		p ++;
		while ((p < end) && (*p == ' ')) p ++;
		if ((p < end) && VCF_RaiseError)
			throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
		return -1;
	}

	int val = 0;
	while (p < end)
	{
		char ch = *p ++;
		if ('0' <= ch && ch <= '9')
		{
			val = val*10 + (ch - '0');
			if (VCF_RaiseError && (val > 32767))
				throw ErrSeqArray(ERR_INT_OUT_RANGE, SHORT(start, end).c_str());
		} else {
			if (ch == ' ')
			{
				while ((p < end) && (*p == ' '))
					p ++;
				if (p >= end) break;
			}
			if (VCF_RaiseError)
				throw ErrSeqArray(ERR_GENO_CONV, SHORT(start, end).c_str());
			return -1;
		}
	}

	if (val >= num_allele)
	{
		if (VCF_RaiseError)
			throw ErrSeqArray(ERR_GENO_OUT_RANGE, SHORT(start, end).c_str());
		val = -1;
	}

	return val;
}


/// get a real number from a string
inline static double getFloat(char *p, char *end)
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
			if (VCF_RaiseError)
				throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
			val = R_NaN;
		} else {
			p = endptr;
			while ((p < end) && (*p == ' ')) p ++;
			if (p < end)
			{
				val = R_NaN;
				if (VCF_RaiseError)
					throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
			}
		}

		return val;
	} else {
		return R_NaN;
	}
}

/// get multiple real numbers from a string
inline static void getFloatArray(char *p, char *end, vector<double> &F64s)
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
			char *s = p + 1;
			while ((s < end) && (*s == ' ')) s ++;
			is_dot = (s >= end) || (*s == ',');
			if (is_dot) p = s;
		}

		if (!is_dot)
		{
			char *endptr = (char*)p;
			val = strtod(p, &endptr);

			if (endptr == p)
			{
				if (VCF_RaiseError)
					throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
				val = R_NaN;
				while ((p < end) && (*p != ',')) p ++;
			} else {
				p = endptr;
				while ((p < end) && (*p == ' ')) p ++;
				if ((p < end) && (*p != ','))
				{
					if (VCF_RaiseError)
						throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
					val = R_NaN;
					while ((p < end) && (*p != ',')) p ++;
				}
			}
		} else
			val = R_NaN;

		F64s.push_back(val);
		if ((p < end) && (*p == ',')) p ++;
	}
}


/// get an integer from a string
inline static void getStringArray(char *p, char *end, vector<string> &UTF8s)
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



// ===========================================================
// VCF structure
// ===========================================================

static string BlankString;

const static int FIELD_TYPE_INT      = 1;
const static int FIELD_TYPE_FLOAT    = 2;
const static int FIELD_TYPE_FLAG     = 3;
const static int FIELD_TYPE_STRING   = 4;


/// the structure of INFO field
struct COREARRAY_DLL_LOCAL TVCF_Info
{
	string name;       //< INFO ID
	int type;          //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;  //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;  //< according to 'Number=' field
	             // -1: variable-length (.), -2: # of alternate alleles (A)
	             // -3: # of possible genotypes (G), -4: # of alleles (R)
	bool used;   //< if TRUE, it has been parsed for the current line

	TVCF_Info()
	{
		type = 0;
		data_obj = len_obj = NULL;
		number = 0;
		used = false;
	}

	template<typename TYPE> inline void Index(vector<TYPE> &array,
		int num_allele, TYPE missing)
	{
		C_Int32 I32 = array.size();
		C_Int32 N;
		switch (number)
		{
		case -1:  // variable-length, .
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
			break;

		case -2:  // # of alternate alleles, A
			N = num_allele - 1;
			if (I32 > N)
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=A) should have %d value(s), but receives %d.",
					name.c_str(), N, I32);
			} else if (I32 < N)
				array.resize(N, missing);
			GDS_Array_AppendData(len_obj, 1, &N, svInt32);
			break;

		case -3:  // # of all possible genotypes, G
			N = (num_allele + 1) * num_allele / 2;
			if (I32 > N)
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=G) should have %d value(s), but receives %d.",
					name.c_str(), N, I32);
			} else if (I32 < N)
				array.resize(N, missing);
			GDS_Array_AppendData(len_obj, 1, &N, svInt32);
			break;

		case -4:  // # of alleles, R
			N = num_allele;
			if (I32 > N)
			{
				throw ErrSeqArray(
					"INFO ID '%s' (Number=R) should have %d value(s), but receives %d.",
					name.c_str(), N, I32);
			} else if (I32 < N)
				array.resize(N, missing);
			GDS_Array_AppendData(len_obj, 1, &N, svInt32);
			break;

		default:
			if (number >= 0)
			{
				N = array.size();
				if (N > number)
				{
					throw ErrSeqArray(
						"INFO ID '%s' should have %d value(s), but receives %d.",
						name.c_str(), number, N);
				} else if (N < number)
					array.resize(number, missing);
			} else
				throw ErrSeqArray("Invalid value 'number' in TVCF_Info.");
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
struct COREARRAY_DLL_LOCAL TVCF_Format
{
public:
	string name;       //< FORMAT ID
	int type;          //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;  //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;  //< according to 'Number=' field
	             // -1: variable-length (.), -2: # of alternate alleles (A)
	             // -3: # of possible genotypes (G), -4: # of alleles (R)
	bool used;   //< if TRUE, it has been parsed for the current line

private:

	vector<C_Int32> I32ss;
	vector<C_Float64> F64ss;
	vector<string> S8ss;
	size_t MaxCellNum;  //< the maximum number of elements per sample
	size_t CellNum;     //< 

	inline void Push_I32(int val, size_t samp_idx)
	{
		if (CellNum >= MaxCellNum)
		{
			MaxCellNum = CellNum + 1;
			I32ss.resize(MaxCellNum * SampleNum, NA_INTEGER);
		}
		I32ss[(CellNum++) * SampleNum + samp_idx] = val;
	}

	inline void Push_F64(double val, size_t samp_idx)
	{
		if (CellNum >= MaxCellNum)
		{
			MaxCellNum = CellNum + 1;
			F64ss.resize(MaxCellNum * SampleNum, R_NaN);
		}
		F64ss[(CellNum++) * SampleNum + samp_idx] = val;
	}

	inline void Push_S8(const string &val, size_t samp_idx)
	{
		if (CellNum >= MaxCellNum)
		{
			MaxCellNum = CellNum + 1;
			S8ss.resize(MaxCellNum * SampleNum, BlankString);
		}
		S8ss[(CellNum++) * SampleNum + samp_idx] = val;
	}

public:

	TVCF_Format()
	{
		type = 0;
		import_flag = used = false;
		data_obj = len_obj = NULL;
		number = 0;
		MaxCellNum = CellNum = 0;
	}

	void Init()
	{
		I32ss.clear(); F64ss.clear(); S8ss.clear();
		MaxCellNum = CellNum = 0;
		used = false;
	}


	/// get multiple integers from a string
	void GetInt32s(const char *p, const char *end, size_t samp_idx)
	{
		CellNum = 0;
		while (p < end)
		{
			while ((p < end) && (*p == ' ')) p ++;
			const char *start = p;

			if ((p < end) && (*p == '.'))
			{
				p ++;
				while ((p < end) && (*p == ' ')) p ++;
				if ((p < end) && (*p != ','))
				{
					if (VCF_RaiseError)
						throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
					while ((p < end) && (*p != ',')) p ++;
				}
				Push_I32(NA_INTEGER, samp_idx);
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
					if (VCF_RaiseError && (val > INT_MAX))
						throw ErrSeqArray(ERR_INT_OUT_RANGE, SHORT(start, end).c_str());
					p ++;
				} else {
					if (ch == ' ')
					{
						while ((p < end) && (*p == ' ')) p ++;
						if (p >= end) break;
						ch = *p;
					}

					if (ch == ',')
					{
						p ++;
						break;
					} else {
						if (VCF_RaiseError)
							throw ErrSeqArray(ERR_INT_CONV, SHORT(start, end).c_str());
						else
							while ((p < end) && (*p != ',')) p ++;
					}

					val = NA_INTEGER; sign = false;
					if ((p < end) && (*p == ',')) p ++;
					break;
				}
			}

			Push_I32(sign ? -val : val, samp_idx);
		}
	}

	/// get multiple real numbers from a string
	void GetFloats(char *p, char *end, size_t samp_idx)
	{
		while ((p < end) && (*(end-1) == ' ')) end --;
		*end = 0;  // no worry, see VCF_BUFFER_SIZE_PLUS
		CellNum = 0;

		while (p < end)
		{
			while ((p < end) && (*p == ' ')) p ++;
			const char *start = p;
			double val;

			bool is_dot = false;
			if ((p < end) && (*p == '.'))
			{
				char *s = p + 1;
				while ((s < end) && (*s == ' ')) s ++;
				is_dot = (s >= end) || (*s == ',');
				if (is_dot) p = s;
			}

			if (!is_dot)
			{
				char *endptr = (char*)p;
				val = strtod(p, &endptr);
				if (endptr == p)
				{
					if (VCF_RaiseError)
						throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
					val = R_NaN;
					while ((p < end) && (*p != ',')) p ++;
				} else {
					p = endptr;
					while ((p < end) && (*p == ' ')) p ++;
					if ((p < end) && (*p != ','))
					{
						if (VCF_RaiseError)
							throw ErrSeqArray(ERR_FLOAT_CONV, SHORT(start, end).c_str());
						val = R_NaN;
						while ((p < end) && (*p != ',')) p ++;
					}
				}
			} else
				val = R_NaN;

			Push_F64(val, samp_idx);
			if ((p < end) && (*p == ',')) p ++;
		}
	}

	/// get an integer from a string
	void GetStrings(char *p, char *end, size_t samp_idx)
	{
		CellNum = 0;
		while (p < end)
		{
			while ((p < end) && (*p == ' ')) p ++;
			const char *s = p;

			while ((p < end) && (*p != ',')) p ++;
			const char *e = p;
			while ((s < e) && (*(e-1) == ' ')) e --;

			Push_S8(string(s, e), samp_idx);
			if ((p < end) && (*p == ',')) p ++;
		}
	}

	/// checking
	inline void Check(size_t num_allele)
	{
		switch (number)
		{
		case -1:
			// Number=.
			break;
		case -2:
			// # of alternate alleles, Number=A
			num_allele --;
			if (CellNum > num_allele)
			{
				throw ErrSeqArray(
					"FORMAT ID '%s' (Number=A) should have %d value(s), but receives %d.",
					name.c_str(), (int)num_allele, (int)CellNum);
			}
			break;
		case -3:
			// # of all possible genotypes, Number=G
			num_allele = (num_allele+1) * num_allele / 2;
			if (CellNum > num_allele)
			{
				throw ErrSeqArray(
					"FORMAT ID '%s' (Number=G) should have %d value(s), but receives %d.",
					name.c_str(), (int)num_allele, (int)CellNum);
			}
			break;
		case -4:
			// # of alleles, Number=R
			if (CellNum > num_allele)
			{
				throw ErrSeqArray(
					"FORMAT ID '%s' (Number=R) should have %d value(s), but receives %d.",
					name.c_str(), (int)num_allele, (int)CellNum);
			}
			break;
		default:
			if (number >= 0)
			{
				if (CellNum > (size_t)number)
				{
					throw ErrSeqArray(
						"FORMAT ID '%s' should have %d value(s), but receives %d.",
						name.c_str(), number, (int)CellNum);
				}
			} else
				throw ErrSeqArray("Invalid value 'number' in TVCF_Format.");
		}
	}

	/// Save
	inline void SaveToGDS()
	{
		if (used)
		{
			size_t n = MaxCellNum * SampleNum;
			switch (type)
			{
			case FIELD_TYPE_INT:
				GDS_Array_AppendData(data_obj, n, &I32ss[0], svInt32);
				break;
			case FIELD_TYPE_FLOAT:
				GDS_Array_AppendData(data_obj, n, &F64ss[0], svFloat64);
				break;
			case FIELD_TYPE_STRING:
				GDS_Array_AppendData(data_obj, n, &S8ss[0], svStrUTF8);
				break;
			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
			}

			C_Int32 I32 = MaxCellNum;
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);

		} else {
			C_Int32 I32 = 0;
			GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
		}
	}
};

}


extern "C"
{
using namespace SeqArray;

// ===========================================================
// Get the number of lines in a VCF file
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_VCF_NumLines(SEXP File, SEXP SkipHead)
{
	Init_VCF_Buffer(File);

	if (Rf_asLogical(SkipHead) == TRUE)
	{
		InitText();
		// get the starting line
		while (!VCF_EOF())
		{
			GetText(NA_INTEGER);
			if (strncmp(Text_pBegin, "#CHROM", 6) == 0)
			{
				SkipLine();
				break;
			}
		}
		DoneText();
	}

	// get the number of left lines
	C_Int64 n = 0;
	while (!VCF_EOF())
	{
		n ++;
		SkipLine();
	}

	Done_VCF_Buffer();
	return ScalarReal(n);
}



// ===========================================================
// Split VCF files
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_VCF_Split(SEXP start, SEXP count, SEXP pnum,
	SEXP avoid_odd)
{
	int num = Rf_asInteger(pnum);
	bool no_odd = Rf_asLogical(avoid_odd)==TRUE;
	SEXP ans = PROTECT(NEW_LIST(2));
	SEXP start_array = PROTECT(NEW_NUMERIC(num));
	SEXP count_array = PROTECT(NEW_NUMERIC(num));
	SET_ELEMENT(ans, 0, start_array);
	SET_ELEMENT(ans, 1, count_array);

	double cnt = Rf_asReal(count);
	double scale = cnt / num;
	double st = Rf_asReal(start);
	for (int i=0; i < num; i++)
	{
		double st1 = round(st);
		REAL(start_array)[i] = st1;
		st += scale;

		C_Int64 m = (C_Int64)(round(st) - REAL(start_array)[i]);
		if (m & 0x01 && no_odd) // avoid odd number
			{ m ++; st ++; }
		if ((st1 + m) > (cnt + 1))
			m = round(cnt + 1 - st1);

		REAL(count_array)[i] = m;
	}

	UNPROTECT(3);
	return ans;
}



// ===========================================================
// Conversion: VCF --> GDS
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
COREARRAY_DLL_EXPORT SEXP SEQ_VCF_Parse(SEXP vcf_fn, SEXP header,
	SEXP gds_root, SEXP param, SEXP line_cnt, SEXP rho)
{
	const char *fn = CHAR(STRING_ELT(vcf_fn, 0));

	COREARRAY_TRY

	#if (GDS_TIMING > 0)
		_timing_ = 0;
		clock_t _start_time = clock();
	#endif

		// the number of calling PROTECT
		int nProtected = 0;

		// =========================================================
		// initialize variables		

		// the total number of samples
		SampleNum = Rf_asInteger(RGetListElement(param, "sample.num"));
		// the variable name for genotypic data
		string geno_id = CHAR(STRING_ELT(RGetListElement(param, "genotype.var.name"), 0));
		// raise an error
		VCF_RaiseError = (Rf_asLogical(RGetListElement(param, "raise.error")) == TRUE);
		// variant start
		C_Int64 variant_start = (C_Int64)Rf_asReal(RGetListElement(param, "start"));
		// variant count
		C_Int64 variant_count = (C_Int64)Rf_asReal(RGetListElement(param, "count"));
		// input file
		Init_VCF_Buffer(RGetListElement(param, "infile"));
		// chromosome prefix
		SEXP ChrPrefix = RGetListElement(param, "chr.prefix");
		// progress file
		SEXP progfile = RGetListElement(param, "progfile");
		// verbose
		// bool Verbose = (LOGICAL(RGetListElement(param, "verbose"))[0] == TRUE);

		// the number of ploidy
		size_t num_ploidy = Rf_asInteger(RGetListElement(header, "ploidy"));
		if (num_ploidy <= 0)
			throw ErrSeqArray("Invalid header$ploidy: %d.", (int)num_ploidy);

		size_t num_ploidy_less = num_ploidy - 1;
		size_t num_samp_ploidy = SampleNum * num_ploidy;
		size_t num_samp_ploidy_less = SampleNum * num_ploidy_less;

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

		PdAbstractArray varGeno = GDS_Node_Path(Root, "genotype/data", FALSE);
		PdAbstractArray varGenoLen = GDS_Node_Path(Root, "genotype/@data", TRUE);
		PdAbstractArray varGenoExtraIdx = GDS_Node_Path(Root, "genotype/extra.index", TRUE);
		PdAbstractArray varGenoExtra = GDS_Node_Path(Root, "genotype/extra", TRUE);

		const int GenoNumBits = varGeno ? GDS_Array_GetBitOf(varGeno) : 2;
		if (GenoNumBits != 2)
			throw ErrSeqArray("Invalid data type in genotype/data, it should be bit2.");

		PdAbstractArray varPhase, varPhaseExtraIdx, varPhaseExtra;
		if (num_ploidy > 1)
		{
			varPhase = GDS_Node_Path(Root, "phase/data", FALSE);
			varPhaseExtraIdx = GDS_Node_Path(Root, "phase/extra.index", FALSE);
			varPhaseExtra = GDS_Node_Path(Root, "phase/extra", FALSE);
		} else {
			varPhase = varPhaseExtraIdx = varPhaseExtra = NULL;
		}

		// INFO
		vector<TVCF_Info> info_list;
		set<string> info_missing;
		{
			SEXP info = RGetListElement(header, "info");
			SEXP info_ID = RGetListElement(info, "ID");
			SEXP info_inttype = RGetListElement(info, "int_type");
			SEXP info_intnum = RGetListElement(info, "int_num");
			SEXP info_flag = RGetListElement(info, "import.flag");
			TVCF_Info val;

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
		vector<TVCF_Format> format_list;
		set<string> format_missing;
		{
			SEXP fmt = RGetListElement(header, "format");
			SEXP fmt_ID = RGetListElement(fmt, "ID");
			SEXP fmt_inttype = RGetListElement(fmt, "int_type");
			SEXP fmt_intnum = RGetListElement(fmt, "int_num");
			SEXP fmt_flag = RGetListElement(fmt, "import.flag");
			TVCF_Format val;

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
			}
		}

		// variant id (integer)
		C_Int64 variant_index = (C_Int64)Rf_asReal(line_cnt);
		C_Int64 variant_start_index = variant_index;

		// the string buffer
		string cell;
		cell.reserve(1024);

		// the numeric buffer
		vector<C_Int8> I8s;
		I8s.reserve(SampleNum);

		C_Int32 I32;
		vector<C_Int32> I32s;
		I32s.reserve(SampleNum);

		C_Float64 F64;
		vector<C_Float64> F64s;
		F64s.reserve(SampleNum);

		// the string buffer
		vector<string> S8s;
		S8s.reserve(SampleNum);

		// genotypes
		vector<C_Int16> Genotypes;
		Genotypes.resize(num_samp_ploidy);

		// phase
		vector<C_Int8> Phases;
		Phases.resize(num_samp_ploidy_less);

		vector< TVCF_Info >::iterator pI;
		vector< TVCF_Format* >::iterator pF;
		vector< TVCF_Format* > fmt_ptr;
		fmt_ptr.reserve(format_list.size());

		// chromosome prefix
		vector<const char *> ChrPref;
		for (size_t i=0; i < RLength(ChrPrefix); i++)
			ChrPref.push_back(CHAR(STRING_ELT(ChrPrefix, i)));


		// =========================================================
		// skip the header and data rows

		InitText();

		if (!Rf_isNull(progfile))
		{
			while (!VCF_EOF())
			{
				GetText(NA_INTEGER);
				if (strncmp(Text_pBegin, "#CHROM", 6) == 0)
				{
					SkipLine();
					break;
				}
			}
		}

		while (!VCF_EOF() && (variant_index+1 < variant_start))
		{
			variant_index ++;
			SkipLine();
		}


		// =========================================================
		// parse the context

		// progress information
		CProgress Progress(variant_index - variant_start + 1, variant_count,
			progfile, true);

		while (!VCF_EOF())
		{
			// -----------------------------------------------------
			// column 1: CHROM
			GetText(FALSE);

			// -----------------------------------------------------
			// variant id
			variant_index ++;
			if (variant_count >= 0)
			{
				if (variant_index >= variant_start+variant_count)
				{
					variant_index --;
					break;
				}
			}
			I32 = variant_index;
			GDS_Array_AppendData(varIdx, 1, &I32, svInt32);

			// column 1: CHROM
			for (vector<const char *>::iterator p=ChrPref.begin(); p != ChrPref.end(); p++)
			{
				if (StrCaseCmp(*p, Text_pBegin, Text_pEnd-Text_pBegin))
				{
					Text_pBegin += strlen(*p);
					break;
				}
			}
			GDS_Array_AppendStrLen(varChr, Text_pBegin, Text_pEnd-Text_pBegin);


			// -----------------------------------------------------
			// column 2: POS
			GetText(FALSE);
			I32 = getInt32(Text_pBegin, Text_pEnd);
			GDS_Array_AppendData(varPos, 1, &I32, svInt32);


			// -----------------------------------------------------
			// column 3: ID
			GetText(FALSE);
			SkipTextWithDot();
			GDS_Array_AppendStrLen(varRSID, Text_pBegin, Text_pEnd-Text_pBegin);


			// -----------------------------------------------------
			// column 4 & 5: REF + ALT 
			GetText(FALSE);  // REF
			SkipWhiteSpace();
			cell.assign(Text_pBegin, Text_pEnd);

			GetText(FALSE);  // ALT
			SkipTextWithDot();
			if (Text_pEnd > Text_pBegin)
			{
				cell.push_back(',');
				cell.append(Text_pBegin, Text_pEnd);
			}

			GDS_Array_AppendData(varAllele, 1, &cell, svStrUTF8);

			// determine how many alleles
			int num_allele;
			if (cell != "." && cell != ".,.")
			{
				num_allele = 0;
				for (const char *p = cell.c_str(); *p; )
				{
					num_allele ++;
					while (*p && (*p != ',')) p ++;
					if (*p == ',') p ++;
				}
			} else {
				num_allele = INT_MAX;
			}


			// -----------------------------------------------------
			// column 6: QUAL
			GetText(FALSE);
			F64 = getFloat(Text_pBegin, Text_pEnd);
			GDS_Array_AppendData(varQual, 1, &F64, svFloat64);


			// -----------------------------------------------------
			// column 7: FILTER
			GetText(FALSE);
			SkipTextWithDot();
			if (Text_pBegin < Text_pEnd)
			{
				cell.assign(Text_pBegin, Text_pEnd);
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

		#if (GDS_TIMING == 3)
			start_timing();
		#endif

			// initialize
			for (pI = info_list.begin(); pI != info_list.end(); pI++)
				pI->used = false;
			GetText(SampleNum<=0);
			SkipTextWithDot();

			// parse
			while (Text_pBegin < Text_pEnd)
			{
				// format: name=val | name
				char *s, *p;
				s = p = Text_pBegin;
				while ((p < Text_pEnd) && (*p != ';') && (*p != '='))
					p ++;
				Text_pBegin = p;

				// variable name
				while ((s < p) && (*(p-1) == ' ')) p --;
				cell.assign(s, p);

				// variable value
				char *ValBegin, *ValEnd;
				ValBegin = ValEnd = p = Text_pBegin;
				if (p < Text_pEnd)
				{
					if (*p == '=')
					{
						p ++;
						while ((p < Text_pEnd) && (*p == ' ')) p ++;
						ValBegin = p;
						while ((p < Text_pEnd) && (*p != ';')) p ++;
						Text_pBegin = p;
						if (p < Text_pEnd) Text_pBegin ++;
						while ((ValBegin < p) && (*(p-1) == ' ')) p --;
						ValEnd = p;
					} else if (*p == ';')
						Text_pBegin = p + 1;
					else
						Text_pBegin = p;
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
							getInt32Array(ValBegin, ValEnd, I32s);
							pI->Index(I32s, num_allele, NA_INTEGER);
							GDS_Array_AppendData(pI->data_obj, I32s.size(),
								&(I32s[0]), svInt32);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(ValBegin, ValEnd, F64s);
							pI->Index(F64s, num_allele, R_NaN);
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
							getStringArray(ValBegin, ValEnd, S8s);
							pI->Index(S8s, num_allele, BlankString);
							GDS_Array_AppendData(pI->data_obj, S8s.size(),
								&(S8s[0]), svStrUTF8);
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
						Rf_warning("Unknown INFO ID '%s' is ignored (it should be defined in the meta-information lines).",
							cell.c_str());
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
						pI->Fill(S8s, BlankString);
						break;
					default:
						throw ErrSeqArray("Invalid INFO Type.");
					}
					pI->used = true;
				}
			}

		#if (GDS_TIMING == 3)
			end_timing();
		#endif


			// -----------------------------------------------------
			// column 9: FORMAT

			if (SampleNum <= 0) continue;

			// initialize
			for (pF = fmt_ptr.begin(); pF != fmt_ptr.end(); pF++)
				(*pF)->Init();
			GetText(FALSE);

			bool first_fmt_id_flag = true;
			bool first_fmt_id_is_geno = false;
			fmt_ptr.clear();

			// parse
			while (Text_pBegin < Text_pEnd)
			{
				while ((Text_pBegin<Text_pEnd) && (*Text_pBegin==' '))
					Text_pBegin ++;

				const char *start = Text_pBegin;
				while ((Text_pBegin<Text_pEnd) && (*Text_pBegin!=':'))
					Text_pBegin ++;

				const char *end = Text_pBegin;
				while ((start < end) && (*(end-1) == ' '))
					end --;
				cell.assign(start, end);

				if ((Text_pBegin<Text_pEnd) && (*Text_pBegin==':'))
					Text_pBegin ++;

				if (first_fmt_id_flag)
				{
					first_fmt_id_flag = false;
					// genotype ID
					first_fmt_id_is_geno = (cell == geno_id);
					if (first_fmt_id_is_geno) continue;
				}

				// find ID
				vector<TVCF_Format>::iterator p;
				for (p=format_list.begin(); p != format_list.end(); p++)
				{
					if (p->name == cell)
						{ p->used = true; break; }
				}
				if (p == format_list.end())
				{
					set<string>::iterator it = format_missing.find(cell);
					if (it == format_missing.end())
					{
						format_missing.insert(cell);
						warning("Unknown FORMAT ID '%s' is ignored (it should be defined in the meta-information lines).",
							cell.c_str());
					}
				} else {
					// push
					fmt_ptr.push_back(&(*p));
				}
			}


			// -----------------------------------------------------
			// Columns for samples

			// pointer to genotype buffer
			C_Int16 *pGeno = &Genotypes[0];
			// pointer to phase buffer
			C_Int8 *pPhase = &Phases[0];

			// for-loop
			for (size_t si=0; si < SampleNum; si ++)
			{
				// read
				GetText(si >= (SampleNum-1));

				// skip whitespace
				while ((Text_pBegin<Text_pEnd) && (*Text_pBegin==' '))
					Text_pBegin ++;

				if (first_fmt_id_is_geno)
				{
					// -------------------------------------------------
					// the first field -- genotypes (GT)

					const char *p = Text_pBegin;
					while ((Text_pBegin<Text_pEnd) && (*Text_pBegin!=':'))
						Text_pBegin ++;
					const char *end = Text_pBegin;

					if ((Text_pBegin<Text_pEnd) && (*Text_pBegin==':'))
						Text_pBegin ++;

					I32s.clear(); // genotype extra data
					I8s.clear(); // phase extra data
					size_t tmp_num_ploidy = 0;

					while (p < end)
					{
						const char *start = p;
						while ((p<end) && (*p!='|') && (*p!='/'))
							p ++;
						C_Int16 g = getGeno(start, p, num_allele);

						tmp_num_ploidy ++;
						if (tmp_num_ploidy <= num_ploidy)
							*pGeno ++ = g;
						else
							I32s.push_back(g);

						if (p < end)
						{
							C_Int8 v;
							if (*p == '|')
							{
								v = 1; p ++;
							} else if (*p == '/')
							{
								v = 0; p ++;
							}
							if (tmp_num_ploidy <= num_ploidy_less)
								*pPhase ++ = v;
							else
								I8s.push_back(v);
						}
					}

					for (size_t m=tmp_num_ploidy; m < num_ploidy; m++)
						*pGeno ++ = -1;
					for (size_t m=tmp_num_ploidy; m < num_ploidy_less; m++)
						*pPhase ++ = 0;

					// check "genotype/extra", e.g., triploid call: 0/0/1
					if (!I32s.empty())
					{
						GDS_Array_AppendData(varGenoExtra, I32s.size(), &(I32s[0]), svInt32);
						I32 = si + 1;
						GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
						I32 = variant_index - variant_start_index;
						GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
						I32 = I32s.size();
						GDS_Array_AppendData(varGenoExtraIdx, 1, &I32, svInt32);
					}

					// check "phase/extra", e.g., triploid call: 0/0/1
					if (varPhase && !I8s.empty())
					{
						GDS_Array_AppendData(varPhaseExtra, I8s.size(), &(I8s[0]), svInt8);
						I32 = si + 1;
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
						I32 = variant_index - variant_start_index;
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
						I32 = I8s.size();
						GDS_Array_AppendData(varPhaseExtraIdx, 1, &I32, svInt32);
					}
				}

				// -------------------------------------------------
				// the other field -- format id
				for (size_t i=0; i < fmt_ptr.size(); i++)
				{
					TVCF_Format *pFmt = fmt_ptr[i];

					char *start = Text_pBegin;
					while ((Text_pBegin<Text_pEnd) && (*Text_pBegin!=':'))
						Text_pBegin ++;
					char *end = Text_pBegin;

					if ((Text_pBegin<Text_pEnd) && (*Text_pBegin==':'))
						Text_pBegin ++;

					// parse the field
					if (pFmt && pFmt->import_flag)
					{
					#if (GDS_TIMING == 1)
						start_timing();
					#endif

						switch (pFmt->type)
						{
						case FIELD_TYPE_INT:
							pFmt->GetInt32s(start, end, si); break;
						case FIELD_TYPE_FLOAT:
							pFmt->GetFloats(start, end, si); break;
						case FIELD_TYPE_STRING:
							pFmt->GetStrings(start, end, si); break;
						default:
							throw ErrSeqArray("Invalid FORMAT Type.");
						}
						pFmt->Check(num_allele);

					#if (GDS_TIMING == 1)
						end_timing();
					#endif
					}
				}
			}


			if (first_fmt_id_is_geno)
			{
				// -------------------------------------------------
				// write genotypes

				// need to identify num_allele if missing
				if (num_allele == INT_MAX)
				{
					num_allele = 0;
					C_Int16 *p = &Genotypes[0];
					for (size_t n=num_samp_ploidy; n > 0; n--, p++)
					{
						if ((*p >= 0) && (*p > num_allele))
							num_allele = *p;
					}
					num_allele ++;
				}

				// determine how many bits
				int num_bits = 2;
				// plus ONE for missing value
				while ((num_allele + 1) > (1 << num_bits))
					num_bits += 2;
				I32 = num_bits >> 1;
				GDS_Array_AppendData(varGenoLen, 1, &I32, svInt32);

				for (int bits=0; bits < num_bits; )
				{
					GDS_Array_AppendData(varGeno, num_samp_ploidy, &Genotypes[0],
						svInt16);
					bits += 2;
					if (bits < num_bits)
						vec_i16_shr_b2(&Genotypes[0], num_samp_ploidy);
				}

				// -------------------------------------------------
				// write phase information

				if (varPhase)
				{
					GDS_Array_AppendData(varPhase, num_samp_ploidy_less,
						&(Phases[0]), svInt8);
				}
			}


			// -------------------------------------------------
			// for-loop all format IDs: write

		#if (GDS_TIMING == 2)
			start_timing();
		#endif
			for (vector<TVCF_Format>::iterator p = format_list.begin();
				p != format_list.end(); p ++)
			{
				if (p->import_flag) p->SaveToGDS();
			}
		#if (GDS_TIMING == 2)
			end_timing();
		#endif

			// update progress
			Progress.Forward();
		}

		// set returned value: levels(filter)
		PROTECT(rv_ans = NEW_CHARACTER(filter_list.size()));
		for (int i=0; i < (int)filter_list.size(); i++)
			SET_STRING_ELT(rv_ans, i, mkChar(filter_list[i].c_str()));
		nProtected ++;

		REAL(line_cnt)[0] = variant_index;

		UNPROTECT(nProtected);

		Done_VCF_Buffer();
		DoneText();


	#if (GDS_TIMING > 0)
		Rprintf("It took %0.2f seconds, in %0.2f%%.\n",
			((double)_timing_)/CLOCKS_PER_SEC,
			((double)_timing_) / (clock() - _start_time) * 100.0);
	#endif

	CORE_CATCH({
		char buf[4096];
		if ((VCF_ColumnNum > 0) && (save_pBegin < save_pEnd))
		{
			snprintf(buf, sizeof(buf),
				"%s\nFILE: %s\nLINE: %lld, COLUMN: %d, %s\n",
				GDS_GetError(), fn, (long long int)VCF_LineNum, VCF_ColumnNum,
				string(save_pBegin, save_pEnd).c_str());
		} else {
			snprintf(buf, sizeof(buf), "%s\nFILE: %s\nLINE: %lld\n",
				GDS_GetError(), fn, (long long int)VCF_LineNum);
		}
		GDS_SetError(buf);
		has_error = true;
	});
	if (has_error) error(GDS_GetError());

	// output
	return rv_ans;
}

} // extern "C"
