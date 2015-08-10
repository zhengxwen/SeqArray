// ===========================================================
//
// ConvVCF2GDS.cpp: the C++ code for the conversion from VCF to GDS
//
// Copyright (C) 2013-2015    Xiuwen Zheng
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



// ===========================================================
// define 
// ===========================================================

static const string BlackString;


// ===========================================================
// the structure of read line
// ===========================================================

/// a class of parsing text
class CReadLine
{
public:
	/// constructor
	CReadLine()
	{
		_ReadFun = _Rho = R_NilValue;
		_ptr_line = _lines.end();
		_ifend = false; _line_no = _column_no = 0;
		_cur_char = NULL;
		nProt = 0;
	}
	/// constructor
	CReadLine(SEXP vFun, SEXP vRho)
	{
		Init(vFun, vRho);
	}

	~CReadLine()
	{
		if (nProt > 0)
			UNPROTECT(nProt);
	}

	/// initialize R call
	void Init(SEXP vFun, SEXP vRho)
	{
		_ReadFun = vFun; _Rho = vRho;	
		_lines.clear(); _ptr_line = _lines.end();
		_ifend = false; _line_no = _column_no = 0;
		_cur_char = NULL;
		nProt = 0;
	}

	/// read a line
	const char *ReadLine()
	{
		if (_ifend) return NULL;
		if (_ptr_line == _lines.end())
		{
			if (_PrepareBuffer())
			{
				const char *rv = *_ptr_line;
				_ptr_line ++; _line_no ++;
				return rv;
			} else
				return NULL;
		} else {
			const char *rv = *_ptr_line;
			_ptr_line ++; _line_no ++;
			return rv;
		}
	}

	/// get a string with a seperator '\t'
	void GetCell(string &buffer, bool last_column)
	{
		if (_ifend)
			throw ErrSeqArray("It is the end.");
		if (!_cur_char)
		{
			_cur_char = ReadLine();
			if (!_cur_char)
				throw ErrSeqArray("It is the end.");
			_column_no = 0;
		}

		const char *str_begin = _cur_char;
		while ((*_cur_char != '\t') && (*_cur_char != 0))
			_cur_char ++;
		const char *str_end = _cur_char;
		_column_no ++;

		// check
		if ((str_begin == str_end) && (*_cur_char == 0))
			throw ErrSeqArray("fewer columns than what expected.");
		if (last_column)
		{
			if (*_cur_char != 0)
				throw ErrSeqArray("more columns than what expected.");
			_cur_char = NULL;
		} else {
			if (*_cur_char == '\t') _cur_char ++;
		}

		if (str_end > str_begin+1)
		{
			if ((str_begin[0] == '\"') && (str_end[-1] == '\"'))
			{
				str_begin ++;
				str_end --;
			} else if ((str_begin[0] == '\'') && (str_end[-1] == '\''))
			{
				str_begin ++;
				str_end --;
			}
		}

		buffer.assign(str_begin, str_end);
	}

	/// return true, if it is of the end
	bool IfEnd()
	{
		if (!_ifend)
		{
			if (_ptr_line == _lines.end())
				_PrepareBuffer();
		}
		return _ifend;
	}

	/// return line number
	COREARRAY_INLINE int LineNo() { return _line_no; }
	/// return column number
	COREARRAY_INLINE int ColumnNo() { return _column_no; }

protected:
	SEXP _ReadFun;  //< R call function
	SEXP _Rho;      //< R environment
	vector<const char *> _lines;               //< store returned string(s)
	vector<const char *>::iterator _ptr_line;  //< the pointer to _lines
	bool _ifend;     //< true for the end of reading
	int _line_no;    //< the index of current line
	int _column_no;  //< the index of current column
	const char *_cur_char;  //< 
	int nProt;

	bool _PrepareBuffer()
	{
		if (nProt > 0)
		{
			UNPROTECT(nProt);
			nProt = 0;
		}

		// call ReadLine R function
		SEXP val = eval(_ReadFun, _Rho);
		PROTECT(val);
		nProt ++;

		// check the returned value
		int n = Rf_length(val);
		if (n > 0)
		{
			_ifend = false;
			_lines.resize(n);
			for (int i=0; i < n; i++)
				_lines[i] = CHAR(STRING_ELT(val, i));
			_ptr_line = _lines.begin();
			return true;
		} else {
			_ifend = true;
			return false;
		}
	}
};




// ===========================================================
// VCF strcture
// ===========================================================

const static int FIELD_TYPE_INT      = 1;
const static int FIELD_TYPE_FLOAT    = 2;
const static int FIELD_TYPE_FLAG     = 3;
const static int FIELD_TYPE_STRING   = 4;


/// the structure of INFO field
struct TVCF_Field_Info
{
	string name;         //< INFO ID
	int type;            //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;    //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;          //< according to 'Number' field, if -1: variable-length, -2: # of alleles, -3: # of genotypes
	bool used;           //< if TRUE, it has been parsed for the current line

	TVCF_Field_Info()
	{
		type = 0;
		data_obj = len_obj = NULL;
		number = 0;
		used = false;
	}

	// INFO field

	template<typename TYPE> void Check(vector<TYPE> &array, string &name, int num_allele)
	{
		C_Int32 I32;
		switch (number)
		{
			case -1:
				// variable-length
				I32 = array.size();
				GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
				break;		

			case -2:
				// # of alleles
				I32 = array.size();
				if (I32 != (num_allele-1))
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), num_allele-1);
				}
				GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
				break;		

			case -3:
				// # of genotypes
				I32 = array.size();
				if (I32 != (num_allele+1)*num_allele/2)
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), (num_allele+1)*num_allele/2);
				}
				GDS_Array_AppendData(len_obj, 1, &I32, svInt32);
				break;		

			default:
				if (number >= 0)
				{
					if (number != (int)array.size())
					{
						throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
							name.c_str(), number);
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
struct TVCF_Field_Format
{
	string name;         //< FORMAT ID
	int type;            //< 1: integer, 2: float, 3: flag, 4: character,
	bool import_flag;    //< true: import, false: not import
	PdAbstractArray data_obj;  //< the pointer to data object
	PdAbstractArray len_obj;   //< can be NULL if variable-length object
	int number;          //< according to 'Number' field, if -1: variable-length, -2: # of alleles, -3: # of genotypes
	bool used;           //< if TRUE, it has been parsed for the current line

	/// data -- Int32
	vector< vector<C_Int32> > I32ss;
	/// data -- C_Float32
	vector< vector<C_Float32> > F32ss;
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

	template<typename TYPE> void Check(vector<TYPE> &array, string &name,
		int num_allele, const TYPE &missing)
	{
		switch (number)
		{
			case -1:
				break;

			case -2:
				// # of alleles
				if ((int)array.size() > (num_allele-1))
				{
					throw ErrSeqArray("FORMAT ID '%s' should have %d value(s).",
						name.c_str(), num_allele-1);
				} else {
					array.resize(num_allele-1, missing);
				}
				break;		

			case -3:
				// # of genotypes
				if ((int)array.size() > (num_allele+1)*num_allele/2)
				{
					throw ErrSeqArray("INFO ID '%s' should have %d value(s).",
						name.c_str(), (num_allele+1)*num_allele/2);
				} else {
					array.resize((num_allele+1)*num_allele/2, missing);
				}
				break;		

			default:
				if (number >= 0)
				{
					if ((int)array.size() > number)
					{
						throw ErrSeqArray("FORMAT ID '%s' should have %d value(s).",
							name.c_str(), number);
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
			throw ErrSeqArray("Wrong call 'WriteFixedLength' in TVCF_Field_Format.");
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
				for (vector< vector<float> >::iterator it = F32ss.begin();
					it != F32ss.end(); it ++)
				{
					GDS_Array_AppendData(data_obj, number, &((*it)[0]), svFloat32);
				}
				break;

			case FIELD_TYPE_STRING:
				for (vector< vector<string> >::iterator it = UTF8ss.begin();
					it != UTF8ss.end(); it ++)
				{
					for (int j=0; j < (int)(*it).size(); j ++)
						GDS_Array_AppendString(data_obj, (*it)[j].c_str());
				}
				break;

			default:
				throw ErrSeqArray("Invalid FORMAT Type.");
		}
	}

	int WriteVariableLength(int nTotalSample, vector<C_Int32> &I32s,
		vector<float> &F32s)
	{
		if (number >= 0)
			throw ErrSeqArray("Wrong call 'WriteVariableLength' in TVCF_Field_Format.");

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
					if (nMax < (int)F32ss[j].size())
						nMax = F32ss[j].size();
				}
				F32s.resize(nTotalSample);
				for (int i=0; i < nMax; i++)
				{
					for (int j=0; j < nTotalSample; j++)
					{
						vector<float> &B = F32ss[j];
						F32s[j] = (i < (int)B.size()) ? B[i] : (float)R_NaN;
					}
					GDS_Array_AppendData(data_obj, nTotalSample,
						&(F32s[0]), svFloat32);
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


/// trim blank characters
static void _Trim_(string &val)
{
	const char *st = val.c_str();
	const char *p = st;
	while ((*p == ' ') || (*p == '\t')) p ++;
	if (p != st) val.erase(0, p - st);

	int L = val.size();
	p = val.c_str() + L - 1;
	while (L > 0)
	{
		if ((*p != ' ') && (*p != '\t')) break;
		L --; p --;
	}
	val.resize(L);
}

/// get an integer from a string
static C_Int32 getInt32(const string &txt, bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	char *endptr = (char*)p;
	long int val = strtol(p, &endptr, 10);

	if (endptr == p)
	{
		if ((*p != '.') && RaiseError)
		{
			throw ErrSeqArray("Invalid integer conversion \"%s\".",
				SHORT_TEXT(p).c_str());
		}
		val = NA_INTEGER;
	} else {
		if ((val < INT_MIN) || (val > INT_MAX))
		{
			val = NA_INTEGER;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
		p = SKIP(endptr);
		if (*p != 0)
		{
			val = NA_INTEGER;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
	}
	return val;
}

/// get multiple integers from a string
static void getInt32Array(const string &txt, vector<C_Int32> &I32,
	bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	I32.clear();
	while (*p)
	{
		char *endptr = (char*)p;
		long int val = strtol(p, &endptr, 10);
		
		if (endptr == p)
		{
			if ((*p != '.') && RaiseError)
			{
				throw ErrSeqArray("Invalid integer conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
			val = NA_INTEGER;
		} else {
			if ((val < INT_MIN) || (val > INT_MAX))
			{
				val = NA_INTEGER;
				if (RaiseError)
				{
					throw ErrSeqArray("Invalid integer conversion \"%s\".",
						SHORT_TEXT(p).c_str());
				}
			}
			p = endptr;
		}

		I32.push_back(val);
		while ((*p != 0) && (*p != ',')) p ++;
		if (*p == ',') p ++;
	}
}


/// get a float number from a string
static float getFloat(string &txt, bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	char *endptr = (char*)p;
	float val = strtof(p, &endptr);
	if (endptr == p)
	{
		if ((*p != '.') && RaiseError)
		{
			throw ErrSeqArray("Invalid float conversion \"%s\".",
				SHORT_TEXT(p).c_str());
		}
		val = R_NaN;
	} else {
		p = SKIP(endptr);
		if (*p != 0)
		{
			val = R_NaN;
			if (RaiseError)
			{
				throw ErrSeqArray("Invalid float conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
		}
	}
	return val;
}

/// get an integer from  a string
static void getFloatArray(const string &txt, vector<float> &F32,
	bool RaiseError)
{
	const char *p = SKIP(txt.c_str());
	F32.clear();
	while (*p)
	{
		char *endptr = (char*)p;
		float val = strtof(p, &endptr);
		if (endptr == p)
		{
			if ((*p != '.') && RaiseError)
			{
				throw ErrSeqArray("Invalid float conversion \"%s\".",
					SHORT_TEXT(p).c_str());
			}
			val = R_NaN;
		} else
			p = endptr;

		F32.push_back(val);
		while ((*p != 0) && (*p != ',')) p ++;
		if (*p == ',') p ++;
	}
}

/// get an integer from  a string
static void getStringArray(string &txt, vector<string> &UTF8s)
{
	string val;
	const char *p = txt.c_str();
	while ((*p == ' ') || (*p == '\t')) p ++;

	UTF8s.clear();
	while (*p)
	{
		val.clear();
		while ((*p != 0) && (*p != ','))
			{ val.push_back(*p); p ++; }
		_Trim_(val);
		UTF8s.push_back(val);
		if (*p == ',') p ++;
	}
}

/// get name and value
static const char *_GetNameValue(const char *p, string &name, string &val)
{
	// name = val
	name.clear();
	while ((*p != 0) && (*p != ';') && (*p != '='))
	{
		name.push_back(*p);
		p ++;
	}

	val.clear();
	if (*p == '=') p ++;
	while ((*p != 0) && (*p != ';'))
	{
		val.push_back(*p);
		p ++;
	}

	if (*p == ';') p ++;
	return p;
}


extern "C"
{
// ===========================================================
// Convert from VCF4: VCF4 -> GDS
// ===========================================================

/// return true, if matching
inline static bool StrCaseCmp(const char *prefix, const char *txt)
{
	while (*prefix && *txt)
	{
		if (toupper(*prefix) != toupper(*txt))
			return false;
		prefix ++; txt ++;
	}
	return (*prefix == 0);
}

/// VCF4 --> GDS
COREARRAY_DLL_EXPORT SEXP SEQ_Parse_VCF4(SEXP vcf_fn, SEXP header,
	SEXP gds_root, SEXP param, SEXP ReadLineFun, SEXP ReadLine_Param,
	SEXP ReadLine_N, SEXP ChrPrefix, SEXP rho)
{
	const char *fn = CHAR(STRING_ELT(vcf_fn, 0));

	// define a variable for reading lines
	CReadLine RL;

	// cell buffer
	string cell;
	cell.reserve(4096);

	COREARRAY_TRY

		// the number of calling PROTECT
		int nProtected = 0;


		// =========================================================
		// initialize variables		

		// the total number of samples
		int nTotalSamp = Rf_asInteger(GetListElement(param, "sample.num"));
		// the variable name for genotypic data
		string geno_id = CHAR(STRING_ELT(GetListElement(param, "genotype.var.name"), 0));
		// raise an error
		bool RaiseError = (Rf_asLogical(GetListElement(param, "raise.error")) == TRUE);
		// verbose
		// bool Verbose = (LOGICAL(GetListElement(param, "verbose"))[0] == TRUE);

		// the number of ploidy
		int num_ploidy = Rf_asInteger(GetListElement(header, "num.ploidy"));
		if (num_ploidy <= 0)
			throw ErrSeqArray("Invalid header$num.ploidy: %d.", num_ploidy);


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
		if ((GenoNumBits!=2) && (GenoNumBits!=8))
			throw ErrSeqArray("Invalid data type in genotype/data.");
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

		// GetListElement: info
		vector<TVCF_Field_Info> info_list;
		set<string> info_missing;
		{
			SEXP info = GetListElement(header, "info");
			SEXP info_ID = GetListElement(info, "ID");
			SEXP info_inttype = GetListElement(info, "int_type");
			SEXP info_intnum = GetListElement(info, "int_num");
			SEXP info_flag = GetListElement(info, "import.flag");
			TVCF_Field_Info val;

			for (size_t i=0; i < GetLength(info_ID); i++)
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

		// GetListElement: format
		vector<TVCF_Field_Format> format_list;
		set<string> format_missing;
		{
			SEXP fmt = GetListElement(header, "format");
			SEXP fmt_ID = GetListElement(fmt, "ID");
			SEXP fmt_inttype = GetListElement(fmt, "int_type");
			SEXP fmt_intnum = GetListElement(fmt, "int_num");
			SEXP fmt_flag = GetListElement(fmt, "import.flag");
			TVCF_Field_Format val;

			for (size_t i=0; i < GetLength(fmt_ID); i++)
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
						format_list.back().F32ss.resize(nTotalSamp);
						break;
					case FIELD_TYPE_STRING:
						format_list.back().UTF8ss.resize(nTotalSamp);
						break;
					default:
						throw ErrSeqArray("Invalid FORMAT Type.");
				}
			}
		}

		// filter level list
		vector<string> filter_list;

		// variant id (integer)
		C_Int32 variant_index = GDS_Array_GetTotalCount(varIdx);

		// the string buffer
		string name, value;
		name.reserve(256);
		value.reserve(4096);

		// the numeric buffer
		C_Int32 I32;
		vector<C_Int32> I32s;
		I32s.reserve(nTotalSamp);

		C_Float32 F32;
		vector<C_Float32> F32s;
		F32s.reserve(nTotalSamp);

		// the string buffer
		vector<string> StrList;
		StrList.reserve(nTotalSamp);

		// genotypes
		vector< vector<C_Int16> > Geno;
		Geno.resize(nTotalSamp);
		vector<C_Int8> I8s;
		I8s.reserve(nTotalSamp * num_ploidy);

		const char *pCh;
		vector<TVCF_Field_Info>::iterator pInfo;
		vector< TVCF_Field_Format* >::iterator pFormat;
		vector< TVCF_Field_Format* > fmt_ptr;
		fmt_ptr.reserve(format_list.size());

		// chr prefix
		vector<string> ChrPref;
		for (R_xlen_t i=0; i < XLENGTH(ChrPrefix); i++)
			ChrPref.push_back(CHAR(STRING_ELT(ChrPrefix, i)));

		// =========================================================
		// initialize calling
		SEXP R_Read_Call;
		PROTECT(R_Read_Call =
			LCONS(ReadLineFun, LCONS(ReadLine_Param,
			LCONS(ReadLine_N, R_NilValue))));
		nProtected ++;
		RL.Init(R_Read_Call, rho);


		// =========================================================
		// skip the header

		while (!RL.IfEnd())
		{
			const char *p = RL.ReadLine();
			if (strncmp(p, "#CHROM", 6) == 0)
				break;
		}


		// the number of alleles in total at a specified site
		int num_allele;

		// =========================================================
		// parse the context

		while (!RL.IfEnd())
		{
			// =====================================================
			// scan line by line

			// -----------------------------------------------------
			// variant id
			variant_index ++;
			GDS_Array_AppendData(varIdx, 1, &variant_index, svInt32);


			// -----------------------------------------------------
			// column 1: CHROM
			RL.GetCell(cell, false);
			{
				const char *s = cell.c_str();
				vector<string>::iterator it = ChrPref.begin();
				for (; it != ChrPref.end(); it++)
				{
					if (StrCaseCmp(it->c_str(), s))
					{
						cell.erase(0, it->size());
						break;
					}
				}
			}
			GDS_Array_AppendString(varChr, cell.c_str());


			// -----------------------------------------------------
			// column 2: POS
			RL.GetCell(cell, false);
			I32 = getInt32(cell, RaiseError);
			GDS_Array_AppendData(varPos, 1, &I32, svInt32);


			// -----------------------------------------------------
			// column 3: ID
			RL.GetCell(cell, false);
			if (cell == ".") cell.clear();
			GDS_Array_AppendString(varRSID, cell.c_str());


			// -----------------------------------------------------
			// column 4 & 5: REF + ALT 
			RL.GetCell(value, false);
			RL.GetCell(cell, false);
			if (!cell.empty() && (cell != "."))
			{
				value.push_back(',');
				value.append(cell);
			}
			GDS_Array_AppendString(varAllele, value.c_str());
			// determine how many alleles
			num_allele = 0;
			pCh = value.c_str();
			while (*pCh != 0)
			{
				num_allele ++;
				while ((*pCh != 0) && (*pCh != ',')) pCh ++;
				if (*pCh == ',') pCh ++;
			}


			// -----------------------------------------------------
			// column 6: QUAL
			RL.GetCell(cell, false);
			F32 = getFloat(cell, RaiseError);
			GDS_Array_AppendData(varQual, 1, &F32, svFloat32);


			// -----------------------------------------------------
			// column 7: FILTER
			RL.GetCell(cell, false);
			if (!cell.empty() && (cell != "."))
			{
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
			for (pInfo = info_list.begin(); pInfo != info_list.end(); pInfo++)
				pInfo->used = false;
			RL.GetCell(cell, false);

			// parse
			pCh = cell.c_str();
			while (*pCh)
			{
				pCh = _GetNameValue(pCh, name, value);
				for (pInfo=info_list.begin(); pInfo != info_list.end(); pInfo++)
				{
					if (pInfo->name == name)
						break;
				}

				if (pInfo != info_list.end())
				{
					if (pInfo->used)
						throw ErrSeqArray("Duplicate INFO ID: %s.", name.c_str());

					if (pInfo->import_flag)
					{
						switch (pInfo->type)
						{
						case FIELD_TYPE_INT:
							getInt32Array(value, I32s, RaiseError);
							pInfo->Check(I32s, name, num_allele);
							GDS_Array_AppendData(pInfo->data_obj, I32s.size(),
								&(I32s[0]), svInt32);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(value, F32s, RaiseError);
							pInfo->Check(F32s, name, num_allele);
							GDS_Array_AppendData(pInfo->data_obj, F32s.size(),
								&(F32s[0]), svFloat32);
							break;

						case FIELD_TYPE_FLAG:
							if (!value.empty())
							{
								throw ErrSeqArray("INFO ID '%s' should be a flag without values.",
									name.c_str());
							}
							I32 = 1;
							GDS_Array_AppendData(pInfo->data_obj, 1,
								&I32, svInt32);
							break;

						case FIELD_TYPE_STRING:
							getStringArray(value, StrList);
							pInfo->Check(StrList, name, num_allele);
							for (int k=0; k < (int)StrList.size(); k++)
							{
								GDS_Array_AppendString(pInfo->data_obj,
									StrList[k].c_str());
							}
							break;

						default:
							throw ErrSeqArray("Invalid INFO Type.");
						}
					}

					pInfo->used = true;
				} else {
					set<string>::iterator it = info_missing.find(name);
					if (it == info_missing.end())
					{
						info_missing.insert(name);
						warning("Unknown INFO ID '%s' is ignored.", name.c_str());
					}
				}
			}

			// for which does not exist
			for (pInfo = info_list.begin(); pInfo != info_list.end(); pInfo++)
			{
				if (!pInfo->used && pInfo->import_flag)
				{
					switch (pInfo->type)
					{
					case FIELD_TYPE_INT:
						pInfo->Fill(I32s, NA_INTEGER);
						break;
					case FIELD_TYPE_FLOAT:
						pInfo->Fill(F32s, (float)R_NaN);
						break;
					case FIELD_TYPE_FLAG:
						I32 = 0;
						GDS_Array_AppendData(pInfo->data_obj, 1, &I32, svInt32);
						break;
					case FIELD_TYPE_STRING:
						pInfo->Fill(StrList, string());
						break;
					default:
						throw ErrSeqArray("Invalid INFO Type.");
					}
					pInfo->used = true;
				}
			}


			// -----------------------------------------------------
			// column 9: FORMAT

			// initialize
			for (pFormat = fmt_ptr.begin(); pFormat != fmt_ptr.end(); pFormat++)
				(*pFormat)->used = false;
			RL.GetCell(cell, false);

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
				RL.GetCell(cell, samp_idx >= (nTotalSamp-1));

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
							pFmt->Check(pFmt->I32ss[samp_idx], name, num_allele, NA_INTEGER);
							break;

						case FIELD_TYPE_FLOAT:
							getFloatArray(value, pFmt->F32ss[samp_idx], RaiseError);
							pFmt->Check(pFmt->F32ss[samp_idx], name, num_allele, (float)R_NaN);
							break;

						case FIELD_TYPE_STRING:
							getStringArray(value, pFmt->UTF8ss[samp_idx]);
							pFmt->Check(pFmt->UTF8ss[samp_idx], name, num_allele, BlackString);
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

			// determine how many bits (GenoNumBits == 2, 4 or 8)
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
							I32 = it->WriteVariableLength(nTotalSamp, I32s, F32s);
						} else
							throw ErrSeqArray("Invalid FORMAT Number.");
						GDS_Array_AppendData(it->len_obj, 1, &I32, svInt32);
					} else {
						I32 = 0;
						GDS_Array_AppendData(it->len_obj, 1, &I32, svInt32);
					}
				}
			}
		}

		// set returned value: levels(filter)
		PROTECT(rv_ans = NEW_CHARACTER(filter_list.size()));
		for (int i=0; i < (int)filter_list.size(); i++)
			SET_STRING_ELT(rv_ans, i, mkChar(filter_list[i].c_str()));
		nProtected ++;

		UNPROTECT(nProtected);

	CORE_CATCH({
		char buf[4096];
		if (RL.ColumnNo() > 0)
		{
			snprintf(buf, sizeof(buf),
				"\nFILE: %s\n\tLINE: %d, COLUMN: %d, %s\n\t%s",
				fn, RL.LineNo(), RL.ColumnNo(), cell.c_str(),
				GDS_GetError());
		} else {
			snprintf(buf, sizeof(buf), "\nFILE: %s\n\tLINE: %d\n\t%s",
				fn, RL.LineNo(), GDS_GetError());
		}
		GDS_SetError(buf);
		has_error = true;
	});
	if (has_error) error(GDS_GetError());

	// output
	return(rv_ans);
}

} // extern "C"
