// ===========================================================
//
// Index.h: Indexing Objects
//
// Copyright (C) 2016-2017    Xiuwen Zheng
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


#ifndef _HEADER_SEQ_INDEX_
#define _HEADER_SEQ_INDEX_

#include <R_GDS_CPP.h>
#include <dTrait.h>

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <ctime>

#include <cctype>
#include <cstring>
#include "vectorization.h"


extern "C"
{
#define class xclass
#define private xprivate
#include <R_ext/Connections.h>
#undef class
#undef private


// NO_R_v3_3 can be defined in Makevars
#ifdef NO_R_v3_3
COREARRAY_DLL_LOCAL Rconnection My_R_GetConnection(SEXP x);
#define R_GetConnection My_R_GetConnection
#endif


/// define missing value of RAW
#define NA_RAW     0xFF
}



namespace SeqArray
{

using namespace std;
using namespace CoreArray;


class ErrSeqArray;


// ===========================================================
// Run-length encoding (RLE) object
// ===========================================================

/// object with run-length encoding
template<typename TYPE> class COREARRAY_DLL_LOCAL C_RLE
{
public:
	/// constructor
	C_RLE()
	{
		TotalLength = 0;
		Position = AccIndex = AccOffset = 0;
	}

	void Init()
	{
		TotalLength = 0;
		vector<C_UInt32>::iterator p;
		for (p=Lengths.begin(); p != Lengths.end(); p++)
			TotalLength += *p;
		Position = AccIndex = AccOffset = 0;
	}

	void Add(TYPE &val, C_UInt32 len)
	{
		Values.push_back(val);
		Lengths.push_back(len);
	}

	void Clear()
	{
		Values.clear(); Lengths.clear();
		TotalLength = 0;
		Position = AccIndex = AccOffset = 0;
	}

	const TYPE &operator [](size_t pos)
	{
		if (pos >= TotalLength)
			throw "Invalid position in C_RLE.";
		if (pos < Position)
			Position = AccIndex = AccOffset = 0;
		for (; Position < pos; )
		{
			size_t L = Lengths[AccIndex];
			size_t n = L - AccOffset;
			if ((Position + n) <= pos)
			{
				AccIndex ++; AccOffset = 0;
			} else {
				n = pos - Position; AccOffset += n;
			}
			Position += n;
		}
		return Values[AccIndex];
	}

	inline bool Empty() const { return (TotalLength <= 0); }

protected:
	/// values according to Lengths, used in run-length encoding
	vector<TYPE> Values;
	/// lengths according to Values, used in run-length encoding
	vector<C_UInt32> Lengths;
	/// total number, = sum(Lengths)
	size_t TotalLength;
	/// the position relative to the total length
	size_t Position;
	/// the index in Lengths according to Position
	size_t AccIndex;
	/// the offset according the value of Lengths[AccIndex]
	size_t AccOffset;
};


// ===========================================================
// Indexing object
// ===========================================================

/// Indexing object with run-length encoding
class COREARRAY_DLL_LOCAL CIndex
{
public:
	/// values according to Lengths, used in run-length encoding
	vector<int> Values;
	/// lengths according to Values, used in run-length encoding
	vector<C_UInt32> Lengths;

	/// constructor
	CIndex();

	/// load data and represent as run-length encoding
	void Init(PdContainer Obj);
	/// load data and represent as run-length encoding
	void InitOne(int num);
	/// return the accumulated sum of values and current value in Lengths and Values given by a position
	void GetInfo(size_t pos, C_Int64 &Sum, int &Value);
	/// get lengths with selection
	SEXP GetLen_Sel(const C_BOOL sel[]);
	/// return true if empty
	inline bool Empty() const { return (TotalLength <= 0); }

protected:
	/// total number, = sum(Lengths)
	size_t TotalLength;
	/// the position relative to the total length
	size_t Position;
	/// the accumulated sum of values in Lengths and Values according to Position
	C_Int64 AccSum;
	/// the index in Lengths according to Position
	size_t AccIndex;
	/// the offset according the value of Lengths[AccIndex]
	size_t AccOffset;
};


/// Indexing object with run-length encoding for genotype indexing
class COREARRAY_DLL_LOCAL CGenoIndex
{
public:
	/// values according to Lengths, used in run-length encoding
	vector<C_UInt16> Values;
	/// lengths according to Values, used in run-length encoding
	vector<C_UInt32> Lengths;

	/// constructor
	CGenoIndex();

	/// load data and represent as run-length encoding
	void Init(PdContainer Obj);
	/// return the accumulated sum of values and current value in Lengths and Values given by a position
	void GetInfo(size_t pos, C_Int64 &Sum, C_UInt8 &Value);
	/// return true if empty
	inline bool Empty() const { return (TotalLength <= 0); }

protected:
	/// total number, = sum(Lengths)
	size_t TotalLength;
	/// the position relative to the total length
	size_t Position;
	/// the accumulated sum of values in Lengths and Values according to Position
	C_Int64 AccSum;
	/// the index in Lengths according to Position
	size_t AccIndex;
	/// the offset according the value of Lengths[AccIndex]
	size_t AccOffset;
};



// ===========================================================
// Chromosome indexing
// ===========================================================

/// Chromosome indexing object
class COREARRAY_DLL_LOCAL CChromIndex
{
public:
	/// range object
	struct TRange
	{
		int Start;   ///< the starting position
		int Length;  ///< the length
	};

	typedef vector<TRange> TRangeList;

	/// constructor
	CChromIndex();

	/// clear
	void Clear();

	/// represent chromosome codes as a RLE object in Map
	void AddChrom(PdGDSFolder Root);

	/// the total length of a TRangeList object
	size_t RangeTotalLength(const TRangeList &RngList);

	/// whether it is empty
	inline bool Empty() const { return Map.empty(); }

	inline const string &operator [](size_t pos) { return PosToChr[pos]; }

	/// map to TRangeList from chromosome coding
	map<string, TRangeList> Map;

protected:
	/// position to chromosome
	C_RLE<string> PosToChr;
};



// ===========================================================
// Genomic Range Sets
// ===========================================================

/// Genomic Range Set Object
class COREARRAY_DLL_LOCAL CRangeSet
{
public:
	/// range object
	struct TRange
	{
		int Start;     ///< the starting position
		int End;       ///< the ending (always, End >= Start)
	};

	void Clear();
	void AddRange(int start, int end);
	bool IsIncluded(int point);

protected:
	/// strict weak ordering for non-overlapping, == when overlapping
	struct less_range
	{
    	bool operator()(const TRange &lhs, const TRange &rhs) const;
	};

	/// 
	set<TRange, less_range> _RangeSet;
};





// ===========================================================
// SeqArray GDS file information
// ===========================================================

/// selection object used in GDS file
struct COREARRAY_DLL_LOCAL TSelection
{
	vector<C_BOOL> Sample;   ///< sample selection
	vector<C_BOOL> Variant;  ///< variant selection

	inline C_BOOL *pSample()
		{ return Sample.empty() ? NULL : &Sample[0]; }
	inline C_BOOL *pVariant()
		{ return Variant.empty() ? NULL : &Variant[0]; }
};


/// GDS file object
class COREARRAY_DLL_LOCAL CFileInfo
{
public:
	list<TSelection> SelList;  ///< a list of sample and variant selections

	/// constructor
	CFileInfo(PdGDSFolder root=NULL);
	/// destructor
	~CFileInfo();

	/// reset the root of GDS file
	void ResetRoot(PdGDSFolder root);
	/// get selection
	TSelection &Selection();

	/// return _Chrom which has been initialized
	CChromIndex &Chromosome();
	/// return _Position which has been initialized
	vector<C_Int32> &Position();

	/// return _GenoIndex which has been initialized
	CGenoIndex &GenoIndex();

	/// return the indexing object according to variable name
	CIndex &VarIndex(const string &varname);

	/// get gds object
	PdAbstractArray GetObj(const char *name, C_BOOL MustExist);

	/// the root of gds file
	inline PdGDSFolder Root() { return _Root; }
	/// the total number of samples
	inline int SampleNum() const { return _SampleNum; }
	/// the total number of variants
	inline int VariantNum() const { return _VariantNum; }
	/// ploidy
	inline int Ploidy() const { return _Ploidy; }

	int SampleSelNum();
	int VariantSelNum();

protected:
	PdGDSFolder _Root;  ///< the root of GDS file
	int _SampleNum;     ///< the total number of samples
	int _VariantNum;    ///< the total number of variants
	int _Ploidy;        ///< ploidy

	CChromIndex _Chrom;  ///< chromosome indexing
	vector<C_Int32> _Position;  ///< position
	CGenoIndex _GenoIndex;  ///< the indexing object for genotypes
	map<string, CIndex> _VarIndex;  ///< the indexing objects for INFO/FORMAT variables
};


extern std::map<int, CFileInfo> COREARRAY_DLL_LOCAL GDSFile_ID_Info;

/// get the associated CFileInfo
COREARRAY_DLL_LOCAL CFileInfo &GetFileInfo(SEXP gdsfile);




// ===========================================================
// GDS Variable Type
// ===========================================================

class COREARRAY_DLL_LOCAL CVariable
{
public:
	enum TVarType
	{
		ctNone,
		ctBasic,       ///< sample.id, variant.id, etc
		ctGenotype,    ///< genotypes or alleles
		ctDosage,      ///< dosage of reference or specified allele
		ctPhase,       ///< phase information
		ctInfo,        ///< variant annotation info field
		ctFormat,      ///< variant annotation format field
		ctSampleAnnot  ///< sample annotation
	};
};


/// The abstract class for applying functions marginally
class COREARRAY_DLL_LOCAL CVarApply: public CVariable
{
protected:
	TVarType fVarType;       ///< VCF data type
	ssize_t MarginalSize;    ///< the size in MarginalSelect
	C_BOOL *MarginalSelect;  ///< pointer to variant selection

public:
	PdAbstractArray Node;  ///< the GDS variable
	C_Int32 Position;  ///< the index of variant/sample, starting from ZERO

	/// constructor
	CVarApply();
	/// destructor
	virtual ~CVarApply();

	/// reset
	virtual void Reset();
	/// move to the next element
	virtual bool Next();

	/// return an R object for the next call 'ReadData()'
	virtual SEXP NeedRData(int &nProtected) = 0;
	/// read data to R object
	virtual void ReadData(SEXP val) = 0;

	/// variable type
	inline TVarType VarType() const { return fVarType; }

	/// need a pointer to size of TRUEs
	C_BOOL *NeedTRUEs(size_t size);

private:
	vector<C_BOOL> _TRUE;
};


/// The abstract class for applying functions by variant
class COREARRAY_DLL_LOCAL CApply_Variant: public CVarApply
{
public:
	/// constructor
	CApply_Variant();
	/// constructor with file information
	CApply_Variant(CFileInfo &File);
};


class COREARRAY_DLL_LOCAL CVarApplyList: public vector<CVarApply*>
{
public:
	~CVarApplyList();

	/// return false if any return false, otherwise return true
	bool CallNext();
};



// ===========================================================
// Progress object
// ===========================================================

class COREARRAY_DLL_LOCAL CProgress
{
public:
	CProgress(C_Int64 start, C_Int64 count, SEXP conn, bool newline);
	virtual ~CProgress();

	void Forward();
	virtual void ShowProgress();

protected:
	C_Int64 TotalCount;  ///< the total number
	C_Int64 Counter;  ///< the current counter
	Rconnection File;  ///< R connection
	bool NewLine;
	double _start, _step;
	C_Int64 _hit;
	vector< pair<double, time_t> > _timer;
};

class COREARRAY_DLL_LOCAL CProgressStdOut: public CProgress
{
public:
	CProgressStdOut(C_Int64 count, bool verbose);
	virtual void ShowProgress();

protected:
	time_t _last_time;
	bool Verbose;
};



// ===========================================================
// Pre-defined R objects
// ===========================================================

extern SEXP R_Geno_Dim2_Name;
extern SEXP R_Geno_Dim3_Name;
extern SEXP R_Dosage_Name;
extern SEXP R_Data_Name;
extern SEXP R_Data_Dim2_Name;
extern SEXP R_Data_Dim3_Name;



// ===========================================================
// Define Functions
// ===========================================================

/// Get the number of TRUEs
#define GetNumOfTRUE(ptr, n)    vec_i8_cnt_nonzero((C_Int8*)(ptr), n)

/// return the length in val, it is safe to call when val=R_NilValue
COREARRAY_DLL_LOCAL size_t RLength(SEXP val);

/// get the list element named str, or return R_NilValue
COREARRAY_DLL_LOCAL SEXP RGetListElement(SEXP list, const char *name);

/// allocate R object given by SVType
COREARRAY_DLL_LOCAL SEXP RObject_GDS(PdAbstractArray Node, size_t n,
	int &nProtected, bool bit1_is_logical);

/// append data to a GDS node
COREARRAY_DLL_LOCAL void RAppendGDS(PdAbstractArray Node, SEXP Val);


/// requires a vector of TRUEs
COREARRAY_DLL_LOCAL C_BOOL *NeedArrayTRUEs(size_t len);

/// Get pretty text for an integer with comma
COREARRAY_DLL_LOCAL const char *PrettyInt(int val);

/// Text matching, return -1 when no maching
COREARRAY_DLL_LOCAL int MatchText(const char *txt, const char *list[]);

/// Get the number of alleles
COREARRAY_DLL_LOCAL int GetNumOfAllele(const char *allele_list);

/// Get the index in an allele list
COREARRAY_DLL_LOCAL int GetIndexOfAllele(const char *allele, const char *allele_list);

/// Get strings split by comma
COREARRAY_DLL_LOCAL void GetAlleles(const char *alleles, vector<string> &out);


/// get PdGDSObj from a SEXP object
COREARRAY_DLL_LOCAL void GDS_PATH_PREFIX_CHECK(const char *path);

/// check variable name
COREARRAY_DLL_LOCAL void GDS_VARIABLE_NAME_CHECK(const char *p);

/// get PdGDSObj from a SEXP object
COREARRAY_DLL_LOCAL string GDS_PATH_PREFIX(const string &path, char prefix);


/// output to a connection
COREARRAY_DLL_LOCAL void ConnPutText(Rconnection file, const char *fmt, ...);



// ===========================================================
// Define Exception
// ===========================================================

class ErrSeqArray: public ErrCoreArray
{
public:
	ErrSeqArray(): ErrCoreArray()
		{ }
	ErrSeqArray(const char *fmt, ...): ErrCoreArray()
		{ _COREARRAY_ERRMACRO_(fmt); }
	ErrSeqArray(const std::string &msg): ErrCoreArray()
		{ fMessage = msg; }
};

}

#endif /* _HEADER_SEQ_INDEX_ */
