// ===========================================================
//
// Index.h: Indexing Objects
//
// Copyright (C) 2016-2025    Xiuwen Zheng
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

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <ctime>
#include <cctype>
#include <cstring>
#include <dTrait.h>
#include <R_GDS_CPP.h>
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
	friend class CChromIndex;

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
	void Init(PdContainer Obj, const char *varname);
	/// load data and represent as run-length encoding
	void InitOne(int num);
	/// return the accumulated sum of values and current value in Lengths and Values given by a position
	void GetInfo(size_t pos, C_Int64 &Sum, int &Value);
	/// get lengths with selection
	SEXP GetLen_Sel(const C_BOOL sel[]);
	/// get lengths and bool selection from a set of selected variants
	SEXP GetLen_Sel(const C_BOOL sel[], int &out_var_start, int &out_var_count,
		vector<C_BOOL> &out_var_sel);
	/// return true if empty
	inline bool Empty() const { return (TotalLength <= 0); }
	/// return true if there is an index stored in GDS
	inline bool HasIndex() const { return has_index; }
	/// return true if fixed length is one
	inline bool IsFixedOne() const { return Values.size()==1 && Values[0]==1; }
	/// the maximum value in Values
	inline int ValLenMax() const { return val_max; }

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
	/// true if there is an index stored in GDS
	bool has_index;
	/// the maximum value in Values
	int val_max;
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
	void Init(PdContainer Obj, const char *varname);
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

	/// map to TRangeList from chromosome coding
	map<string, TRangeList> Map;

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
	/// return chromosome coding with the index 
	inline const string &operator [](size_t idx) { return RleChr[idx]; }
	/// return the RLE representation of chromosome: value
	inline vector<string> &RLE_Values() { return RleChr.Values; }
	/// return the RLE representation of chromosome: length
	inline vector<C_UInt32> &RLE_Lengths() { return RleChr.Lengths; }

protected:
	/// indexing chromosome
	C_RLE<string> RleChr;
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
	void GetRanges(int Start[], int End[]);
	inline size_t Size() const { return _RangeSet.size(); }

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

class CFileInfo;

/// selection object used in GDS file
struct COREARRAY_DLL_LOCAL TSelection
{
	struct COREARRAY_DLL_LOCAL TSampStruct {
		ssize_t length;
		ssize_t offset;
		C_BOOL *sel;
	};

	TSelection *Link;  ///< the pointer to the last one
	C_BOOL *pSample;   ///< sample selection
	C_BOOL *pVariant;  ///< variant selection

	ssize_t varTrueNum;  ///< the number of TRUEs in pVariant, -1 for requiring initialization
	ssize_t varStart;    ///< the start position of the first TRUE in pVariant
	ssize_t varEnd;      ///< the next position of the last TRUE in pVariant

	/// constructor
	TSelection(CFileInfo &File, bool init);
	/// destructor
	~TSelection();

	/// get the pointer to the sample reading structure
	TSampStruct *GetStructSample();
	/// clear the structure of selected samples for resetting the sample filter
	void ClearStructSample();

	/// get the structure of selected variants
	void GetStructVariant();
	/// clear selected varaints
	void ClearSelectVariant();
	/// clear the structure of selected variants for resetting the variant filter
	void ClearStructVariant();

private:
	size_t numSamp;    ///< the total number of samples
	size_t numVar;     ///< the total number of variants
	size_t numPloidy;  ///< the ploidy
	C_BOOL *pFlagGenoSel;  ///< the pointer to the genotype selection according to the selected samples
	vector<TSampStruct> pSampList;
};


/// GDS variable structure used in SeqArray::seqGetData()
struct COREARRAY_DLL_LOCAL TVarMap
{
public:
	/// function call according to a GDS node
	typedef SEXP (*TFunction)(CFileInfo &File, TVarMap &Var, void *param);

	string Name;     ///< the variable name
	PdGDSObj Obj;    ///< GDS node
	int ObjID;       ///< GDS node ID stored in gdsfmt
	int NDim;        ///< the number of dimensions
	C_Int32 Dim[4];  ///< the size of each dimension
	TFunction Func;  ///< function pointer to get an R object
	bool IsBit1;     ///< ture if it is dBit1
	CIndex Index;    ///< indexing

	TVarMap();       ///< constructor
	/// initialize assuming no indexing
	void Init(CFileInfo &file, const string &varnm, TFunction fc);
	/// initialize with possible indexing
	void InitWtIndex(CFileInfo &file, const string &varnm, TFunction fc);

private:
	void get_obj(CFileInfo &file, const string &varnm);
};


/// GDS file object
class COREARRAY_DLL_LOCAL CFileInfo
{
public:
	friend struct TVarMap;

	/// constructor
	CFileInfo(PdGDSFolder root=NULL);
	/// destructor
	~CFileInfo();

	/// reset the root of GDS file
	void ResetRoot(PdGDSFolder root);

	/// get the current selection
	TSelection &Selection();
	/// push a new selection to the list of selection
	TSelection &Push_Selection(bool init_samp, bool init_var);
	/// pop back a selection
	void Pop_Selection();

	/// return _Chrom which has been initialized
	CChromIndex &Chromosome();
	/// reload chromosome coding when it is changed
	void ResetChromosome();
	/// return _Position which has been initialized
	vector<C_Int32> &Position();
	/// clear the buffer for variant positions
	void ClearPosition();

	/// return _GenoIndex which has been initialized
	CGenoIndex &GenoIndex();

	/// return variable structure with possible indexing
	map<string, TVarMap> &VarMap() { return _VarMap; }

	/// get gds object
	PdAbstractArray GetObj(const char *name, C_BOOL MustExist);

	/// the gds file
	inline PdGDSFile File() { return _File; }
	/// the root of gds file
	inline PdGDSFolder Root() { return _Root; }
	/// the total number of samples
	inline int SampleNum() const { return _SampleNum; }
	/// the total number of variants
	inline int VariantNum() const { return _VariantNum; }
	/// ploidy
	inline int Ploidy() const { return _Ploidy; }

	/// get the number of selected samples
	int SampleSelNum();
	/// get the number of selected variants
	int VariantSelNum();

protected:
	PdGDSFile _File;       ///< the GDS file
	PdGDSFolder _Root;     ///< the root of GDS file
	TSelection *_SelList;  ///< the pointer to the sample and variant selections
	int _SampleNum;   ///< the total number of samples
	int _VariantNum;  ///< the total number of variants
	int _Ploidy;      ///< ploidy

	CChromIndex _Chrom;  ///< chromosome indexing
	vector<C_Int32> _Position;  ///< position
	CGenoIndex _GenoIndex;  ///< the indexing object for genotypes
	map<string, TVarMap> _VarMap;  ///< the indexing objects for seqGetData()

private:
	inline void clear_selection();
};

/// GDS files with specified IDs
extern std::map<int, CFileInfo> COREARRAY_DLL_LOCAL GDSFile_ID_Info;

/// get the associated CFileInfo
COREARRAY_DLL_LOCAL CFileInfo &GetFileInfo(SEXP gdsfile);

/// get TVarMap from a variable name
COREARRAY_DLL_LOCAL TVarMap &VarGetStruct(CFileInfo &File, const string &name);


/// Vector read from a GDS node
template<typename TYPE> class COREARRAY_DLL_LOCAL CVectorRead
{
public:
	/// constructor
	CVectorRead(PdContainer obj, C_BOOL *selbase, C_Int32 st, C_Int32 nTRUE)
	{
		Obj = obj; SelBase = selbase;
		Start = st; NumTRUE = nTRUE;
	}
	/// read
	template<typename OUTTYPE> int Read(OUTTYPE out[], int nmax)
	{
		if (nmax > NumTRUE) nmax = NumTRUE;
		if (nmax > 0)
		{
			C_BOOL *ss = SelBase + Start, *p = ss;
			for (int m=nmax; m > 0;) if (*p++) m--;
			C_Int32 Len = p - ss;
			GDS_Array_ReadDataEx(Obj, &Start, &Len, &ss, out, TdTraits<OUTTYPE>::SVType);
			Start += Len;
			NumTRUE -= nmax;
		}
		return nmax;
	}

private:
	PdContainer Obj;    ///< a GDS node
	C_BOOL *SelBase;    ///< the pointer to the start of selection
	C_Int32 Start;      ///< the start position of selection
	C_Int32 NumTRUE;    ///< the total number of selected variants after Start
};



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
	ssize_t MarginalStart;   ///< the start position in MarginalSelect
	ssize_t MarginalEnd;     ///< the ending position in MarginalSelect
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

protected:
	void InitMarginal(CFileInfo &File);
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
	CProgress(C_Int64 count, SEXP Rconn, bool verbose);
	~CProgress();

	void Forward(C_Int64 Inc=1);
	void ShowProgress();

	inline C_Int64 Counter() const { return vCounter; }
	inline C_Int64 TotalCount() const { return vTotalCount; }

protected:
	C_Int64 vTotalCount;  ///< the total number
	C_Int64 vCounter;     ///< the current counter
	Rconnection OutFile;  ///< R connection output
	bool Verbose;         ///< whether display on stdout via Rprintf()
	C_Int64 FwdCnt;       ///< the number of calling Forward()
	time_t _start_time;   ///< the starting time
	time_t _last_time;    ///< last saved time for calculating the interval
	double _start, _step;
	C_Int64 _hit;
	vector< pair<double, time_t> > _timer;
};



// ===========================================================
// Pre-defined R objects
// ===========================================================

extern SEXP R_Geno_Dim2_Name;
extern SEXP R_Geno_Dim3_Name;
extern SEXP R_Dosage_Name;
extern SEXP R_Data_Name;
extern SEXP R_Data_Dim2_Name;
extern SEXP R_Data_ListClass;

// the index of child process
extern int* R_Process_Index;
// the number of child processes
extern int* R_Process_Count;
// list of file names for child processes
extern vector<string> R_Process_StatusFName;
// the index of block being used currently
extern int R_Block_Index;
// the number of blocks, 0 when it is not used
extern int R_Block_Count;


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
