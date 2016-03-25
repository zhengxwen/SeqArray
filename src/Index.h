// ===========================================================
//
// Index.h: Indexing Objects
//
// Copyright (C) 2016    Xiuwen Zheng
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

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>

#include <cctype>
#include <cstring>
#include "vectorization.h"


namespace SeqArray
{

using namespace std;
using namespace CoreArray;

// ===========================================================
// Indexing object
// ===========================================================

/// Indexing object
template<typename TYPE> class COREARRAY_DLL_LOCAL CIndex
{
public:
	/// Run-length encoding
	struct TRLE
	{
		int Count;     ///< the number of the value
		TYPE Value;    ///< value
	};

	/// represent chromosome codes as a RLE object in Map
	// void Init(PdGDSObj GDSObj);

	///
	void Seek(size_t pos) {}


protected:

	/// 
	vector<TRLE> _List;
	///
	size_t ListIdx;
	///
	size_t TotalSize, Position;
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

	/// map to TRangeList from chromosome coding
	map<string, TRangeList> Map;
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

	inline C_BOOL *pSample() { return &Sample[0]; }
	inline C_BOOL *pVariant() { return &Variant[0]; }
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

	/// get gds object
	PdAbstractArray GetObj(const char *name, C_BOOL MustExist);

	/// the root of gds file
	inline PdGDSFolder Root() { return _Root; }
	/// the total number of samples
	inline int SampleNum() { return _SampleNum; }
	/// the total number of variants
	inline int VariantNum() { return _VariantNum; }

	int SampleSelNum();
	int VariantSelNum();

protected:
	PdGDSFolder _Root;  ///< the root of GDS file
	int _SampleNum;   ///< the total number of samples
	int _VariantNum;  ///< the total number of variants

	CChromIndex _Chrom;  ///< chromosome indexing
	vector<C_Int32> _Position;  ///< position
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
	enum TType
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


class COREARRAY_DLL_LOCAL CVarApply: public CVariable
{
public:
	C_BOOL *NeedTRUEs(size_t size);

	CVarApply() {}
	virtual ~CVarApply() {}

private:
	vector<C_BOOL> _TRUE;
};



// ===========================================================
// Define Functions
// ===========================================================

/// return the length in val, it is safe to call when val=R_NilValue
COREARRAY_DLL_LOCAL size_t RLength(SEXP val);

/// get the list element named str, or return R_NilValue
COREARRAY_DLL_LOCAL SEXP RGetListElement(SEXP list, const char *name);

/// requires a vector of TRUEs
COREARRAY_DLL_LOCAL C_BOOL *NeedArrayTRUEs(size_t len);

/// Get pretty text for an integer with comma
COREARRAY_DLL_LOCAL const char *PrettyInt(int val);

/// Text matching, return -1 when no maching
COREARRAY_DLL_LOCAL int MatchText(const char *txt, const char *list[]);



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
