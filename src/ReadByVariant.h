// ===========================================================
//
// ReadByVariant.h: Read data variant by variant
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


namespace SeqArray
{

using namespace Vectorization;


// =====================================================================

/// Object for reading basic variabls variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Basic: public CVarApply
{
protected:
	C_SVType SVType;
	SEXP VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Basic(CFileInfo &File, const char *varname);
	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


// =====================================================================

/// Object for reading genotypes variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Geno: public CVarApply
{
protected:
	CIndex<C_UInt8> *GenoIndex;  ///< indexing genotypes
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	ssize_t _SampNum;   ///< the number of selected samples
	int _Ploidy;  ///< ploidy
	bool UseRaw;  ///< whether use RAW type
	vector<C_BOOL> Selection;  ///< the buffer of selection
	AUTO_PTR ExtPtr;           ///< a pointer to the additional buffer
	SEXP VarGeno;    ///< genotype R object

	inline int _ReadGenoData(int *Base);
	inline C_UInt8 _ReadGenoData(C_UInt8 *Base);

public:
	/// constructor
	CApply_Variant_Geno();
	CApply_Variant_Geno(CFileInfo &File, bool use_raw);

	void Init(CFileInfo &File, bool use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);

	/// read genotypes in 32-bit integer
	void ReadGenoData(int *Base);
	/// read genotypes in unsigned 8-bit intetger
	void ReadGenoData(C_UInt8 *Base);

	inline int SampNum() const { return _SampNum; }
	inline int Ploidy() const { return _Ploidy; }
};


// =====================================================================

/// Object for reading genotypes (dosages) variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Dosage: public CApply_Variant_Geno
{
public:
	/// constructor
	CApply_Variant_Dosage(CFileInfo &File, bool use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);

	/// read dosages in 32-bit integer
	void ReadDosage(int *Base);
	/// read dosages in unsigned 8-bit intetger
	void ReadDosage(C_UInt8 *Base);
};


// =====================================================================

/// Object for reading phasing information variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Phase: public CVarApply
{
protected:
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	ssize_t _SampNum;   ///< the number of selected samples
	int _Ploidy;  ///< ploidy
	bool UseRaw;  ///< whether use RAW type
	vector<C_BOOL> Selection;  ///< the buffer of selection
	SEXP VarPhase;  ///< genotype R object

public:
	/// constructor
	CApply_Variant_Phase();
	CApply_Variant_Phase(CFileInfo &File, bool use_raw);

	void Init(CFileInfo &File, bool use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);

	inline int SampNum() const { return _SampNum; }
	inline int Ploidy() const { return _Ploidy; }
};


// =====================================================================
/*
/// Object for reading format variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Format: public CVarApply
{
protected:
	CIndex<C_UInt8> *VarIndex;  ///< indexing the format variable
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	ssize_t _SampNum;   ///< the number of selected samples
	int _Ploidy;  ///< ploidy

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[3];      ///< pointers to selection
	bool UseRaw;  ///< whether use RAW type
	map<size_t, SEXP> VarList; ///< a list of SEXP variables

public:
	/// constructor
	CApply_Variant_Format();
	CApply_Variant_Format(CFileInfo &File, const char *var_name, bool use_raw);

	void Init(CFileInfo &File, const char *var_name, bool use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);

	inline int SampNum() const { return _SampNum; }
};
*/


// =====================================================================

/// Object for reading a variable variant by variant
class COREARRAY_DLL_LOCAL CApplyByVariant: public CVarApply
{
protected:
	PdAbstractArray IndexNode;  ///< the corresponding index variable

	C_Int32 IndexRaw;          ///< the index according to the raw data set
	C_Int32 NumIndexRaw;       ///< the increment of raw index
	int NumOfBits;             ///< the number of bits
	size_t CellCount;          ///< the number of entries for the current variant
	map<size_t, SEXP> VarList; ///< a list of SEXP variables

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[3];      ///< pointers to selection
	bool UseRaw;            ///< whether use RAW type

	vector<C_BOOL> Selection;  ///< the buffer of selection
	AUTO_PTR ExtPtr;           ///< a pointer to the additional buffer

private:
	inline int _ReadGenoData(int *Base);
	inline C_UInt8 _ReadGenoData(C_UInt8 *Base);

public:
	int TotalNum_Variant;   ///< the total number of variants
	int _SampNum;         ///< the number of selected samples

	int DimCnt;             ///< the number of dimensions
	C_Int32 DLen[4];        ///< the dimension size

	CApplyByVariant();

	void InitObject(TVarType Type, const char *Path, CFileInfo &FileInfo,
		bool _UseRaw);

	virtual void Reset();
	virtual bool Next();

	/// read data to R object
	virtual void ReadData(SEXP val);
	/// return an R object for the next call 'ReadData()'
	virtual SEXP NeedRData(int &nProtected);
};

}


extern "C"
{

/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho);

} // extern "C"
