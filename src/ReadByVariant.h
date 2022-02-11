// ===========================================================
//
// ReadByVariant.h: Read data variant by variant
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

#include "Index.h"


namespace SeqArray
{

using namespace Vectorization;


// =====================================================================

/// Object for reading basic variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Basic: public CApply_Variant
{
protected:
	C_SVType SVType;
	SEXP VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Basic(CFileInfo &File, const char *var_name);
	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


/// Object for reading positions variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Pos: public CApply_Variant
{
protected:
	int *PtrPos;
	SEXP VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Pos(CFileInfo &File);
	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


/// Object for reading chromosomes variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Chrom: public CApply_Variant
{
protected:
	CChromIndex *ChromIndex;
	SEXP VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_Chrom(CFileInfo &File);
	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


// =====================================================================

/// Object for reading genotypes variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Geno: public CApply_Variant
{
protected:
	CGenoIndex *GenoIndex;  ///< indexing genotypes
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	int UseRaw;  ///< whether use RAW type: FALSE, int; TRUE, raw; NA: auto
	TSelection::TSampStruct *pSampSel;   ///< the structure for selected samples
	VEC_AUTO_PTR ExtPtr;  ///< a pointer to the additional buffer
	SEXP VarIntGeno;    ///< genotype R integer object
	SEXP VarRawGeno;    ///< genotype R RAW object

	inline int _ReadGenoData(int *Base);
	inline C_UInt8 _ReadGenoData(C_UInt8 *Base);

public:
	ssize_t SampNum;  ///< the number of selected samples
	int Ploidy;       ///< ploidy

	/// constructor
	CApply_Variant_Geno();
	CApply_Variant_Geno(CFileInfo &File, int use_raw);

	void Init(CFileInfo &File, int use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
	bool NeedIntType();

	/// read genotypes in 32-bit integer
	void ReadGenoData(int *Base);
	/// read genotypes in unsigned 8-bit intetger
	void ReadGenoData(C_UInt8 *Base);
};


// =====================================================================

/// Object for reading genotypes (dosages) variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Dosage: public CApply_Variant_Geno
{
protected:
	SEXP VarDosage;        ///< dosage R object
	VEC_AUTO_PTR ExtPtr2;  ///< a pointer to the additional buffer for dosages
	bool IsAlt;            ///< if true, ReadData() returns the dosage of alternative alleles
public:
	/// constructor
	CApply_Variant_Dosage(CFileInfo &File, int use_raw, bool alt);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);

	/// read dosages in 32-bit integer
	void ReadDosage(int *Base);
	/// read dosages of alternative alleles in 32-bit integer
	void ReadDosageAlt(int *Base);
	/// read dosages in unsigned 8-bit intetger
	void ReadDosage(C_UInt8 *Base);
	/// read dosages of alternative alleles in unsigned 8-bit intetger
	void ReadDosageAlt(C_UInt8 *Base);
};


// =====================================================================

/// Object for reading phasing information variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Phase: public CApply_Variant
{
protected:
	ssize_t SiteCount;  ///< the total number of entries at a site
	ssize_t CellCount;  ///< the selected number of entries at a site
	bool UseRaw;  ///< whether use RAW type
	vector<C_BOOL> Selection;  ///< the buffer of selection
	SEXP VarPhase;  ///< genotype R object

public:
	ssize_t SampNum;  ///< the number of selected samples
	int Ploidy;       ///< ploidy

	/// constructor
	CApply_Variant_Phase();
	CApply_Variant_Phase(CFileInfo &File, bool use_raw);

	void Init(CFileInfo &File, bool use_raw);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


// =====================================================================

/// Object for reading info variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Info: public CApply_Variant
{
protected:
	CIndex *VarIndex;  ///< indexing the format variable
	C_SVType SVType;        ///< data type for GDS reading
	C_Int32 BaseNum;        ///< if 2-dim, the size of the first dimension
	map<int, SEXP> VarList;  ///< a list of SEXP variables

public:
	/// constructor
	CApply_Variant_Info(CFileInfo &File, const char *var_name);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


// =====================================================================

/// Object for reading format variables variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_Format: public CApply_Variant
{
protected:
	CIndex *VarIndex;  ///< indexing the format variable
	ssize_t _TotalSampNum;  ///< the total number of samples

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[2];      ///< pointers to selection
	map<int, SEXP> VarList;  ///< a list of SEXP variables

public:
	ssize_t SampNum;  ///< the number of selected samples

	/// constructor
	CApply_Variant_Format();
	CApply_Variant_Format(CFileInfo &File, const char *var_name);

	void Init(CFileInfo &File, const char *var_name);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
};


// =====================================================================

/// Object for calculating the number of distinct alleles variant by variant
class COREARRAY_DLL_LOCAL CApply_Variant_NumAllele: public CApply_Variant
{
private:
	string strbuf;
protected:
	SEXP VarNode;  ///< R object
public:
	/// constructor
	CApply_Variant_NumAllele(CFileInfo &File);

	virtual void ReadData(SEXP val);
	virtual SEXP NeedRData(int &nProtected);
	int GetNumAllele();
};

}


extern "C"
{

/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho);

} // extern "C"
