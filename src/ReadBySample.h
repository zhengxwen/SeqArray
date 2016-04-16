// ===========================================================
//
// ReadBySample.h: Read data sample by sample
//
// Copyright (C) 2015-2016    Xiuwen Zheng
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

/// 
class COREARRAY_DLL_LOCAL CVarApplyBySample: public CVarApply
{
protected:
	C_Int32 VariantStart;   ///< start index according to the variants
	C_Int32 VariantCount;   ///< the length according to the variants
	int NumOfBits;             ///< the number of bits
	size_t CellCount;       ///< the number of entries for the current sample
	vector<C_UInt8> GenoCellCnt;  ///< 
	map<size_t, SEXP> VarList;    ///< a list of SEXP variables

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[3];      ///< pointers to selection
	bool UseRaw;            ///< whether use RAW type

	vector<C_BOOL> Selection;  ///< the buffer of selection

public:
	int TotalNum_Sample;    ///< the total number of samples
	int Num_Variant;        ///< the number of selected variants

	int DimCnt;             ///< the number of dimensions
	C_Int32 DLen[4];        ///< the dimension size


	CVarApplyBySample();
	virtual ~CVarApplyBySample() {}

	void InitObject(TVarType Type, const char *Path, PdGDSObj Root,
		int nVariant, C_BOOL *VariantSel, int nSample, C_BOOL *SampleSel,
		bool _UseRaw);

	virtual void Reset();
	virtual bool Next();

	/// read genotypes in 32-bit integer
	void ReadGenoData(int *Base);
	/// read genotypes in unsigned 8-bit intetger
	void ReadGenoData(C_UInt8 *Base);

	virtual void ReadData(SEXP Val);

	virtual SEXP NeedRData(int &nProtected);
};

}


extern "C"
{

COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Sample(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP use_raw, SEXP rho);

} // extern "C"
