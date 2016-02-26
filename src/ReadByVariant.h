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


// ===================================================================== //

/// Object for reading a variable variant by variant
class COREARRAY_DLL_LOCAL CVarApplyByVariant: public CVarApply
{
protected:
	PdAbstractArray Node;       ///< the GDS variable
	PdAbstractArray IndexNode;  ///< the corresponding index variable

	C_Int32 IndexRaw;          ///< the index according to the raw data set
	C_Int32 NumIndexRaw;       ///< the increment of raw index
	size_t CellCount;          ///< the number of entries for the current variant
	map<size_t, SEXP> VarList; ///< a list of SEXP variables

	C_SVType SVType;        ///< data type for GDS reading
	C_BOOL *SelPtr[3];      ///< pointers to selection
	C_BOOL *VariantSelect;  ///< pointer to variant selection
	bool UseRaw;            ///< whether use RAW type

	vector<C_BOOL> Selection;  ///< the buffer of selection
	int NumOfBits;             ///< the number of bits

public:
	TType VarType;          ///< VCF data type
	int TotalNum_Variant;   ///< the total number of variants
	int Num_Sample;         ///< the number of selected samples

	int DimCnt;             ///< the number of dimensions
	C_Int32 DLen[4];        ///< the dimension size
	C_Int32 CurIndex;       ///< the index of variant, starting from ZERO

	CVarApplyByVariant();
	virtual ~CVarApplyByVariant() {}

	void InitObject(TType Type, const char *Path, PdGDSObj Root,
		int nVariant, C_BOOL *VariantSel, int nSample, C_BOOL *SampleSel,
		bool _UseRaw);

	void ResetObject();
	bool NextCell();

	/// read genotypes in 32-bit integer
	void ReadGenoData(int *Base);
	/// read genotypes in unsigned 8-bit intetger
	void ReadGenoData(C_UInt8 *Base);

	void ReadData(SEXP Val);

	SEXP NeedRData(int &nProtected);
};



extern "C"
{
/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP use_raw, SEXP list_duplicate,
	SEXP rho);
} // extern "C"
