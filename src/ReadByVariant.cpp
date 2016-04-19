// ===========================================================
//
// ReadByVariant.cpp: Read data variant by variant
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

#include "ReadByVariant.h"


namespace SeqArray
{

using namespace Vectorization;

static const char *ERR_DIM = "Invalid dimension of '%s'.";


// =====================================================================
// Object for reading basic variabls variant by variant

CApply_Variant_Basic::CApply_Variant_Basic(CFileInfo &File,
	const char *varname): CVarApply()
{
	fVarType = ctBasic;
	MarginalSize = File.VariantNum();
	MarginalSelect = File.Selection().pVariant();
	Node = File.GetObj(varname, TRUE);
	SVType = GDS_Array_GetSVType(Node);
	VarNode = NULL;
	Reset();
}

void CApply_Variant_Basic::ReadData(SEXP val)
{
	C_Int32 st = Position;
	C_Int32 one = 1;
	if (COREARRAY_SV_INTEGER(SVType))
	{
		GDS_Array_ReadData(Node, &st, &one, INTEGER(val), svInt32);
	} else if (COREARRAY_SV_FLOAT(SVType))
	{
		GDS_Array_ReadData(Node, &st, &one, REAL(val), svFloat64);
	} else if (COREARRAY_SV_STRING(SVType))
	{
		string s;
		GDS_Array_ReadData(Node, &st, &one, &s, svStrUTF8);
		SET_STRING_ELT(val, 0, mkChar(s.c_str()));
	}
}

SEXP CApply_Variant_Basic::NeedRData(int &nProtected)
{
	if (VarNode == NULL)
		VarNode = RObject_GDS(Node, 1, nProtected, false);
	return VarNode;
}



// =====================================================================
// Object for reading genotypes variant by variant

CApply_Variant_Geno::CApply_Variant_Geno():
	CVarApply()
{
	fVarType = ctGenotype;
	SiteCount = CellCount = 0;
	_SampNum = 0; _Ploidy = 0;
	UseRaw = FALSE;
	VarGeno = NULL;
}

CApply_Variant_Geno::CApply_Variant_Geno(CFileInfo &File, bool use_raw):
	CVarApply()
{
	fVarType = ctGenotype;
	Init(File, use_raw);
}

void CApply_Variant_Geno::Init(CFileInfo &File, bool use_raw)
{
	static const char *VAR_NAME = "genotype/data";

	// initialize
	Node = File.GetObj(VAR_NAME, TRUE);

	// check
	if (GDS_Array_DimCnt(Node) != 3)
		throw ErrSeqArray(ERR_DIM, VAR_NAME);
	C_Int32 DLen[3];
	GDS_Array_GetDim(Node, DLen, 3);
	if ((DLen[0] < File.VariantNum()) || (DLen[1] != File.SampleNum()))
		throw ErrSeqArray(ERR_DIM, VAR_NAME);

	// initialize
	MarginalSize = File.VariantNum();
	MarginalSelect = File.Selection().pVariant();
	GenoIndex = &File.GenoIndex();
	SiteCount = ssize_t(DLen[1]) * DLen[2];
	_SampNum = File.SampleSelNum();
	CellCount = _SampNum * DLen[2];
	_Ploidy = File.Ploidy();
	UseRaw = use_raw;

	// initialize selection
	Selection.resize(SiteCount);
	C_BOOL *p = &Selection[0];
	memset(p, TRUE, SiteCount);
	C_BOOL *s = File.Selection().pSample();
	for (int n=DLen[1]; n > 0; n--)
	{
		if (*s++ == FALSE)
		{
			for (int m=DLen[2]; m > 0; m--) *p ++ = FALSE;
		} else {
			p += DLen[2];
		}
	}

	ExtPtr.reset(SiteCount);
	VarGeno = NULL;
	Reset();
}

int CApply_Variant_Geno::_ReadGenoData(int *Base)
{
	C_UInt8 NumIndexRaw;
	C_Int64 Index;
	GenoIndex->GetInfo(Position, Index, NumIndexRaw);

	if (NumIndexRaw >= 1)
	{
		CdIterator it;
		GDS_Iter_Position(Node, &it, Index*SiteCount);
		GDS_Iter_RDataEx(&it, Base, SiteCount, svInt32, &Selection[0]);

		const int bit_mask = 0x03;
		int missing = bit_mask;
		for (C_UInt8 i=1; i < NumIndexRaw; i++)
		{
			GDS_Iter_RDataEx(&it, ExtPtr.get(), SiteCount, svUInt8, &Selection[0]);

			C_UInt8 shift = i * 2;
			C_UInt8 *s = (C_UInt8*)ExtPtr.get();
			int *p = Base;
			for (ssize_t n=CellCount; n > 0; n--)
				*p++ |= int(*s++) << shift;

			missing = (missing << 2) | bit_mask;
		}

		return missing;
	} else {
		memset(Base, 0, sizeof(int)*CellCount);
		return 0;
	}
}

C_UInt8 CApply_Variant_Geno::_ReadGenoData(C_UInt8 *Base)
{
	C_UInt8 NumIndexRaw;
	C_Int64 Index;
	GenoIndex->GetInfo(Position, Index, NumIndexRaw);

	if (NumIndexRaw >= 1)
	{
		CdIterator it;
		GDS_Iter_Position(Node, &it, Index*SiteCount);
		GDS_Iter_RDataEx(&it, Base, SiteCount, svUInt8, &Selection[0]);

		const C_UInt8 bit_mask = 0x03;
		C_UInt8 missing = bit_mask;
		if (NumIndexRaw > 4)
		{
			NumIndexRaw = 4;
			warning("RAW type may not be sufficient to store genotypes.");
		}

		for (C_UInt8 i=1; i < NumIndexRaw; i++)
		{
			GDS_Iter_RDataEx(&it, ExtPtr.get(), SiteCount, svUInt8, &Selection[0]);

			C_UInt8 shift = i * 2;
			C_UInt8 *s = (C_UInt8*)ExtPtr.get();
			C_UInt8 *p = Base;
			for (ssize_t n=CellCount; n > 0; n--)
				*p++ |= (*s++) << shift;

			missing = (missing << 2) | bit_mask;
		}

		return missing;
	} else {
		memset(Base, 0, CellCount);
		return 0;
	}
}

void CApply_Variant_Geno::ReadData(SEXP val)
{
	if (UseRaw)
		ReadGenoData(RAW(val));
	else
		ReadGenoData(INTEGER(val));
}

SEXP CApply_Variant_Geno::NeedRData(int &nProtected)
{
	if (VarGeno == NULL)
	{
		VarGeno = UseRaw ? NEW_RAW(CellCount) : NEW_INTEGER(CellCount);
		PROTECT(VarGeno);
		nProtected ++;

		SEXP dim = NEW_INTEGER(2);
		int *p = INTEGER(dim);
		p[0] = _Ploidy; p[1] = _SampNum;
		SET_DIM(VarGeno, dim);

		SEXP name_list = PROTECT(NEW_LIST(2));
		SEXP tmp = PROTECT(NEW_CHARACTER(2));
		SET_STRING_ELT(tmp, 0, mkChar("allele"));
		SET_STRING_ELT(tmp, 1, mkChar("sample"));
		SET_NAMES(name_list, tmp);
		SET_DIMNAMES(VarGeno, name_list);
		UNPROTECT(2);
	}
	return VarGeno;
}

void CApply_Variant_Geno::ReadGenoData(int *Base)
{
	int missing = _ReadGenoData(Base);
	vec_i32_replace(Base, CellCount, missing, NA_INTEGER);
}

void CApply_Variant_Geno::ReadGenoData(C_UInt8 *Base)
{
	C_UInt8 missing = _ReadGenoData(Base);
	vec_i8_replace((C_Int8*)Base, CellCount, missing, NA_RAW);
}



// =====================================================================
// Object for reading genotypes variant by variant

CApply_Variant_Dosage::CApply_Variant_Dosage(CFileInfo &File, bool use_raw):
	CApply_Variant_Geno(File, use_raw)
{
	fVarType = ctDosage;
	ExtPtr.reset(sizeof(int)*CellCount);
}

void CApply_Variant_Dosage::ReadData(SEXP val)
{
	if (UseRaw)
		ReadDosage(RAW(val));
	else
		ReadDosage(INTEGER(val));
}

SEXP CApply_Variant_Dosage::NeedRData(int &nProtected)
{
	if (VarGeno == NULL)
	{
		VarGeno = UseRaw ? NEW_RAW(_SampNum) : NEW_INTEGER(_SampNum);
		PROTECT(VarGeno);
		nProtected ++;
	}
	return VarGeno;
}

void CApply_Variant_Dosage::ReadDosage(int *Base)
{
	int *p = (int *)ExtPtr.get();
	int missing = _ReadGenoData(p);

	// count the number of reference allele
	if (_Ploidy == 2) // diploid
	{
		vec_i32_cnt_dosage2(p, Base, _SampNum, 0, missing, NA_INTEGER);
	} else {
		for (int n=_SampNum; n > 0; n--)
		{
			int cnt = 0;
			for (int m=_Ploidy; m > 0; m--, p++)
			{
				if (*p == 0)
				{
					if (cnt != NA_INTEGER)
						cnt ++;
				} else if (*p == missing)
					cnt = NA_INTEGER;
			}
			*Base ++ = cnt;
		}
	}
}

void CApply_Variant_Dosage::ReadDosage(C_UInt8 *Base)
{
	C_UInt8 *p = (C_UInt8 *)ExtPtr.get();
	C_UInt8 missing = _ReadGenoData(p);

	// count the number of reference allele
	if (_Ploidy == 2) // diploid
	{
		vec_i8_cnt_dosage2((int8_t *)p, (int8_t *)Base, _SampNum, 0,
			missing, NA_RAW);
	} else {
		C_UInt8 *p = (C_UInt8 *)ExtPtr.get();
		for (int n=_SampNum; n > 0; n--)
		{
			C_UInt8 cnt = 0;
			for (int m=_Ploidy; m > 0; m--, p++)
			{
				if (*p == 0)
				{
					if (cnt != NA_RAW)
						cnt ++;
				} else if (*p == missing)
					cnt = NA_RAW;
			}
			*Base ++ = cnt;
		}
	}
}



// =====================================================================
// Object for reading phasing information variant by variant

CApply_Variant_Phase::CApply_Variant_Phase():
	CVarApply()
{
	fVarType = ctPhase;
	SiteCount = CellCount = 0;
	_SampNum = 0; _Ploidy = 0;
	UseRaw = FALSE;
	VarPhase = NULL;
}

CApply_Variant_Phase::CApply_Variant_Phase(CFileInfo &File, bool use_raw):
	CVarApply()
{
	fVarType = ctPhase;
	Init(File, use_raw);
}

void CApply_Variant_Phase::Init(CFileInfo &File, bool use_raw)
{
	static const char *VAR_NAME = "phase/data";

	// initialize
	Node = File.GetObj(VAR_NAME, TRUE);

	// check
	int DimCnt = GDS_Array_DimCnt(Node);
	if ((DimCnt != 2) && (DimCnt != 3))
		throw ErrSeqArray(ERR_DIM, VAR_NAME);
	C_Int32 DLen[3] = { 0, 0, 1 };
	GDS_Array_GetDim(Node, DLen, 3);
	if ((DLen[0] != File.VariantNum()) || (DLen[1] != File.SampleNum()))
		throw ErrSeqArray(ERR_DIM, VAR_NAME);

	// initialize
	MarginalSize = File.VariantNum();
	MarginalSelect = File.Selection().pVariant();
	SiteCount = ssize_t(DLen[1]) * DLen[2];
	_SampNum = File.SampleSelNum();
	CellCount = _SampNum * DLen[2];
	_Ploidy = File.Ploidy();
	UseRaw = use_raw;

	// initialize selection
	Selection.resize(SiteCount);
	C_BOOL *p = &Selection[0];
	memset(p, TRUE, SiteCount);
	C_BOOL *s = File.Selection().pSample();
	for (int n=DLen[1]; n > 0; n--)
	{
		if (*s++ == FALSE)
		{
			for (int m=DLen[2]; m > 0; m--) *p ++ = FALSE;
		} else {
			p += DLen[2];
		}
	}

	VarPhase = NULL;
	Reset();
}

void CApply_Variant_Phase::ReadData(SEXP val)
{
	CdIterator it;
	GDS_Iter_Position(Node, &it, ssize_t(Position)*SiteCount);
	if (UseRaw)
		GDS_Iter_RDataEx(&it, RAW(val), SiteCount, svInt8, &Selection[0]);
	else
		GDS_Iter_RDataEx(&it, INTEGER(val), SiteCount, svInt32, &Selection[0]);
}

SEXP CApply_Variant_Phase::NeedRData(int &nProtected)
{
	if (VarPhase == NULL)
	{
		VarPhase = UseRaw ? NEW_RAW(CellCount) : NEW_INTEGER(CellCount);
		PROTECT(VarPhase);
		nProtected ++;
		if (_Ploidy > 2)
		{
			SEXP dim = NEW_INTEGER(2);
			int *p = INTEGER(dim);
			p[0] = _Ploidy-1; p[1] = _SampNum;
			SET_DIM(VarPhase, dim);
		}
	}
	return VarPhase;
}



// =====================================================================
// Object for reading format variables variant by variant

/*
CApply_Variant_Format::CApply_Variant_Format()
{
	fVarType = ctFormat;
	
}

CApply_Variant_Format::CApply_Variant_Format(CFileInfo &File,
	const char *var_name, bool use_raw)
{
	fVarType = ctFormat;
	Init(File, var_name, use_raw);
}

void CApply_Variant_Format::Init(CFileInfo &File, const char *var_name,
	bool use_raw)
{


			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = File.GetObj(Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ERR_DIM, Path2.c_str());
			} else
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());

			SelPtr[1] = Sel.pSample();
			if (DimCnt > 2)
				SelPtr[2] = NeedTRUEs(DLen[2]);
			break;
}

void CApply_Variant_Format::ReadData(SEXP val)
{

}

SEXP CApply_Variant_Format::NeedRData(int &nProtected)
{

}

*/



// ===========================================================
// Object for reading a variable variant by variant

/// 
CApplyByVariant::CApplyByVariant(): CVarApply()
{
	Node = IndexNode = NULL;
	MarginalSelect = NULL;
	UseRaw = false;
}

void CApplyByVariant::InitObject(TVarType Type, const char *Path,
	CFileInfo &File, bool _UseRaw)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	// initialize
	GDS_PATH_PREFIX_CHECK(Path);
	fVarType = Type;
	Node = File.GetObj(Path, TRUE);
	SVType = GDS_Array_GetSVType(Node);
	DimCnt = GDS_Array_DimCnt(Node);

	int nVariant = File.VariantNum();
	// int nSample  = File.SampleNum();

	TSelection &Sel = File.Selection();
	TotalNum_Variant = File.VariantNum();
	MarginalSelect = &Sel.Variant[0];
	_SampNum = GetNumOfTRUE(&Sel.Sample[0], File.SampleNum());
	UseRaw = _UseRaw;
	NumOfBits = GDS_Array_GetBitOf(Node);

	string Path2; // the path with '@'
	switch (Type)
	{
		case ctBasic:
			// ===========================================================
			// VARIABLE: variant.id, position, allele
			if ((DimCnt != 1) || (GDS_Array_GetTotalCount(Node) != nVariant))
				throw ErrSeqArray(ERR_DIM, Path);
			break;

		case ctInfo:
			// ===========================================================
			// VARIABLE: info/...
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 2);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = File.GetObj(Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ERR_DIM, Path2.c_str());
			} else {
				if (DLen[0] != nVariant)
					throw ErrSeqArray(ERR_DIM, Path);
			}

			if (DimCnt > 1)
				SelPtr[1] = NeedTRUEs(DLen[1]);
			break;

		case ctFormat:
			// ===========================================================
			// VARIABLE: format/...
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = File.GetObj(Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ERR_DIM, Path2.c_str());
			} else
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());

			SelPtr[1] = Sel.pSample();
			if (DimCnt > 2)
				SelPtr[2] = NeedTRUEs(DLen[2]);
			break;

		default:
			throw ErrSeqArray("Internal Error in 'CApplyByVariant::InitObject'.");
	}

	Reset();
}


void CApplyByVariant::Reset()
{
	Position = 0;
	IndexRaw = 0;
	if (IndexNode)
	{
		C_Int32 Cnt=1;
		GDS_Array_ReadData(IndexNode, &Position, &Cnt, &NumIndexRaw, svInt32);
		if (NumIndexRaw < 0) NumIndexRaw = 0;
	} else
		NumIndexRaw = 1;

	if (!MarginalSelect[0]) Next();
}


bool CApplyByVariant::Next()
{
	Position ++;

	if (IndexNode)
	{
		IndexRaw += NumIndexRaw;
		C_Int32 Cnt=1, L;
		while ((Position<TotalNum_Variant) && !MarginalSelect[Position])
		{
			GDS_Array_ReadData(IndexNode, &Position, &Cnt, &L, svInt32);
			if (L > 0)
				IndexRaw += L;
			Position ++;
		}
		if (Position < TotalNum_Variant)
		{
			GDS_Array_ReadData(IndexNode, &Position, &Cnt, &NumIndexRaw,
				svInt32);
			if (NumIndexRaw < 0) NumIndexRaw = 0;
		} else
			NumIndexRaw = 0;
	} else {
		while ((Position<TotalNum_Variant) && !MarginalSelect[Position])
			Position ++;
		IndexRaw = Position;
		NumIndexRaw = 1;
	}

	return (Position < TotalNum_Variant);
}

void CApplyByVariant::ReadData(SEXP val)
{
	if (NumIndexRaw <= 0) return;

	switch (fVarType)
	{
	case ctGenotype:
	case ctDosage:
		throw "ERROR";

	default:
		C_Int32 st[3] = { IndexRaw, 0, 0 };
		DLen[0] = NumIndexRaw;
		SelPtr[0] = NeedTRUEs(NumIndexRaw);

		if (COREARRAY_SV_INTEGER(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, INTEGER(val), svInt32);
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, REAL(val), svFloat64);
		} else if (COREARRAY_SV_STRING(SVType))
		{
			vector<string> buffer(CellCount);
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, &buffer[0], svStrUTF8);
			for (int i=0; i < (int)buffer.size(); i++)
				SET_STRING_ELT(val, i, mkChar(buffer[i].c_str()));
		}
	}
}


SEXP CApplyByVariant::NeedRData(int &nProtected)
{
	if (NumIndexRaw <= 0) return R_NilValue;

	map<size_t, SEXP>::iterator it = VarList.find(NumIndexRaw);
	if (it == VarList.end())
	{
		switch (fVarType)
		{
		case ctBasic:
			CellCount = 1; break;
		case ctGenotype:
			CellCount = _SampNum * DLen[2]; break;
		case ctDosage:
			CellCount = _SampNum; break;
		case ctPhase:
			CellCount = (DimCnt>2) ? _SampNum*DLen[2] : _SampNum;
			break;
		case ctInfo:
			CellCount = ((DimCnt>1) ? DLen[1] : 1) * NumIndexRaw;
			break;
		case ctFormat:
			CellCount = ((DimCnt>2) ? _SampNum*DLen[2] : _SampNum) *
						NumIndexRaw;
			break;
		default:
			CellCount = 0;
		}

		SEXP ans=R_NilValue, dim;
		if (COREARRAY_SV_INTEGER(SVType))
		{
			if ((fVarType == ctGenotype) || (fVarType == ctDosage))
			{
				if (UseRaw)
					PROTECT(ans = NEW_RAW(CellCount));
				else
					PROTECT(ans = NEW_INTEGER(CellCount));
				nProtected ++;
			} else {
				char classname[64];
				classname[0] = 0;
				GDS_Node_GetClassName(Node, classname, sizeof(classname));
				if (strcmp(classname, "dBit1") == 0)
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else if (GDS_R_Is_Logical(Node))
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else {
					PROTECT(ans = NEW_INTEGER(CellCount));
					nProtected += GDS_R_Set_IfFactor(Node, ans);
				}
				nProtected ++;
			}
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			PROTECT(ans = NEW_NUMERIC(CellCount));
			nProtected ++;
		} else if (COREARRAY_SV_STRING(SVType))
		{
			PROTECT(ans = NEW_CHARACTER(CellCount));
			nProtected ++;
		}

		int *p;
		switch (fVarType)
		{
		case ctPhase:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = DLen[2]; p[1] = _SampNum;
				SET_DIM(ans, dim);
			}
			break;

		case ctFormat:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				p = INTEGER(dim = NEW_INTEGER(3));
				p[0] = DLen[2]; p[1] = _SampNum; p[2] = NumIndexRaw;
				SET_DIM(ans, dim);
			} else if (DimCnt == 2)
			{
				p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = _SampNum; p[1] = NumIndexRaw;
				SET_DIM(ans, dim);
			}
			break;

		default:
			break;
		}

		VarList.insert(pair<int, SEXP>(NumIndexRaw, ans));
		return ans;
	} else
		return it->second;
}

}


extern "C"
{
using namespace SeqArray;

// ===========================================================
// Apply functions over margins on a working space
// ===========================================================

COREARRAY_DLL_LOCAL const char *Txt_Apply_AsIs[] =
{
	"none", "list", "integer", "double", "character", "logical",
	"raw", NULL
};

COREARRAY_DLL_LOCAL const char *Txt_Apply_VarIdx[] =
{
	"none", "relative", "absolute", NULL
};


/// output to a connection
inline static void put_text(Rconnection file, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*file->vfprintf)(file, fmt, args);
	va_end(args);
}


/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho)
{
	int use_raw_flag = Rf_asLogical(RGetListElement(param, "useraw"));
	if (use_raw_flag == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");

	int prog_flag = Rf_asLogical(RGetListElement(param, "progress"));
	if (prog_flag == NA_LOGICAL)
		error("'.progress' must be TRUE or FALSE.");

	int dup_flag = Rf_asLogical(RGetListElement(param, "list_dup"));
	if (dup_flag == NA_LOGICAL)
		error("'.list_dup' must be TRUE or FALSE.");

	COREARRAY_TRY

		// the selection
		CFileInfo &File = GetFileInfo(gdsfile);

		// the number of calling PROTECT
		int nProtected = 0;

		// the number of selected variants
		int nVariant = File.VariantSelNum();
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");


		// ===========================================================
		// initialize the GDS Node list

		CVarApplyList NodeList;

		// for-loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));

			if ( s=="variant.id" || s=="position" || s=="chromosome" ||
				s=="allele" || s=="annotation/id" || s=="annotation/qual" ||
				s=="annotation/filter" )
			{
				NodeList.push_back(
					new CApply_Variant_Basic(File, s.c_str()));
			} else if (s == "genotype")
			{
				NodeList.push_back(
					new CApply_Variant_Geno(File, use_raw_flag!=FALSE));
			} else if (s == "phase")
			{
				NodeList.push_back(
					new CApply_Variant_Phase(File, use_raw_flag!=FALSE));
			} else if (strncmp(s.c_str(), "annotation/info/", 16) == 0)
			{
				CApplyByVariant *Obj = new CApplyByVariant;
				Obj->InitObject(CApplyByVariant::ctInfo, s.c_str(), File, use_raw_flag!=FALSE);
				NodeList.push_back(Obj);
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				s.append("/data");
				CApplyByVariant *Obj = new CApplyByVariant;
				Obj->InitObject(CApplyByVariant::ctFormat, s.c_str(), File, use_raw_flag!=FALSE);
				NodeList.push_back(Obj);
			} else if (s == "$dosage")
			{
				NodeList.push_back(
					new CApply_Variant_Dosage(File, use_raw_flag!=FALSE));
			} else {
				throw ErrSeqArray(
					"'%s' is not a standard variable name, and the standard format:\n"
					"\tvariant.id, position, chromosome, allele, "
					"annotation/id, annotation/qual, annotation/filter\n"
					"\tannotation/info/VARIABLE_NAME', annotation/format/VARIABLE_NAME",
					s.c_str());
			}
		}


		// ===========================================================
		// as.is

		Rconnection OutputConn = NULL;
		int DatType;
		if (!Rf_inherits(as_is, "connection"))
		{
			DatType = MatchText(CHAR(STRING_ELT(as_is, 0)), Txt_Apply_AsIs);
			if (DatType < 0)
				throw ErrSeqArray("'as.is' is not valid!");
		} else {
			OutputConn = R_GetConnection(as_is);
			DatType = 7;
		}

		C_Int8 *R_rv_ptr = NULL;
		switch (DatType)
		{
		case 1:
			rv_ans = PROTECT(NEW_LIST(nVariant)); nProtected ++;
			break;
		case 2:
			rv_ans = PROTECT(NEW_INTEGER(nVariant)); nProtected ++;
			R_rv_ptr = (C_Int8 *)INTEGER(rv_ans);
			break;
		case 3:
			rv_ans = PROTECT(NEW_NUMERIC(nVariant)); nProtected ++;
			R_rv_ptr = (C_Int8 *)REAL(rv_ans);
			break;
		case 4:
			rv_ans = PROTECT(NEW_CHARACTER(nVariant)); nProtected ++;
			break;
		case 5:
			rv_ans = PROTECT(NEW_LOGICAL(nVariant)); nProtected ++;
			R_rv_ptr = (C_Int8 *)LOGICAL(rv_ans);
			break;
		case 6:
			rv_ans = PROTECT(NEW_RAW(nVariant)); nProtected ++;
			R_rv_ptr = (C_Int8 *)RAW(rv_ans);
			break;
		}

		// ===========================================================
		// rho
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");


		// ===========================================================
		// initialize calling

		SEXP R_call_param = R_NilValue;
		if (NodeList.size() > 1)
		{
			PROTECT(R_call_param = NEW_LIST(NodeList.size()));
			nProtected ++;
			// set name to R_call_param
			SET_NAMES(R_call_param, GET_NAMES(var_name));
		}

		// ===============================================================
		// var.index
		int VarIdx = MatchText(CHAR(STRING_ELT(var_index, 0)), Txt_Apply_VarIdx);
		if (VarIdx < 0)
			throw ErrSeqArray("'var.index' is not valid!");

		SEXP R_fcall;
		SEXP R_Index = NULL;
		if (VarIdx > 0)
		{
			PROTECT(R_Index = NEW_INTEGER(1));
			nProtected ++;
			PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
			nProtected ++;
		} else {
			PROTECT(R_fcall = LCONS(FUN,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
			nProtected ++;
		}

		map<SEXP, SEXP> R_fcall_map;
		R_fcall_map[R_call_param] = R_fcall;


		// ===========================================================
		// for-loop calling

		CProgressStdOut progress(nVariant, prog_flag);

		int ans_index = 0;
		do {
			switch (VarIdx)
			{
			case 1:  // relative
				INTEGER(R_Index)[0] = ans_index + 1; break;
			case 2:  // absolute
				INTEGER(R_Index)[0] = NodeList[0]->Position + 1; break;
			}

			if (NodeList.size() <= 1)
			{
				R_call_param = NodeList[0]->NeedRData(nProtected);
				map<SEXP, SEXP>::iterator it = R_fcall_map.find(R_call_param);
				if (it == R_fcall_map.end())
				{
					if (VarIdx > 0)
					{
						PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
					} else {
						PROTECT(R_fcall = LCONS(FUN,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
					}
					nProtected ++;
					R_fcall_map[R_call_param] = R_fcall;
				} else
					R_fcall = it->second;

				NodeList[0]->ReadData(R_call_param);

			} else {
				CVarApply **p = &NodeList[0];
				size_t n = NodeList.size();
				for (size_t i=0; i < n; i++, p++)
				{
					SEXP tmp = (*p)->NeedRData(nProtected);
					(*p)->ReadData(tmp);
					SET_ELEMENT(R_call_param, i, tmp);
				}
			}

			// call R function
			SEXP val = eval(R_fcall, rho);

			// store data
			switch (DatType)
			{
			case 1:  // list
				if (dup_flag) val = duplicate(val);
				SET_ELEMENT(rv_ans, ans_index, val);
				break;
			case 2:  // integer
				*((int*)R_rv_ptr) = Rf_asInteger(val);
				R_rv_ptr += sizeof(int);
				break;
			case 3:  // double
				*((double*)R_rv_ptr) = Rf_asReal(val);
				R_rv_ptr += sizeof(double);
				break;
			case 4:  // character
				SET_STRING_ELT(rv_ans, ans_index, Rf_asChar(val));
				break;
			case 5:  // logical
				*((int*)R_rv_ptr) = Rf_asLogical(val);
				R_rv_ptr += sizeof(int);
				break;
			case 6:  // raw
				*R_rv_ptr = Rf_asInteger(val);
				R_rv_ptr ++;
				break;
			case 7:  // connection
				if (OutputConn->text)
				{
					if (Rf_isList(val))
					{
						throw ErrSeqArray("the user-defined function should return a character vector.");
					} else if (!Rf_isString(val))
					{
						val = AS_CHARACTER(val);
					}
					size_t n = XLENGTH(val);
					for (size_t i=0; i < n; i++)
					{
						put_text(OutputConn, "%s\n", CHAR(STRING_ELT(val, i)));
					}
				} else {
					if (TYPEOF(val) != RAWSXP)
						throw ErrSeqArray("the user-defined function should return a RAW vector.");
					size_t n = XLENGTH(val);
					size_t m = R_WriteConnection(OutputConn, RAW(val), n);
					if (n != m)
						throw ErrSeqArray("error in writing to a connection.");
				}
				break;
			}
			ans_index ++;

			progress.Forward();

		// check the end
		} while (NodeList.CallNext());

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

} // extern "C"
