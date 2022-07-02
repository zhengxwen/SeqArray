// ===========================================================
//
// ReadByVariant.cpp: Read data variant by variant
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

#include "ReadByVariant.h"


namespace SeqArray
{

using namespace Vectorization;

static const char *ERR_DIM = "Invalid dimension of '%s'.";


// =====================================================================
// Object for reading basic variables variant by variant

CApply_Variant_Basic::CApply_Variant_Basic(CFileInfo &File,
	const char *var_name): CApply_Variant(File)
{
	fVarType = ctBasic;
	Node = File.GetObj(var_name, TRUE);
	SVType = GDS_Array_GetSVType(Node);
	VarNode = NULL;
	Reset();
}

void CApply_Variant_Basic::ReadData(SEXP val)
{
	C_Int32 st = Position, one = 1;
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

// ====

CApply_Variant_Pos::CApply_Variant_Pos(CFileInfo &File):
	CApply_Variant(File)
{
	fVarType = ctBasic;
	Node = File.GetObj("position", TRUE);
	PtrPos = &File.Position()[0];
	VarNode = NULL;
	Reset();
}

void CApply_Variant_Pos::ReadData(SEXP val)
{
	INTEGER(val)[0] = PtrPos[Position];
}

SEXP CApply_Variant_Pos::NeedRData(int &nProtected)
{
	if (VarNode == NULL)
	{
		VarNode = PROTECT(NEW_INTEGER(1));
		nProtected ++;
	}
	return VarNode;
}

// ====

CApply_Variant_Chrom::CApply_Variant_Chrom(CFileInfo &File):
	CApply_Variant(File)
{
	fVarType = ctBasic;
	Node = File.GetObj("chromosome", TRUE);
	ChromIndex = &File.Chromosome();
	VarNode = NULL;
	Reset();
}

void CApply_Variant_Chrom::ReadData(SEXP val)
{
	const string &s1 = (*ChromIndex)[Position];
	const char *s2 = CHAR(STRING_ELT(val, 0));
	if (s1 != s2)
		SET_STRING_ELT(val, 0, mkChar(s1.c_str()));
}

SEXP CApply_Variant_Chrom::NeedRData(int &nProtected)
{
	if (VarNode == NULL)
	{
		VarNode = PROTECT(mkString(""));
		nProtected ++;
	}
	return VarNode;
}



// =====================================================================
// Object for reading genotypes variant by variant

static inline void read_geno(CdIterator &it, int *out, TSelection::TSampStruct *p)
{
	C_Int64 base = it.Ptr;
	for (; p->length > 0; p++)
	{
		it.Ptr = base + p->offset;
		if (!p->sel)
			out = (int*)GDS_Iter_RData(&it, out, p->length, svInt32);
		else
			out = (int*)GDS_Iter_RDataEx(&it, out, p->length, svInt32, p->sel);
	}
}

static inline void read_geno(CdIterator &it, C_UInt8 *out, TSelection::TSampStruct *p)
{
	C_Int64 base = it.Ptr;
	for (; p->length > 0; p++)
	{
		it.Ptr = base + p->offset;
		if (!p->sel)
			out = (C_UInt8*)GDS_Iter_RData(&it, out, p->length, svUInt8);
		else
			out = (C_UInt8*)GDS_Iter_RDataEx(&it, out, p->length, svUInt8, p->sel);
	}
}


CApply_Variant_Geno::CApply_Variant_Geno(): CApply_Variant()
{
	fVarType = ctGenotype;
	SiteCount = CellCount = 0;
	SampNum = 0; Ploidy = 0;
	UseRaw = FALSE;
	VarIntGeno = VarRawGeno = NULL;
}

CApply_Variant_Geno::CApply_Variant_Geno(CFileInfo &File, int use_raw):
	CApply_Variant()
{
	fVarType = ctGenotype;
	Init(File, use_raw);
}

void CApply_Variant_Geno::Init(CFileInfo &File, int use_raw)
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
	InitMarginal(File);
	GenoIndex = &File.GenoIndex();
	SiteCount = ssize_t(DLen[1]) * DLen[2];
	SampNum = File.SampleSelNum();
	CellCount = SampNum * DLen[2];
	Ploidy = File.Ploidy();
	UseRaw = use_raw;

	// initialize selection
	pSampSel = File.Selection().GetStructSample();

	ExtPtr.reset(SiteCount);
	VarIntGeno = VarRawGeno = NULL;
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
		read_geno(it, Base, pSampSel);

		const int bit_mask = 0x03;
		int missing = bit_mask;
		for (C_UInt8 i=1; i < NumIndexRaw; i++)
		{
			GDS_Iter_Position(Node, &it, (Index+i)*SiteCount);
			read_geno(it, (C_UInt8*)ExtPtr.get(), pSampSel);

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
		read_geno(it, Base, pSampSel);

		const C_UInt8 bit_mask = 0x03;
		C_UInt8 missing = bit_mask;
		if (NumIndexRaw > 4)
		{
			NumIndexRaw = 4;
			warning("RAW type may not be sufficient to store genotypes.");
		}

		for (C_UInt8 i=1; i < NumIndexRaw; i++)
		{
			GDS_Iter_Position(Node, &it, (Index+i)*SiteCount);
			read_geno(it, (C_UInt8*)ExtPtr.get(), pSampSel);

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
	switch (TYPEOF(val))
	{
		case INTSXP:
			ReadGenoData(INTEGER(val)); break;
		case RAWSXP:
			ReadGenoData(RAW(val)); break;
		default:
			throw ErrSeqArray("Invalid type (%d) in CApply_Variant_Geno::ReadData()",
				(int)TYPEOF(val));
	}
}

SEXP CApply_Variant_Geno::NeedRData(int &nProtected)
{
	bool int_type;
	if (UseRaw == NA_INTEGER)
	{
		C_UInt8 NumIndexRaw;
		C_Int64 Index;
		GenoIndex->GetInfo(Position, Index, NumIndexRaw);
		int_type = (NumIndexRaw > 4);
	} else if (UseRaw == FALSE)
		int_type = true;
	else
		int_type = false;

	if (int_type)
	{
		if (VarIntGeno == NULL)
		{
			VarIntGeno = PROTECT(allocMatrix(INTSXP, Ploidy, SampNum));
			nProtected ++;
			SET_DIMNAMES(VarIntGeno, R_Geno_Dim2_Name);
		}
		return VarIntGeno;
	} else {
		if (VarRawGeno == NULL)
		{
			VarRawGeno = PROTECT(allocMatrix(RAWSXP, Ploidy, SampNum));
			nProtected ++;
			SET_DIMNAMES(VarRawGeno, R_Geno_Dim2_Name);
		}
		return VarRawGeno;
	}
}

bool CApply_Variant_Geno::NeedIntType()
{
	const C_BOOL *b = MarginalSelect;
	C_Int32 P = Position;
	while (P < MarginalEnd)
	{
		C_UInt8 NumIndexRaw;
		C_Int64 Index;
		GenoIndex->GetInfo(P, Index, NumIndexRaw);
		if (NumIndexRaw > 4) return true;
		// find next
		P = VEC_BOOL_FIND_TRUE(b+P+1, b+MarginalEnd) - b;
	}
	return false;
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

CApply_Variant_Dosage::CApply_Variant_Dosage(CFileInfo &File, int use_raw, bool alt):
	CApply_Variant_Geno(File, use_raw)
{
	fVarType = ctDosage;
	IsAlt = alt;
	ExtPtr2.reset(sizeof(int)*CellCount);
	VarDosage = NULL;
}

void CApply_Variant_Dosage::ReadData(SEXP val)
{
	switch (TYPEOF(val))
	{
		case INTSXP:
			if (IsAlt) ReadDosageAlt(INTEGER(val)); else ReadDosage(INTEGER(val));
			break;
		case RAWSXP:
			if (IsAlt) ReadDosageAlt(RAW(val)); else ReadDosage(RAW(val));
			break;
		default:
			throw ErrSeqArray("Invalid type (%d) in CApply_Variant_Dosage::ReadData()",
				(int)TYPEOF(val));
	}
}

SEXP CApply_Variant_Dosage::NeedRData(int &nProtected)
{
	if (VarDosage == NULL)
	{
		VarDosage = UseRaw ? NEW_RAW(SampNum) : NEW_INTEGER(SampNum);
		PROTECT(VarDosage);
		nProtected ++;
	}
	return VarDosage;
}

void CApply_Variant_Dosage::ReadDosage(int *Base)
{
	int *p = (int *)ExtPtr2.get();
	int missing = _ReadGenoData(p);

	// count the number of reference allele
	if (Ploidy == 2) // diploid
	{
		vec_i32_cnt_dosage2(p, Base, SampNum, 0, missing, NA_INTEGER);
	} else {
		for (int n=SampNum; n > 0; n--)
		{
			int cnt = 0;
			for (int m=Ploidy; m > 0; m--, p++)
			{
				if (*p == 0)
				{
					if (cnt != NA_INTEGER) cnt ++;
				} else if (*p == missing)
					cnt = NA_INTEGER;
			}
			*Base ++ = cnt;
		}
	}
}

void CApply_Variant_Dosage::ReadDosageAlt(int *Base)
{
	int *p = (int *)ExtPtr2.get();
	int missing = _ReadGenoData(p);

	// count the number of reference allele
	if (Ploidy == 2) // diploid
	{
		vec_i32_cnt_dosage_alt2(p, Base, SampNum, 0, missing, NA_INTEGER);
	} else {
		for (int n=SampNum; n > 0; n--)
		{
			int cnt = 0;
			for (int m=Ploidy; m > 0; m--, p++)
			{
				if (*p != 0)
				{
					if (cnt != NA_INTEGER) cnt ++;
				} else if (*p == missing)
					cnt = NA_INTEGER;
			}
			*Base ++ = cnt;
		}
	}
}

void CApply_Variant_Dosage::ReadDosage(C_UInt8 *Base)
{
	C_UInt8 *p = (C_UInt8 *)ExtPtr2.get();
	C_UInt8 missing = _ReadGenoData(p);

	// count the number of reference allele
	if (Ploidy == 2) // diploid
	{
		vec_i8_cnt_dosage2((int8_t *)p, (int8_t *)Base, SampNum, 0,
			missing, NA_RAW);
	} else {
		C_UInt8 *p = (C_UInt8 *)ExtPtr.get();
		for (int n=SampNum; n > 0; n--)
		{
			C_UInt8 cnt = 0;
			for (int m=Ploidy; m > 0; m--, p++)
			{
				if (*p == 0)
				{
					if (cnt != NA_RAW) cnt ++;
				} else if (*p == missing)
					cnt = NA_RAW;
			}
			*Base ++ = cnt;
		}
	}
}

void CApply_Variant_Dosage::ReadDosageAlt(C_UInt8 *Base)
{
	C_UInt8 *p = (C_UInt8 *)ExtPtr2.get();
	C_UInt8 missing = _ReadGenoData(p);

	// count the number of reference allele
	if (Ploidy == 2) // diploid
	{
		vec_i8_cnt_dosage_alt2((int8_t *)p, (int8_t *)Base, SampNum, 0,
			missing, NA_RAW);
	} else {
		C_UInt8 *p = (C_UInt8 *)ExtPtr.get();
		for (int n=SampNum; n > 0; n--)
		{
			C_UInt8 cnt = 0;
			for (int m=Ploidy; m > 0; m--, p++)
			{
				if (*p != 0)
				{
					if (cnt != NA_RAW) cnt ++;
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
	CApply_Variant()
{
	fVarType = ctPhase;
	SiteCount = CellCount = 0;
	SampNum = 0; Ploidy = 0;
	UseRaw = FALSE;
	VarPhase = NULL;
}

CApply_Variant_Phase::CApply_Variant_Phase(CFileInfo &File, bool use_raw):
	CApply_Variant()
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
	InitMarginal(File);
	SiteCount = ssize_t(DLen[1]) * DLen[2];
	SampNum = File.SampleSelNum();
	CellCount = SampNum * DLen[2];
	Ploidy = File.Ploidy();
	UseRaw = use_raw;

	// initialize selection
	Selection.resize(SiteCount);
	C_BOOL *p = &Selection[0];
	memset(p, TRUE, SiteCount);
	C_BOOL *s = File.Selection().pSample;
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
		if (Ploidy > 2)
		{
			SEXP dim = NEW_INTEGER(2);
			int *p = INTEGER(dim);
			p[0] = Ploidy-1; p[1] = SampNum;
			SET_DIM(VarPhase, dim);
		}
	}
	return VarPhase;
}



// =====================================================================
// Object for reading info variables variant by variant

CApply_Variant_Info::CApply_Variant_Info(CFileInfo &File,
	const char *var_name): CApply_Variant(File)
{
	// initialize
	fVarType = ctInfo;
	Node = File.GetObj(var_name, TRUE);

	// check
	int DimCnt = GDS_Array_DimCnt(Node);
	if ((DimCnt != 1) && (DimCnt != 2))
		throw ErrSeqArray(ERR_DIM, var_name);

	// initialize
	C_Int32 DLen[2];
	GDS_Array_GetDim(Node, DLen, 2);
	BaseNum = (DimCnt == 2) ? DLen[1] : 1;
	VarIndex = &VarGetStruct(File, var_name).Index;
	SVType = GDS_Array_GetSVType(Node);

	Reset();
}

void CApply_Variant_Info::ReadData(SEXP val)
{
	C_Int64 IndexRaw;
	int NumIndexRaw;
	VarIndex->GetInfo(Position, IndexRaw, NumIndexRaw);

	if (NumIndexRaw > 0)
	{
		C_Int32 st[2]  = { (C_Int32)IndexRaw, 0 };
		C_Int32 cnt[2] = { NumIndexRaw, BaseNum };

		if (COREARRAY_SV_INTEGER(SVType))
		{
			GDS_Array_ReadData(Node, st, cnt, INTEGER(val), svInt32);
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			GDS_Array_ReadData(Node, st, cnt, REAL(val), svFloat64);
		} else if (COREARRAY_SV_STRING(SVType))
		{
			vector<string> buffer(XLENGTH(val));
			GDS_Array_ReadData(Node, st, cnt, &buffer[0], svStrUTF8);
			for (size_t i=0; i < buffer.size(); i++)
				SET_STRING_ELT(val, i, mkChar(buffer[i].c_str()));
		}
	}
}

SEXP CApply_Variant_Info::NeedRData(int &nProtected)
{
	C_Int64 IndexRaw;
	int NumIndexRaw;
	VarIndex->GetInfo(Position, IndexRaw, NumIndexRaw);
	if (NumIndexRaw <= 0) return R_NilValue;

	map<int, SEXP>::iterator it = VarList.find(NumIndexRaw);
	if (it == VarList.end())
	{
		SEXP ans = RObject_GDS(Node, BaseNum*NumIndexRaw, nProtected, true);
		if (BaseNum > 1)
		{
			SEXP dim = NEW_INTEGER(2);
			int *p = INTEGER(dim);
			p[0] = BaseNum; p[1] = NumIndexRaw;
			SET_DIM(ans, dim);
		}

		VarList.insert(pair<int, SEXP>(NumIndexRaw, ans));
		return ans;
	} else
		return it->second;
}



// =====================================================================
// Object for reading format variables variant by variant

CApply_Variant_Format::CApply_Variant_Format(): CApply_Variant()
{
	fVarType = ctFormat;
}

CApply_Variant_Format::CApply_Variant_Format(CFileInfo &File,
	const char *var_name): CApply_Variant()
{
	fVarType = ctFormat;
	Init(File, var_name);
}

void CApply_Variant_Format::Init(CFileInfo &File, const char *var_name)
{
	// initialize
	Node = File.GetObj((string(var_name)+"/data").c_str(), TRUE);

	// check
	int DimCnt = GDS_Array_DimCnt(Node);
	if (DimCnt != 2)
		throw ErrSeqArray(ERR_DIM, var_name);
	C_Int32 DLen[2];
	GDS_Array_GetDim(Node, DLen, 2);
	if (DLen[1] != File.SampleNum())
		throw ErrSeqArray(ERR_DIM, var_name);

	// initialize
	InitMarginal(File);
	SVType = GDS_Array_GetSVType(Node);
	VarIndex = &VarGetStruct(File, var_name).Index;
	SampNum = File.SampleSelNum();
	_TotalSampNum = File.SampleNum();

	// initialize selection
	SelPtr[0] = NULL;
	SelPtr[1] = File.Selection().pSample;

	Reset();
}

void CApply_Variant_Format::ReadData(SEXP val)
{
	C_Int64 IndexRaw;
	int NumIndexRaw;
	VarIndex->GetInfo(Position, IndexRaw, NumIndexRaw);

	if (NumIndexRaw > 0)
	{
		C_Int32 st[2]  = { (C_Int32)IndexRaw, 0 };
		C_Int32 cnt[2] = { NumIndexRaw, (C_Int32)_TotalSampNum };
		SelPtr[0] = NeedTRUEs(NumIndexRaw);

		if (COREARRAY_SV_INTEGER(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, cnt, SelPtr, INTEGER(val), svInt32);
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, cnt, SelPtr, REAL(val), svFloat64);
		} else if (COREARRAY_SV_STRING(SVType))
		{
			vector<string> buffer(XLENGTH(val));
			GDS_Array_ReadDataEx(Node, st, cnt, SelPtr, &buffer[0], svStrUTF8);
			for (size_t i=0; i < buffer.size(); i++)
				SET_STRING_ELT(val, i, mkChar(buffer[i].c_str()));
		}
	}
}

SEXP CApply_Variant_Format::NeedRData(int &nProtected)
{
	C_Int64 IndexRaw;
	int NumIndexRaw;
	VarIndex->GetInfo(Position, IndexRaw, NumIndexRaw);
	if (NumIndexRaw <= 0) return R_NilValue;

	map<int, SEXP>::iterator it = VarList.find(NumIndexRaw);
	if (it == VarList.end())
	{
		SEXP ans = RObject_GDS(Node, SampNum*NumIndexRaw, nProtected, false);
		SEXP dim = NEW_INTEGER(2);
		int *p = INTEGER(dim);
		p[0] = SampNum; p[1] = NumIndexRaw;
		SET_DIM(ans, dim);

		SEXP name_list = PROTECT(NEW_LIST(2));
		SEXP tmp = PROTECT(NEW_CHARACTER(2));
		SET_STRING_ELT(tmp, 0, mkChar("sample"));
		SET_STRING_ELT(tmp, 1, mkChar("index"));
		SET_NAMES(name_list, tmp);
		SET_DIMNAMES(ans, name_list);
		UNPROTECT(2);

		VarList.insert(pair<int, SEXP>(NumIndexRaw, ans));
		return ans;
	} else
		return it->second;
}



// =====================================================================
// Object for reading format variables variant by variant

CApply_Variant_NumAllele::CApply_Variant_NumAllele(CFileInfo &File):
	CApply_Variant(File)
{
	strbuf.reserve(128);
	fVarType = ctBasic;
	Node = File.GetObj("allele", TRUE);
	VarNode = NULL;
	Reset();
}

void CApply_Variant_NumAllele::ReadData(SEXP val)
{
	INTEGER(val)[0] = GetNumAllele();
}

SEXP CApply_Variant_NumAllele::NeedRData(int &nProtected)
{
	if (VarNode == NULL)
	{
		VarNode = PROTECT(NEW_INTEGER(1));
		nProtected ++;
	}
	return VarNode;
}

int CApply_Variant_NumAllele::GetNumAllele()
{
	C_Int32 st = Position, one = 1;
	GDS_Array_ReadData(Node, &st, &one, &strbuf, svStrUTF8);
	return GetNumOfAllele(strbuf.c_str());
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


/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho)
{
	SEXP pam_use_raw = RGetListElement(param, "useraw");
	if (!Rf_isLogical(pam_use_raw))
		error("'.useraw' must be TRUE, FALSE or NA.");
	int use_raw_flag = Rf_asLogical(pam_use_raw);

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
		{
			if (Rf_inherits(as_is, "connection") || Rf_inherits(as_is, "gdsn.class"))
			{
				return R_NilValue;
			} else {
				int dt = MatchText(CHAR(STRING_ELT(as_is, 0)), Txt_Apply_AsIs);
				switch (dt)
				{
				case 0:
					return R_NilValue;
				case 1:
					return NEW_LIST(nVariant);
				case 2:
					return NEW_INTEGER(nVariant);
				case 3:
					return NEW_NUMERIC(nVariant);
				case 4:
					return NEW_CHARACTER(nVariant);
				case 5:
					return NEW_LOGICAL(nVariant);
				case 6:
					return NEW_RAW(nVariant);
				default:
					throw ErrSeqArray("'as.is' is not valid!");
				}
			}
		}


		// ===========================================================
		// initialize the GDS Node list

		CVarApplyList NodeList;

		// for-loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));

			if (s=="variant.id" || s=="allele" || s=="annotation/id" ||
				s=="annotation/qual" || s=="annotation/filter")
			{
				NodeList.push_back(
					new CApply_Variant_Basic(File, s.c_str()));
			} else if (s == "position")
			{
				NodeList.push_back(new CApply_Variant_Pos(File));
			} else if (s == "chromosome")
			{
				NodeList.push_back(new CApply_Variant_Chrom(File));
			} else if (s == "genotype")
			{
				NodeList.push_back(
					new CApply_Variant_Geno(File, use_raw_flag));
			} else if (s == "phase")
			{
				NodeList.push_back(
					new CApply_Variant_Phase(File, use_raw_flag!=FALSE));
			} else if (strncmp(s.c_str(), "annotation/info/", 16) == 0)
			{
				NodeList.push_back(
					new CApply_Variant_Info(File, s.c_str()));
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				NodeList.push_back(
					new CApply_Variant_Format(File, s.c_str()));
			} else if (s == "$dosage")
			{
				NodeList.push_back(
					new CApply_Variant_Dosage(File, use_raw_flag, false));
			} else if (s == "$dosage_alt")
			{
				NodeList.push_back(
					new CApply_Variant_Dosage(File, use_raw_flag, true));
			} else if (s == "$num_allele")
			{
				NodeList.push_back(new CApply_Variant_NumAllele(File));
			} else {
				throw ErrSeqArray(
					"'%s' is not a standard variable name, and the standard format:\n"
					"    variant.id, position, chromosome, allele, annotation/id, annotation/qual, annotation/filter\n"
					"    annotation/info/VARIABLE_NAME', annotation/format/VARIABLE_NAME, etc",
					s.c_str());
			}
		}


		// ===========================================================
		// as.is

		Rconnection OutputConn = NULL;
		PdGDSObj OutputGDS = NULL;
		int DatType;
		if (Rf_inherits(as_is, "connection"))
		{
			OutputConn = R_GetConnection(as_is);
			DatType = 7;
		} else if (Rf_inherits(as_is, "gdsn.class"))
		{
			OutputGDS = GDS_R_SEXP2Obj(as_is, FALSE);
			DatType = 8;
		} else {
			DatType = MatchText(CHAR(STRING_ELT(as_is, 0)), Txt_Apply_AsIs);
			if (DatType < 0)
				throw ErrSeqArray("'as.is' is not valid!");
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


		// ===========================================================
		// for-loop calling

		CProgressStdOut progress(nVariant, 1, prog_flag);

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
				SEXP pm = NodeList[0]->NeedRData(nProtected);
				if (pm != R_call_param)
				{
					if (VarIdx > 0)
						SETCADDR(R_fcall, pm);
					else
						SETCADR(R_fcall, pm);
					R_call_param = pm;
				}
				NodeList[0]->ReadData(pm);
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
						ConnPutText(OutputConn, "%s\n", CHAR(STRING_ELT(val, i)));
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
			case 8:  // gdsn.class
				RAppendGDS(OutputGDS, val);
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
