// ===========================================================
//
// LinkSNPRelate.cpp: C interface for the SNPRelate package
//
// Copyright (C) 2015    Xiuwen Zheng
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
#include "ReadBySample.h"

#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C"
{
// TTypeGenoDim and TParam are also defined in "SNPRelate/src/dGenGWAS.h"

enum TTypeGenoDim
{
	RDim_Sample_X_SNP = 0,  ///< genotype matrix: sample X snp
	RDim_SNP_X_Sample = 1   ///< genotype matrix: snp X sample
};

typedef struct
{
	TTypeGenoDim *pGenoDimType;
	C_Int32 *pTotalSampleNum;
	C_Int32 *pTotalSNPNum;
	C_Int32 *pSampleNum;
	C_Int32 *pSNPNum;

	SEXP SeqGDSFile;
	CVarApply *Object;
	int *GenoBuffer;
	int Index;
} TParam;

inline static void Done_Object(TParam *Param)
{
	if (Param->Object)
	{
		delete Param->Object;
		Param->Object = NULL;
	}
	if (Param->GenoBuffer)
	{
		delete [](Param->GenoBuffer);
		Param->GenoBuffer = NULL;
	}
}

static void SNPRelate_InitSeqArray(TParam *Param)
{
	Done_Object(Param);

	// the GDS root node
	PdGDSObj Root = GDS_R_SEXP2FileRoot(Param->SeqGDSFile);
	PdAbstractArray N;

	N = GDS_Node_Path(Root, "sample.id", TRUE);
	*Param->pTotalSampleNum = GDS_Array_GetTotalCount(N);
	N = GDS_Node_Path(Root, "variant.id", TRUE);
	*Param->pTotalSNPNum = GDS_Array_GetTotalCount(N);

	*Param->pSampleNum = *Param->pSNPNum = 0;
	*Param->pGenoDimType = RDim_Sample_X_SNP;
}

static void SNPRelate_DoneSeqArray(TParam *Param)
{
	Done_Object(Param);
}

static void SNPRelate_InitSelSampOnly(C_BOOL *Sel, TParam *Param)
{
	int sum = 0;
	C_BOOL *p = Sel;
	for (int i=0; i < *Param->pTotalSampleNum; i++)
		if (*p++) sum ++;
	*Param->pSampleNum = sum;

	TInitObject::TSelection &s = Init.Selection(Param->SeqGDSFile);
	s.Sample.resize(*Param->pTotalSampleNum);
	memcpy(&s.Sample[0], Sel, *Param->pTotalSampleNum);

	Done_Object(Param);
}

static void SNPRelate_InitSelSNPOnly(C_BOOL *Sel, TParam *Param)
{
	int sum = 0;
	C_BOOL *p = Sel;
	for (int i=0; i < *Param->pTotalSNPNum; i++)
		if (*p++) sum ++;
	*Param->pSNPNum = sum;

	TInitObject::TSelection &s = Init.Selection(Param->SeqGDSFile);
	s.Variant.resize(*Param->pTotalSNPNum);
	memcpy(&s.Variant[0], Sel, *Param->pTotalSNPNum);

	Done_Object(Param);
}

static void SNPRelate_SnpRead(C_Int32 SnpStart, C_Int32 SnpCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim, TParam *Param)
{
	CVarApplyByVariant *Obj =
		(CVarApplyByVariant*)(Param->Object);

	if (!Obj)
	{
		Obj = new CVarApplyByVariant;
		Param->Object = Obj;

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(Param->SeqGDSFile);
		TInitObject::TSelection &Sel = Init.Selection(Param->SeqGDSFile);
		Obj->InitObject(CVariable::ctGenotype,
			"genotype/data", Root, Sel.Variant.size(),
			&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);

		size_t SIZE = (Obj->Num_Sample) * (Obj->DLen[2]);
		Param->GenoBuffer = new int[SIZE];
		Param->Index = 0;
	}

	if (Param->Index > SnpStart)
	{
		Obj->ResetObject();
		Param->Index = 0;
	}
	while (Param->Index < SnpStart)
	{
		Obj->NextCell();
		Param->Index ++;
	}

	for (int sn = SnpCount; sn > 0; sn--)
	{
		Obj->ReadGenoData(Param->GenoBuffer);
		Obj->NextCell();
		Param->Index ++;

		if (OutDim == RDim_Sample_X_SNP)
		{
			int *p = Param->GenoBuffer;
			for (size_t n=Obj->Num_Sample; n > 0; n--)
			{
				int *pp = p;
				p += Obj->DLen[2];
				C_UInt8 val = 0;
				for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
				{
					if (*pp == 0)
					{
						val ++;
						if (val > 2) val = 2;
					} else if (*pp == NA_INTEGER)
					{
						val = 3; break;
					}
				}
				*OutBuf ++ = val;
			}
		} else {
			C_UInt8 *g = (OutBuf ++);
			int *p = Param->GenoBuffer;
			for (size_t n=Obj->Num_Sample; n > 0; n--)
			{
				int *pp = p;
				p += Obj->DLen[2];
				C_UInt8 val = 0;
				for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
				{
					if (*pp == 0)
					{
						val ++;
						if (val > 2) val = 2;
					} else if (*pp == NA_INTEGER)
					{
						val = 3; break;
					}
				}
				*g = val; g += SnpCount;
			}
		}
	}
}

static void SNPRelate_SampleRead(C_Int32 SampStart, C_Int32 SampCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim, TParam *Param)
{
	CVarApplyBySample *Obj =
		(CVarApplyBySample*)(Param->Object);

	if (!Obj)
	{
		Obj = new CVarApplyBySample;
		Param->Object = Obj;

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(Param->SeqGDSFile);
		TInitObject::TSelection &Sel = Init.Selection(Param->SeqGDSFile);
		Obj->InitObject(CVariable::ctGenotype,
			"genotype/data", Root, Sel.Variant.size(),
			&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);

		size_t SIZE = (Obj->Num_Variant) * (Obj->DLen[2]);
		Param->GenoBuffer = new int[SIZE];
		Param->Index = 0;
	}

	if (Param->Index > SampStart)
	{
		Obj->ResetObject();
		Param->Index = 0;
	}
	while (Param->Index < SampStart)
	{
		Obj->NextCell();
		Param->Index ++;
	}

	for (int sn = SampCount; sn > 0; sn--)
	{
		Obj->ReadGenoData(Param->GenoBuffer);
		Obj->NextCell();
		Param->Index ++;

		if (OutDim == RDim_SNP_X_Sample)
		{
			int *p = Param->GenoBuffer;
			for (size_t n=Obj->Num_Variant; n > 0; n--)
			{
				int *pp = p;
				p += Obj->DLen[2];
				C_UInt8 val = 0;
				for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
				{
					if (*pp == 0)
					{
						val ++;
						if (val > 2) val = 2;
					} else if (*pp == NA_INTEGER)
					{
						val = 3; break;
					}
				}
				*OutBuf ++ = val;
			}

		} else {
			C_UInt8 *g = (OutBuf ++);
			int *p = Param->GenoBuffer;
			for (size_t n=Obj->Num_Variant; n > 0; n--)
			{
				int *pp = p;
				p += Obj->DLen[2];
				C_UInt8 val = 0;
				for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
				{
					if (*pp == 0)
					{
						val ++;
						if (val > 2) val = 2;
					} else if (*pp == NA_INTEGER)
					{
						val = 3; break;
					}
				}
				*g = val; g += SampCount;
			}
		}
	}
}

static void SNPRelate_SetSnpSelection(C_BOOL *sel, TParam *Param)
{
	TInitObject::TSelection &s = Init.Selection(Param->SeqGDSFile);
	C_BOOL *p = &s.Variant[0];

	int sum = 0;
	for (int i=0; i < *Param->pTotalSNPNum; i++, p++)
	{
		if (*p)
		{
			if (*sel ++)
				sum ++;
			else
				*p = FALSE;
		}
	}
	*Param->pSNPNum = sum;

	Done_Object(Param);
}

static void SNPRelate_SetSampSelection(C_BOOL *sel, TParam *Param)
{
	TInitObject::TSelection &s = Init.Selection(Param->SeqGDSFile);
	C_BOOL *p = &s.Sample[0];

	int sum = 0;
	for (int i=0; i < *Param->pTotalSampleNum; i++, p++)
	{
		if (*p)
		{
			if (*sel ++)
				sum ++;
			else
				*p = FALSE;
		}
	}
	*Param->pSampleNum = sum;

	Done_Object(Param);
}

COREARRAY_DLL_LOCAL void Register_SNPRelate_Functions()
{
	static const char *pkg_name = "SeqArray";

	#define REG(nm)    R_RegisterCCallable(pkg_name, #nm, (DL_FUNC)&nm)

	REG(SNPRelate_InitSeqArray);
	REG(SNPRelate_DoneSeqArray);
	REG(SNPRelate_InitSelSampOnly);
	REG(SNPRelate_InitSelSNPOnly);
	REG(SNPRelate_SnpRead);
	REG(SNPRelate_SampleRead);
	REG(SNPRelate_SetSnpSelection);
	REG(SNPRelate_SetSampSelection);
}

} // extern "C"
