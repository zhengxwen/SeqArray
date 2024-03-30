// ===========================================================
//
// LinkSNPRelate.cpp: C interface for the SNPRelate package
//
// Copyright (C) 2015-2018    Xiuwen Zheng
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
using namespace SeqArray;

// TTypeGenoDim and TParam are also defined in "SNPRelate/src/dGenGWAS.h"

enum TTypeGenoDim
{
	RDim_Sample_X_SNP = 0,  ///< genotype matrix: sample X snp
	RDim_SNP_X_Sample = 1   ///< genotype matrix: snp X sample
};

typedef struct COREARRAY_DLL_LOCAL SParam
{
	TTypeGenoDim *pGenoDimType;
	C_Int32 *pTotalSampleNum;
	C_Int32 *pTotalSNPNum;
	C_Int32 *pSampleNum;
	C_Int32 *pSNPNum;

	SEXP SeqGDSFile;
	CVarApply *Object;
	C_UInt8 *GenoBuffer;
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
	*Param->pSampleNum = GetNumOfTRUE(Sel, *Param->pTotalSampleNum);

	CFileInfo &File = GetFileInfo(Param->SeqGDSFile);
	TSelection &s = File.Selection();
	s.ClearStructSample();
	memcpy(s.pSample, Sel, *Param->pTotalSampleNum);

	Done_Object(Param);
}


static void SNPRelate_InitSelSNPOnly(C_BOOL *Sel, TParam *Param)
{
	*Param->pSNPNum = GetNumOfTRUE(Sel, *Param->pTotalSNPNum);

	CFileInfo &File = GetFileInfo(Param->SeqGDSFile);
	TSelection &s = File.Selection();
	s.ClearStructVariant();
	memcpy(s.pVariant, Sel, *Param->pTotalSNPNum);

	Done_Object(Param);
}


static void SNPRelate_SnpRead(C_Int32 SnpStart, C_Int32 SnpCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim, TParam *Param)
{
	if (Param->Object && !dynamic_cast<CApply_Variant_Dosage*>(Param->Object))
	{
		delete Param->Object;
		Param->Object = NULL;
	}

	CApply_Variant_Dosage *Obj = (CApply_Variant_Dosage*)(Param->Object);
	if (!Obj)
	{
		Obj = new CApply_Variant_Dosage(GetFileInfo(Param->SeqGDSFile), true, false);
		Param->Object = Obj;
		Param->GenoBuffer = new C_UInt8[Obj->SampNum];
		Param->Index = 0;
	}

	if (Param->Index > SnpStart)
	{
		Obj->Reset();
		Param->Index = 0;
	}
	while (Param->Index < SnpStart)
	{
		Obj->Next();
		Param->Index ++;
	}

	if (OutDim == RDim_Sample_X_SNP)
	{
		for (int sn = SnpCount; sn > 0; sn--)
		{
			Obj->ReadDosage(OutBuf);
			Obj->Next();
			OutBuf += Obj->SampNum;
			Param->Index ++;
		}
	} else {
		size_t size = SnpCount;
		for (size_t sn = size; sn > 0; sn--)
		{
			Obj->ReadDosage(Param->GenoBuffer);
			Obj->Next();
			Param->Index ++;

			C_UInt8 *g = (OutBuf ++);
			C_UInt8 *p = Param->GenoBuffer;
			for (size_t n=Obj->SampNum; n > 0; n--)
			{
				*g = *p ++;
				g += size;
			}
		}
	}
}


static void SNPRelate_SampleRead(C_Int32 SampStart, C_Int32 SampCount,
	C_UInt8 *OutBuf, TTypeGenoDim OutDim, TParam *Param)
{
	if (dynamic_cast<CApply_Variant_Dosage*>(Param->Object))
	{
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(Param->SeqGDSFile);
		if (GDS_Node_Path(Root, "genotype/~data", FALSE))
		{
			delete Param->Object;
			Param->Object = NULL;
		}
	}

	if (Param->Object == NULL)
	{
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(Param->SeqGDSFile);
		if (GDS_Node_Path(Root, "genotype/~data", FALSE))
		{
			CVarApplyBySample *Obj = new CVarApplyBySample;
			Param->Object = Obj;

			CFileInfo &File = GetFileInfo(Param->SeqGDSFile);
			TSelection &Sel = File.Selection();
			Obj->InitObject(CVariable::ctGenotype,
				"genotype/data", Root, File.VariantNum(),
				Sel.pVariant, File.SampleNum(), Sel.pSample, false);

			size_t SIZE = (Obj->Num_Variant) * (Obj->DLen[2]);
			Param->GenoBuffer = new C_UInt8[SIZE];
		} else {
			CApply_Variant_Dosage *Obj = new CApply_Variant_Dosage(
				GetFileInfo(Param->SeqGDSFile), true, false);
			Param->Object = Obj;
			size_t SIZE = (Obj->SampNum) * (Obj->Ploidy);
			Param->GenoBuffer = new C_UInt8[SIZE];
		}
		Param->Index = 0;
	}	

	// reading
	if (dynamic_cast<CApply_Variant_Dosage*>(Param->Object))
	{
		CApply_Variant_Dosage *Obj = (CApply_Variant_Dosage*)(Param->Object);
		Obj->Reset();

		if (OutDim == RDim_Sample_X_SNP)
		{
			do {
				Obj->ReadGenoData(Param->GenoBuffer);
				C_UInt8 *p = Param->GenoBuffer + (SampStart * Obj->Ploidy);

				for (size_t n=SampCount; n > 0; n--)
				{
					C_UInt8 *pp = p;
					p += Obj->Ploidy;
					C_UInt8 val = 0;
					for (size_t m=Obj->Ploidy; m > 0; m--, pp++)
					{
						if (*pp == 0)
						{
							val ++;
							if (val > 2) val = 2;
						} else if (*pp == NA_RAW)
						{
							val = 3; break;
						}
					}
					*OutBuf ++ = val;
				}
			} while (Obj->Next());

		} else {
			const size_t SNPNum = *Param->pSNPNum;
			do {
				Obj->ReadGenoData(Param->GenoBuffer);
				C_UInt8 *p = Param->GenoBuffer + (SampStart * Obj->Ploidy);
				C_UInt8 *g = (OutBuf ++);

				for (size_t n=SampCount; n > 0; n--)
				{
					C_UInt8 *pp = p;
					p += Obj->Ploidy;
					C_UInt8 val = 0;
					for (size_t m=Obj->Ploidy; m > 0; m--, pp++)
					{
						if (*pp == 0)
						{
							val ++;
							if (val > 2) val = 2;
						} else if (*pp == NA_RAW)
						{
							val = 3; break;
						}
					}
					*g = val; g += SNPNum;
				}
			} while (Obj->Next());
		}

	} else {
		CVarApplyBySample *Obj = (CVarApplyBySample*)(Param->Object);

		// indexing
		if (Param->Index > SampStart)
		{
			Obj->Reset();
			Param->Index = 0;
		}
		while (Param->Index < SampStart)
		{
			Obj->Next();
			Param->Index ++;
		}

		for (int sn = SampCount; sn > 0; sn--)
		{
			Obj->ReadGenoData(Param->GenoBuffer);
			Obj->Next();
			Param->Index ++;

			if (OutDim == RDim_SNP_X_Sample)
			{
				C_UInt8 *p = Param->GenoBuffer;
				for (size_t n=Obj->Num_Variant; n > 0; n--)
				{
					C_UInt8 *pp = p;
					p += Obj->DLen[2];
					C_UInt8 val = 0;
					for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
					{
						if (*pp == 0)
						{
							val ++;
							if (val > 2) val = 2;
						} else if (*pp == NA_RAW)
						{
							val = 3; break;
						}
					}
					*OutBuf ++ = val;
				}

			} else {
				C_UInt8 *g = (OutBuf ++);
				C_UInt8 *p = Param->GenoBuffer;
				for (size_t n=Obj->Num_Variant; n > 0; n--)
				{
					C_UInt8 *pp = p;
					p += Obj->DLen[2];
					C_UInt8 val = 0;
					for (size_t m=Obj->DLen[2]; m > 0; m--, pp++)
					{
						if (*pp == 0)
						{
							val ++;
							if (val > 2) val = 2;
						} else if (*pp == NA_RAW)
						{
							val = 3; break;
						}
					}
					*g = val; g += SampCount;
				}
			}
		}
	}
}


static void SNPRelate_SetSnpSelection(C_BOOL *sel, TParam *Param)
{
	CFileInfo &File = GetFileInfo(Param->SeqGDSFile);
	TSelection &s = File.Selection();
	s.ClearStructVariant();
	C_BOOL *p = s.pVariant;

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
	CFileInfo &File = GetFileInfo(Param->SeqGDSFile);
	TSelection &s = File.Selection();
	s.ClearStructSample();
	C_BOOL *p = s.pSample;

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
