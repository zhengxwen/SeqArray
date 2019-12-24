// ===========================================================
//
// FileMerge.cpp: GDS file merging
//
// Copyright (C) 2016-2019    Xiuwen Zheng
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

#include <set>
#include <algorithm>

#include "ReadByVariant.h"
#include "ReadBySample.h"



// ===========================================================
// File Merging
// ===========================================================

extern "C"
{
using namespace SeqArray;


#define MERGE_VAR_DEF    \
	const int TotalNum = Rf_asInteger(num); \
	const int FileCnt  = Rf_length(varidx); \
	vector<int*> pIdx(FileCnt); \
	vector<C_Int32> pI(FileCnt); \
	for (int i=0; i < FileCnt; i++) \
	{ \
		pIdx[i] = INTEGER(VECTOR_ELT(varidx, i)); \
		pI[i] = 0; \
	}


static const C_Int32 ZERO = 0;
static const C_Int32 ONE  = 1;


/// merge alleles from multiple files
COREARRAY_DLL_EXPORT SEXP SEQ_MergeAllele(SEXP num, SEXP varidx, SEXP files,
	SEXP export_var)
{
	COREARRAY_TRY

		MERGE_VAR_DEF

		vector<PdAbstractArray> pVar(FileCnt);
		for (int i=0; i < FileCnt; i++)
		{
			PdGDSFolder Root = GDS_R_SEXP2FileRoot(VECTOR_ELT(files, i));
			pVar[i] = GDS_Node_Path(Root, "allele", TRUE);
		}

		PdAbstractArray exp_var = GDS_R_SEXP2Obj(export_var, FALSE);

		// for-loop
		vector<string> vec, vv;
		string ss, val;
		for (int i=1; i <= TotalNum; i++)
		{
			vec.clear();
			for (int j=0; j < FileCnt; j++)
			{
				if (*pIdx[j] == i)  // deal with this variant?
				{
					++ pIdx[j];
					GDS_Array_ReadData(pVar[j], &pI[j], &ONE, &val, svStrUTF8);
					++ pI[j];
					// parse alleles
					GetAlleles(val.c_str(), vv);
					for (int k=0; k < (int)vv.size(); k++)
					{
						vector<string>::iterator it = find(vec.begin(), vec.end(), vv[k]);
						if (it == vec.end())
							vec.push_back(vv[k]);
					}
				}
			}
			// save
			ss.clear();
			for (int j=0; j < (int)vec.size(); j++)
			{
				if (j > 0) ss.push_back(',');
				ss.append(vec[j]);
			}
			GDS_Array_AppendString(exp_var, ss.c_str());
		}

	COREARRAY_CATCH
}


/// merge genotypes from multiple files
COREARRAY_DLL_EXPORT SEXP SEQ_MergeGeno(SEXP num, SEXP varidx, SEXP files,
	SEXP export_file, SEXP param)
{
	COREARRAY_TRY

		MERGE_VAR_DEF

		vector<CApply_Variant_Geno> Files(FileCnt);
		for (int i=0; i < FileCnt; i++)
			Files[i].Init(GetFileInfo(VECTOR_ELT(files, i)), false);

		vector<PdAbstractArray> pAllele(FileCnt);
		for (int i=0; i < FileCnt; i++)
		{
			PdGDSFolder Root = GDS_R_SEXP2FileRoot(VECTOR_ELT(files, i));
			pAllele[i] = GDS_Node_Path(Root, "allele", TRUE);
		}

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(export_file);
		PdAbstractArray allele   = GDS_Node_Path(Root, "allele", TRUE);
		PdAbstractArray geno_var = GDS_Node_Path(Root, "genotype/data", TRUE);
		PdAbstractArray geno_idx = GDS_Node_Path(Root, "genotype/@data", TRUE);

		// the genotype buffer
		const int nsamp = INTEGER(num)[1];
		const int ploidy = INTEGER(num)[2];
		const int geno_cnt = nsamp * ploidy;

		vector<int> geno_buffer(geno_cnt);
		vector<C_Int8> I8s(geno_cnt);
		vector<string> ss;
		vector<int> allele_map;
		string allele_list, s;

		int div = TotalNum / 25;
		if (div <= 0) div = 1;
		bool Verbose = (Rf_asLogical(RGetListElement(param, "verbose")) == TRUE);

		// for-loop
		for (int i=1; i <= TotalNum; i++)
		{
			C_Int32 st = i - 1;
			GDS_Array_ReadData(allele, &st, &ONE, &allele_list, svStrUTF8);
			int *pGeno = &(geno_buffer[0]);

			for (int j=0; j < FileCnt; j++)
			{
				CApply_Variant_Geno &FILE = Files[j];
				const size_t size = (size_t)FILE.SampNum * ploidy;

				if (*pIdx[j] == i)  // deal with this variant?
				{
					++ pIdx[j];
					GDS_Array_ReadData(pAllele[j], &pI[j], &ONE, &s, svStrUTF8);
					++ pI[j];
					// parse alleles
					GetAlleles(s.c_str(), ss);
					const int nAllele = ss.size();
					allele_map.resize(nAllele);
					for (int k=0; k < nAllele; k++)
					{
						int x = GetIndexOfAllele(ss[k].c_str(), allele_list.c_str());
						if (x < 0)
							throw ErrSeqArray("internal error in SEQ_MergeGeno");
						allele_map[k] = x;
					}
					FILE.ReadGenoData(pGeno);
					FILE.Next();
					// replace
					size_t m = size;
					const int *map = &allele_map[0];
					for (int *p = pGeno; m > 0; m--, p++)
					{
						int v = *p;
						if ((0 <= v) && (v < nAllele))
							*p = map[v];
						else if (v != NA_INTEGER)
							warning("Genotype in File(%d), out of range.", j+1);
					}
				} else
					vec_int32_set(pGeno, size, NA_INTEGER);
				pGeno += size;
			}

			// determine how many bits
			const int GenoNumBits = 2;
			const int GenoBitMask = 0x03;
			int num_allele = GetNumOfAllele(allele_list.c_str());
			int num_bits = GenoNumBits;
			while ((num_allele + 1) > (1 << num_bits))
				num_bits += GenoNumBits;
			C_Int32 I32 = num_bits / GenoNumBits;
			GDS_Array_AppendData(geno_idx, 1, &I32, svInt32);

			// write to the variable "genotype"
			for (int bits=0; bits < num_bits; bits += GenoNumBits)
			{
				int *p = &(geno_buffer[0]);
				C_Int8 *s = &(I8s[0]);
				for (int k=0; k < geno_cnt; k++)
				{
					int v = *p++;
					*s++ = (v == NA_INTEGER) ? GenoBitMask :
						((v >> bits) & GenoBitMask);
				}
				GDS_Array_AppendData(geno_var, geno_cnt, &(I8s[0]), svInt8);
			}

			if (Verbose)
				if (i % div == 0) Rprintf("<");
		}

		if (Verbose) Rprintf("]");

	COREARRAY_CATCH
}


/// merge phasing status from multiple files
COREARRAY_DLL_EXPORT SEXP SEQ_MergePhase(SEXP num, SEXP varidx, SEXP files,
	SEXP export_file, SEXP param)
{
	COREARRAY_TRY

		MERGE_VAR_DEF

		int nProtected = 0;

		vector<CApply_Variant_Phase> Files(FileCnt);
		for (int i=0; i < FileCnt; i++)
			Files[i].Init(GetFileInfo(VECTOR_ELT(files, i)), false);

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(export_file);
		PdAbstractArray phase_var = GDS_Node_Path(Root, "phase/data", TRUE);

		// the phase buffer
		const int nsamp = INTEGER(num)[1];
		const int ploidy = INTEGER(num)[2];
		const int pcnt = nsamp * (ploidy - 1);

		int div = TotalNum / 25;
		if (div <= 0) div = 1;
		bool Verbose = (Rf_asLogical(RGetListElement(param, "verbose"))==TRUE);
		vector<int> phase_buf(pcnt);

		// for-loop
		for (int i=1; i <= TotalNum; i++)
		{
			int *pp = &(phase_buf[0]);

			for (int j=0; j < FileCnt; j++)
			{
				CApply_Variant_Phase &FILE = Files[j];
				const size_t size = (size_t)FILE.SampNum * (ploidy-1);

				if (*pIdx[j] == i)  // deal with this variant?
				{
					++ pIdx[j];
					SEXP RD = FILE.NeedRData(nProtected);
					FILE.ReadData(RD);
					FILE.Next();
					memcpy(pp, INTEGER(RD), sizeof(int)*size);
				} else
					vec_int32_set(pp, size, 0);
				pp += size;
			}

			// write to the variable "phase"
			GDS_Array_AppendData(phase_var, pcnt, &(phase_buf[0]), svInt32);

			if (Verbose)
				if (i % div == 0) Rprintf("<");
		}

		// finally
		if (Verbose) Rprintf("]");
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// merge INFO variables from multiple files
COREARRAY_DLL_EXPORT SEXP SEQ_MergeInfo(SEXP num, SEXP varidx, SEXP files,
	SEXP varname, SEXP export_file, SEXP param)
{
	COREARRAY_TRY

		MERGE_VAR_DEF

		int nProtected = 0;
		string VarName  = CHAR(STRING_ELT(varname, 0));
		string VarName2 = GDS_PATH_PREFIX(VarName, '@');

		CVarApplyList NodeList;
		for (int i=0; i < FileCnt; i++)
		{
			CFileInfo &File = GetFileInfo(VECTOR_ELT(files, i));
			if (VarName=="annotation/id" || VarName=="annotation/qual" ||
					VarName=="annotation/filter")
			{
				NodeList.push_back(
					new CApply_Variant_Basic(File, VarName.c_str()));
			} else {
				NodeList.push_back(
					new CApply_Variant_Info(File, VarName.c_str()));
			}
		}

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(export_file);
		PdAbstractArray info_var = GDS_Node_Path(Root, VarName.c_str(), TRUE);
		PdAbstractArray info_idx = GDS_Node_Path(Root, VarName2.c_str(), FALSE);

		// for-loop
		for (int i=1; i <= TotalNum; i++)
		{
			bool has = false;
			for (int j=0; j < FileCnt; j++)
			{
				CVarApply &FILE = *NodeList[j];
				if (*pIdx[j] == i)  // deal with this variant?
				{
					++ pIdx[j];
					SEXP RD = FILE.NeedRData(nProtected);
					FILE.ReadData(RD);
					FILE.Next();
					if (!Rf_isNull(RD))
						GDS_R_Append(info_var, RD);
					if (info_idx)
					{
						C_Int32 I32 = RLength(RD);
						GDS_Array_AppendData(info_idx, 1, &I32, svInt32);
					}
					has = true;
					break;
				}
			}
			if (!has)
			{
				if (info_idx)
					GDS_Array_AppendData(info_idx, 1, &ZERO, svInt32);
				else
					GDS_R_Append(info_var, ScalarInteger(NA_INTEGER));
			}
		}

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}



/// merge FORMAT variable from multiple files
COREARRAY_DLL_EXPORT SEXP SEQ_MergeFormat(SEXP num, SEXP varidx, SEXP files,
	SEXP varname, SEXP export_file, SEXP param)
{
	COREARRAY_TRY

		MERGE_VAR_DEF

		int nProtected = 0;
		string VarName  = CHAR(STRING_ELT(varname, 0));
		string VarName2 = GDS_PATH_PREFIX(VarName, '@');

		vector<CApply_Variant_Format> Files(FileCnt);
		for (int i=0; i < FileCnt; i++)
		{
			SEXP file = VECTOR_ELT(files, i);
			Files[i].Init(GetFileInfo(file), VarName.c_str());
		}

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(export_file);
		PdAbstractArray fmt_var = GDS_Node_Path(Root, VarName.c_str(), TRUE);
		PdAbstractArray fmt_idx = GDS_Node_Path(Root, VarName2.c_str(), TRUE);

		int div = TotalNum / 25;
		if (div <= 0) div = 1;
		SEXP NAs = RGetListElement(param, "na");
		bool Verbose = (Rf_asLogical(RGetListElement(param, "verbose"))==TRUE);
		vector<SEXP> RDList(FileCnt);

		// for-loop
		for (int i=1; i <= TotalNum; i++)
		{
			for (int j=0; j < FileCnt; j++)
			{
				SEXP RD = R_NilValue;
				if (*pIdx[j] == i)  // deal with this variant?
				{
					++ pIdx[j];
					CApply_Variant_Format &FILE = Files[j];
					RD = FILE.NeedRData(nProtected);
					FILE.ReadData(RD);
					FILE.Next();
				}
				RDList[j] = RD;
			}

			// 2-dim FORMAT variables
			int step = 0;
			for (int j=0; j < FileCnt; j++)
			{
				if (!Rf_isNull(RDList[j]))
				{
					size_t len = XLENGTH(RDList[j]);
					int m = len / Files[j].SampNum;
					if (m > step) step = m;
				}
			}

			// write to the variable
			for (int k=0; k < step; k++)
			{
				for (int j=0; j < FileCnt; j++)
				{
					CApply_Variant_Format &FILE = Files[j];
					if (!Rf_isNull(RDList[j]))
					{
						size_t len = XLENGTH(RDList[j]);
						int m = len / FILE.SampNum;
						if (k < m)
						{
							GDS_R_AppendEx(fmt_var, RDList[j],
								k * FILE.SampNum, FILE.SampNum);
						} else
							GDS_R_AppendEx(fmt_var, NAs, 0, FILE.SampNum);
					} else
						GDS_R_AppendEx(fmt_var, NAs, 0, FILE.SampNum);
				}
			}

			// write to index variable
			GDS_Array_AppendData(fmt_idx, 1, &step, svInt32);

			if (Verbose)
				if (i % div == 0) Rprintf("<");
		}

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


} // extern "C"
