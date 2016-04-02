// ===========================================================
//
// SeqArray.cpp: the C/C++ codes for the SeqArray package
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

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <set>
#include <algorithm>

#include "ReadByVariant.h"
#include "ReadBySample.h"



// ===========================================================
// The Initialized Object
// ===========================================================

TInitObject::TInitObject(): GENO_BUFFER(1024)
{
	memset(TRUE_ARRAY, TRUE, sizeof(TRUE_ARRAY));
}

void TInitObject::Need_GenoBuffer(size_t size)
{
	if (size > GENO_BUFFER.size())
		GENO_BUFFER.resize(size);
}

TInitObject Init;



// ===========================================================
// Library Functions
// ===========================================================

/// Get the total count requiring the number of dimension is one
COREARRAY_DLL_LOCAL int GetGDSObjCount(PdAbstractArray Obj, const char *varname)
{
	if (GDS_Array_DimCnt(Obj) != 1)
		throw SeqArray::ErrSeqArray("Invalid dimension of '%s'!", varname);
	return GDS_Array_GetTotalCount(Obj);
}

/// Get the number of alleles
COREARRAY_DLL_LOCAL int GetNumOfAllele(const char *allele_list)
{
	int n = 0;
	while (*allele_list)
	{
		if (*allele_list != ',')
		{
			n ++;
			while ((*allele_list != ',') && (*allele_list != 0))
				allele_list ++;
			if (*allele_list == ',')
			{
				allele_list ++;
				if (*allele_list == 0)
				{
					n ++;
					break;
				}
			}
		}
	}
	return n;
}

/// Get the index in an allele list
COREARRAY_DLL_LOCAL int GetIndexOfAllele(const char *allele, const char *allele_list)
{
	const size_t len = strlen(allele);
	const char *st = allele_list;
	int idx = 0;
	while (*allele_list)
	{
		while ((*allele_list != ',') && (*allele_list != 0))
			allele_list ++;
		size_t n = allele_list - st;
		if ((len==n) && (strncmp(allele, st, n)==0))
			return idx;
		if (*allele_list == ',')
		{
			idx ++;
			allele_list ++;
			st = allele_list;
		}
	}
	return -1;
}

/// Get strings split by comma
COREARRAY_DLL_LOCAL void GetAlleles(const char *alleles, vector<string> &out)
{
	out.clear();
	const char *p, *s;
	p = s = alleles;
	do {
		if ((*p == 0) || (*p == ','))
		{
			out.push_back(string(s, p));
			if (*p == ',') p ++;
			s = p;
			if (*p == 0) break;
		}
		p ++;
	} while (1);
}



extern "C"
{
using namespace SeqArray;

// ===========================================================
// Open a GDS file
// ===========================================================

/// initialize a SeqArray file
COREARRAY_DLL_EXPORT SEXP SEQ_File_Init(SEXP gdsfile)
{
	COREARRAY_TRY
		GetFileInfo(gdsfile);
	COREARRAY_CATCH
}

/// finalize a SeqArray file
COREARRAY_DLL_EXPORT SEXP SEQ_File_Done(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
			GDSFile_ID_Info.erase(it);
	COREARRAY_CATCH
}



// ===========================================================
// Set a working space
// ===========================================================

/// push the current filter to the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushEmpty(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it =
			GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			it->second.SelList.push_back(TSelection());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// push the current filter to the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushLast(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it =
			GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			if (!it->second.SelList.empty())
				it->second.SelList.push_back(it->second.SelList.back());
			else
				it->second.SelList.push_back(TSelection());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// pop up the previous filter from the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPop(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it =
			GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			if (it->second.SelList.size() <= 1)
				throw ErrSeqArray("No filter can be pop up.");
			it->second.SelList.pop_back();
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// set a working space with selected sample id
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceSample(SEXP gdsfile, SEXP samp_id,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pSample();
		int Count = File.SampleNum();
		PdAbstractArray varSamp = File.GetObj("sample.id", TRUE);

		if (Rf_isInteger(samp_id))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(samp_id), INTEGER(samp_id) + XLENGTH(samp_id));
			// sample id
			vector<int> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varSamp, &_st, &_cnt, &sample_id[0], svInt32);

			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(sample_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(sample_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isReal(samp_id))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(samp_id), REAL(samp_id) + XLENGTH(samp_id));
			// sample id
			vector<double> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varSamp, &_st, &_cnt, &sample_id[0], svFloat64);

			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(sample_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(sample_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isString(samp_id))
		{
			// initialize
			set<string> set_id;
			R_xlen_t m = XLENGTH(samp_id);
			for (R_xlen_t i=0; i < m; i++)
				set_id.insert(string(CHAR(STRING_ELT(samp_id, i))));
			// sample id
			vector<string> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varSamp, &_st, &_cnt, &sample_id[0], svStrUTF8);

			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(sample_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(sample_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isNull(samp_id))
		{
			memset(pArray, TRUE, Count);
		} else
			throw ErrSeqArray("Invalid type of 'sample.id'.");

		int n = File.SampleSelNum();
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected samples: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


/// set a working space with selected sample id (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceSample2(SEXP gdsfile, SEXP samp_sel,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pSample();
		int Count = File.SampleNum();

		if (Rf_isLogical(samp_sel) || IS_RAW(samp_sel))
		{
			// a logical vector for selected samples
			if (!intersect_flag)
			{
				if (XLENGTH(samp_sel) != Count)
					throw ErrSeqArray("Invalid length of 'sample.sel'.");
				// set selection
				if (Rf_isLogical(samp_sel))
				{
					int *base = LOGICAL(samp_sel);
					for (int i=0; i < Count; i++)
						*pArray++ = ((*base++) == TRUE);
				} else {
					Rbyte *base = RAW(samp_sel);
					for (int i=0; i < Count; i++)
						*pArray++ = ((*base++) != 0);
				}
			} else {
				if (XLENGTH(samp_sel) != File.SampleSelNum())
				{
					throw ErrSeqArray(
						"Invalid length of 'sample.sel' "
						"(should be equal to the number of selected samples).");
				}
				// set selection
				if (Rf_isLogical(samp_sel))
				{
					int *base = LOGICAL(samp_sel);
					for (int i=0; i < Count; i++, pArray++)
					{
						if (*pArray)
							*pArray = ((*base++) == TRUE);
					}
				} else {
					Rbyte *base = RAW(samp_sel);
					for (int i=0; i < Count; i++)
					{
						if (*pArray)
							*pArray = ((*base++) != 0);
					}
				}
			}
		} else if (Rf_isInteger(samp_sel) || Rf_isReal(samp_sel))
		{
			if (Rf_isReal(samp_sel))
				samp_sel = AS_INTEGER(samp_sel);

			if (!intersect_flag)
			{
				int *pI = INTEGER(samp_sel);
				R_xlen_t N = XLENGTH(samp_sel);
				// check
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if ((I != NA_INTEGER) && ((I < 1) || (I > Count)))
						throw ErrSeqArray("Out of range 'sample.sel'.");
				}
				// set values
				memset((void*)pArray, 0, Count);
				pI = INTEGER(samp_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I != NA_INTEGER) pArray[I-1] = TRUE;
				}
			} else {
				int Cnt = File.SampleSelNum();
				int *pI = INTEGER(samp_sel);
				R_xlen_t N = XLENGTH(samp_sel);
				// check
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if ((I != NA_INTEGER) && ((I < 1) || (I > Cnt)))
						throw ErrSeqArray("Out of range 'sample.sel'.");
				}
				// get the current index
				vector<int> Idx;
				Idx.reserve(Cnt);
				for (int i=0; i < Count; i++)
				{
					if (pArray[i]) Idx.push_back(i);
				}
				// set values
				memset((void*)pArray, 0, Count);
				pI = INTEGER(samp_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I != NA_INTEGER) pArray[Idx[I-1]] = TRUE;
				}
			}
		} else if (Rf_isNull(samp_sel))
		{
			memset(pArray, TRUE, Count);
		} else
			throw ErrSeqArray("Invalid type of 'sample.sel'.");

		int n = File.SampleSelNum();
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected samples: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant(SEXP gdsfile, SEXP var_id,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pVariant();
		int Count = File.VariantNum();
		PdAbstractArray varVariant = File.GetObj("variant.id", TRUE);

		if (Rf_isInteger(var_id))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(var_id), INTEGER(var_id) + XLENGTH(var_id));
			// sample id
			vector<int> var_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varVariant, &_st, &_cnt, &var_id[0], svInt32);

			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(var_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(var_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isReal(var_id))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(var_id), REAL(var_id) + XLENGTH(var_id));
			// variant id
			vector<double> var_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varVariant, &_st, &_cnt, &var_id[0], svFloat64);

			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(var_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(var_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isString(var_id))
		{
			// initialize
			set<string> set_id;
			R_xlen_t m = XLENGTH(var_id);
			for (R_xlen_t i=0; i < m; i++)
				set_id.insert(string(CHAR(STRING_ELT(var_id, i))));
			// sample id
			vector<string> var_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varVariant, &_st, &_cnt, &var_id[0], svStrUTF8);

			// set selection
			// set selection
			if (!intersect_flag)
			{
				for (int i=0; i < Count; i++)
					*pArray++ = (set_id.find(var_id[i]) != set_id.end());
			} else {
				for (int i=0; i < Count; i++, pArray++)
				{
					if (*pArray)
						*pArray = (set_id.find(var_id[i]) != set_id.end());
				}
			}
		} else if (Rf_isNull(var_id))
		{
			memset(pArray, TRUE, Count);
		} else
			throw ErrSeqArray("Invalid type of 'variant.id'.");

		int n = File.VariantSelNum();
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected variants: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


/// set a working space with selected variant id (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant2(SEXP gdsfile, SEXP var_sel,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pVariant();
		int Count = File.VariantNum();

		if (Rf_isLogical(var_sel) || IS_RAW(var_sel))
		{
			// a logical vector for selected samples
			if (!intersect_flag)
			{
				if (XLENGTH(var_sel) != Count)
					throw ErrSeqArray("Invalid length of 'variant.sel'.");
				// set selection
				if (Rf_isLogical(var_sel))
				{
					int *base = LOGICAL(var_sel);
					for (int i=0; i < Count; i++)
						*pArray++ = ((*base++) == TRUE);
				} else {
					Rbyte *base = RAW(var_sel);
					for (int i=0; i < Count; i++)
						*pArray++ = ((*base++) != 0);
				}
			} else {
				if (XLENGTH(var_sel) != File.VariantSelNum())
				{
					throw ErrSeqArray(
						"Invalid length of 'variant.sel' "
						"(should be equal to the number of selected variants).");
				}
				// set selection
				if (Rf_isLogical(var_sel))
				{
					int *base = LOGICAL(var_sel);
					for (int i=0; i < Count; i++, pArray++)
					{
						if (*pArray)
							*pArray = ((*base++) == TRUE);
					}
				} else {
					Rbyte *base = RAW(var_sel);
					for (int i=0; i < Count; i++)
					{
						if (*pArray)
							*pArray = ((*base++) != 0);
					}
				}
			}
		} else if (Rf_isInteger(var_sel) || Rf_isReal(var_sel))
		{
			if (Rf_isReal(var_sel))
				var_sel = AS_INTEGER(var_sel);

			if (!intersect_flag)
			{
				int *pI = INTEGER(var_sel);
				R_xlen_t N = XLENGTH(var_sel);
				// check
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if ((I != NA_INTEGER) && ((I < 1) || (I > Count)))
						throw ErrSeqArray("Out of range 'variant.sel'.");
				}
				// set values
				memset((void*)pArray, 0, Count);
				pI = INTEGER(var_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I != NA_INTEGER) pArray[I-1] = TRUE;
				}
			} else {
				int Cnt = File.VariantSelNum();
				int *pI = INTEGER(var_sel);
				R_xlen_t N = XLENGTH(var_sel);
				// check
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if ((I != NA_INTEGER) && ((I < 1) || (I > Cnt)))
						throw ErrSeqArray("Out of range 'variant.sel'.");
				}
				// get the current index
				vector<int> Idx;
				Idx.reserve(Cnt);
				for (int i=0; i < Count; i++)
				{
					if (pArray[i]) Idx.push_back(i);
				}
				// set values
				memset((void*)pArray, 0, Count);
				pI = INTEGER(var_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I != NA_INTEGER) pArray[Idx[I-1]] = TRUE;
				}
			}
		} else if (Rf_isNull(var_sel))
		{
			memset(pArray, TRUE, Count);
		} else
			throw ErrSeqArray("Invalid type of 'variant.sel'.");

		int n = File.VariantSelNum();
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected variants: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


static bool is_numeric(const string &txt)
{
	char *endptr = (char*)(txt.c_str());
	strtol(txt.c_str(), &endptr, 10);
	return (endptr != txt.c_str()) && (*endptr == 0);
}

/// set a working space flag with selected chromosome(s)
COREARRAY_DLL_EXPORT SEXP SEQ_SetChrom(SEXP gdsfile, SEXP include,
	SEXP is_num, SEXP frombp, SEXP tobp, SEXP verbose)
{
	int nProtected = 0;
	int *pFrom=NULL, *pTo=NULL;

	int IsNum = Rf_asLogical(is_num);

	if (Rf_isNull(include))
	{
		if (!Rf_isNull(frombp))
			error("'from.bp' should be NULL.");
		if (!Rf_isNull(tobp))
			error("'to.bp' should be NULL.");
	} else {
		include = PROTECT(AS_CHARACTER(include));
		nProtected ++;
		if (!Rf_isNull(frombp) || !Rf_isNull(tobp))
		{
			if (RLength(include) != RLength(frombp))
				error("'from.bp' should have the same length as 'include'.");
			if (RLength(include) != RLength(tobp))
				error("'to.bp' should have the same length as 'include'.");
			frombp = PROTECT(AS_INTEGER(frombp));
			tobp = PROTECT(AS_INTEGER(tobp));
			pFrom = INTEGER(frombp); pTo = INTEGER(tobp);
			nProtected += 2;
		}
	}

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();

		vector<C_BOOL> &array = Sel.Variant;
		memset(&array[0], FALSE, array.size());

		if (Rf_isNull(include))
		{
			// include = NULL
			if (IsNum == NA_INTEGER)
			{
				memset(&array[0], TRUE, array.size());
			} else {
				CChromIndex &Chrom = File.Chromosome();
				map<string, CChromIndex::TRangeList>::iterator it;
				for (it=Chrom.Map.begin(); it != Chrom.Map.end(); it++)
				{
					bool flag = is_numeric(it->first);
					if (((IsNum==TRUE) && flag) || ((IsNum==FALSE) && !flag))
					{
						CChromIndex::TRangeList &rng = it->second;
						vector<CChromIndex::TRange>::iterator it;
						for (it=rng.begin(); it != rng.end(); it++)
						{
							memset(&array[it->Start], TRUE, it->Length);
						}
					}
				}
			}

		} else {
			// include != NULL
			vector<C_Int32> *varPos = NULL;
			if (pFrom && pTo)
				varPos = &File.Position();

			CChromIndex &Chrom = File.Chromosome();
			map<string, CRangeSet> RngSets;

			R_xlen_t n = XLENGTH(include);
			for (R_xlen_t idx=0; idx < n; idx++)
			{
				string s = CHAR(STRING_ELT(include, idx));

				if (IsNum == TRUE)
				{
					if (!is_numeric(s)) continue;
				} else if (IsNum == FALSE)
				{
					if (is_numeric(s)) continue;
				}

				map<string, CChromIndex::TRangeList>::iterator it =
					Chrom.Map.find(s);
				if (it != Chrom.Map.end())
				{
					if (varPos)
					{
						// if from.bp and to.bp
						int from = pFrom[idx], to = pTo[idx];
						if (from == NA_INTEGER) from = 0;
						if (to == NA_INTEGER) to = 2147483647;
						RngSets[s].AddRange(from, to);
					} else {
						// no from.bp and to.bp
						CChromIndex::TRangeList &rng = it->second;
						vector<CChromIndex::TRange>::iterator p;
						for (p=rng.begin(); p != rng.end(); p++)
						{
							memset(&array[p->Start], TRUE, p->Length);
						}
					}
				}
			}

			if (varPos)
			{
				map<string, CRangeSet>::iterator it;
				for (it=RngSets.begin(); it != RngSets.end(); it++)
				{
					CChromIndex::TRangeList &rng = Chrom.Map[it->first];
					CRangeSet &RngSet = it->second;
					vector<CChromIndex::TRange>::const_iterator p;
					for (p=rng.begin(); p != rng.end(); p++)
					{
						size_t i = p->Start;
						size_t n = p->Length;
						C_Int32 *s = &((*varPos)[0]) + i;
						for (; n > 0; n--)
						{
							if (RngSet.IsIncluded(*s++))
								array[i] = TRUE;
							i ++;
						}
					}
				}
			}
		}

		if (Rf_asLogical(verbose) == TRUE)
		{
			int n = GetNumOfTRUE(&array[0], array.size());
			Rprintf("# of selected variants: %s\n", PrettyInt(n));
		}

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// set a working space flag with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_GetSpace(SEXP gdsfile, SEXP UseRaw)
{
	int use_raw_flag = Rf_asLogical(UseRaw);
	if (use_raw_flag == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SEXP tmp;

		// sample selection
		size_t n = File.SampleNum();
		if (use_raw_flag)
		{
			PROTECT(tmp = NEW_RAW(n));
			memcpy(RAW(tmp), Sel.pSample(), n);
		} else {
			PROTECT(tmp = NEW_LOGICAL(n));
			int *p = LOGICAL(tmp);
			C_BOOL *s = Sel.pSample();
			for (; n > 0; n--) *p++ = *s++;
		}
		SET_ELEMENT(rv_ans, 0, tmp);

		// variant selection
		n = File.VariantNum();
		if (use_raw_flag)
		{
			PROTECT(tmp = NEW_RAW(n));
			memcpy(RAW(tmp), Sel.pVariant(), n);
		} else {
			PROTECT(tmp = NEW_LOGICAL(n));
			int *p = LOGICAL(tmp);
			C_BOOL *s = Sel.pVariant();
			for (; n > 0; n--) *p++ = *s++;
		}
		SET_ELEMENT(rv_ans, 1, tmp);

		PROTECT(tmp = NEW_CHARACTER(2));
			SET_STRING_ELT(tmp, 0, mkChar("sample.sel"));
			SET_STRING_ELT(tmp, 1, mkChar("variant.sel"));
		SET_NAMES(rv_ans, tmp);

		UNPROTECT(4);

	COREARRAY_CATCH
}


// ===========================================================

inline static C_BOOL *CLEAR_SELECTION(size_t num, C_BOOL *p)
{
	while (num > 0)
	{
		if (*p != FALSE) { num--; *p = FALSE; }
		p ++;
	}
	return p;
}
inline static C_BOOL *SKIP_SELECTION(size_t num, C_BOOL *p)
{
	while (num > 0)
	{
		if (*p != FALSE) num--;
		p ++;
	}
	return p;
}

/// split the selected variants according to multiple processes
COREARRAY_DLL_EXPORT SEXP SEQ_SplitSelection(SEXP gdsfile, SEXP split,
	SEXP index, SEXP n_process, SEXP selection_flag)
{
	const char *split_str = CHAR(STRING_ELT(split, 0));
	int Process_Index = Rf_asInteger(index) - 1;  // starting from 0
	int Num_Process = Rf_asInteger(n_process);
	int SelFlag = Rf_asLogical(selection_flag);

	COREARRAY_TRY

		// selection object
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &s = File.Selection();

		// the total number of selected elements
		int SelectCount;
		C_BOOL *sel;
		if (strcmp(split_str, "by.variant") == 0)
		{
			if (s.Variant.empty())
			{
				s.Variant.resize(
					GDS_Array_GetTotalCount(GDS_Node_Path(
					GDS_R_SEXP2FileRoot(gdsfile), "variant.id", TRUE)), TRUE);
			}
			sel = &s.Variant[0];
			SelectCount = GetNumOfTRUE(sel, s.Variant.size());
		} else if (strcmp(split_str, "by.sample") == 0)
		{
			if (s.Sample.empty())
			{
				s.Sample.resize(
					GDS_Array_GetTotalCount(GDS_Node_Path(
					GDS_R_SEXP2FileRoot(gdsfile), "sample.id", TRUE)), TRUE);
			}
			sel = &s.Sample[0];
			SelectCount = GetNumOfTRUE(sel, s.Sample.size());
		} else {
			return rv_ans;
		}

		// split a list
		vector<int> split(Num_Process);
		double avg = (double)SelectCount / Num_Process;
		double start = 0;
		for (int i=0; i < Num_Process; i++)
		{
			start += avg;
			split[i] = (int)(start + 0.5);
		}

		// ---------------------------------------------------
		int st = 0;
		for (int i=0; i < Process_Index; i++)
		{
			sel = CLEAR_SELECTION(split[i] - st, sel);
			st = split[i];
		}
		int ans_n = split[Process_Index] - st;
		sel = SKIP_SELECTION(ans_n, sel);
		st = split[Process_Index];
		for (int i=Process_Index+1; i < Num_Process; i++)
		{
			sel = CLEAR_SELECTION(split[i] - st, sel);
			st = split[i];
		}

		// ---------------------------------------------------
		// output
		if (SelFlag == TRUE)
		{
			rv_ans = NEW_LOGICAL(SelectCount);
			int *p = INTEGER(rv_ans);
			memset((void*)p, 0, sizeof(int) * size_t(SelectCount));
			if (Process_Index > 0)
				p += split[Process_Index-1];
			for (; ans_n > 0; ans_n--) *p++ = TRUE;
		} else {
			rv_ans = ScalarInteger(ans_n);
		}

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_Summary(SEXP gdsfile, SEXP varname)
{
	COREARRAY_TRY

		// the selection
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		// the variable name
		string vn = CHAR(STRING_ELT(varname, 0));

		if ((vn=="genotype") || (vn=="phase"))
		{
			PdGDSObj vGeno = GDS_Node_Path(Root, "genotype/data", TRUE);
			if (vGeno == NULL)
			{
				vGeno = GDS_Node_Path(Root, "genotype/~data", FALSE);
				if (vGeno == NULL)
				{
					throw ErrSeqArray(
						"There is no 'genotype/data' or 'genotype/~data'.");
				}
			}

			PROTECT(rv_ans = NEW_LIST(2));

				SEXP I32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 0, I32);
				C_Int32 Buf[4];
				GDS_Array_GetDim(vGeno, Buf, 3);
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = Sel.Sample.size();
				INTEGER(I32)[2] = Sel.Variant.size();

				SEXP S32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 1, S32);
				INTEGER(S32)[0] = Buf[2];
				INTEGER(S32)[1] = GetNumOfTRUE(&Sel.Sample[0], Sel.Sample.size());
				INTEGER(S32)[2] = GetNumOfTRUE(&Sel.Variant[0], Sel.Variant.size());

			SEXP tmp = PROTECT(NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("dim"));
				SET_STRING_ELT(tmp, 1, mkChar("seldim"));
				SET_NAMES(rv_ans, tmp);
			UNPROTECT(4);

		} else {
			PdGDSObj var = GDS_Node_Path(Root, vn.c_str(), TRUE);
			rv_ans = ScalarInteger(GDS_Array_GetTotalCount(var));
		}

	COREARRAY_CATCH
}


/// get a logical vector with selection
COREARRAY_DLL_EXPORT SEXP SEQ_SelectFlag(SEXP select, SEXP len)
{
	R_len_t n = XLENGTH(select);
	if (XLENGTH(len) != n)
		error("Index variable error.");

	int *p = INTEGER(len);
	R_len_t m = 0;
	for (R_len_t k=n; k > 0; k--, p++)
	{
		if (*p > 0) m += *p;
	}

	SEXP rv_ans = NEW_LOGICAL(m);
	int *r = INTEGER(rv_ans), *s = INTEGER(select);
	p = INTEGER(len);
	for (; n > 0; n--, s++, p++)
	{
		for (int k=*p; k > 0; k--)
			*r++ = *s;
	}

	return rv_ans;
}



// ===========================================================
// get system configuration
// ===========================================================

/// the number of alleles per site
COREARRAY_DLL_EXPORT SEXP SEQ_System()
{
	COREARRAY_TRY

		int nProtect = 0;
		rv_ans = PROTECT(NEW_LIST(2));
		SEXP nm = PROTECT(NEW_CHARACTER(2));
		nProtect += 2;
		SET_NAMES(rv_ans, nm);

		// the number of logical cores
		SET_ELEMENT(rv_ans, 0, ScalarInteger(GDS_Mach_GetNumOfCores()));
		SET_STRING_ELT(nm, 0, mkChar("num.logical.core"));

		// compiler flags
		vector<string> ss;

	#ifdef __SSE__
		ss.push_back("SSE");
	#endif
	#ifdef __SSE2__
		ss.push_back("SSE2");
	#endif
	#ifdef __SSE3__
		ss.push_back("SSE3");
	#endif
	#ifdef __SSSE3__
		ss.push_back("SSSE3");
	#endif
	#ifdef __SSE4_1__
		ss.push_back("SSE4.1");
	#endif
	#ifdef __SSE4_2__
		ss.push_back("SSE4.2");
	#endif
	#ifdef __AVX__
		ss.push_back("AVX");
	#endif
	#ifdef __AVX2__
		ss.push_back("AVX2");
	#endif
	#ifdef __FMA__
		ss.push_back("FMA");
	#endif
	#ifdef __FMA4__
		ss.push_back("FMA4");
	#endif
		SEXP SIMD = PROTECT(NEW_CHARACTER(ss.size()));
		nProtect ++;
		SET_ELEMENT(rv_ans, 1, SIMD);
		SET_STRING_ELT(nm, 1, mkChar("compiler.flag"));
		for (int i=0; i < (int)ss.size(); i++)
			SET_STRING_ELT(SIMD, i, mkChar(ss[i].c_str()));

		UNPROTECT(nProtect);

	COREARRAY_CATCH
}


// ===========================================================
// the initial function when the package is loaded
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName0()
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName1(SEXP x)
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName2(SEXP x, SEXP y)
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName3(SEXP x, SEXP y, SEXP z)
{
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName4(SEXP w, SEXP x, SEXP y, SEXP z)
{
	return R_NilValue;
}


/// initialize the package
COREARRAY_DLL_EXPORT void R_init_SeqArray(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	extern void Register_SNPRelate_Functions();

	extern SEXP SEQ_GetData(SEXP, SEXP);
	extern SEXP SEQ_ConvBEDFlag(SEXP, SEXP, SEXP);
	extern SEXP SEQ_ConvBED2GDS(SEXP, SEXP, SEXP, SEXP, SEXP);

	extern SEXP SEQ_MergeAllele(SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeGeno(SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergePhase(SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeInfo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeFormat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_ExternalName0, 0),         CALL(SEQ_ExternalName1, 1),
		CALL(SEQ_ExternalName2, 2),         CALL(SEQ_ExternalName3, 3),
		CALL(SEQ_ExternalName4, 4),

		CALL(SEQ_File_Init, 1),             CALL(SEQ_File_Done, 1),
		CALL(SEQ_FilterPushEmpty, 1),       CALL(SEQ_FilterPushLast, 1),
		CALL(SEQ_FilterPop, 1),

		CALL(SEQ_MergeAllele, 4),           CALL(SEQ_MergeGeno, 5),
		CALL(SEQ_MergePhase, 5),            CALL(SEQ_MergeInfo, 6),
		CALL(SEQ_MergeFormat, 6),

		CALL(SEQ_SetSpaceSample, 4),        CALL(SEQ_SetSpaceSample2, 4),
		CALL(SEQ_SetSpaceVariant, 4),       CALL(SEQ_SetSpaceVariant2, 4),
		CALL(SEQ_SplitSelection, 5),        CALL(SEQ_SetChrom, 6),
		CALL(SEQ_GetSpace, 2),

		CALL(SEQ_Summary, 2),               CALL(SEQ_System, 0),

		CALL(SEQ_GetData, 3),
		CALL(SEQ_Apply_Sample, 7),          CALL(SEQ_Apply_Variant, 7),

		CALL(SEQ_ConvBEDFlag, 3),           CALL(SEQ_ConvBED2GDS, 5),
		CALL(SEQ_SelectFlag, 2),

		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Register_SNPRelate_Functions();
	Init_GDS_Routines();
}

} // extern "C"
