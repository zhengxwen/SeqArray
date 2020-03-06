// ===========================================================
//
// SeqArray.cpp: the C++ codes for the SeqArray package
//
// Copyright (C) 2013-2020    Xiuwen Zheng
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

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <set>
#include <algorithm>

#include "ReadByVariant.h"
#include "ReadBySample.h"
#include <ctype.h>



// ===========================================================
// Library Functions
// ===========================================================

extern "C"
{

static const char *INFO_SEL_NUM_SAMPLE  = "# of selected samples: %s\n";
static const char *INFO_SEL_NUM_VARIANT = "# of selected variants: %s\n";


using namespace SeqArray;

// ===========================================================
// Open a GDS file
// ===========================================================

/// initialize a SeqArray file
COREARRAY_DLL_EXPORT SEXP SEQ_File_Init(SEXP gdsfile)
{
	COREARRAY_TRY
		CFileInfo &file = GetFileInfo(gdsfile);
		file.Selection();  // force to initialize selection
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

/// push the current filter to the stack and reset the filter
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushEmpty(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			CFileInfo &f = it->second;
			TSelection &s = f.Push_Selection(false, false);
			memset(s.pSample, TRUE, f.SampleNum());
			memset(s.pVariant, TRUE, f.VariantNum());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// push the current filter to the stack and not reset the filter
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushLast(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			CFileInfo &f = it->second;
			f.Push_Selection(true, true);
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// pop up the previous filter from the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPop(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(RGetListElement(gdsfile, "id"));
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(id);
		if (it != GDSFile_ID_Info.end())
		{
			CFileInfo &f = it->second;
			f.Pop_Selection();
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

		int nProtected = 0;
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		Sel.ClearStructSample();

		C_BOOL *pArray = Sel.pSample;
		int Count = File.SampleNum();
		PdAbstractArray varSamp = File.GetObj("sample.id", TRUE);
		C_SVType sv = GDS_Array_GetSVType(varSamp);
		if (COREARRAY_SV_STRING(sv) && !Rf_isNull(samp_id) && !Rf_isString(samp_id))
		{
			samp_id = PROTECT(AS_CHARACTER(samp_id));
			nProtected ++;
		}

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

		if (Rf_asLogical(verbose) == TRUE)
			Rprintf(INFO_SEL_NUM_SAMPLE, PrettyInt(File.SampleSelNum()));
		if (nProtected > 0) UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// set a working space with selected samples (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceSample2(SEXP gdsfile, SEXP samp_sel,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		Sel.ClearStructSample();

		C_BOOL *pArray = Sel.pSample;
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
				Sel.ClearStructSample();
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

		if (Rf_asLogical(verbose) == TRUE)
			Rprintf(INFO_SEL_NUM_SAMPLE, PrettyInt(File.SampleSelNum()));

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant(SEXP gdsfile, SEXP var_id,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		int nProtected = 0;
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		Sel.ClearStructVariant();

		C_BOOL *pArray = Sel.pVariant;
		int Count = File.VariantNum();
		PdAbstractArray varVariant = File.GetObj("variant.id", TRUE);
		C_SVType sv = GDS_Array_GetSVType(varVariant);
		if (COREARRAY_SV_STRING(sv) && !Rf_isNull(var_id) && !Rf_isString(var_id))
		{
			var_id = PROTECT(AS_CHARACTER(var_id));
			nProtected ++;
		}

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
			// variant id
			vector<string> var_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Array_ReadData(varVariant, &_st, &_cnt, &var_id[0], svStrUTF8);

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

		if (Rf_asLogical(verbose) == TRUE)
			Rprintf(INFO_SEL_NUM_VARIANT, PrettyInt(File.VariantSelNum()));
		if (nProtected > 0) UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// set a working space with selected variants (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant2(SEXP gdsfile, SEXP var_sel,
	SEXP intersect, SEXP verbose)
{
	static const char *ERR_OUT_RANGE = "Out of range 'variant.sel'.";
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();

		C_BOOL *pArray = Sel.pVariant;
		int Count = File.VariantNum();

		if (Rf_isLogical(var_sel) || IS_RAW(var_sel))
		{
			// a logical vector for selected samples
			if (!intersect_flag)
			{
				Sel.ClearStructVariant();
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
				Count = XLENGTH(var_sel);
				if (Count != File.VariantSelNum())
				{
					throw ErrSeqArray(
						"Invalid length of 'variant.sel' "
						"(should be equal to the number of selected variants).");
				}
				// call File.SampleSelNum() making varTrueNum, varStart and varEnd available
				// set selection
				pArray += Sel.varStart;
				ssize_t num = 0;
				if (Rf_isLogical(var_sel))
				{
					int *base = LOGICAL(var_sel);
					for (; Count > 0; pArray++)
					{
						if (*pArray)
						{
							Count--;
							if ((*base++) == TRUE) num++; else *pArray = FALSE;
						}
					}
				} else {
					Rbyte *base = RAW(var_sel);
					for (; Count > 0; pArray++)
					{
						if (*pArray)
						{
							Count--;
							if ((*base++) != 0) num++; else *pArray = FALSE;
						}
					}
				}
				// set TSelection structure
				if (num > 0)
				{
					C_BOOL *p = VEC_BOOL_FIND_TRUE(Sel.pVariant + Sel.varStart,
						Sel.pVariant + Sel.varEnd);
					Sel.varTrueNum = num;
					Sel.varStart = p - Sel.pVariant;
					for (; num > 0; ) if (*p++) num--;
					Sel.varEnd = p - Sel.pVariant;
				} else {
					Sel.varTrueNum = Sel.varStart = Sel.varEnd = 0;
				}
			}
		} else if (Rf_isInteger(var_sel) || Rf_isReal(var_sel))
		{
			if (Rf_isReal(var_sel))
				var_sel = AS_INTEGER(var_sel);
			if (!intersect_flag)
			{
				R_xlen_t N = XLENGTH(var_sel);
				// check
				if (!vec_i32_bound_check(INTEGER(var_sel), N, Count))
					throw ErrSeqArray(ERR_OUT_RANGE);
				// clear
				Sel.ClearSelectVariant();
				// set values
				ssize_t num=0, st=Count, ed=0;
				int *pI = INTEGER(var_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I!=NA_INTEGER && !pArray[I-1])
					{
						if (I > ed) ed = I;
						ssize_t ii = I - 1;
						if (ii < st) st = ii;
						pArray[ii] = TRUE;
						num ++;
					}
				}
				// set the structure of selected variants
				Sel.varTrueNum = num;
				Sel.varStart = st;
				Sel.varEnd = (ed < st) ? st : ed;
			} else {
				int Cnt = File.VariantSelNum();
				R_xlen_t N = XLENGTH(var_sel);
				// check
				if (!vec_i32_bound_check(INTEGER(var_sel), N, Cnt))
					throw ErrSeqArray(ERR_OUT_RANGE);
				// get the current index
				vector<int> Idx;
				Idx.reserve(Cnt);
				for (int i=0; i < Count; i++)
				{
					if (pArray[i]) Idx.push_back(i);
				}
				// set values
				memset((void*)pArray, 0, Count);
				int *pI = INTEGER(var_sel);
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if (I != NA_INTEGER) pArray[Idx[I-1]] = TRUE;
				}
				Sel.ClearStructVariant();
			}
		} else if (Rf_isNull(var_sel))
		{
			memset(pArray, TRUE, Count);
			Sel.varStart = 0;
			Sel.varEnd = Sel.varTrueNum = Count;
		} else
			throw ErrSeqArray("Invalid type of 'variant.sel'.");

		if (Rf_asLogical(verbose) == TRUE)
			Rprintf(INFO_SEL_NUM_VARIANT, PrettyInt(File.VariantSelNum()));

	COREARRAY_CATCH
}


// ================================================================

static bool is_numeric(const string &txt)
{
	char *endptr = (char*)(txt.c_str());
	strtol(txt.c_str(), &endptr, 10);
	return (endptr != txt.c_str()) && (*endptr == 0);
}

/// set a working space flag with selected chromosome(s)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceChrom(SEXP gdsfile, SEXP include,
	SEXP is_num, SEXP frombp, SEXP tobp, SEXP intersect, SEXP verbose)
{
	int nProtected = 0;
	int *pFrom=NULL, *pTo=NULL;

	int IsNum = Rf_asLogical(is_num);
	int IsIntersect = Rf_asLogical(intersect);
	if (IsIntersect == NA_INTEGER)
		error("'intersect' should be either FALSE or TRUE.");

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
		Sel.ClearStructVariant();

		const size_t array_size = File.VariantNum();
		C_BOOL *sel_array = Sel.pVariant;
		vector<C_BOOL> tmp_array;
		if (IsIntersect) tmp_array.resize(array_size);

		C_BOOL *array = IsIntersect ? &tmp_array[0] : sel_array;
		memset(array, FALSE, array_size);

		if (Rf_isNull(include))
		{
			// include = NULL
			if (IsNum == NA_INTEGER)
			{
				memset(array, TRUE, array_size);
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
			map<string, CRangeSet> RngSets;  // Chromosome ==> CRangeSet

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
						// if specify from.bp and to.bp
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
				// Chromosome ==> CRangeSet
				map<string, CRangeSet>::iterator it;
				for (it=RngSets.begin(); it != RngSets.end(); it++)
				{
					CChromIndex::TRangeList &rng = Chrom.Map[it->first];
					CRangeSet &RngSet = it->second;
					vector<CChromIndex::TRange>::const_iterator p;
					for (p=rng.begin(); p != rng.end(); p++)
					{
						size_t i=p->Start, n=p->Length;
						C_Int32 *s = &((*varPos)[0]) + i;
						if (RngSet.Size() == 1)
						{
							// there is only a range, optimized for this situation
							int st, ed;
							RngSet.GetRanges(&st, &ed);
							if (!IsIntersect)
							{
								for (; n > 0; n--, i++, s++)
									if (st<=*s && *s<=ed) array[i] = TRUE;
							} else {
								C_BOOL *b = &sel_array[i];
								for (; n > 0; n--, i++, s++)
									if (*b++ && st<=*s && *s<=ed) array[i] = TRUE;
							}
						} else {
							if (!IsIntersect)
							{
								for (; n > 0; n--, i++)
									if (RngSet.IsIncluded(*s++)) array[i] = TRUE;
							} else {
								C_BOOL *b = &sel_array[i];
								for (; n > 0; n--, i++, s++)
								{
									if (*b++)
										if (RngSet.IsIncluded(*s)) array[i] = TRUE;
								}
							}
						}
					}
				}
			}
		}

		if (IsIntersect)
		{
			// TODO: optimized by SIMD
			C_BOOL *p = sel_array, *s = array;
			for (size_t n=array_size; n > 0; n--)
				(*p++) &= (*s++);
		}
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf(INFO_SEL_NUM_VARIANT, PrettyInt(File.VariantSelNum()));

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


// ================================================================

/// set a working space flag with selected annotation id
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceAnnotID(SEXP gdsfile, SEXP ID, SEXP Verbose)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";
	static const char *VarName = "annotation/id";

	int verbose = Rf_asLogical(Verbose);
	if (verbose == NA_LOGICAL)
		error("'verbose' must be TRUE or FALSE.");

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);

		// check
		PdAbstractArray N = File.GetObj(VarName, TRUE);
		int ndim = GDS_Array_DimCnt(N);
		if (ndim != 1)
			throw ErrSeqArray(ERR_DIM, VarName);
		C_Int32 len;
		GDS_Array_GetDim(N, &len, 1);
		if (len != File.VariantNum())
			throw ErrSeqArray(ERR_DIM, VarName);

		TSelection &Sel = File.Selection();
		set<string> id_set;
		const size_t n = XLENGTH(ID);
		for (size_t i=0; i < n; i++)
		{
			SEXP s = STRING_ELT(ID, i);
			if ((s != NA_STRING) && (CHAR(s) != 0))
				id_set.insert(CHAR(s));
		}

		const int SIZE = 4096;
		C_BOOL *p = Sel.pVariant;
		vector<string> buffer(SIZE);
		for (C_Int32 st=0; len > 0; )
		{
			C_Int32 m = (len <= SIZE) ? len : SIZE;
			GDS_Array_ReadData(N, &st, &m, &buffer[0], svStrUTF8);
			for (C_Int32 i=0; i < m; i++)
				*p++ = (id_set.find(buffer[i]) != id_set.end());
			st += m; len -= m;
		}

		Sel.varTrueNum = -1;
		if (verbose)
			Rprintf(INFO_SEL_NUM_VARIANT, PrettyInt(File.VariantSelNum()));

	COREARRAY_CATCH
}


// ================================================================

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
			memcpy(RAW(tmp), Sel.pSample, n);
		} else {
			PROTECT(tmp = NEW_LOGICAL(n));
			int *p = LOGICAL(tmp);
			C_BOOL *s = Sel.pSample;
			for (; n > 0; n--) *p++ = *s++;
		}
		SET_ELEMENT(rv_ans, 0, tmp);

		// variant selection
		n = File.VariantNum();
		if (use_raw_flag)
		{
			PROTECT(tmp = NEW_RAW(n));
			memcpy(RAW(tmp), Sel.pVariant, n);
		} else {
			PROTECT(tmp = NEW_LOGICAL(n));
			int *p = LOGICAL(tmp);
			C_BOOL *s = Sel.pVariant;
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
			sel = s.pVariant;
			SelectCount = File.VariantSelNum();
			s.ClearStructVariant();
		} else if (strcmp(split_str, "by.sample") == 0)
		{
			sel = s.pSample;
			SelectCount = File.SampleSelNum();
			s.ClearStructSample();
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


/// split the selected variants according to multiple processes
COREARRAY_DLL_EXPORT SEXP SEQ_SplitSelectionX(SEXP gdsfile, SEXP index, SEXP split,
	SEXP sel_idx, SEXP sel_variant, SEXP sel_sample, SEXP bl_size, SEXP selection_flag,
	SEXP totlen)
{
	int job_idx = Rf_asInteger(index) - 1;  // starting from 0
	const bool split_by_variant = Rf_asLogical(split)==TRUE;
	const bool sel_flag = Rf_asLogical(selection_flag)==TRUE;
	const int *p_sel_idx = INTEGER(sel_idx);
	const int blsize = Rf_asInteger(bl_size);
	const int tlen = Rf_asInteger(totlen);

	COREARRAY_TRY

		// selection object
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &s = File.Selection();

		// clear and initialize
		int ntot;
		C_BOOL *base_sel, *p_sel;
		if (split_by_variant)
		{
			ntot = File.VariantNum();
			base_sel = (C_BOOL*)RAW(sel_variant);
			p_sel = s.pVariant;
			s.ClearSelectVariant();
		} else {
			ntot = File.SampleNum();
			base_sel = (C_BOOL*)RAW(sel_sample);
			p_sel = s.pSample;
			memset(p_sel, 0, ntot);
		}

		// set selection
		const int st = p_sel_idx[job_idx] - 1;
		int n = 0, i = st;
		for (; n<blsize && i<ntot; i++)
		{
			if (base_sel[i])
				{ p_sel[i] = TRUE; n++; }
		}

		// finalize
		if (split_by_variant)
		{
			s.varTrueNum = n;
			s.varStart = st; s.varEnd = i;
		} else {
			s.ClearStructSample();
		}

		// ---------------------------------------------------
		// output
		if (sel_flag)
		{
			rv_ans = NEW_LOGICAL(tlen);
			int *p = INTEGER(rv_ans);
			memset((void*)p, 0, sizeof(int) * size_t(tlen));
			p += blsize * job_idx;
			for (; n > 0; n--) *p++ = TRUE;
		} else {
			rv_ans = ScalarInteger(n);
		}

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_Summary(SEXP gdsfile, SEXP varname)
{
	COREARRAY_TRY

		// the selection
		CFileInfo &File = GetFileInfo(gdsfile);
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		// the variable name
		string vn = CHAR(STRING_ELT(varname, 0));

		if ((vn=="genotype") || (vn=="phase"))
		{
			PdGDSObj vGeno = GDS_Node_Path(Root, "genotype/data", FALSE);
			if (vGeno == NULL)
			{
				vGeno = GDS_Node_Path(Root, "genotype/~data", FALSE);
			}

			PROTECT(rv_ans = NEW_LIST(2));

				SEXP I32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 0, I32);
				C_Int32 Buf[4];
				if (vGeno)
					GDS_Array_GetDim(vGeno, Buf, 3);
				else
					Buf[2] = NA_INTEGER;
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = File.SampleNum();
				INTEGER(I32)[2] = File.VariantNum();

				SEXP S32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 1, S32);
				INTEGER(S32)[0] = Buf[2];
				INTEGER(S32)[1] = File.SampleSelNum();
				INTEGER(S32)[2] = File.VariantSelNum();

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


/// require reloading chromosom coding
COREARRAY_DLL_EXPORT SEXP SEQ_ResetChrom(SEXP gdsfile)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(gdsfile);
		File.ResetChromosome();
	COREARRAY_CATCH
}



// ===========================================================
// Get system configuration
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_IntAssign(SEXP Dst, SEXP Src)
{
	INTEGER(Dst)[0] = Rf_asInteger(Src);
//	void *p = INTEGER(Dst);
//	Rprintf("addr: %p, val: %d\n", p, Rf_asInteger(Src));
	return R_NilValue;
}



// ===========================================================
// Get system configuration
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_AppendFill(SEXP gdsnode, SEXP val, SEXP cnt)
{
	COREARRAY_TRY
		PdGDSObj obj = GDS_R_SEXP2Obj(gdsnode, FALSE);
		C_Int64 n = (C_Int64)Rf_asReal(cnt);
		const ssize_t SIZE = 65536;
		switch(TYPEOF(val))
		{
		case RAWSXP:
			{
				vector<C_Int8> buf(SIZE, C_Int8(RAW(val)[0]));
				for (ssize_t m; n > 0; n -= m)
				{
					m = (n <= SIZE) ? n : SIZE;
					GDS_Array_AppendData(obj, m, &buf[0], svInt8);
				}
				break;
			}
		case INTSXP:
			{
				vector<C_Int32> buf(SIZE, Rf_asInteger(val));
				for (ssize_t m; n > 0; n -= m)
				{
					m = (n <= SIZE) ? n : SIZE;
					GDS_Array_AppendData(obj, m, &buf[0], svInt32);
				}
				break;
			}
		case REALSXP:
			{
				vector<double> buf(SIZE, Rf_asReal(val));
				for (ssize_t m; n > 0; n -= m)
				{
					m = (n <= SIZE) ? n : SIZE;
					GDS_Array_AppendData(obj, m, &buf[0], svFloat64);
				}
				break;
			}
		default:
			throw ErrSeqArray(
				"Invalid type of 'elm', it should be raw, int or real");
		}
	COREARRAY_CATCH
}



// ===========================================================
// Clear VarMap in a GDS file
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_ClearVarMap(SEXP gdsfile)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(gdsfile);
		File.VarMap().clear();
	COREARRAY_CATCH
}



// ===========================================================
// Get system configuration
// ===========================================================

/// the number of alleles per site
COREARRAY_DLL_EXPORT SEXP SEQ_System()
{
	COREARRAY_TRY

		int nProtect = 0;
		rv_ans = PROTECT(NEW_LIST(3));
		SEXP nm = PROTECT(NEW_CHARACTER(3));
		nProtect += 2;
		SET_NAMES(rv_ans, nm);

		// the number of logical cores
		SET_ELEMENT(rv_ans, 0, ScalarInteger(GDS_Mach_GetNumOfCores()));
		SET_STRING_ELT(nm, 0, mkChar("num.logical.core"));

		// compiler
		SEXP Compiler = PROTECT(NEW_CHARACTER(2));
		nProtect ++;
		SET_ELEMENT(rv_ans, 1, Compiler);
		SET_STRING_ELT(nm, 1, mkChar("compiler"));
	#ifdef __VERSION__
		SET_STRING_ELT(Compiler, 0, mkChar(__VERSION__));
	#endif
	#ifdef __GNUC__
		char buf_compiler[128] = { 0 };
		#ifndef __GNUC_PATCHLEVEL__
		#   define __GNUC_PATCHLEVEL__    0
		#endif
		sprintf(buf_compiler, "GNUG_v%d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		SET_STRING_ELT(Compiler, 1, mkChar(buf_compiler));
	#endif

		// compiler flags
		vector<string> ss;
	#ifdef COREARRAY_SIMD_SSE
		ss.push_back("SSE");
	#endif
	#ifdef COREARRAY_SIMD_SSE2
		ss.push_back("SSE2");
	#endif
	#ifdef COREARRAY_SIMD_SSE3
		ss.push_back("SSE3");
	#endif
	#ifdef COREARRAY_SIMD_SSSE3
		ss.push_back("SSSE3");
	#endif
	#ifdef COREARRAY_SIMD_SSE4_1
		ss.push_back("SSE4.1");
	#endif
	#ifdef COREARRAY_SIMD_SSE4_2
		ss.push_back("SSE4.2");
	#endif
	#ifdef COREARRAY_SIMD_AVX
		ss.push_back("AVX");
	#endif
	#ifdef COREARRAY_SIMD_AVX2
		ss.push_back("AVX2");
	#endif
	#ifdef COREARRAY_SIMD_AVX512F
		ss.push_back("AVX512F");
	#endif
	#ifdef COREARRAY_SIMD_AVX512BW
		ss.push_back("AVX512BW");
	#endif
	#ifdef COREARRAY_SIMD_AVX512CD
		ss.push_back("AVX512CD");
	#endif
	#ifdef COREARRAY_SIMD_AVX512DQ
		ss.push_back("AVX512DQ");
	#endif
	#ifdef COREARRAY_SIMD_AVX512VL
		ss.push_back("AVX512VL");
	#endif
	#ifdef COREARRAY_SIMD_FMA
		ss.push_back("FMA");
	#endif
	#ifdef COREARRAY_SIMD_FMA4
		ss.push_back("FMA4");
	#endif
	#ifdef COREARRAY_POPCNT
		ss.push_back("POPCNT");
	#endif
		SEXP SIMD = PROTECT(NEW_CHARACTER(ss.size()));
		nProtect ++;
		SET_ELEMENT(rv_ans, 2, SIMD);
		SET_STRING_ELT(nm, 2, mkChar("compiler.flag"));
		for (int i=0; i < (int)ss.size(); i++)
			SET_STRING_ELT(SIMD, i, mkChar(ss[i].c_str()));

		UNPROTECT(nProtect);

	COREARRAY_CATCH
}



// ===========================================================
// Debug information
// ===========================================================

/// the number of alleles per site
COREARRAY_DLL_EXPORT SEXP SEQ_Debug(SEXP gdsfile)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(gdsfile);
		int ploidy = File.Ploidy();
		TSelection &Sel = File.Selection();

		Rprintf("Selected samples:\n");
		TSelection::TSampStruct *p = Sel.GetStructSample();
		while (p->length > 0)
		{
			Rprintf("    start: %d, length: %d, sel: %p\n", p->offset/ploidy,
				p->length/ploidy, p->sel);
			p ++;
		}

		Rprintf("Selected variants:\n");
		Sel.GetStructVariant();
		Rprintf("    start: %d, end: %d, num: %d\n", (int)Sel.varStart,
			(int)Sel.varEnd, (int)Sel.varTrueNum);

		return R_NilValue;
	COREARRAY_CATCH
}



// ===========================================================
// Progress object
// ===========================================================

static void free_progress(SEXP ref)
{
    CProgressStdOut *obj = (CProgressStdOut*)R_ExternalPtrAddr(ref);
    if (obj) delete obj;
}

/// Get a progress bar object
COREARRAY_DLL_EXPORT SEXP SEQ_Progress(SEXP Count, SEXP NProc)
{
	C_Int64 TotalCount = (C_Int64)Rf_asReal(Count);
	if (TotalCount < 0)
		error(".seqProgress(): the total number should be >= 0.");
	int nproc = Rf_asInteger(NProc);
	if (nproc <= 0)
		error(".seqProgress(): the number of processes should be > 0.");
	COREARRAY_TRY
		CProgressStdOut *obj = new CProgressStdOut(TotalCount, nproc, true);
		rv_ans = PROTECT(R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
		R_RegisterCFinalizerEx(rv_ans, free_progress, TRUE);
		Rf_setAttrib(rv_ans, R_ClassSymbol, mkString("SeqClass_Progress"));
		UNPROTECT(1);
	COREARRAY_CATCH
}

/// Get a progress bar object
COREARRAY_DLL_EXPORT SEXP SEQ_ProgressAdd(SEXP ref, SEXP inc)
{
	COREARRAY_TRY
		C_Int64 v = (C_Int64)Rf_asReal(inc);
		CProgressStdOut *obj = (CProgressStdOut*)R_ExternalPtrAddr(ref);
		if (obj) obj->Forward(v);
	COREARRAY_CATCH
}



// ===========================================================
// Initialize R objects when the package is loaded
// ===========================================================

COREARRAY_DLL_EXPORT SEXP SEQ_Pkg_Init(SEXP dim_name, SEXP proc_cnt, SEXP proc_idx)
{
	R_Geno_Dim2_Name = VECTOR_ELT(dim_name, 0);
	R_Geno_Dim3_Name = VECTOR_ELT(dim_name, 1);
	R_Dosage_Name = VECTOR_ELT(dim_name, 2);
	R_Data_Name = VECTOR_ELT(dim_name, 3);
	R_Data_Dim2_Name = VECTOR_ELT(dim_name, 4);
	R_Process_Count = INTEGER(proc_cnt);
	R_Process_Index = INTEGER(proc_idx);
	return R_NilValue;
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

COREARRAY_DLL_EXPORT SEXP SEQ_ExternalName5(SEXP v, SEXP w, SEXP x, SEXP y, SEXP z)
{
	return R_NilValue;
}


/// initialize the package
COREARRAY_DLL_EXPORT void R_init_SeqArray(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	extern void Register_SNPRelate_Functions();

	extern SEXP SEQ_GetData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_ConvBED2GDS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	extern SEXP SEQ_MergeAllele(SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeGeno(SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergePhase(SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeInfo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_MergeFormat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	extern SEXP SEQ_BApply_Variant(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	extern SEXP SEQ_Unit_SlidingWindows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	extern SEXP SEQ_bgzip_create(SEXP);

	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_Pkg_Init, 3),
		CALL(SEQ_ExternalName0, 0),         CALL(SEQ_ExternalName1, 1),
		CALL(SEQ_ExternalName2, 2),         CALL(SEQ_ExternalName3, 3),
		CALL(SEQ_ExternalName4, 4),         CALL(SEQ_ExternalName5, 5),

		CALL(SEQ_File_Init, 1),             CALL(SEQ_File_Done, 1),
		CALL(SEQ_FilterPushEmpty, 1),       CALL(SEQ_FilterPushLast, 1),
		CALL(SEQ_FilterPop, 1),

		CALL(SEQ_MergeAllele, 4),           CALL(SEQ_MergeGeno, 5),
		CALL(SEQ_MergePhase, 5),            CALL(SEQ_MergeInfo, 6),
		CALL(SEQ_MergeFormat, 6),

		CALL(SEQ_SetSpaceSample, 4),        CALL(SEQ_SetSpaceSample2, 4),
		CALL(SEQ_SetSpaceVariant, 4),       CALL(SEQ_SetSpaceVariant2, 4),
		CALL(SEQ_SetSpaceChrom, 7),         CALL(SEQ_SetSpaceAnnotID, 3),

		CALL(SEQ_SplitSelection, 5),        CALL(SEQ_SplitSelectionX, 9),
		CALL(SEQ_GetSpace, 2),

		CALL(SEQ_Summary, 2),               CALL(SEQ_System, 0),

		CALL(SEQ_GetData, 6),
		CALL(SEQ_Apply_Sample, 7),          CALL(SEQ_Apply_Variant, 7),
		CALL(SEQ_BApply_Variant, 7),        CALL(SEQ_Unit_SlidingWindows, 7),

		CALL(SEQ_ConvBED2GDS, 6),
		CALL(SEQ_SelectFlag, 2),            CALL(SEQ_ResetChrom, 1),

		CALL(SEQ_IntAssign, 2),             CALL(SEQ_AppendFill, 3),
		CALL(SEQ_ClearVarMap, 1),

		CALL(SEQ_bgzip_create, 1),

		CALL(SEQ_Progress, 2),              CALL(SEQ_ProgressAdd, 2),

		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Register_SNPRelate_Functions();
	Init_GDS_Routines();
}

} // extern "C"
