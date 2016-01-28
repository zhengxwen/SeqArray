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

/// the initial data
TInitObject::TInitObject(): GENO_BUFFER(1024)
{
	memset(TRUE_ARRAY, TRUE, sizeof(TRUE_ARRAY));
}

TInitObject::TSelection &TInitObject::Selection(SEXP gdsfile, bool alloc)
{
	// TODO: check whether handle is valid
	int id = INTEGER(GetListElement(gdsfile, "id"))[0];
	TSelList &m = _Map[id];
	if (m.empty()) m.push_back(TSelection());
	TSelection &s = m.back();
	if (alloc && (s.Sample.empty() | s.Variant.empty()))
	{
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		// selection
		if (s.Sample.empty())
		{
			PdAbstractArray N = GDS_Node_Path(Root, "sample.id", TRUE);
			s.Sample.resize(GDS_Array_GetTotalCount(N), TRUE);
		}
		if (s.Variant.empty())
		{
			PdAbstractArray N = GDS_Node_Path(Root, "variant.id", TRUE);
			s.Variant.resize(GDS_Array_GetTotalCount(N), TRUE);
		}
	}
	return s;
}

void TInitObject::Need_GenoBuffer(size_t size)
{
	if (size > GENO_BUFFER.size())
		GENO_BUFFER.resize(size);
}

TInitObject Init;



// ===========================================================
// GDS Variable Type
// ===========================================================

C_BOOL *CVarApply::NeedTRUE(size_t size)
{
	if (size <= sizeof(Init.TRUE_ARRAY))
	{
		return Init.TRUE_ARRAY;
	} else {
		_TRUE.resize(size, TRUE);
		return &_TRUE[0];
	}
}



// ===========================================================
// Library Functions
// ===========================================================

/// Get the list element named str, or return NULL
COREARRAY_DLL_LOCAL int MatchElement(const char *txt, const char *list[],
	size_t nlist)
{
	for (size_t i=0; i < nlist; i++)
	{
		if (strcmp(txt, *list) == 0)
			return i;
		list ++;
	}
	return -1;
}

/// Get the list element named str, or return NULL
COREARRAY_DLL_LOCAL SEXP GetListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	const size_t n = GetLength(names);
	for (size_t i = 0; i < n; i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

/// Get the total count requiring the number of dimension is one
COREARRAY_DLL_LOCAL int GetGDSObjCount(PdAbstractArray Obj, const char *varname)
{
	if (GDS_Array_DimCnt(Obj) != 1)
		throw ErrSeqArray("Invalid dimension of '%s'!", varname);
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
// ===========================================================
// Open a GDS file
// ===========================================================

/// initialize a SeqArray file
COREARRAY_DLL_EXPORT SEXP SEQ_File_Init(SEXP gdsfile)
{
	COREARRAY_TRY
		TInitObject::TSelection &s = Init.Selection(gdsfile);
		s.Sample.clear();
		s.Variant.clear();
	COREARRAY_CATCH
}

/// finalize a SeqArray file
COREARRAY_DLL_EXPORT SEXP SEQ_File_Done(SEXP gdsfile)
{
	COREARRAY_TRY
		int gds_file_id = Rf_asInteger(GetListElement(gdsfile, "id"));
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(gds_file_id);
		if (it != Init._Map.end())
			Init._Map.erase(it);
	COREARRAY_CATCH
}



// ===========================================================
// Set a working space
// ===========================================================

/// push the current filter to the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushEmpty(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(GetListElement(gdsfile, "id"));
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(id);
		if (it != Init._Map.end())
		{
			it->second.push_back(TInitObject::TSelection());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// push the current filter to the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPushLast(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(GetListElement(gdsfile, "id"));
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(id);
		if (it != Init._Map.end())
		{
			if (!it->second.empty())
				it->second.push_back(it->second.back());
			else
				it->second.push_back(TInitObject::TSelection());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// pop up the previous filter from the stack
COREARRAY_DLL_EXPORT SEXP SEQ_FilterPop(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = Rf_asInteger(GetListElement(gdsfile, "id"));
		map<int, TInitObject::TSelList>::iterator it =
			Init._Map.find(id);
		if (it != Init._Map.end())
		{
			if (it->second.size() <= 1)
				throw ErrSeqArray("No filter can be pop up.");
			it->second.pop_back();
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

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		PdAbstractArray varSamp = GDS_Node_Path(Root, "sample.id", TRUE);
		int Count = GetGDSObjCount(varSamp, "sample.id");

		vector<C_BOOL> &flag_array = Init.Selection(gdsfile).Sample;
		if (flag_array.empty())
			flag_array.resize(Count, TRUE);
		C_BOOL *pArray = &flag_array[0];

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
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'sample.id'.");

		int n = GetNumOfTRUE(&flag_array[0], flag_array.size());
		if (Rf_isNull(samp_id)) n = Count;
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected samples: %d\n", n);

	COREARRAY_CATCH
}


/// set a working space with selected sample id (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceSample2(SEXP gdsfile, SEXP samp_sel,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		PdAbstractArray varSamp = GDS_Node_Path(Root, "sample.id", TRUE);
		int Count = GetGDSObjCount(varSamp, "sample.id");

		vector<C_BOOL> &flag_array = Init.Selection(gdsfile).Sample;
		if (flag_array.empty())
			flag_array.resize(Count, TRUE);
		C_BOOL *pArray = &flag_array[0];

		if (Rf_isLogical(samp_sel) || IS_RAW(samp_sel))
		{
			// a logical vector for selected samples
			if (!intersect_flag)
			{
				if (XLENGTH(samp_sel) != Count)
					throw ErrSeqArray("Invalid length of 'samp.sel'.");
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
				if ((size_t)XLENGTH(samp_sel) != GetNumOfTRUE(pArray, flag_array.size()))
				{
					throw ErrSeqArray(
						"Invalid length of 'samp.sel' "
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
						throw ErrSeqArray("Out of range 'samp.sel'.");
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
				int Cnt = GetNumOfTRUE(&flag_array[0], flag_array.size());
				int *pI = INTEGER(samp_sel);
				R_xlen_t N = XLENGTH(samp_sel);
				// check
				for (R_xlen_t i=0; i < N; i++)
				{
					int I = *pI ++;
					if ((I != NA_INTEGER) && ((I < 1) || (I > Cnt)))
						throw ErrSeqArray("Out of range 'samp.sel'.");
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
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'samp.sel'.");

		int n = GetNumOfTRUE(&flag_array[0], flag_array.size());
		if (Rf_isNull(samp_sel)) n = Count;
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected samples: %d\n", n);

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant(SEXP gdsfile, SEXP var_id,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		PdAbstractArray varVariant = GDS_Node_Path(Root, "variant.id", TRUE);
		int Count = GetGDSObjCount(varVariant, "variant.id");

		vector<C_BOOL> &flag_array = Init.Selection(gdsfile).Variant;
		if (flag_array.empty())
			flag_array.resize(Count, TRUE);
		C_BOOL *pArray = &flag_array[0];

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
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'variant.id'.");

		int n = GetNumOfTRUE(&flag_array[0], flag_array.size());
		if (Rf_isNull(var_id)) n = Count;
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected variants: %d\n", n);

	COREARRAY_CATCH
}


/// set a working space with selected variant id (logical/raw vector, or index)
COREARRAY_DLL_EXPORT SEXP SEQ_SetSpaceVariant2(SEXP gdsfile, SEXP var_sel,
	SEXP intersect, SEXP verbose)
{
	int intersect_flag = Rf_asLogical(intersect);

	COREARRAY_TRY

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		PdAbstractArray varVariant = GDS_Node_Path(Root, "variant.id", TRUE);
		int Count = GetGDSObjCount(varVariant, "variant.id");

		vector<C_BOOL> &flag_array = Init.Selection(gdsfile).Variant;
		if (flag_array.empty())
			flag_array.resize(Count, TRUE);
		C_BOOL *pArray = &flag_array[0];

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
				if ((size_t)XLENGTH(var_sel) != GetNumOfTRUE(pArray, flag_array.size()))
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
				int Cnt = GetNumOfTRUE(&flag_array[0], flag_array.size());
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
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'variant.sel'.");

		int n = GetNumOfTRUE(&flag_array[0], flag_array.size());
		if (Rf_isNull(var_sel)) n = Count;
		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected variants: %d\n", n);

	COREARRAY_CATCH
}


/// set a working space flag with selected chromosome(s)
COREARRAY_DLL_EXPORT SEXP SEQ_SetChrom(SEXP gdsfile, SEXP include,
	SEXP is_num, SEXP verbose)
{
	int IsNum = Rf_asLogical(is_num);
	if (!Rf_isNull(include))
		include = AS_CHARACTER(include);

	COREARRAY_TRY

		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		PdAbstractArray varVariant = GDS_Node_Path(Root, "variant.id", TRUE);
		int nVariant = GDS_Array_GetTotalCount(varVariant);

		PdAbstractArray varChrom = GDS_Node_Path(Root, "chromosome", TRUE);
		if ((GDS_Array_DimCnt(varChrom) != 1) ||
				(nVariant != GDS_Array_GetTotalCount(varChrom)))
			throw ErrSeqArray("Invalid 'chromosome'.");

		set<string> Inc;
		bool IncFlag = (Rf_isNull(include) != TRUE);
		if (IncFlag)
		{
			R_xlen_t n = XLENGTH(include);
			for (R_xlen_t i=0; i < n; i++)
				Inc.insert(CHAR(STRING_ELT(include, i)));
		}

		vector<C_BOOL> &array = Init.Selection(gdsfile).Variant;
		array.resize(nVariant);
		string txt;
		int nSum = 0;

		for (C_Int32 i=0; i < nVariant; i++)
		{
			static const C_Int32 ONE = 1;
			GDS_Array_ReadData(varChrom, &i, &ONE, &txt, svStrUTF8);

			bool flag = true;
			if (IsNum != NA_INTEGER)
			{
				// whether val is numeric
				char *endptr = (char*)(txt.c_str());
				strtol(txt.c_str(), &endptr, 10);
				flag = (endptr != txt.c_str());
				if (IsNum == FALSE) flag = !flag;
			}
			if (IncFlag && flag)
				flag = (Inc.find(txt) != Inc.end());

			array[i] = flag;
			if (flag) nSum ++;
		}

		if (Rf_asLogical(verbose) == TRUE)
			Rprintf("# of selected variants: %d\n", nSum);

	COREARRAY_CATCH
}


/// set a working space flag with selected variant id
COREARRAY_DLL_EXPORT SEXP SEQ_GetSpace(SEXP gdsfile, SEXP UseRaw)
{
	int use_raw_flag = Rf_asLogical(UseRaw);
	if (use_raw_flag == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");

	COREARRAY_TRY

		// the GDS root node
		TInitObject::TSelection &sel = Init.Selection(gdsfile);
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SEXP tmp;

		if (sel.Sample.empty())
		{
			PdAbstractArray var = GDS_Node_Path(Root, "sample.id", TRUE);
			size_t n = GDS_Array_GetTotalCount(var);
			if (use_raw_flag)
			{
				PROTECT(tmp = NEW_RAW(n));
				memset(RAW(tmp), TRUE, n);
			} else {
				PROTECT(tmp = NEW_LOGICAL(n));
				int *p = LOGICAL(tmp);
				for (; n > 0; n--) *p++ = TRUE;
			}
		} else {
			size_t n = sel.Sample.size();
			if (use_raw_flag)
			{
				PROTECT(tmp = NEW_RAW(n));
				memcpy(RAW(tmp), &sel.Sample[0], n);
			} else {
				PROTECT(tmp = NEW_LOGICAL(n));
				int *p = LOGICAL(tmp);
				C_BOOL *s = &(sel.Sample[0]);
				for (; n > 0; n--) *p++ = *s++;
			}
		}
		SET_ELEMENT(rv_ans, 0, tmp);

		if (sel.Variant.empty())
		{
			PdAbstractArray var = GDS_Node_Path(Root, "variant.id", TRUE);
			size_t n = GDS_Array_GetTotalCount(var);
			if (use_raw_flag)
			{
				PROTECT(tmp = NEW_RAW(n));
				memset(RAW(tmp), TRUE, n);
			} else {
				PROTECT(tmp = NEW_LOGICAL(n));
				int *p = LOGICAL(tmp);
				for (; n > 0; n--) *p++ = TRUE;
			}
		} else {
			size_t n = sel.Variant.size();
			if (use_raw_flag)
			{
				PROTECT(tmp = NEW_RAW(n));
				memcpy(RAW(tmp), &sel.Variant[0], n);
			} else {
				PROTECT(tmp = NEW_LOGICAL(n));
				int *p = LOGICAL(tmp);
				C_BOOL *s = &(sel.Variant[0]);
				for (; n > 0; n--) *p++ = *s++;
			}
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
		TInitObject::TSelection &s = Init.Selection(gdsfile);

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
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		// the variable name
		string vn = CHAR(STRING_ELT(varname, 0));

		if ((vn=="genotype") || (vn=="phase"))
		{
			PdGDSObj vSample  = GDS_Node_Path(Root, "sample.id", TRUE);
			PdGDSObj vVariant = GDS_Node_Path(Root, "variant.id", TRUE);
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
				SEXP I32, S32;

				PROTECT(I32 = NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 0, I32);
				int Buf[256];
				GDS_Array_GetDim(vGeno, Buf, 3);
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = GDS_Array_GetTotalCount(vSample);
				INTEGER(I32)[2] = GDS_Array_GetTotalCount(vVariant);

				PROTECT(S32 = NEW_INTEGER(2));
				SET_ELEMENT(rv_ans, 1, S32);
				if (!Sel.Sample.empty())
				{
					int &n = INTEGER(S32)[0]; n = 0;
					vector<C_BOOL>::iterator it;
					for (it=Sel.Sample.begin(); it != Sel.Sample.end(); it ++)
						if (*it) n ++;
				} else
					INTEGER(S32)[0] = INTEGER(I32)[1];
				if (!Sel.Variant.empty())
				{
					int &n = INTEGER(S32)[1]; n = 0;
					vector<C_BOOL>::iterator it;
					for (it=Sel.Variant.begin(); it != Sel.Variant.end(); it ++)
						if (*it) n ++;
				} else
					INTEGER(S32)[1] = INTEGER(I32)[2];

			SEXP tmp;
			PROTECT(tmp = NEW_CHARACTER(2));
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
// analysis
// ===========================================================

/// the number of alleles per site
COREARRAY_DLL_EXPORT SEXP SEQ_NumOfAllele(SEXP allele_node)
{
	COREARRAY_TRY

		// GDS nodes
		PdAbstractArray N = GDS_R_SEXP2Obj(allele_node, TRUE);
		if (GDS_Array_DimCnt(N) != 1)
			throw ErrSeqArray("Invalid dimension of 'allele'!");
		int Count = GetGDSObjCount(N, "allele");

		// allocate integers
		PROTECT(rv_ans = NEW_INTEGER(Count));
		int *base = INTEGER(rv_ans);
		string s;

		for (C_Int32 i=0; i < Count; i ++)
		{
			static const C_Int32 ONE = 1;
			GDS_Array_ReadData(N, &i, &ONE, &s, svStrUTF8);
			base[i] = GetNumOfAllele(s.c_str());
		}

		UNPROTECT(1);

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
		CALL(SEQ_SplitSelection, 5),        CALL(SEQ_SetChrom, 4),
		CALL(SEQ_GetSpace, 2),

		CALL(SEQ_Summary, 2),

		CALL(SEQ_GetData, 2),
		CALL(SEQ_Apply_Sample, 7),          CALL(SEQ_Apply_Variant, 8),

		CALL(SEQ_ConvBEDFlag, 3),           CALL(SEQ_ConvBED2GDS, 5),
		CALL(SEQ_SelectFlag, 2),

		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Register_SNPRelate_Functions();
	Init_GDS_Routines();
}

} // extern "C"
