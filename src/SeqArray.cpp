// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// SeqArray.cpp: the C/C++ codes for the SeqArray package
//
// Copyright (C) 2013-2015    Xiuwen Zheng [zhengx@u.washington.edu]
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



// ===========================================================
// The initialized object
// ===========================================================

/// the initial data
TInitObject::TInitObject(): TRUE_ARRAY(256, TRUE), GENO_BUFFER(1024)
{ }

TInitObject::TSelection &TInitObject::Selection(SEXP gds)
{
	// TODO: check whether handle is valid
	int id = INTEGER(GetListElement(gds, "id"))[0];
	TSelList &m = _Map[id];
	if (m.empty()) m.push_back(TSelection());
	return m.back();
}

void TInitObject::Check_TrueArray(int Cnt)
{
	if (Cnt > (int)TRUE_ARRAY.size())
		TRUE_ARRAY.resize(Cnt, TRUE);
}

TInitObject Init;



extern "C"
{
// ===========================================================
// the initial and final functions
// ===========================================================

/// initialize the package
COREARRAY_DLL_EXPORT void R_init_SeqArray(DllInfo *info)
{
	Init_GDS_Routines();
}



// ===========================================================
// Open a GDS file
// ===========================================================

/// initialize a SeqArray file
COREARRAY_DLL_EXPORT SEXP sqa_Open_Init(SEXP gdsfile)
{
	COREARRAY_TRY
		TInitObject::TSelection &s = Init.Selection(gdsfile);
		s.Sample.clear();
		s.Variant.clear();
	COREARRAY_CATCH
}

/// finalize a SeqArray file
COREARRAY_DLL_EXPORT SEXP sqa_File_Done(SEXP gdsfile)
{
	COREARRAY_TRY
		int gds_file_id = INTEGER(GetListElement(gdsfile, "id"))[0];
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
COREARRAY_DLL_EXPORT SEXP sqa_FilterPushEmpty(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = INTEGER(GetListElement(gdsfile, "id"))[0];
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
COREARRAY_DLL_EXPORT SEXP sqa_FilterPushLast(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = INTEGER(GetListElement(gdsfile, "id"))[0];
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
COREARRAY_DLL_EXPORT SEXP sqa_FilterPop(SEXP gdsfile)
{
	COREARRAY_TRY
		int id = INTEGER(GetListElement(gdsfile, "id"))[0];
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
COREARRAY_DLL_EXPORT SEXP sqa_SetSpaceSample(SEXP gds, SEXP samp_sel,
	SEXP verbose)
{
	COREARRAY_TRY

		TInitObject::TSelection &s = Init.Selection(gds);

		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gds, "root"));
		PdSequenceX varSamp = GDS_Node_Path(Root, "sample.id", TRUE);

		if (GDS_Seq_DimCnt(varSamp) != 1)
			throw ErrSeqArray("Invalid dimension of 'sample.id'!");
		int Count = GDS_Seq_GetTotalCount(varSamp);

		vector<C_BOOL> flag_array(Count);

		if (IS_LOGICAL(samp_sel))
		{
			// a logical vector for selected samples
			if (Rf_length(samp_sel) != Count)
				throw ErrSeqArray("Invalid length of 'samp.sel'.");
			// set selection
			int *base = LOGICAL(samp_sel);
			for (int i=0; i < Count; i++)
				flag_array[i] = (base[i] == TRUE);

		} else if (IS_INTEGER(samp_sel))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(samp_sel), &INTEGER(samp_sel)[Rf_length(samp_sel)]);
			// sample id
			vector<int> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varSamp, &_st, &_cnt, &sample_id[0], svInt32);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<int>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_NUMERIC(samp_sel))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(samp_sel), &REAL(samp_sel)[Rf_length(samp_sel)]);
			// sample id
			vector<double> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varSamp, &_st, &_cnt, &sample_id[0], svFloat64);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<double>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_CHARACTER(samp_sel))
		{
			// initialize
			set<string> set_id;
			for (int i=0; i < Rf_length(samp_sel); i++)
				set_id.insert(string(CHAR(STRING_ELT(samp_sel, i))));
			// sample id
			vector<string> sample_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varSamp, &_st, &_cnt, &sample_id[0], svStrUTF8);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<string>::iterator it = set_id.find(sample_id[i]);
				flag_array[i] = (it != set_id.end());
			}
		
		} else if (isNull(samp_sel))
		{
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'samp_sel'.");

		int n = 0;
		for (vector<C_BOOL>::iterator it=flag_array.begin();
			it != flag_array.end(); it ++)
		{
			if (*it != 0) n ++;
		}
		if (isNull(samp_sel)) n = Count;
		if (n > 0)
		{
			s.Sample = flag_array;
			if (LOGICAL(verbose)[0] == TRUE)
				Rprintf("# of selected samples: %d\n", n);
		} else
			throw ErrSeqArray("No sample selected!");

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP sqa_SetSpaceVariant(SEXP gds, SEXP var_sel,
	SEXP verbose)
{
	COREARRAY_TRY

		TInitObject::TSelection &s = Init.Selection(gds);

		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gds, "root"));
		PdSequenceX varVariant = GDS_Node_Path(Root, "variant.id", TRUE);

		if (GDS_Seq_DimCnt(varVariant) != 1)
			throw ErrSeqArray("Invalid dimension of 'variant.id'!");
		int Count = GDS_Seq_GetTotalCount(varVariant);

		vector<C_BOOL> flag_array(Count);

		if (IS_LOGICAL(var_sel))
		{
			// a logical vector for selected samples
			if (Rf_length(var_sel) != Count)
				throw ErrSeqArray("Invalid length of 'variant.sel'.");
			// set selection
			int *base = LOGICAL(var_sel);
			for (int i=0; i < Count; i++)
				flag_array[i] = (base[i] == TRUE);

		} else if (IS_INTEGER(var_sel))
		{
			// initialize
			set<int> set_id;
			set_id.insert(INTEGER(var_sel), &INTEGER(var_sel)[Rf_length(var_sel)]);
			// sample id
			vector<int> var_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varVariant, &_st, &_cnt, &var_id[0], svInt32);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<int>::iterator it = set_id.find(var_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_NUMERIC(var_sel))
		{
			// initialize
			set<double> set_id;
			set_id.insert(REAL(var_sel), &REAL(var_sel)[Rf_length(var_sel)]);
			// variant id
			vector<double> variant_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varVariant, &_st, &_cnt, &variant_id[0], svFloat64);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<double>::iterator it = set_id.find(variant_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (IS_CHARACTER(var_sel))
		{
			// initialize
			set<string> set_id;
			for (int i=0; i < Rf_length(var_sel); i++)
				set_id.insert(string(CHAR(STRING_ELT(var_sel, i))));
			// sample id
			vector<string> variant_id(Count);
			C_Int32 _st=0, _cnt=Count;
			GDS_Seq_rData(varVariant, &_st, &_cnt, &variant_id[0], svStrUTF8);
			// set selection
			for (int i=0; i < Count; i++)
			{
				set<string>::iterator it = set_id.find(variant_id[i]);
				flag_array[i] = (it != set_id.end());
			}

		} else if (isNull(var_sel))
		{
			flag_array.clear();
		} else
			throw ErrSeqArray("Invalid type of 'samp_sel'.");

		int n = 0;
		for (vector<C_BOOL>::iterator it=flag_array.begin();
			it != flag_array.end(); it ++)
		{
			if (*it != 0) n ++;
		}
		if (isNull(var_sel)) n = Count;
		if (n > 0)
		{
			s.Variant = flag_array;
			if (LOGICAL(verbose)[0] == TRUE)
				Rprintf("# of selected variants: %d\n", n);
		} else
			throw ErrSeqArray("No variant selected!");

	COREARRAY_CATCH
}


static void CLEAR_SEL_VALUE(int num, vector<C_BOOL>::iterator &it)
{
	while (num > 0)
	{
		if (*it != 0) { num --; *it = FALSE; }
		it ++;
	}
}
static void SKIP_SEL(int num, vector<C_BOOL>::iterator &it)
{
	while (num > 0)
	{
		if (*it != 0) num --;
		it ++;
	}
}

/// split the selected variants according to multiple processes
COREARRAY_DLL_EXPORT SEXP sqa_SplitSelectedVariant(SEXP gdsfile, SEXP Index,
	SEXP n_process)
{
	COREARRAY_TRY

		// selection object
		TInitObject::TSelection &s = Init.Selection(gdsfile);

		// the index process starting from 1
		int Process_Index = INTEGER(AS_INTEGER(Index))[0] - 1;
		int Num_Process = INTEGER(AS_INTEGER(n_process))[0];

		// the total number of selected variants
		vector<C_BOOL>::iterator it;
		int N_Variant = 0;
		for (it=s.Variant.begin(); it != s.Variant.end();)
		{
			if (*it != 0) N_Variant ++;
			it ++;
		}
		if (N_Variant <= 0) error("No variant!");

		// split a list
		vector<int> split(Num_Process);
		double avg = (double)N_Variant / Num_Process;
		double start = 0;
		for (int i=0; i < Num_Process; i++)
		{
			start += avg;
			split[i] = (int)(start + 0.5);
		}

		// ---------------------------------------------------
		it = s.Variant.begin();
		int st = 0;
		for (int i=0; i < Process_Index; i++)
		{
			CLEAR_SEL_VALUE(split[i] - st, it);
			st = split[i];
		}
		int ans_n = split[Process_Index] - st;
		SKIP_SEL(ans_n, it);
		st = split[Process_Index];
		for (int i=Process_Index+1; i < Num_Process; i++)
		{
			CLEAR_SEL_VALUE(split[i] - st, it);
			st = split[i];
		}

		// ---------------------------------------------------
		// output
		rv_ans = ScalarInteger(ans_n);

	COREARRAY_CATCH
}


/// split the selected samples according to multiple processes
COREARRAY_DLL_EXPORT SEXP sqa_SplitSelectedSample(SEXP gdsfile, SEXP Index,
	SEXP n_process)
{
	COREARRAY_TRY

		// selection object
		TInitObject::TSelection &s = Init.Selection(gdsfile);

		// the index process starting from 1
		int Process_Index = INTEGER(AS_INTEGER(Index))[0] - 1;
		int Num_Process = INTEGER(AS_INTEGER(n_process))[0];

		// the total number of selected samples
		vector<C_BOOL>::iterator it;
		int N_Sample = 0;
		for (it=s.Sample.begin(); it != s.Sample.end();)
		{
			if (*it != 0) N_Sample ++;
			it ++;
		}
		if (N_Sample <= 0) error("No sample!");

		// split a list
		vector<int> split(Num_Process);
		double avg = (double)N_Sample / Num_Process;
		double start = 0;
		for (int i=0; i < Num_Process; i++)
		{
			start += avg;
			split[i] = (int)(start + 0.5);
		}

		// ---------------------------------------------------
		it = s.Sample.begin();
		int st = 0;
		for (int i=0; i < Process_Index; i++)
		{
			CLEAR_SEL_VALUE(split[i] - st, it);
			st = split[i];
		}
		int ans_n = split[Process_Index] - st;
		SKIP_SEL(ans_n, it);
		st = split[Process_Index];
		for (int i=Process_Index+1; i < Num_Process; i++)
		{
			CLEAR_SEL_VALUE(split[i] - st, it);
			st = split[i];
		}

		// ---------------------------------------------------
		// output
		rv_ans = ScalarInteger(ans_n);

	COREARRAY_CATCH
}


/// set a working space flag with selected variant id
COREARRAY_DLL_EXPORT SEXP sqa_GetSpace(SEXP gdsfile)
{
	COREARRAY_TRY

		TInitObject::TSelection &s = Init.Selection(gdsfile);

		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gdsfile, "root"));
		PdSequenceX varSample = GDS_Node_Path(Root, "sample.id", TRUE);
		PdSequenceX varVariant = GDS_Node_Path(Root, "variant.id", TRUE);

		int nProtected = 0;
		SEXP tmp;

		PROTECT(rv_ans = NEW_LIST(2));
		nProtected ++;

		if (s.Sample.empty())
		{
			int L = GDS_Seq_GetTotalCount(varSample);
			PROTECT(tmp = NEW_LOGICAL(L));
			nProtected ++;
			for (int i=0; i < L; i++) LOGICAL(tmp)[i] = TRUE;
		} else {
			PROTECT(tmp = NEW_LOGICAL(s.Sample.size()));
			nProtected ++;
			for (int i=0; i < (int)s.Sample.size(); i++)
				LOGICAL(tmp)[i] = (s.Sample[i] != 0);
		}
		SET_ELEMENT(rv_ans, 0, tmp);

		if (s.Variant.empty())
		{
			int L = GDS_Seq_GetTotalCount(varVariant);
			PROTECT(tmp = NEW_LOGICAL(L));
			nProtected ++;
			for (int i=0; i < L; i++) LOGICAL(tmp)[i] = TRUE;
		} else {
			PROTECT(tmp = NEW_LOGICAL(s.Variant.size()));
			nProtected ++;
			for (int i=0; i < (int)s.Variant.size(); i++)
				LOGICAL(tmp)[i] = (s.Variant[i] != 0);
		}
		SET_ELEMENT(rv_ans, 1, tmp);

		PROTECT(tmp = NEW_CHARACTER(2));
		nProtected ++;
			SET_STRING_ELT(tmp, 0, mkChar("sample.sel"));
			SET_STRING_ELT(tmp, 1, mkChar("variant.sel"));
		SET_NAMES(rv_ans, tmp);

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// set a working space with selected variant id
COREARRAY_DLL_EXPORT SEXP sqa_VarSummary(SEXP gdsfile, SEXP varname)
{
	COREARRAY_TRY

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gdsfile, "root"));
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
				GDS_Seq_GetDim(vGeno, Buf, 3);
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = GDS_Seq_GetTotalCount(vSample);
				INTEGER(I32)[2] = GDS_Seq_GetTotalCount(vVariant);

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
		}

	COREARRAY_CATCH
}



// ===========================================================
// analysis
// ===========================================================

/// the number of alleles per site
COREARRAY_DLL_EXPORT SEXP sqa_NumOfAllele(SEXP allele_node)
{
	COREARRAY_TRY

		// GDS nodes
		PdSequenceX N;
		memcpy(&N, INTEGER(allele_node), sizeof(N));

		if (GDS_Seq_DimCnt(N) != 1)
			throw ErrSeqArray("Invalid dimension!");
		int Count = GDS_Seq_GetTotalCount(N);

		// allocate integers
		PROTECT(rv_ans = NEW_INTEGER(Count));

		int *base = INTEGER(rv_ans);
		string s;

		for (int i=0; i < Count; i ++)
		{
			C_Int32 _st = i;
			C_Int32 _cnt = 1;
			GDS_Seq_rData(N, &_st, &_cnt, &s, svStrUTF8);

			// determine how many alleles
			int num_allele = 0;
			const char *p = s.c_str();
			while (*p != 0)
			{
				num_allele ++;
				while ((*p != 0) && (*p != ',')) p ++;
				if (*p == ',') p ++;
			}
			base[i] = num_allele;
		}

		UNPROTECT(1);

	COREARRAY_CATCH
}




// ===========================================================
// analysis
// ===========================================================

COREARRAY_DLL_EXPORT SEXP seq_Merge_Pos(SEXP opfile, SEXP outgds_root)
{
/*
	// GDS nodes
	PdSequenceX Root;
	memcpy(&Root, INTEGER(outgds_root), sizeof(Root));
	PdSequenceX varPos = GDS_Node_Path(Root, "position");

	// for - loop
	for (int i=0; i < Rf_length(opfile); i++)
	{
		// get variable
		PdSequenceX _R_;
		memcpy(&_R_, INTEGER(VECTOR_ELT(opfile, i)), sizeof(_R_));
		PdSequenceX sPos = GDS_Node_Path(_R_, "position");

		// read and write
		
	}
*/
	return R_NilValue;
}



COREARRAY_DLL_EXPORT SEXP seq_missing_snp(SEXP geno)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int miss_cnt = 0;

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (C_UInt32(*pp) > 2)
				{ miss_cnt ++; break; }
		}
		p += num_ploidy;
	}

	SEXP rv;
	PROTECT(rv = NEW_NUMERIC(1));
	REAL(rv)[0] = (double)miss_cnt / num_sample;
	UNPROTECT(1);

	return rv;
}


COREARRAY_DLL_EXPORT SEXP seq_missing_samp(SEXP geno, SEXP miss_cnt)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int *miss = INTEGER(miss_cnt);

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (C_UInt32(*pp) > 2)
				{ miss[i] ++; break; }
		}
		p += num_ploidy;
	}

	return R_NilValue;
}


COREARRAY_DLL_EXPORT SEXP seq_allele_freq(SEXP geno)
{
	SEXP dim = getAttrib(geno, R_DimSymbol);
	int num_ploidy = INTEGER(dim)[0];
	int num_sample = INTEGER(dim)[1];
	int ref_cnt=0, valid_cnt=0;

	int *p = INTEGER(geno);
	for (int i=0; i < num_sample; i++)
	{
		int *pp = p;
		for (int j=0; j < num_ploidy; j++, pp++)
		{
			if (C_UInt32(*pp) == 0) ref_cnt ++;
			if (C_UInt32(*pp) <= 2) valid_cnt ++;
		}
		p += num_ploidy;
	}

	SEXP rv;
	PROTECT(rv = NEW_NUMERIC(1));
	REAL(rv)[0] = (double)ref_cnt / valid_cnt;
	UNPROTECT(1);

	return rv;
}

} // extern "C"
