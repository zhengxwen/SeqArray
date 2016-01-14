// ===========================================================
//
// GetData.cpp: Get data from the GDS file
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

#include "ReadBySample.h"
#include "ReadByVariant.h"
#include <stdio.h>


extern "C"
{
// ===========================================================
// Get data from a working space
// ===========================================================

static SEXP VAR_LOGICAL(PdGDSObj Node, SEXP Array)
{
	char classname[32];
	classname[0] = 0;
	GDS_Node_GetClassName(Node, classname, sizeof(classname));
	if (strcmp(classname, "dBit1") == 0)
	{
		PROTECT(Array);
		Array = AS_LOGICAL(Array);
		UNPROTECT(1);
	}
	return Array;
}

static void MAP_INDEX(PdAbstractArray Node, const vector<C_BOOL> &sel,
	vector<int> &out_len, vector<C_BOOL> &out_var_sel,
	C_Int32 &out_var_start, C_Int32 &out_var_count)
{
	if (GDS_Array_DimCnt(Node) != 1)
		throw ErrSeqArray("Invalid dimension.");
	C_Int64 Cnt = GDS_Array_GetTotalCount(Node);

	if (sel.empty())
	{
		out_len.resize(Cnt);
		C_Int32 _st=0, _cnt=Cnt;
		GDS_Array_ReadData(Node, &_st, &_cnt, &out_len[0], svInt32);

		out_var_start = 0;
		out_var_count = 0;
		for (vector<int>::iterator it=out_len.begin();
			it != out_len.end(); it++)
		{
			if (*it > 0) out_var_count += *it;
		}
		out_var_sel.clear();
		out_var_sel.resize(out_var_count, TRUE);

	} else {
		// check
		if ((int)sel.size() != Cnt)
			throw ErrSeqArray("Invalid dimension.");

		// find the start
		int _start = 0;
		for (; _start < (int)sel.size(); _start++)
			if (sel[_start]) break;
		// find the end
		int _end = sel.size()-1;
		for (; _end >= 0; _end --)
			if (sel[_end]) break;

		if (_end >= 0)
		{
			const int N_MAX = 16384;
			C_Int32 buffer[N_MAX];

			out_var_start = 0;
			int pos = 0;
			while (pos < _start)
			{
				int L = _start - pos;
				if (L > N_MAX) L = N_MAX;
				GDS_Array_ReadData(Node, &pos, &L, buffer, svInt32);
				pos += L;
				for (int i=0; i < L; i++)
				{
					if (buffer[i] > 0)
						out_var_start += buffer[i];
				}
			}

			out_len.clear();
			out_var_sel.clear();
			while (pos <= _end)
			{
				int L = _end - pos + 1;
				if (L > N_MAX) L = N_MAX;
				GDS_Array_ReadData(Node, &pos, &L, buffer, svInt32);
				for (int i=0; i < L; i++)
				{
					int LL = (buffer[i] > 0) ? buffer[i] : 0;
					if (sel[pos+i])
					{
						out_len.push_back(LL);
						for (int j=0; j < LL; j++)
							out_var_sel.push_back(TRUE);
					} else {
						for (int j=0; j < LL; j++)
							out_var_sel.push_back(FALSE);
					}
				}
				pos += L;
			}
			
			out_var_count = out_var_sel.size();
		} else {
			out_len.clear(); out_var_sel.clear();
			out_var_start = out_var_count = 0;
		}
	}
}


/// Get data from a working space
COREARRAY_DLL_EXPORT SEXP SEQ_GetData(SEXP gdsfile, SEXP var_name)
{
	COREARRAY_TRY

		SEXP tmp;

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gdsfile, "root"), TRUE);

		// 
		C_BOOL *SelPtr[GDS_MAX_NUM_DIMENSION];
		int DStart[GDS_MAX_NUM_DIMENSION], DLen[GDS_MAX_NUM_DIMENSION];
		int DimCnt;

		// the path of GDS variable
		const char *s = CHAR(STRING_ELT(var_name, 0));
		if (strcmp(s, "sample.id") == 0)
		{
			// ===========================================================
			// sample.id

			PdAbstractArray N = GDS_Node_Path(Root, s, TRUE);
			DimCnt = GDS_Array_DimCnt(N);
			if (DimCnt != 1)
				throw ErrSeqArray("Invalid dimension of 'sample.id'.");
			if (Sel.Sample.empty())
			{
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, NULL, 0);
			} else {
				GDS_Array_GetDim(N, DLen, 1);
				if ((int)Sel.Sample.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of 'sample.id'.");
				SelPtr[0] = &Sel.Sample[0];
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, &SelPtr[0], 0);
			}

		} else if ( (strcmp(s, "variant.id")==0) || (strcmp(s, "position")==0) ||
			(strcmp(s, "chromosome")==0) || (strcmp(s, "allele")==0) ||
			(strcmp(s, "annotation/id")==0) || (strcmp(s, "annotation/qual")==0) ||
			(strcmp(s, "annotation/filter")==0) )
		{
			// ===========================================================
			// variant.id, position, chromosome, allele, annotation/id
			// annotation/qual, annotation/filter

			PdAbstractArray N = GDS_Node_Path(Root, s, TRUE);
			DimCnt = GDS_Array_DimCnt(N);
			if (DimCnt != 1)
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			if (Sel.Variant.empty())
			{
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, NULL, 0);
			} else {
				GDS_Array_GetDim(N, DLen, 1);
				if ((int)Sel.Variant.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);
				SelPtr[0] = &Sel.Variant[0];
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, &SelPtr[0], 0);
			}

		} else if (strcmp(s, "phase") == 0)
		{
			// ===========================================================
			// phase/

			PdAbstractArray N = GDS_Node_Path(Root, "phase/data", TRUE);
			DimCnt = GDS_Array_DimCnt(N);
			if ((DimCnt != 2) && (DimCnt != 3))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			if (!Sel.Sample.empty() || !Sel.Variant.empty())
			{
				GDS_Array_GetDim(N, DLen, 3);

				if (Sel.Variant.empty())
					Sel.Variant.resize(DLen[0], TRUE);
				else if ((int)Sel.Variant.size() != DLen[0])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);

				if (Sel.Sample.empty())
					Sel.Sample.resize(DLen[1], TRUE);
				else if ((int)Sel.Sample.size() != DLen[1])
					throw ErrSeqArray("Invalid dimension of '%s'.", s);

				CVarApply Var;

				SelPtr[0] = &Sel.Variant[0];
				SelPtr[1] = &Sel.Sample[0];
				if (DimCnt == 3)
					SelPtr[2] = Var.NeedTRUE(DLen[2]);

				rv_ans = GDS_R_Array_Read(N, NULL, NULL, &SelPtr[0], 0);
			} else {
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, NULL, 0);
			}

		} else if (strcmp(s, "genotype") == 0)
		{
			// ===========================================================
			// genotypic data

			// init selection
			if (Sel.Sample.empty())
			{
				PdAbstractArray N = GDS_Node_Path(Root, "sample.id", TRUE);
				int Cnt = GDS_Array_GetTotalCount(N);
				if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'sample.id'.");
				Sel.Sample.resize(Cnt, TRUE);
			}
			if (Sel.Variant.empty())
			{
				PdAbstractArray N = GDS_Node_Path(Root, "variant.id", TRUE);
				int Cnt = GDS_Array_GetTotalCount(N);
				if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'variant.id'.");
				Sel.Variant.resize(Cnt, TRUE);
			}

			// the number of selected variants
			int nVariant = 0;
			for (vector<C_BOOL>::iterator it = Sel.Variant.begin();
				it != Sel.Variant.end(); it ++)
			{
				if (*it) nVariant ++;
			}
			if (nVariant > 0)
			{
				// initialize the GDS Node list
				CVarApplyByVariant NodeVar;
				NodeVar.InitObject(CVariable::ctGenotype,
					"genotype/data", Root, Sel.Variant.size(),
					&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0], false);

				// the number of calling PROTECT
				int SIZE = NodeVar.Num_Sample * NodeVar.DLen[2];
				PROTECT(rv_ans = NEW_INTEGER(nVariant * SIZE));
				PROTECT(tmp = NEW_INTEGER(3));
					INTEGER(tmp)[0] = NodeVar.DLen[2];
					INTEGER(tmp)[1] = NodeVar.Num_Sample;
					INTEGER(tmp)[2] = nVariant;
				SET_DIM(rv_ans, tmp);
				SEXP name_list;
				PROTECT(name_list = NEW_LIST(3));
				PROTECT(tmp = NEW_CHARACTER(3));
					SET_STRING_ELT(tmp, 0, mkChar("allele"));
					SET_STRING_ELT(tmp, 1, mkChar("sample"));
					SET_STRING_ELT(tmp, 2, mkChar("variant"));
					SET_NAMES(name_list, tmp);
				SET_DIMNAMES(rv_ans, name_list);

				int *base = INTEGER(rv_ans);
				do {
					NodeVar.ReadGenoData(base);
					base += SIZE;
				} while (NodeVar.NextCell());

				// finally
				UNPROTECT(4);
			}

		} else if (strcmp(s, "@genotype") == 0)
		{
			static const char *VarName = "genotype/@data";
			PdAbstractArray N = GDS_Node_Path(Root, VarName, TRUE);
			C_Int32 st  = 0;
			C_Int32 cnt = GetGDSObjCount(N, VarName);
			if (Sel.Variant.empty())
			{
				rv_ans = GDS_R_Array_Read(N, &st, &cnt, NULL, 0);
			} else {
				C_BOOL *SelList = &Sel.Variant[0];
				rv_ans = GDS_R_Array_Read(N, &st, &cnt, &SelList, 0);
			}

		} else if (strncmp(s, "annotation/info/@", 17) == 0)
		{
			PdAbstractArray N = GDS_Node_Path(Root, s, FALSE);
			if (N != NULL)
			{
				C_Int32 st  = 0;
				C_Int32 cnt = GetGDSObjCount(N, s);
				if (Sel.Variant.empty())
				{
					rv_ans = GDS_R_Array_Read(N, &st, &cnt, NULL, 0);
				} else {
					C_BOOL *SelList = &Sel.Variant[0];
					rv_ans = GDS_R_Array_Read(N, &st, &cnt, &SelList, 0);
				}
			}
		} else if (strncmp(s, "annotation/info/", 16) == 0)
		{
			GDS_PATH_PREFIX_CHECK(s);
			PdAbstractArray N = GDS_Node_Path(Root, s, TRUE);
			DimCnt = GDS_Array_DimCnt(N);
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);

			string path_ex = GDS_PATH_PREFIX(s, '@');
			PdAbstractArray N_idx = GDS_Node_Path(Root, path_ex.c_str(), FALSE);
			if (N_idx == NULL)
			{
				// no index
				if (!Sel.Variant.empty())
				{
					GDS_Array_GetDim(N, DLen, 2);
					CVarApply Var;
					SelPtr[0] = &Sel.Variant[0];
					if (DimCnt == 2)
						SelPtr[1] = Var.NeedTRUE(DLen[1]);
					rv_ans = GDS_R_Array_Read(N, NULL, NULL, &SelPtr[0], 0);
				} else
					rv_ans = GDS_R_Array_Read(N, NULL, NULL, NULL, 0);

				rv_ans = VAR_LOGICAL(N, rv_ans);

			} else {
				// with index
				if (!Sel.Variant.empty())
				{
					memset(DStart, 0, sizeof(DStart));
					GDS_Array_GetDim(N, DLen, 2);

					vector<int> len;
					vector<C_BOOL> var_sel;
					MAP_INDEX(N_idx, Sel.Variant, len, var_sel, DStart[0], DLen[0]);

					CVarApply Var;
					SelPtr[0] = &var_sel[0];
					if (DimCnt == 2)
						SelPtr[1] = Var.NeedTRUE(DLen[1]);

					PROTECT(rv_ans = NEW_LIST(2));
						SEXP I32;
						PROTECT(I32 = NEW_INTEGER(len.size()));
						int *base = INTEGER(I32);
						for (int i=0; i < (int)len.size(); i++)
							base[i] = len[i];
						SET_ELEMENT(rv_ans, 0, I32);
						SET_ELEMENT(rv_ans, 1,
							VAR_LOGICAL(N, GDS_R_Array_Read(N, DStart, DLen, &SelPtr[0], 0)));
					PROTECT(tmp = NEW_CHARACTER(2));
						SET_STRING_ELT(tmp, 0, mkChar("length"));
						SET_STRING_ELT(tmp, 1, mkChar("data"));
						SET_NAMES(rv_ans, tmp);
					UNPROTECT(3);

				} else {
					PROTECT(rv_ans = NEW_LIST(2));
						SET_ELEMENT(rv_ans, 0,
							GDS_R_Array_Read(N_idx, NULL, NULL, NULL, 0));
						SET_ELEMENT(rv_ans, 1,
							VAR_LOGICAL(N, GDS_R_Array_Read(N, NULL, NULL, NULL, 0)));
					PROTECT(tmp = NEW_CHARACTER(2));
						SET_STRING_ELT(tmp, 0, mkChar("length"));
						SET_STRING_ELT(tmp, 1, mkChar("data"));
						SET_NAMES(rv_ans, tmp);
					UNPROTECT(2);
				}
			}

		} else if (strncmp(s, "annotation/format/@", 19) == 0)
		{
			string name(s);
			name.erase(18, 1).append("/@data");
			PdAbstractArray N = GDS_Node_Path(Root, name.c_str(), FALSE);
			if (N != NULL)
			{
				C_Int32 st  = 0;
				C_Int32 cnt = GetGDSObjCount(N, name.c_str());
				if (Sel.Variant.empty())
				{
					rv_ans = GDS_R_Array_Read(N, &st, &cnt, NULL, 0);
				} else {
					C_BOOL *SelList = &Sel.Variant[0];
					rv_ans = GDS_R_Array_Read(N, &st, &cnt, &SelList, 0);
				}
			}
		} else if (strncmp(s, "annotation/format/", 18) == 0)
		{
			GDS_PATH_PREFIX_CHECK(s);
			PdAbstractArray N =
				GDS_Node_Path(Root, string(string(s)+"/data").c_str(), TRUE);
			PdAbstractArray N_idx =
				GDS_Node_Path(Root, string(string(s)+"/@data").c_str(), TRUE);

			DimCnt = GDS_Array_DimCnt(N);
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			memset(DStart, 0, sizeof(DStart));
			GDS_Array_GetDim(N, DLen, 3);

			if (Sel.Sample.empty())
				Sel.Sample.resize(DLen[1], TRUE);
			if (Sel.Variant.empty())
				Sel.Variant.resize(GDS_Array_GetTotalCount(N_idx), TRUE);

			vector<int> len;
			vector<C_BOOL> var_sel;
			MAP_INDEX(N_idx, Sel.Variant, len, var_sel, DStart[0], DLen[0]);

			CVarApply Var;
			SelPtr[0] = &var_sel[0];
			SelPtr[1] = &Sel.Sample[0];
			if (DimCnt == 3)
				SelPtr[2] = Var.NeedTRUE(DLen[2]);

			PROTECT(rv_ans = NEW_LIST(2));
				SEXP I32;
				PROTECT(I32 = NEW_INTEGER(len.size()));
				int *base = INTEGER(I32);
				for (int i=0; i < (int)len.size(); i++)
					base[i] = len[i];
				SET_ELEMENT(rv_ans, 0, I32);
				SEXP DAT = GDS_R_Array_Read(N, DStart, DLen, &SelPtr[0], 0);
				SET_ELEMENT(rv_ans, 1, DAT);
			PROTECT(tmp = NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("length"));
				SET_STRING_ELT(tmp, 1, mkChar("data"));
				SET_NAMES(rv_ans, tmp);

				if (Rf_length(DAT) > 0)
				{
					SEXP name_list;
					PROTECT(name_list = NEW_LIST(DimCnt));
					PROTECT(tmp = NEW_CHARACTER(DimCnt));
						SET_STRING_ELT(tmp, 0, mkChar("sample"));
						SET_STRING_ELT(tmp, 1, mkChar("variant"));
						SET_NAMES(name_list, tmp);
					SET_DIMNAMES(VECTOR_ELT(rv_ans, 1), name_list);
					UNPROTECT(5);
				} else {
					UNPROTECT(3);
				}

		} else if (strncmp(s, "sample.annotation/", 18) == 0)
		{
			GDS_PATH_PREFIX_CHECK(s);
			PdAbstractArray N = GDS_Node_Path(Root, s, TRUE);
			int nSamp = GDS_Array_GetTotalCount(
				GDS_Node_Path(Root, "sample.id", TRUE));

			DimCnt = GDS_Array_DimCnt(N);
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray("Invalid dimension of '%s'.", s);
			GDS_Array_GetDim(N, DLen, 2);
			if (DLen[0] != nSamp)
				throw ErrSeqArray("Invalid dimension of '%s'.", s);

			CVarApply Var;
			if (Sel.Sample.empty())
				Sel.Sample.resize(nSamp, TRUE);
			SelPtr[0] = &Sel.Sample[0];
			if (DimCnt == 2)
				SelPtr[1] = Var.NeedTRUE(DLen[1]);

			memset(DStart, 0, sizeof(DStart));
			rv_ans = GDS_R_Array_Read(N, DStart, DLen, &SelPtr[0], 0);

		} else if (strcmp(s, "chrom-pos") == 0)
		{
			// ===========================================================
			// chromosome-position
			static const char *ERR_CHR_POS =
				"Invalid dimension of 'chromosome' and 'position'.";

			PdAbstractArray N1 = GDS_Node_Path(Root, "chromosome", TRUE);
			PdAbstractArray N2 = GDS_Node_Path(Root, "position", TRUE);
			C_Int64 n1 = GDS_Array_GetTotalCount(N1),
				n2 = GDS_Array_GetTotalCount(N2);
			if (n1 != n2)
				throw ErrSeqArray(ERR_CHR_POS);

			vector<string> chr;
			vector<C_Int32> pos;
			C_Int32 st=0, cnt=n1;

			if (Sel.Variant.empty())
			{
				chr.resize(n1); pos.resize(n1);
				GDS_Array_ReadData(N1, &st, &cnt, &chr[0], svStrUTF8);
				GDS_Array_ReadData(N2, &st, &cnt, &pos[0], svInt32);
			} else {
				if (Sel.Variant.size() != (size_t)n1)
					throw ErrSeqArray(ERR_CHR_POS);
				n1 = GetNumOfTRUE(&Sel.Variant[0], Sel.Variant.size());
				chr.resize(n1); pos.resize(n1);
				SelPtr[0] = &Sel.Variant[0];
				GDS_Array_ReadDataEx(N1, &st, &cnt, &SelPtr[0], &chr[0], svStrUTF8);
				GDS_Array_ReadDataEx(N2, &st, &cnt, &SelPtr[0], &pos[0], svInt32);
			}

			char buf1[1024] = { 0 };
			char buf2[1024] = { 0 };
			char *p1 = buf1, *p2 = buf2;
			int dup = 0;
			rv_ans = PROTECT(NEW_CHARACTER(n1));
			for (size_t i=0; i < (size_t)n1; i++)
			{
				snprintf(p1, sizeof(buf1), "%s-%d", chr[i].c_str(), pos[i]);
				if (strcmp(p1, p2) == 0)
				{
					dup ++;
					snprintf(p1, sizeof(buf1), "%s-%d.%d", chr[i].c_str(),
						pos[i], dup);
					SET_STRING_ELT(rv_ans, i, mkChar(p1));
				} else {
					char *tmp;
					tmp = p1; p1 = p2; p2 = tmp;
					SET_STRING_ELT(rv_ans, i, mkChar(p2));
					dup = 0;
				}
			}
			UNPROTECT(1);

		} else {
			throw ErrSeqArray(
				"'%s' is not a standard variable name, and the standard format:\n"
				"\tsample.id, variant.id, position, chromosome, allele, genotype\n"
				"\tannotation/id, annotation/qual, annotation/filter\n"
				"\tannotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME\n"
				"\tsample.annotation/VARIABLE_NAME", s);
		}

	COREARRAY_CATCH
}

} // extern "C"
