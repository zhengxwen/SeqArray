// ===========================================================
//
// GetData.cpp: Get data from the GDS file
//
// Copyright (C) 2015-2016    Xiuwen Zheng
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
using namespace SeqArray;

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
COREARRAY_DLL_EXPORT SEXP SEQ_GetData(SEXP gdsfile, SEXP var_name, SEXP UseRaw)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	int use_raw = Rf_asLogical(UseRaw);
	if (use_raw == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");
	const C_UInt32 UseMode = GDS_R_READ_DEFAULT_MODE |
		(use_raw ? GDS_R_READ_ALLOW_RAW_TYPE : 0);

	COREARRAY_TRY

		// File information
		CFileInfo &File = GetFileInfo(gdsfile);
		// the selection
		TSelection &Sel = File.Selection();

		// the path of GDS variable
		const char *name = CHAR(STRING_ELT(var_name, 0));
		if (strcmp(name, "sample.id") == 0)
		{
			// ===========================================================
			// sample.id

			PdAbstractArray N = File.GetObj(name, TRUE);
			// check
			if ((GDS_Array_DimCnt(N) != 1) ||
					(GDS_Array_GetTotalCount(N) != File.SampleNum()))
				throw ErrSeqArray(ERR_DIM, name);
			// read
			C_BOOL *ss = Sel.pSample();
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

		} else if (strcmp(name, "position") == 0)
		{
			int n = File.VariantSelNum();
			if (n > 0)
			{
				const int *base = &File.Position()[0];
				rv_ans = NEW_INTEGER(n);
				int *p = INTEGER(rv_ans);
				C_BOOL *s = Sel.pVariant();
				for (size_t m=File.VariantNum(); m > 0; m--)
				{
					if (*s++) *p++ = *base;
					base ++;
				}
			} else
				rv_ans = NEW_INTEGER(0);

		} else if (strcmp(name, "chromosome") == 0)
		{
			int n = File.VariantSelNum();
			if (n > 0)
			{
				CChromIndex &Chrom = File.Chromosome();
				rv_ans = PROTECT(NEW_CHARACTER(n));
				C_BOOL *s = Sel.pVariant();
				size_t m = File.VariantNum();
				size_t p = 0;
				SEXP last = mkChar("");
				for (size_t i=0; i < m; i++)
				{
					if (*s++)
					{
						const string &ss = Chrom[i];
						if (ss != CHAR(last))
							last = mkChar(ss.c_str());
						SET_STRING_ELT(rv_ans, p++, last);
					}
				}
				UNPROTECT(1);
			} else
				rv_ans = NEW_CHARACTER(0);
		
		} else if ( (strcmp(name, "variant.id")==0) ||
			(strcmp(name, "allele")==0) ||
			(strcmp(name, "annotation/id")==0) ||
			(strcmp(name, "annotation/qual")==0) ||
			(strcmp(name, "annotation/filter")==0) )
		{
			// ===========================================================
			// variant.id, position, chromosome, allele, annotation/id
			// annotation/qual, annotation/filter

			PdAbstractArray N = File.GetObj(name, TRUE);
			// check
			if ((GDS_Array_DimCnt(N) != 1) ||
					(GDS_Array_GetTotalCount(N) != File.VariantNum()))
				throw ErrSeqArray(ERR_DIM, name);
			// read
			C_BOOL *ss = Sel.pVariant();
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

		} else if (strcmp(name, "phase") == 0)
		{
			// ===========================================================
			// phase/

			PdAbstractArray N = File.GetObj("phase/data", TRUE);
			// check
			int ndim = GDS_Array_DimCnt(N);
			C_Int32 dim[4];
			GDS_Array_GetDim(N, dim, 3);
			if (ndim<2 || ndim>3 || dim[0]!= File.VariantNum() ||
					dim[1]!=File.SampleNum())
				throw ErrSeqArray(ERR_DIM, name);
			// read
			C_BOOL *ss[3] = { Sel.pVariant(), Sel.pSample(), NULL };
			if (ndim == 3)
				ss[2] = NeedArrayTRUEs(dim[2]);
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);

		} else if (strcmp(name, "genotype") == 0)
		{
			// ===========================================================
			// genotypic data

			int nSample  = File.SampleSelNum();
			int nVariant = File.VariantSelNum();

			if ((nSample > 0) && (nVariant > 0))
			{
				// initialize GDS genotype Node
				CApply_Variant_Geno NodeVar(File, use_raw);
				// size to be allocated
				ssize_t SIZE = (ssize_t)nSample * File.Ploidy();
				if (use_raw)
				{
					rv_ans = PROTECT(NEW_RAW(nVariant * SIZE));
					C_UInt8 *base = (C_UInt8 *)RAW(rv_ans);
					do {
						NodeVar.ReadGenoData(base);
						base += SIZE;
					} while (NodeVar.Next());
				} else {
					rv_ans = PROTECT(NEW_INTEGER(nVariant * SIZE));
					int *base = INTEGER(rv_ans);
					do {
						NodeVar.ReadGenoData(base);
						base += SIZE;
					} while (NodeVar.Next());
				}

				SEXP dim = PROTECT(NEW_INTEGER(3));
					int *p = INTEGER(dim);
					p[0] = File.Ploidy(); p[1] = nSample; p[2] = nVariant;
				SET_DIM(rv_ans, dim);

				SEXP name_list = PROTECT(NEW_LIST(3));
				SEXP tmp = PROTECT(NEW_CHARACTER(3));
					SET_STRING_ELT(tmp, 0, mkChar("allele"));
					SET_STRING_ELT(tmp, 1, mkChar("sample"));
					SET_STRING_ELT(tmp, 2, mkChar("variant"));
					SET_NAMES(name_list, tmp);
				SET_DIMNAMES(rv_ans, name_list);

				// finally
				UNPROTECT(4);
			}

		} else if (strcmp(name, "@genotype") == 0)
		{
			static const char *VarName = "genotype/@data";
			PdAbstractArray N = File.GetObj(VarName, TRUE);
			// check
			if ((GDS_Array_DimCnt(N) != 1) ||
					(GDS_Array_GetTotalCount(N) != File.VariantNum()))
				throw ErrSeqArray(ERR_DIM, VarName);
			// read
			C_BOOL *ss = Sel.pVariant();
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

		} else if (strncmp(name, "annotation/info/@", 17) == 0)
		{
			PdAbstractArray N = File.GetObj(name, FALSE);
			if (N != NULL)
			{
				// check
				if ((GDS_Array_DimCnt(N) != 1) ||
						(GDS_Array_GetTotalCount(N) != File.VariantNum()))
					throw ErrSeqArray(ERR_DIM, name);
				// read
				C_BOOL *ss = Sel.pVariant();
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);
			}

		} else if (strncmp(name, "annotation/info/", 16) == 0)
		{
			// ===========================================================
			// annotation/info

			GDS_PATH_PREFIX_CHECK(name);
			PdAbstractArray N = File.GetObj(name, TRUE);
			int ndim = GDS_Array_DimCnt(N);
			if ((ndim!=1) && (ndim!=2))
				throw ErrSeqArray(ERR_DIM, name);

			string name2 = GDS_PATH_PREFIX(name, '@');
			PdAbstractArray N_idx = File.GetObj(name2.c_str(), FALSE);
			if (N_idx == NULL)
			{
				// no index
				C_Int32 dim[4];
				GDS_Array_GetDim(N, dim, 2);
				C_BOOL *ss[2] = { Sel.pVariant(), NULL };
				if (ndim == 2)
					ss[1] = NeedArrayTRUEs(dim[1]);
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);
				rv_ans = VAR_LOGICAL(N, rv_ans);

			} else {
				// with index

				// check
				if ((GDS_Array_DimCnt(N_idx) != 1) ||
						(GDS_Array_GetTotalCount(N_idx) != File.VariantNum()))
					throw ErrSeqArray(ERR_DIM, name2.c_str());

				C_Int32 dim[4], dimst[4];
				memset(dimst, 0, sizeof(dimst));
				GDS_Array_GetDim(N, dim, 2);

				vector<int> len;
				vector<C_BOOL> var_sel;
				MAP_INDEX(N_idx, Sel.Variant, len, var_sel, dimst[0], dim[0]);

				C_BOOL *ss[2] = { &var_sel[0], NULL };
				if (ndim == 2)
					ss[1] = NeedArrayTRUEs(dim[1]);

				PROTECT(rv_ans = NEW_LIST(2));
					SEXP I32;
					PROTECT(I32 = NEW_INTEGER(len.size()));
					int *base = INTEGER(I32);
					for (int i=0; i < (int)len.size(); i++)
						base[i] = len[i];
					SET_ELEMENT(rv_ans, 0, I32);
					SET_ELEMENT(rv_ans, 1,
						VAR_LOGICAL(N, GDS_R_Array_Read(N, dimst, dim, ss,
						UseMode)));
				SEXP tmp = PROTECT(NEW_CHARACTER(2));
					SET_STRING_ELT(tmp, 0, mkChar("length"));
					SET_STRING_ELT(tmp, 1, mkChar("data"));
					SET_NAMES(rv_ans, tmp);
				UNPROTECT(3);
			}

		} else if (strncmp(name, "annotation/format/@", 19) == 0)
		{
			string name2(name);
			name2.erase(18, 1).append("/@data");
			PdAbstractArray N = File.GetObj(name2.c_str(), FALSE);
			if (N != NULL)
			{
				// check
				if ((GDS_Array_DimCnt(N) != 1) ||
						(GDS_Array_GetTotalCount(N) != File.VariantNum()))
					throw ErrSeqArray(ERR_DIM, name2.c_str());
				// read
				C_BOOL *ss = Sel.pVariant();
				rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);
			}

		} else if (strncmp(name, "annotation/format/", 18) == 0)
		{
			// ===========================================================
			// annotation/format

			GDS_PATH_PREFIX_CHECK(name);
			string name1 = string(name) + "/data";
			string name2 = string(name) + "/@data";
			PdAbstractArray N = File.GetObj(name1.c_str(), TRUE);
			PdAbstractArray N_idx = File.GetObj(name2.c_str(), TRUE);

			// check
			int ndim = GDS_Array_DimCnt(N);
			if ((ndim!=2) && (ndim!=3))
				throw ErrSeqArray(ERR_DIM, name1.c_str());
			C_Int32 dim[4];
			GDS_Array_GetDim(N, dim, 3);
			if (dim[1] != File.SampleNum())
				throw ErrSeqArray(ERR_DIM, name1.c_str());
			if ((GDS_Array_DimCnt(N_idx) != 1) ||
					(GDS_Array_GetTotalCount(N_idx) != File.VariantNum()))
				throw ErrSeqArray(ERR_DIM, name2.c_str());

			C_Int32 dimst[4];
			memset(dimst, 0, sizeof(dimst));
			vector<int> len;
			vector<C_BOOL> var_sel;
			MAP_INDEX(N_idx, Sel.Variant, len, var_sel, dimst[0], dim[0]);

			C_BOOL *ss[3] = { &var_sel[0], Sel.pSample(), NULL };
			if (ndim == 3)
				ss[2] = NeedArrayTRUEs(dim[2]);

			PROTECT(rv_ans = NEW_LIST(2));
				SEXP I32 = PROTECT(NEW_INTEGER(len.size()));
				int *base = INTEGER(I32);
				for (int i=0; i < (int)len.size(); i++)
					base[i] = len[i];
				SET_ELEMENT(rv_ans, 0, I32);
				SEXP DAT = GDS_R_Array_Read(N, dimst, dim, ss, UseMode);
				SET_ELEMENT(rv_ans, 1, DAT);
			SEXP tmp = PROTECT(NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("length"));
				SET_STRING_ELT(tmp, 1, mkChar("data"));
				SET_NAMES(rv_ans, tmp);
				if (XLENGTH(DAT) > 0)
				{
					SEXP name_list = PROTECT(NEW_LIST(ndim));
					tmp = PROTECT(NEW_CHARACTER(ndim));
					if (ndim == 2)
					{
						SET_STRING_ELT(tmp, 0, mkChar("sample"));
						SET_STRING_ELT(tmp, 1, mkChar("variant"));
					} else {
						SET_STRING_ELT(tmp, 0, mkChar("n"));
						SET_STRING_ELT(tmp, 1, mkChar("sample"));
						SET_STRING_ELT(tmp, 2, mkChar("variant"));
					}
					SET_NAMES(name_list, tmp);
					SET_DIMNAMES(VECTOR_ELT(rv_ans, 1), name_list);
					UNPROTECT(2);
				}
			UNPROTECT(3);

		} else if (strncmp(name, "sample.annotation/", 18) == 0)
		{
			// ===========================================================
			// sample.annotation

			GDS_PATH_PREFIX_CHECK(name);
			PdAbstractArray N = File.GetObj(name, TRUE);
			// check
			int ndim = GDS_Array_DimCnt(N);
			if ((ndim!=1) && (ndim!=2))
				throw ErrSeqArray(ERR_DIM, name);
			C_Int32 dim[2];
			GDS_Array_GetDim(N, dim, 2);
			if (dim[0] != File.SampleNum())
				throw ErrSeqArray(ERR_DIM, name);

			C_BOOL *ss[2] = { Sel.pSample(), NULL };
			if (ndim == 2)
				ss[1] = NeedArrayTRUEs(dim[1]);
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);

		} else if (strcmp(name, "$chrom_pos") == 0)
		{
			// ===========================================================
			// chromosome-position

			PdAbstractArray N1 = File.GetObj("chromosome", TRUE);
			PdAbstractArray N2 = File.GetObj("position", TRUE);
			C_Int64 n1 = GDS_Array_GetTotalCount(N1);
			C_Int64 n2 = GDS_Array_GetTotalCount(N2);
			if ((n1 != n2) || (n1 != File.VariantNum()))
				throw ErrSeqArray("Invalid dimension of 'chromosome' and 'position'.");

			vector<string> chr;
			vector<C_Int32> pos;

			int n = File.VariantSelNum();
			chr.resize(n);
			pos.resize(n);
			C_BOOL *ss = Sel.pVariant();

			GDS_Array_ReadDataEx(N1, NULL, NULL, &ss, &chr[0], svStrUTF8);
			GDS_Array_ReadDataEx(N2, NULL, NULL, &ss, &pos[0], svInt32);

			char buf1[1024] = { 0 };
			char buf2[1024] = { 0 };
			char *p1 = buf1, *p2 = buf2;
			int dup = 0;
			rv_ans = PROTECT(NEW_CHARACTER(n1));
			for (size_t i=0; i < (size_t)n1; i++)
			{
				snprintf(p1, sizeof(buf1), "%s_%d", chr[i].c_str(), pos[i]);
				if (strcmp(p1, p2) == 0)
				{
					dup ++;
					snprintf(p1, sizeof(buf1), "%s_%d_%d", chr[i].c_str(),
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

		} else if (strcmp(name, "$dosage") == 0)
		{
			// ===========================================================
			// dosage data

			ssize_t nSample  = File.SampleSelNum();
			ssize_t nVariant = File.VariantSelNum();

			if ((nSample > 0) && (nVariant > 0))
			{
				// initialize GDS genotype Node
				CApply_Variant_Dosage NodeVar(File, false);

				if (use_raw)
				{
					rv_ans = PROTECT(allocMatrix(RAWSXP, nSample, nVariant));
					C_UInt8 *base = (C_UInt8 *)RAW(rv_ans);
					do {
						NodeVar.ReadDosage(base);
						base += nSample;
					} while (NodeVar.Next());
				} else {
					rv_ans = PROTECT(allocMatrix(INTSXP, nSample, nVariant));
					int *base = INTEGER(rv_ans);
					do {
						NodeVar.ReadDosage(base);
						base += nSample;
					} while (NodeVar.Next());
				}

				SEXP name_list = PROTECT(NEW_LIST(2));
				SEXP tmp = PROTECT(NEW_CHARACTER(2));
					SET_STRING_ELT(tmp, 0, mkChar("sample"));
					SET_STRING_ELT(tmp, 1, mkChar("variant"));
					SET_NAMES(name_list, tmp);
				SET_DIMNAMES(rv_ans, name_list);
				// finally
				UNPROTECT(3);
			}

		} else if (strcmp(name, "$num_allele") == 0)
		{
			// ===========================================================
			// the number of distinct alleles

			ssize_t nVariant = File.VariantSelNum();
			rv_ans = PROTECT(NEW_INTEGER(nVariant));
			int *p = INTEGER(rv_ans);

			CApply_Variant_NumAllele NodeVar(File);
			for (ssize_t i=0; i < nVariant; i++)
			{
				p[i] = NodeVar.GetNumAllele();
				NodeVar.Next();
			}
			UNPROTECT(1);

		} else {
			throw ErrSeqArray(
				"'%s' is not a standard variable name, and the standard format:\n"
				"\tsample.id, variant.id, position, chromosome, allele, genotype\n"
				"\tannotation/id, annotation/qual, annotation/filter\n"
				"\tannotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME\n"
				"\tsample.annotation/VARIABLE_NAME", name);
		}

	COREARRAY_CATCH
}

} // extern "C"
