// ===========================================================
//
// GetData.cpp: Get data from the GDS file
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

#include "Index.h"
#include "ReadByVariant.h"
#include "ReadBySample.h"


using namespace SeqArray;

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


// get data
static SEXP VarGetData(CFileInfo &File, const char *name, bool use_raw)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	SEXP rv_ans = R_NilValue;
	const C_UInt32 UseMode = GDS_R_READ_DEFAULT_MODE |
		(use_raw ? GDS_R_READ_ALLOW_RAW_TYPE : 0);

	TSelection &Sel = File.Selection();

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
		C_BOOL *ss = Sel.pSample;
		rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

	} else if (strcmp(name, "position") == 0)
	{
		int n = File.VariantSelNum();
		if (n > 0)
		{
			const int *base = &File.Position()[0];
			rv_ans = NEW_INTEGER(n);
			int *p = INTEGER(rv_ans);
			C_BOOL *s = Sel.pVariant;
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
			C_BOOL *s = Sel.pVariant;
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
		// variant.id, allele, annotation/id, annotation/qual, annotation/filter

		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss = Sel.pVariant;
		rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

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
			SET_DIMNAMES(rv_ans, R_Geno_Dim3_Name);

			// finally
			UNPROTECT(2);
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
		C_BOOL *ss = Sel.pVariant;
		rv_ans = GDS_R_Array_Read(N, NULL, NULL, &ss, UseMode);

	} else if (strcmp(name, "$dosage") == 0)
	{
		// ===========================================================
		// dosage data

		ssize_t nSample  = File.SampleSelNum();
		ssize_t nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Dosage NodeVar(File, false, false);

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

			SET_DIMNAMES(rv_ans, R_Dosage_Name);
			// finally
			UNPROTECT(1);
		}

	} else if (strcmp(name, "$dosage_alt") == 0)
	{
		// ===========================================================
		// dosage data

		ssize_t nSample  = File.SampleSelNum();
		ssize_t nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Dosage NodeVar(File, false, true);

			if (use_raw)
			{
				rv_ans = PROTECT(allocMatrix(RAWSXP, nSample, nVariant));
				C_UInt8 *base = (C_UInt8 *)RAW(rv_ans);
				do {
					NodeVar.ReadDosageAlt(base);
					base += nSample;
				} while (NodeVar.Next());
			} else {
				rv_ans = PROTECT(allocMatrix(INTSXP, nSample, nVariant));
				int *base = INTEGER(rv_ans);
				do {
					NodeVar.ReadDosageAlt(base);
					base += nSample;
				} while (NodeVar.Next());
			}

			SET_DIMNAMES(rv_ans, R_Dosage_Name);
			// finally
			UNPROTECT(1);
		}

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
		C_BOOL *ss[3] = { Sel.pVariant, Sel.pSample, NULL };
		if (ndim == 3)
			ss[2] = NeedArrayTRUEs(dim[2]);
		rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);

	} else if (strncmp(name, "annotation/info/@", 17) == 0)
	{
		if (File.GetObj(name, FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name);
			rv_ans = V.GetLen_Sel(Sel.pVariant);
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
			C_BOOL *ss[2] = { Sel.pVariant, NULL };
			if (ndim == 2)
				ss[1] = NeedArrayTRUEs(dim[1]);
			rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);
			rv_ans = VAR_LOGICAL(N, rv_ans);

		} else {
			// with index
			CIndex &V = File.VarIndex(name2);
			int var_start, var_count;
			vector<C_BOOL> var_sel;
			SEXP I32 = PROTECT(V.GetLen_Sel(Sel.pVariant, var_start, var_count, var_sel));

			C_BOOL *ss[2] = { &var_sel[0], NULL };
			C_Int32 dimst[2]  = { var_start, 0 };
			C_Int32 dimcnt[2] = { var_count, 0 };
			if (ndim == 2)
			{
				GDS_Array_GetDim(N, dimcnt, 2);
				dimcnt[0] = var_count;
			}

			PROTECT(rv_ans = NEW_LIST(2));
				SET_ELEMENT(rv_ans, 0, I32);
				SET_ELEMENT(rv_ans, 1,
					VAR_LOGICAL(N, GDS_R_Array_Read(N, dimst, dimcnt, ss,
					UseMode)));
			SET_NAMES(rv_ans, R_Data_Name);
			UNPROTECT(2);
		}

	} else if (strncmp(name, "annotation/format/@", 19) == 0)
	{
		string name2(name);
		name2.erase(18, 1).append("/@data");
		if (File.GetObj(name2.c_str(), FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name2.c_str());
			rv_ans = V.GetLen_Sel(Sel.pVariant);
		}

	} else if (strncmp(name, "annotation/format/", 18) == 0)
	{
		// ===========================================================
		// annotation/format

		GDS_PATH_PREFIX_CHECK(name);
		string name1 = string(name) + "/data";
		string name2 = string(name) + "/@data";
		PdAbstractArray N = File.GetObj(name1.c_str(), TRUE);

		// with index
		CIndex &V = File.VarIndex(name2);
		int var_start, var_count;
		vector<C_BOOL> var_sel;
		SEXP I32 = PROTECT(V.GetLen_Sel(Sel.pVariant, var_start, var_count, var_sel));

		C_BOOL *ss[2] = { &var_sel[0], Sel.pSample };
		C_Int32 dimst[2]  = { var_start, 0 };
		C_Int32 dimcnt[2];
		GDS_Array_GetDim(N, dimcnt, 2);
		dimcnt[0] = var_count;

		PROTECT(rv_ans = NEW_LIST(2));
			SET_ELEMENT(rv_ans, 0, I32);
			SEXP DAT = GDS_R_Array_Read(N, dimst, dimcnt, ss, UseMode);
			SET_ELEMENT(rv_ans, 1, DAT);
			SET_NAMES(rv_ans, R_Data_Name);
			if (XLENGTH(DAT) > 0)
				SET_DIMNAMES(DAT, R_Data_Dim2_Name);
		UNPROTECT(2);

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

		C_BOOL *ss[2] = { Sel.pSample, NULL };
		if (ndim == 2)
			ss[1] = NeedArrayTRUEs(dim[1]);
		rv_ans = GDS_R_Array_Read(N, NULL, NULL, ss, UseMode);

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

	} else if (strcmp(name, "$ref") == 0)
	{
		// ===========================================================
		// the reference allele

		PdAbstractArray N = File.GetObj("allele", TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		size_t n = File.VariantSelNum();
		vector<string> buffer(n);
		C_BOOL *ss = Sel.pVariant;
		GDS_Array_ReadDataEx(N, NULL, NULL, &ss, &buffer[0], svStrUTF8);
		// output
		rv_ans = PROTECT(NEW_CHARACTER(n));
		for (size_t i=0; i < n; i++)
		{
			const char *p = buffer[i].c_str();
			size_t m = 0;
			for (const char *s=p; *s!=',' && *s!=0; s++) m++;
			SET_STRING_ELT(rv_ans, i, mkCharLen(p, m));
		}
		UNPROTECT(1);

	} else if (strcmp(name, "$alt") == 0)
	{
		// ===========================================================
		// the reference allele

		PdAbstractArray N = File.GetObj("allele", TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		size_t n = File.VariantSelNum();
		vector<string> buffer(n);
		C_BOOL *ss = Sel.pVariant;
		GDS_Array_ReadDataEx(N, NULL, NULL, &ss, &buffer[0], svStrUTF8);
		// output
		rv_ans = PROTECT(NEW_CHARACTER(n));
		for (size_t i=0; i < n; i++)
		{
			const char *p = buffer[i].c_str();
			for (; *p!=',' && *p!=0; p++);
			if (*p == ',') p++;
			SET_STRING_ELT(rv_ans, i, mkChar(p));
		}
		UNPROTECT(1);

	} else if (strcmp(name, "$chrom_pos") == 0)
	{
		// ===========================================================
		// chromosome-position

		PdAbstractArray N1 = File.GetObj("chromosome", TRUE);
		PdAbstractArray N2 = File.GetObj("position", TRUE);
		C_Int64 n1 = GDS_Array_GetTotalCount(N1);
		C_Int64 n2 = GDS_Array_GetTotalCount(N2);
		if ((n1 != n2) || (n1 != File.VariantNum()))
			throw ErrSeqArray("Invalid dimension of 'chromosome' or 'position'.");

		int n = File.VariantSelNum();
		vector<string> chr(n);
		vector<C_Int32> pos(n);

		C_BOOL *ss = Sel.pVariant;
		GDS_Array_ReadDataEx(N1, NULL, NULL, &ss, &chr[0], svStrUTF8);
		GDS_Array_ReadDataEx(N2, NULL, NULL, &ss, &pos[0], svInt32);

		char buf1[1024] = { 0 };
		char buf2[1024] = { 0 };
		char *p1 = buf1, *p2 = buf2;
		int dup = 0;
		rv_ans = PROTECT(NEW_CHARACTER(n1));
		for (size_t i=0; i < (size_t)n1; i++)
		{
			snprintf(p1, sizeof(buf1), "%s:%d", chr[i].c_str(), pos[i]);
			if (strcmp(p1, p2) == 0)
			{
				dup ++;
				snprintf(p1, sizeof(buf1), "%s:%d_%d", chr[i].c_str(),
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

	} else if (strcmp(name, "$chrom_pos_allele") == 0)
	{
		// ===========================================================
		// chromosome-position-allele

		PdAbstractArray N1 = File.GetObj("chromosome", TRUE);
		PdAbstractArray N2 = File.GetObj("position", TRUE);
		PdAbstractArray N3 = File.GetObj("allele", TRUE);
		C_Int64 n1 = GDS_Array_GetTotalCount(N1);
		C_Int64 n2 = GDS_Array_GetTotalCount(N2);
		C_Int64 n3 = GDS_Array_GetTotalCount(N3);
		if ((n1 != n2) || (n1 != n3) || (n1 != File.VariantNum()))
			throw ErrSeqArray("Invalid dimension of 'chromosome', 'position' or 'allele'.");

		int n = File.VariantSelNum();
		vector<string> chr(n);
		vector<C_Int32> pos(n);
		vector<string> allele(n);

		C_BOOL *ss = Sel.pVariant;
		GDS_Array_ReadDataEx(N1, NULL, NULL, &ss, &chr[0], svStrUTF8);
		GDS_Array_ReadDataEx(N2, NULL, NULL, &ss, &pos[0], svInt32);
		GDS_Array_ReadDataEx(N3, NULL, NULL, &ss, &allele[0], svStrUTF8);

		char buf[65536] = { 0 };
		rv_ans = PROTECT(NEW_CHARACTER(n1));
		for (size_t i=0; i < (size_t)n1; i++)
		{
			const char *s = allele[i].c_str();
			for (char *p=(char*)s; *p != 0; p++)
				if (*p == ',') *p = '_';
			snprintf(buf, sizeof(buf), "%s:%d_%s", chr[i].c_str(), pos[i], s);
			SET_STRING_ELT(rv_ans, i, mkChar(buf));
		}
		UNPROTECT(1);

	} else {
		throw ErrSeqArray(
			"'%s' is not a standard variable name, and the standard format:\n"
			"    sample.id, variant.id, position, chromosome, allele, genotype\n"
			"    annotation/id, annotation/qual, annotation/filter\n"
			"    annotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME\n"
			"    sample.annotation/VARIABLE_NAME, etc", name);
	}

	return rv_ans;
}


/// Get data from a working space
COREARRAY_DLL_EXPORT SEXP SEQ_GetData(SEXP gdsfile, SEXP var_name, SEXP UseRaw)
{
	if (!Rf_isString(var_name) || RLength(var_name)!=1)
		error("'var.name' should be a string of length one.");
	int use_raw = Rf_asLogical(UseRaw);
	if (use_raw == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");

	COREARRAY_TRY
		// File information
		CFileInfo &File = GetFileInfo(gdsfile);
		// Get data
		rv_ans = VarGetData(File, CHAR(STRING_ELT(var_name, 0)), use_raw);
	COREARRAY_CATCH
}



COREARRAY_DLL_LOCAL extern const char *Txt_Apply_VarIdx[];


/// Apply functions over variants in block
COREARRAY_DLL_EXPORT SEXP SEQ_BApply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho)
{
	int bsize = Rf_asInteger(RGetListElement(param, "bsize"));
	if (bsize < 1)
		error("'bsize' must be >= 1.");

	SEXP pam_use_raw = RGetListElement(param, "useraw");
	if (!Rf_isLogical(pam_use_raw))
		error("'.useraw' must be TRUE, FALSE or NA.");
	int use_raw_flag = Rf_asLogical(pam_use_raw);

	int prog_flag = Rf_asLogical(RGetListElement(param, "progress"));
	if (prog_flag == NA_LOGICAL)
		error("'.progress' must be TRUE or FALSE.");

	COREARRAY_TRY

		// File information
		CFileInfo &File = GetFileInfo(gdsfile);
		// Selection
		TSelection &Selection = File.Selection();

		// the number of selected variants
		int nVariant = File.VariantSelNum();
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");

		// the number of data blocks
		int NumBlock = nVariant / bsize;
		if (nVariant % bsize) NumBlock ++;

		// the number of calling PROTECT
		int nProtected = 0;

		// as.is
		Rconnection OutputConn = NULL;
		PdGDSObj OutputGDS = NULL;
		int DatType;
		if (Rf_inherits(as_is, "connection"))
		{
			OutputConn = R_GetConnection(as_is);
			DatType = 2;
		} else if (Rf_inherits(as_is, "gdsn.class"))
		{
			OutputGDS = GDS_R_SEXP2Obj(as_is, FALSE);
			DatType = 3;
		} else {
			const char *s = CHAR(STRING_ELT(as_is, 0));
			if (strcmp(s, "list")==0 || strcmp(s, "unlist")==0)
			{
				DatType = 1;
				rv_ans = PROTECT(NEW_LIST(NumBlock)); nProtected ++;
			} else {
				DatType = 0;
			}
		}

		// rho environment
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");

		// var.index
		int VarIdx = MatchText(CHAR(STRING_ELT(var_index, 0)), Txt_Apply_VarIdx);
		if (VarIdx < 0)
			throw ErrSeqArray("'var.index' is not valid!");
		SEXP R_Index = R_NilValue;
		if (VarIdx > 0)
		{
			PROTECT(R_Index = NEW_INTEGER(1));
			nProtected ++;
		}

		// calling
		SEXP R_fcall = R_NilValue;
		SEXP R_call_param = R_NilValue;
		int num_var = RLength(var_name);
		if (num_var > 1)
		{
			PROTECT(R_call_param = NEW_LIST(num_var));
			nProtected ++;
			// set name to R_call_param
			SET_NAMES(R_call_param, GET_NAMES(var_name));
			// make a call function
			if (VarIdx > 0)
			{
				PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
					LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
				nProtected ++;
			} else {
				PROTECT(R_fcall = LCONS(FUN,
					LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
				nProtected ++;
			}
		}


		// local selection
		TSelection &Sel = File.Push_Selection(true, false);
		Sel.ClearStructVariant();
		memset(Sel.pVariant, 0, File.VariantNum());

		C_BOOL *pBase, *pSel, *pEnd;
		pBase = pSel = Selection.pVariant;
		pEnd = pBase + File.VariantNum();

		// progress object
		CProgressStdOut progress(NumBlock, prog_flag);

		// for-loop
		for (int idx=0; idx < NumBlock; idx++)
		{
			switch (VarIdx)
			{
			case 1:  // relative
				INTEGER(R_Index)[0] = idx*bsize + 1; break;
			case 2:  // absolute
				while ((pSel < pEnd) && (*pSel == FALSE))
					pSel ++;
				INTEGER(R_Index)[0] = pSel - pBase + 1; break;
			}

			// assign sub-selection
			{
				// clear selection
				Sel.ClearSelectVariant();
				// find the first TRUE
				pSel = VEC_BOOL_FIND_TRUE(pSel, pEnd);
				Sel.varStart = pSel - pBase;
				// for-loop
				C_BOOL *pNewSel = Sel.pVariant;
				int bs = bsize;
				for (; bs > 0; bs--)
				{
					while ((pSel < pEnd) && (*pSel == FALSE))
						pSel ++;
					if (pSel < pEnd)
					{
						pNewSel[pSel - pBase] = TRUE;
						pSel ++;
					} else
						break;
				}
				Sel.varTrueNum = bsize - bs;
				Sel.varEnd = pSel - pBase;
			}

			// load data and call the user-defined function
			SEXP call_val = R_NilValue;
			if (num_var > 1)
			{
				for (int i=0; i < num_var; i++)
				{
					SET_ELEMENT(R_call_param, i,
						VarGetData(File, CHAR(STRING_ELT(var_name, i)),
						use_raw_flag));
				}
				// call R function
				call_val = eval(R_fcall, rho);

			} else {
				R_call_param = VarGetData(File, CHAR(STRING_ELT(var_name, 0)),
					use_raw_flag);
				// make a call function
				if (VarIdx > 0)
				{
					PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
						LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
				} else {
					PROTECT(R_fcall = LCONS(FUN,
						LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
				}

				// call R function
				call_val = eval(R_fcall, rho);
			}

			// store data
			switch (DatType)
			{
			case 1:  // list
				SET_ELEMENT(rv_ans, idx, call_val);
				break;
			case 2:  // connection
				if (OutputConn->text)
				{
					if (Rf_isList(call_val))
					{
						throw ErrSeqArray("the user-defined function should return a character vector.");
					} else if (!Rf_isString(call_val))
					{
						call_val = AS_CHARACTER(call_val);
					}
					size_t n = XLENGTH(call_val);
					for (size_t i=0; i < n; i++)
					{
						ConnPutText(OutputConn, "%s\n", CHAR(STRING_ELT(call_val, i)));
					}
				} else {
					if (TYPEOF(call_val) != RAWSXP)
						throw ErrSeqArray("the user-defined function should return a RAW vector.");
					size_t n = XLENGTH(call_val);
					size_t m = R_WriteConnection(OutputConn, RAW(call_val), n);
					if (n != m)
						throw ErrSeqArray("error in writing to a connection.");
				}
				break;
			case 3:  // gdsn.class
				RAppendGDS(OutputGDS, call_val);
				break;
			}

			// release R_fcall
			if (num_var <= 1) UNPROTECT(1);

			progress.Forward();
		}

		File.Pop_Selection();

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

} // extern "C"
