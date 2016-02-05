// ===========================================================
//
// ReadByVariant.cpp: Read data variant by variant
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

#include "ReadByVariant.h"


// ===================================================================== //

/// 
CVarApplyByVariant::CVarApplyByVariant()
{
	Node = IndexNode = NULL;
	VariantSelect = NULL;
	UseRaw = false;
}

void CVarApplyByVariant::InitObject(TType Type, const char *Path,
	PdGDSObj Root, int nVariant, C_BOOL *VariantSel, int nSample,
	C_BOOL *SampleSel, bool _UseRaw)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	// initialize
	GDS_PATH_PREFIX_CHECK(Path);
	VarType = Type;
	Node = GDS_Node_Path(Root, Path, TRUE);
	SVType = GDS_Array_GetSVType(Node);
	DimCnt = GDS_Array_DimCnt(Node);

	TotalNum_Variant = nVariant;
	VariantSelect = VariantSel;
	Num_Sample = GetNumOfTRUE(SampleSel, nSample);
	UseRaw = _UseRaw;
	NumOfBits = GDS_Array_GetBitOf(Node);

	string Path2; // the path with '@'
	switch (Type)
	{
		case ctBasic:
			// ===========================================================
			// VARIABLE: variant.id, position, allele
			if ((DimCnt != 1) || (GDS_Array_GetTotalCount(Node) != nVariant))
				throw ErrSeqArray(ERR_DIM, Path);
			break;

		case ctGenotype:
			// ===========================================================
			// VARIABLE: genotype/data, genotype/@data
			if (DimCnt != 3)
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] < nVariant) || (DLen[1] != nSample))
				throw ErrSeqArray(ERR_DIM, Path);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode == NULL)
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());
			if ((GDS_Array_DimCnt(IndexNode) != 1) ||
					(GDS_Array_GetTotalCount(IndexNode) != nVariant))
				throw ErrSeqArray(ERR_DIM, Path2.c_str());

			CellCount = Num_Sample * DLen[2];
			Init.Need_GenoBuffer(CellCount);
			{
				Selection.resize(DLen[1] * DLen[2]);
				C_BOOL *p = SelPtr[1] = &Selection[0];
				memset(p, TRUE, Selection.size());
				C_BOOL *s = SampleSel;
				for (int n=DLen[1]; n > 0; n--)
				{
					if (*s++ == FALSE)
					{
						for (int m=DLen[2]; m > 0; m--)
							*p ++ = FALSE;
					} else
						p += DLen[2];
				}
			}
			break;

		case ctPhase:
			// ===========================================================
			// VARIABLE: phase/data
			if ((DimCnt != 2) && (DimCnt != 3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] != nVariant) || (DLen[1] != nSample))
				throw ErrSeqArray(ERR_DIM, Path);

			SelPtr[1] = SampleSel;
			if (DimCnt > 2)
				SelPtr[2] = NeedTRUE(DLen[2]);
			break;

		case ctInfo:
			// ===========================================================
			// VARIABLE: info/...
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 2);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ERR_DIM, Path2.c_str());
			} else {
				if (DLen[0] != nVariant)
					throw ErrSeqArray(ERR_DIM, Path);
			}

			if (DimCnt > 1)
				SelPtr[1] = NeedTRUE(DLen[1]);
			break;

		case ctFormat:
			// ===========================================================
			// VARIABLE: format/...
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ERR_DIM, Path2.c_str());
			} else
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());

			SelPtr[1] = SampleSel;
			if (DimCnt > 2)
				SelPtr[2] = NeedTRUE(DLen[2]);
			break;

		default:
			throw ErrSeqArray("Internal Error in 'CVarApplyByVariant::InitObject'.");
	}

	ResetObject();
}

void CVarApplyByVariant::ResetObject()
{
	CurIndex = 0;
	IndexRaw = 0;
	if (IndexNode)
	{
		C_Int32 Cnt=1;
		GDS_Array_ReadData(IndexNode, &CurIndex, &Cnt, &NumIndexRaw, svInt32);
		if (NumIndexRaw < 0) NumIndexRaw = 0;
	} else
		NumIndexRaw = 1;

	if (!VariantSelect[0])
		NextCell();
}

bool CVarApplyByVariant::NextCell()
{
	CurIndex ++;

	if (IndexNode)
	{
		IndexRaw += NumIndexRaw;
		C_Int32 Cnt=1, L;
		while ((CurIndex<TotalNum_Variant) && !VariantSelect[CurIndex])
		{
			GDS_Array_ReadData(IndexNode, &CurIndex, &Cnt, &L, svInt32);
			if (L > 0)
				IndexRaw += L;
			CurIndex ++;
		}
		if (CurIndex < TotalNum_Variant)
		{
			GDS_Array_ReadData(IndexNode, &CurIndex, &Cnt, &NumIndexRaw,
				svInt32);
			if (NumIndexRaw < 0) NumIndexRaw = 0;
		} else
			NumIndexRaw = 0;
	} else {
		while ((CurIndex<TotalNum_Variant) && !VariantSelect[CurIndex])
			CurIndex ++;
		IndexRaw = CurIndex;
		NumIndexRaw = 1;
	}

	return (CurIndex < TotalNum_Variant);
}

void CVarApplyByVariant::ReadGenoData(int *Base)
{
	// the size of Init.GENO_BUFFER has been checked in 'Init()'
	const ssize_t SlideCnt = ssize_t(DLen[1]) * ssize_t(DLen[2]);

	// NumIndexRaw always >= 1
	CdIterator it;
	GDS_Iter_Position(Node, &it, C_Int64(IndexRaw)*SlideCnt);
	GDS_Iter_RDataEx(&it, Base, SlideCnt, svInt32, SelPtr[1]);

	const int bit_mask = ~((-1) << NumOfBits);
	int missing = bit_mask;
	for (int idx=1; idx < NumIndexRaw; idx ++)
	{
		GDS_Iter_Position(Node, &it, (C_Int64(IndexRaw) + idx)*SlideCnt);
		GDS_Iter_RDataEx(&it, &Init.GENO_BUFFER[0], SlideCnt, svUInt8, SelPtr[1]);

		int shift = idx * NumOfBits;
		C_UInt8 *s = &Init.GENO_BUFFER[0];
		int *p = Base;
		for (int n=Num_Sample; n > 0 ; n--)
		{
			for (int m=DLen[2]; m > 0 ; m--)
				*p++ |= int(*s++) << shift;
		}

		missing = (missing << NumOfBits) | bit_mask;
	}

	// CellCount = Num_Sample * DLen[2] in 'NeedRData'
	for (size_t n=CellCount; n > 0; n--)
	{
		if (*Base == missing) *Base = NA_INTEGER;
		Base ++;
	}	
}

void CVarApplyByVariant::ReadGenoData(C_UInt8 *Base)
{
	// the size of Init.GENO_BUFFER has been checked in 'Init()'
	const ssize_t SlideCnt = ssize_t(DLen[1]) * ssize_t(DLen[2]);

	// NumIndexRaw always >= 1
	CdIterator it;
	GDS_Iter_Position(Node, &it, C_Int64(IndexRaw)*SlideCnt);
	GDS_Iter_RDataEx(&it, Base, SlideCnt, svUInt8, SelPtr[1]);

	const C_UInt8 bit_mask = ~((-1) << NumOfBits);
	C_UInt8 missing = bit_mask;
	int MyNumIndexRaw = NumIndexRaw;
	switch (NumOfBits)
	{
	case 2:
		if (NumIndexRaw > 4)
			warning("RAW type may not be sufficient to store genotypes.");
		break;
	case 4:
		if (NumIndexRaw > 2)
			warning("RAW type may not be sufficient to store genotypes.");
		break;
	case 8:
		MyNumIndexRaw = 1;
		if (NumIndexRaw > 1)
			warning("RAW type may not be sufficient to store genotypes.");
		break;
	}

	for (int idx=1; idx < MyNumIndexRaw; idx ++)
	{
		GDS_Iter_Position(Node, &it, (C_Int64(IndexRaw) + idx)*SlideCnt);
		GDS_Iter_RDataEx(&it, &Init.GENO_BUFFER[0], SlideCnt, svUInt8, SelPtr[1]);

		C_UInt8 shift = idx * NumOfBits;
		C_UInt8 *s = &Init.GENO_BUFFER[0];
		C_UInt8 *p = Base;
		for (int n=Num_Sample; n > 0 ; n--)
		{
			for (int m=DLen[2]; m > 0 ; m--)
				*p++ |= (*s++) << shift;
		}

		missing = (missing << NumOfBits) | bit_mask;
	}

	// CellCount = Num_Sample * DLen[2] in 'NeedRData'
	for (size_t n=CellCount; n > 0; n--)
	{
		if (*Base == missing) *Base = NA_RAW;
		Base ++;
	}	
}

void CVarApplyByVariant::ReadData(SEXP Val)
{
	if (NumIndexRaw <= 0) return;
	if (VarType == ctGenotype)
	{
		if (UseRaw)
			ReadGenoData(RAW(Val));
		else
			ReadGenoData(INTEGER(Val));
	} else {
		C_Int32 st[3] = { IndexRaw, 0, 0 };
		DLen[0] = NumIndexRaw;
		SelPtr[0] = NeedTRUE(NumIndexRaw);

		if (COREARRAY_SV_INTEGER(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, INTEGER(Val), svInt32);
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, REAL(Val), svFloat64);
		} else if (COREARRAY_SV_STRING(SVType))
		{
			vector<string> buffer(CellCount);
			GDS_Array_ReadDataEx(Node, st, DLen, SelPtr, &buffer[0], svStrUTF8);
			for (int i=0; i < (int)buffer.size(); i++)
				SET_STRING_ELT(Val, i, mkChar(buffer[i].c_str()));
		}
	}
}

SEXP CVarApplyByVariant::NeedRData(int &nProtected)
{
	if (NumIndexRaw <= 0) return R_NilValue;

	map<size_t, SEXP>::iterator it = VarList.find(NumIndexRaw);
	if (it == VarList.end())
	{
		switch (VarType)
		{
		case ctBasic:
			CellCount = 1; break;
		case ctGenotype:
			CellCount = Num_Sample * DLen[2]; break;
		case ctPhase:
			CellCount = (DimCnt>2) ? Num_Sample*DLen[2] : Num_Sample;
			break;
		case ctInfo:
			CellCount = ((DimCnt>1) ? DLen[1] : 1) * NumIndexRaw;
			break;
		case ctFormat:
			CellCount = ((DimCnt>2) ? Num_Sample*DLen[2] : Num_Sample) *
						NumIndexRaw;
			break;
		default:
			CellCount = 0;
		}

		SEXP ans = R_NilValue, dim;
		if (COREARRAY_SV_INTEGER(SVType))
		{
			if (VarType == ctGenotype)
			{
				if (UseRaw)
					PROTECT(ans = NEW_RAW(CellCount));
				else
					PROTECT(ans = NEW_INTEGER(CellCount));
				nProtected ++;
			} else {
				char classname[64];
				classname[0] = 0;
				GDS_Node_GetClassName(Node, classname, sizeof(classname));
				if (strcmp(classname, "dBit1") == 0)
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else if (GDS_R_Is_Logical(Node))
				{
					PROTECT(ans = NEW_LOGICAL(CellCount));
				} else {
					PROTECT(ans = NEW_INTEGER(CellCount));
					nProtected += GDS_R_Set_IfFactor(Node, ans);
				}
				nProtected ++;
			}
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			PROTECT(ans = NEW_NUMERIC(CellCount));
			nProtected ++;
		} else if (COREARRAY_SV_STRING(SVType))
		{
			PROTECT(ans = NEW_CHARACTER(CellCount));
			nProtected ++;
		}

		SEXP name_list, tmp;
		switch (VarType)
		{
		case ctGenotype:
			{
				int *p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = DLen[2]; p[1] = Num_Sample;
				SET_DIM(ans, dim);
				PROTECT(name_list = NEW_LIST(2));
				PROTECT(tmp = NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("allele"));
				SET_STRING_ELT(tmp, 1, mkChar("sample"));
				SET_NAMES(name_list, tmp);
				SET_DIMNAMES(ans, name_list);
				UNPROTECT(2);
			}
			break;

		case ctPhase:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				int *p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = DLen[2]; p[1] = Num_Sample;
				SET_DIM(ans, dim);
			}
			break;

		case ctFormat:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				int *p = INTEGER(dim = NEW_INTEGER(3));
				p[0] = DLen[2]; p[1] = Num_Sample; p[2] = NumIndexRaw;
				SET_DIM(ans, dim);
			} else if (DimCnt == 2)
			{
				int *p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = Num_Sample; p[1] = NumIndexRaw;
				SET_DIM(ans, dim);
			}
			break;

		default:
			break;
		}

		VarList.insert(pair<int, SEXP>(NumIndexRaw, ans));
		return ans;
	} else
		return it->second;
}



extern "C"
{
// ===========================================================
// Apply functions over margins on a working space
// ===========================================================

/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Variant(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP use_raw, SEXP list_duplicate,
	SEXP rho)
{
	int use_raw_flag = Rf_asLogical(use_raw);
	if (use_raw_flag == NA_LOGICAL)
		error("'.useraw' must be TRUE or FALSE.");

	int dup_flag = Rf_asLogical(list_duplicate);
	if (dup_flag == NA_LOGICAL)
		error("'.duplicate' must be TRUE or FALSE.");

	COREARRAY_TRY

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);

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

		// the number of calling PROTECT
		int nProtected = 0;

		// the number of selected variants
		int nVariant = GetNumOfTRUE(&Sel.Variant[0], Sel.Variant.size());
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");


		// ===========================================================
		// initialize the GDS Node list

		vector<CVarApplyByVariant> NodeList(Rf_length(var_name));
		vector<CVarApplyByVariant>::iterator it;

		// for - loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));
			CVarApplyByVariant::TType VarType;

			if ( s=="variant.id" || s=="position" || s=="chromosome" ||
				s=="allele" || s=="annotation/id" || s=="annotation/qual" ||
				s=="annotation/filter" )
			{
				// =======================================================
				// variant.id, position, chromosome, allele, annotation/id
				// annotation/qual, annotation/filter
				VarType = CVarApplyByVariant::ctBasic;
			} else if (s == "genotype")
			{
				VarType = CVarApplyByVariant::ctGenotype;
				s.append("/data");
			} else if (s == "phase")
			{
				// =======================================================
				// phase/
				VarType = CVarApplyByVariant::ctPhase;
				s.append("/data");
			} else if (strncmp(s.c_str(), "annotation/info/", 16) == 0)
			{
				VarType = CVarApplyByVariant::ctInfo;
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				VarType = CVarApplyByVariant::ctFormat;
				s.append("/data");
			} else {
				throw ErrSeqArray(
					"'%s' is not a standard variable name, and the standard format:\n"
					"\tvariant.id, position, chromosome, allele, "
					"annotation/id, annotation/qual, annotation/filter\n"
					"\tannotation/info/VARIABLE_NAME', annotation/format/VARIABLE_NAME",
					s.c_str());
			}

			NodeList[i].InitObject(VarType, s.c_str(), Root, Sel.Variant.size(),
				&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0],
				use_raw_flag != FALSE);
		}

		// ===========================================================
		// as.is
		static const char *AsList[] =
		{
			"none", "list", "integer", "double", "character", "logical", "raw"
		};
		int DatType = MatchElement(CHAR(STRING_ELT(as_is, 0)), AsList, 7);
		if (DatType < 0)
			throw ErrSeqArray("'as.is' is not valid!");
		switch (DatType)
		{
		case 1:
			PROTECT(rv_ans = NEW_LIST(nVariant)); nProtected ++; break;
		case 2:
			PROTECT(rv_ans = NEW_INTEGER(nVariant)); nProtected ++; break;
		case 3:
			PROTECT(rv_ans = NEW_NUMERIC(nVariant)); nProtected ++; break;
		case 4:
			PROTECT(rv_ans = NEW_CHARACTER(nVariant)); nProtected ++; break;
		case 5:
			PROTECT(rv_ans = NEW_LOGICAL(nVariant)); nProtected ++; break;
		case 6:
			PROTECT(rv_ans = NEW_RAW(nVariant)); nProtected ++; break;
		default:
			rv_ans = R_NilValue;
		}

		// ===========================================================
		// rho
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");


		// ===========================================================
		// initialize calling

		SEXP R_call_param = R_NilValue;
		if (NodeList.size() > 1)
		{
			PROTECT(R_call_param = NEW_LIST(NodeList.size()));
			nProtected ++;
			// set name to R_call_param
			SET_NAMES(R_call_param, GET_NAMES(var_name));
		}

		// 1 -- none, 2 -- relative, 3 -- absolute
		int VarIdx = INTEGER(var_index)[0];

		SEXP R_fcall;
		SEXP R_Index = NULL;
		if (VarIdx > 1)
		{
			PROTECT(R_Index = NEW_INTEGER(1));
			nProtected ++;
			PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
			nProtected ++;
		} else {
			PROTECT(R_fcall = LCONS(FUN,
				LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
			nProtected ++;
		}

		map<SEXP, SEXP> R_fcall_map;
		R_fcall_map[R_call_param] = R_fcall;


		// ===========================================================
		// for-loop calling

		bool ifend = false;
		int ans_index = 0;
		do {
			switch (VarIdx)
			{
				case 2:
					INTEGER(R_Index)[0] = ans_index + 1;
					break;
				case 3:
					INTEGER(R_Index)[0] = NodeList.begin()->CurIndex + 1;
					break;
			}

			if (NodeList.size() <= 1)
			{
				R_call_param = NodeList[0].NeedRData(nProtected);
				map<SEXP, SEXP>::iterator it = R_fcall_map.find(R_call_param);
				if (it == R_fcall_map.end())
				{
					if (VarIdx > 1)
					{
						PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
					} else {
						PROTECT(R_fcall = LCONS(FUN,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
					}
					nProtected ++;
					R_fcall_map[R_call_param] = R_fcall;
				} else
					R_fcall = it->second;

				NodeList[0].ReadData(R_call_param);

			} else {
				int idx = 0;
				for (it=NodeList.begin(); it != NodeList.end(); it ++)
				{
					SEXP tmp = it->NeedRData(nProtected);
					it->ReadData(tmp);
					SET_ELEMENT(R_call_param, idx, tmp);
					idx ++;
				}
			}

			// call R function
			SEXP val = eval(R_fcall, rho);
			switch (DatType)
			{
			case 1:
				if (dup_flag) val = duplicate(val);
				SET_ELEMENT(rv_ans, ans_index, val);
				break;
			case 2:
				INTEGER(rv_ans)[ans_index] = Rf_asInteger(val); break;
			case 3:
				REAL(rv_ans)[ans_index] = Rf_asReal(val); break;
			case 4:
				SET_STRING_ELT(rv_ans, ans_index, Rf_asChar(val)); break;
			case 5:
				LOGICAL(rv_ans)[ans_index] = Rf_asLogical(val); break;
			case 6:
				RAW(rv_ans)[ans_index] = Rf_asInteger(val); break;
			}
			ans_index ++;

			// check the end
			for (it=NodeList.begin(); it != NodeList.end(); it ++)
			{
				if (!it->NextCell())
					{ ifend = true; break; }
			}

		} while (!ifend);

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

} // extern "C"
