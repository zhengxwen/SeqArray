// ===========================================================
//
// ReadBySample.cpp: Read data sample by sample
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


static void GetFirstAndLength(C_BOOL *sel, size_t n, C_Int32 &st, C_Int32 &len)
{
	st = 0; len = 0;
	for (size_t i=0; i < n; i++)
	{
		if (sel[i]) { st = i; break; }
	}
	for (ssize_t i=n-1; i >= 0; i++)
	{
		if (sel[i]) { len = i - st + 1; break; }
	}
}


CVarApplyBySample::CVarApplyBySample()
{
	Node = NULL;
	SampleSelect = NULL;
}

void CVarApplyBySample::InitObject(TType Type, const char *Path, PdGDSObj Root,
	int nVariant, C_BOOL *VariantSel, int nSample, C_BOOL *SampleSel)
{
	static const char *ErrDim = "Invalid dimension of '%s'.";

	// initialize
	GDS_PATH_PREFIX_CHECK(Path);
	VarType = Type;
	Node = GDS_Node_Path(Root, Path, TRUE);
	SVType = GDS_Array_GetSVType(Node);
	DimCnt = GDS_Array_DimCnt(Node);

	TotalNum_Sample = nSample;
	SampleSelect = SampleSel;
	Num_Variant = NUM_OF_TRUE(VariantSel, nVariant);

	string Path2; // the path with '@'
	PdSequenceX IndexNode;  // the corresponding index variable

	switch (Type)
	{
		case ctBasic:
			// VARIABLE: sample.id
			// TODO
			if ((DimCnt != 1) || (GDS_Array_GetTotalCount(Node) != nVariant))
				throw ErrSeqArray(ErrDim, Path);
			break;

		case ctGenotype:
			// VARIABLE: genotype/~data, genotype/@data
			if (DimCnt != 3)
				throw ErrSeqArray(ErrDim, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] != nSample) || (DLen[1] < nVariant))
				throw ErrSeqArray(ErrDim, Path);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode == NULL)
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());
			if ((GDS_Array_DimCnt(IndexNode) != 1) ||
					(GDS_Array_GetTotalCount(IndexNode) != nVariant))
				throw ErrSeqArray(ErrDim, Path2.c_str());

			{
				C_Int32 I, Cnt;
				GetFirstAndLength(VariantSel, nVariant, I, Cnt);

				C_Int32 II=0, ICnt=I+Cnt;
				vector<C_Int32> ILen(ICnt);

				GDS_Array_ReadData(IndexNode, &II, &ICnt, &ILen[0], svInt32);

				VariantStart = 0;
				for (int i=0; i < I; i++) VariantStart += ILen[i];
				VariantCount = 0;
				for (int i=I; i < ICnt; i++) VariantCount += ILen[i];

				VariantSelectBuffer.resize(VariantCount);
				VariantCellCnt.resize(Num_Variant);

				C_BOOL *p = &VariantSelectBuffer[0];
				C_UInt8 *p8 = &VariantCellCnt[0];
				CellCount = 0;
				for (int i=I; i < ICnt; i++)
				{
					C_BOOL flag = VariantSel[i];
					int m = ILen[i];
					if (flag)
					{
						if ((m <= 0) || (m >255))
						{
							throw ErrSeqArray("Invalid '%s': should be 1..255.",
								Path2.c_str());
						}
						CellCount += m;
						*p8 ++ = m;
					}
					for (; m > 0; m--) *p++ = flag;
				}
			}

			CellCount *= DLen[2];
			Init.Need_GenoBuffer(CellCount);

			SelPtr[0] = NeedTRUE(1);
			SelPtr[1] = &VariantSelectBuffer[0];
			SelPtr[2] = NeedTRUE(DLen[2]);
			break;

/*		case ctPhase:
			// VARIABLE: phase/data
			if ((DimCnt != 2) && (DimCnt != 3))
				throw ErrSeqArray(ErrDim, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] != nVariant) || (DLen[1] != nSample))
				throw ErrSeqArray(ErrDim, Path);

			SelPtr[1] = SampleSel;
			if (DimCnt > 2)
			{
				Check_TRUE_ARRAY(DLen[2]);
				SelPtr[2] = &Init.TRUE_ARRAY[0];  //< ToDo: check
			}
			break;

		case ctInfo:
			// VARIABLE: info/...
			if ((DimCnt!=1) && (DimCnt!=2))
				throw ErrSeqArray(ErrDim, Path);
			GDS_Array_GetDim(Node, DLen, 2);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ErrDim, Path2.c_str());
			} else {
				if (DLen[0] != nVariant)
					throw ErrSeqArray(ErrDim, Path);
			}

			if (DimCnt > 1)
			{
				Check_TRUE_ARRAY(DLen[1]);
				SelPtr[1] = &Init.TRUE_ARRAY[0];
			}
			break;

		case ctFormat:
			// VARIABLE: format/...
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray(ErrDim, Path);
			GDS_Array_GetDim(Node, DLen, 3);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode != NULL)
			{
				if ((GDS_Array_DimCnt(IndexNode) != 1) || (GDS_Array_GetTotalCount(IndexNode) != nVariant))
					throw ErrSeqArray(ErrDim, Path2.c_str());
			} else
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());

			SelPtr[1] = SampleSel;
			if (DimCnt > 2)
			{
				Check_TRUE_ARRAY(DLen[2]);
				SelPtr[2] = &Init.TRUE_ARRAY[0];
			}
			break;
*/
		default:
			throw ErrSeqArray("Internal Error in 'CVarApplyBySample::InitObject'.");
	}

	ResetObject();
}

void CVarApplyBySample::ResetObject()
{
	CurIndex = 0;
	if (!SampleSelect[0]) NextCell();
}

bool CVarApplyBySample::NextCell()
{
	CurIndex ++;
	while ((CurIndex<TotalNum_Sample) && !SampleSelect[CurIndex])
		CurIndex ++;
	return (CurIndex < TotalNum_Sample);
}

void CVarApplyBySample::ReadGenoData(int *Base)
{
	C_Int32 st[3] = { CurIndex, VariantStart, 0 };
	C_Int32 cn[3] = { 1, VariantCount, DLen[2] };

	C_UInt8 *s = &Init.GENO_BUFFER[0];
	GDS_Array_ReadDataEx(Node, st, cn, SelPtr, s, svUInt8);

	for (int i=0; i < Num_Variant; i++)
	{
		int *p;
		int missing = 3;

		// the first 2 bits
		p = Base;
		for (int j=DLen[2]; j > 0; j--)
			*p++ = *s++;

		/// the left bits
		C_UInt8 shift = 2;
		for (int m=VariantCellCnt[i]; m > 1; m--)
		{
			p = Base;
			for (int j=DLen[2]; j > 0; j--)
			{
				*p |= int(*s) << shift;
				p ++; s ++;
			}
			shift += 2;
			missing = (missing << 2) | 0x03;
		}

		for (int j=DLen[2]; j > 0; j--)
		{
			if (*Base == missing) *Base = NA_INTEGER;
			Base ++;
		}
	}
}

void CVarApplyBySample::ReadData(SEXP Val)
{
	if (VarType == ctGenotype)
	{
		ReadGenoData(INTEGER(Val));
	} else {
/*		C_Int32 st[3] = { IndexRaw, 0, 0 };
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
*/	}
}

SEXP CVarApplyBySample::NeedRData(int &nProtected)
{
	map<size_t, SEXP>::iterator it = VarList.find(CellCount);
	if (it == VarList.end())
	{
		SEXP ans = R_NilValue, dim;
		if (COREARRAY_SV_INTEGER(SVType))
		{
			char classname[32];
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
			PROTECT(dim = NEW_INTEGER(2));
			INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Variant;
			SET_DIM(ans, dim);
			PROTECT(name_list = NEW_LIST(2));
			PROTECT(tmp = NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("allele"));
				SET_STRING_ELT(tmp, 1, mkChar("variant"));
				SET_NAMES(name_list, tmp);
			SET_DIMNAMES(ans, name_list);
			nProtected += 3;
			break;

/*		case ctPhase:
			if (DimCnt > 2)
			{
				PROTECT(dim = NEW_INTEGER(2)); nProtected ++;
				INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Sample;
				SET_DIM(ans, dim);
			}
			break;

		case ctFormat:
			if (DimCnt == 2)
			{
				PROTECT(dim = NEW_INTEGER(2)); nProtected ++;
				INTEGER(dim)[0] = Num_Sample; INTEGER(dim)[1] = NumIndexRaw;
				SET_DIM(ans, dim);
			} else if (DimCnt > 2)
			{
				PROTECT(dim = NEW_INTEGER(3)); nProtected ++;
				INTEGER(dim)[0] = DLen[2]; INTEGER(dim)[1] = Num_Sample;
				INTEGER(dim)[2] = NumIndexRaw;
				SET_DIM(ans, dim);
			}
			break;
*/
		default:
			break;
		}

		VarList.insert(pair<size_t, SEXP>(CellCount, ans));
		return ans;
	} else
		return it->second;
}



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


// ===========================================================
// Apply functions over margins on a working space
// ===========================================================

/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP sqa_Apply_Sample(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP rho)
{
	COREARRAY_TRY

		// the selection
		TInitObject::TSelection &Sel = Init.Selection(gdsfile);
		// the GDS root node
		PdGDSObj Root = GDS_R_SEXP2Obj(GetListElement(gdsfile, "root"));

		// init selection
		if (Sel.Sample.empty())
		{
			PdSequenceX N = GDS_Node_Path(Root, "sample.id", TRUE);
			int Cnt = GDS_Array_GetTotalCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'sample.id'.");
			Sel.Sample.resize(Cnt, TRUE);
		}
		if (Sel.Variant.empty())
		{
			PdSequenceX N = GDS_Node_Path(Root, "variant.id", TRUE);
			int Cnt = GDS_Array_GetTotalCount(N);
			if (Cnt < 0) throw ErrSeqArray("Invalid dimension of 'variant.id'.");
			Sel.Variant.resize(Cnt, TRUE);
		}

		// the number of calling PROTECT
		int nProtected = 0;

		// the number of selected variants
		int nSample = NUM_OF_TRUE(&Sel.Sample[0], Sel.Sample.size());
		if (nSample <= 0)
			throw ErrSeqArray("There is no selected sample.");


		// ===============================================================
		// initialize the GDS Node list

		vector<CVarApplyBySample> NodeList(Rf_length(var_name));
		vector<CVarApplyBySample>::iterator it;

		// for - loop
		for (int i=0; i < Rf_length(var_name); i++)
		{
			// the path of GDS variable
			string s = CHAR(STRING_ELT(var_name, i));
			CVarApplyBySample::TType VarType;

			if (s == "sample.id")
			{
				// =======================================================
				// sample.id
				// annotation/qual, annotation/filter
				VarType = CVarApplyBySample::ctBasic;
			} else if (s == "genotype")
			{
				VarType = CVarApplyBySample::ctGenotype;
				s.append("/~data");
			} else if (s == "phase")
			{
				// =======================================================
				// phase/
				VarType = CVarApplyBySample::ctPhase;
				s.append("/~data");
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				VarType = CVarApplyBySample::ctFormat;
				s.append("/~data");
			} else {
				throw ErrSeqArray(
					"'%s' is not a standard variable name, and the standard format:\n"
					"\tsample.id, annotation/format/VARIABLE_NAME",
					s.c_str());
			}

			NodeList[i].InitObject(VarType, s.c_str(), Root, Sel.Variant.size(),
				&Sel.Variant[0], Sel.Sample.size(), &Sel.Sample[0]);
		}

		// ===============================================================
		// as.is
		//     0: integer, 1: double, 2: character, 3: list, other: NULL
		int DatType;
		const char *as = CHAR(STRING_ELT(as_is, 0));
		if (strcmp(as, "integer") == 0)
			DatType = 0;
		else if (strcmp(as, "double") == 0)
			DatType = 1;
		else if (strcmp(as, "character") == 0)
			DatType = 2;
		else if (strcmp(as, "list") == 0)
			DatType = 3;
		else if (strcmp(as, "none") == 0)
			DatType = -1;
		else
			throw ErrSeqArray("'as.is' is not valid!");

		// init return values
		// int DatType;  //< 0: integer, 1: double, 2: character, 3: list, other: NULL
		switch (DatType)
		{
		case 0:
			PROTECT(rv_ans = NEW_INTEGER(nSample)); nProtected ++;
			break;
		case 1:
			PROTECT(rv_ans = NEW_NUMERIC(nSample)); nProtected ++;
			break;
		case 2:
			PROTECT(rv_ans = NEW_CHARACTER(nSample)); nProtected ++;
			break;
		case 3:
			PROTECT(rv_ans = NEW_LIST(nSample)); nProtected ++;
			break;
		default:
			rv_ans = R_NilValue;
		}

		// ===============================================================
		// rho
		if (!isEnvironment(rho))
			throw ErrSeqArray("'rho' should be an environment");


		// ===============================================================
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

		// ===============================================================
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
				// ToDo: optimize this
				SEXP tmp = NodeList[0].NeedRData(nProtected);
				if (tmp != R_call_param)
				{
					R_call_param = tmp;
					if (VarIdx > 1)
					{
						PROTECT(R_fcall = LCONS(FUN, LCONS(R_Index,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue)))));
					} else {
						PROTECT(R_fcall = LCONS(FUN,
							LCONS(R_call_param, LCONS(R_DotsSymbol, R_NilValue))));
					}
					nProtected ++;
				}
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
			case 0:    // integer
				INTEGER(rv_ans)[ans_index] = Rf_asInteger(val);
				break;
			case 1:    // double
				REAL(rv_ans)[ans_index] = Rf_asReal(val);
				break;
			case 2:    // character
				val = AS_CHARACTER(val);
				SET_STRING_ELT(rv_ans, ans_index,
					(LENGTH(val) > 0) ? STRING_ELT(val, 0) : NA_STRING);
				break;
			case 3:    // others
				SET_ELEMENT(rv_ans, ans_index, duplicate(val));
				break;
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
