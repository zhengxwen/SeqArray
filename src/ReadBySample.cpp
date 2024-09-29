// ===========================================================
//
// ReadBySample.cpp: Read data sample by sample
//
// Copyright (C) 2015-2024    Xiuwen Zheng
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


namespace SeqArray
{

static void GetFirstAndLength(C_BOOL *sel, size_t n, C_Int32 &st, C_Int32 &len)
{
	st = 0; len = 0;
	for (size_t i=0; i < n; i++)
	{
		if (sel[i]) { st = i; break; }
	}
	for (ssize_t i=n-1; i >= 0; i--)
	{
		if (sel[i]) { len = i - st + 1; break; }
	}
}


CVarApplyBySample::CVarApplyBySample()
{
	Node = NULL;
	MarginalSelect = NULL;
}

void CVarApplyBySample::InitObject(TVarType Type, const char *Path, PdGDSObj Root,
	int nVariant, C_BOOL *VariantSel, int nSample, C_BOOL *SampleSel,
	bool _UseRaw)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	// initialize
	GDS_PATH_PREFIX_CHECK(Path);

	if ((Type == ctBasic) || (Type == ctSampleAnnot))
	{
		Node = GDS_Node_Path(Root, Path, TRUE);
	} else {
		string s = GDS_PATH_PREFIX(Path, '~'); // the path with '~'
		Node = GDS_Node_Path(Root, s.c_str(), FALSE);
		if (Node == NULL)
		{
			throw ErrSeqArray(
				"'%s' is missing, please call 'seqOptimize(..., target="
				"\"by.sample\")' to prepare the data for scan by sample.",
				s.c_str());
		}
	}

	fVarType = Type;
	SVType = GDS_Array_GetSVType(Node);
	DimCnt = GDS_Array_DimCnt(Node);

	TotalNum_Sample = nSample;
	MarginalSelect = SampleSel;
	Num_Variant = GetNumOfTRUE(VariantSel, nVariant);
	UseRaw = _UseRaw;
	NumOfBits = GDS_Array_GetBitOf(Node);

	string Path2; // the path with '@'
	PdAbstractArray IndexNode = NULL;  // the corresponding index variable

	switch (Type)
	{
		case ctBasic:
			// ===========================================================
			// VARIABLE: sample.id
			if ((DimCnt != 1) || (GDS_Array_GetTotalCount(Node) != nSample))
				throw ErrSeqArray(ERR_DIM, Path);
			CellCount = 1;
			SelPtr[0] = NeedTRUEs(1);
			break;

		case ctGenotype:
			// ===========================================================
			// VARIABLE: genotype/~data, genotype/@data
			if (DimCnt != 3)
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] != nSample) || (DLen[1] < nVariant))
				throw ErrSeqArray(ERR_DIM, Path);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode == NULL)
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());
			if ((GDS_Array_DimCnt(IndexNode) != 1) ||
					(GDS_Array_GetTotalCount(IndexNode) != nVariant))
				throw ErrSeqArray(ERR_DIM, Path2.c_str());

			{
				C_Int32 I, Cnt;
				GetFirstAndLength(VariantSel, nVariant, I, Cnt);
				C_Int32 II=0, ICnt=I+Cnt;
				vector<C_UInt8> ILen(ICnt);

				GDS_Array_ReadData(IndexNode, &II, &ICnt, &ILen[0], svUInt8);

				VariantStart = 0;
				for (int i=0; i < I; i++) VariantStart += ILen[i];
				VariantCount = 0;
				for (int i=I; i < ICnt; i++) VariantCount += ILen[i];

				Selection.resize(VariantCount);
				GenoCellCnt.resize(Num_Variant);

				C_BOOL *p = &Selection[0];
				C_UInt8 *p8 = &GenoCellCnt[0];
				CellCount = 0;
				for (int i=I; i < ICnt; i++)
				{
					C_BOOL flag = VariantSel[i];
					int m = ILen[i];
					if (flag)
					{
						if ((m <= 0) || (m > 255))
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
			GENO_BUFFER.resize(CellCount);

			SelPtr[0] = NeedTRUEs(1);
			SelPtr[1] = &Selection[0];
			SelPtr[2] = NeedTRUEs(DLen[2]);
			break;

		case ctPhase:
			// ===========================================================
			// VARIABLE: phase/~data
			if ((DimCnt != 2) && (DimCnt != 3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);
			if ((DLen[0] != nSample) || (DLen[1] != nVariant))
				throw ErrSeqArray(ERR_DIM, Path);

			CellCount = Num_Variant;
			SelPtr[0] = NeedTRUEs(1);
			GetFirstAndLength(VariantSel, nVariant, VariantStart, VariantCount);
			SelPtr[1] = VariantSel + VariantStart;
			if (DimCnt > 2)
			{
				CellCount *= DLen[2];
				SelPtr[2] = NeedTRUEs(DLen[2]);
			}
			break;

		case ctFormat:
			// ===========================================================
			// VARIABLE: format/...
			if ((DimCnt!=2) && (DimCnt!=3))
				throw ErrSeqArray(ERR_DIM, Path);
			GDS_Array_GetDim(Node, DLen, 3);

			Path2 = GDS_PATH_PREFIX(Path, '@');
			IndexNode = GDS_Node_Path(Root, Path2.c_str(), FALSE);
			if (IndexNode == NULL)
				throw ErrSeqArray("'%s' is missing!", Path2.c_str());
			if ((GDS_Array_DimCnt(IndexNode) != 1) ||
					(GDS_Array_GetTotalCount(IndexNode) != nVariant))
				throw ErrSeqArray(ERR_DIM, Path2.c_str());

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

				Selection.resize(VariantCount);
				C_BOOL *p = &Selection[0];
				CellCount = 0;
				for (int i=I; i < ICnt; i++)
				{
					C_BOOL flag = VariantSel[i];
					int m = ILen[i];
					if (flag)
					{
						if (m < 0)
						{
							throw ErrSeqArray("Invalid '%s': should be >= 0.",
								Path2.c_str());
						}
						CellCount += m;
					}
					for (; m > 0; m--) *p++ = flag;
				}
			}

			SelPtr[0] = NeedTRUEs(1);
			SelPtr[1] = &Selection[0];
			if (DimCnt > 2)
			{
				SelPtr[2] = NeedTRUEs(DLen[2]);
				CellCount *= DLen[2];
			}
			break;

		case ctSampleAnnot:
			// ===========================================================
			// VARIABLE: sample.annotation/
			if ((DimCnt != 1) && (DimCnt != 2))
				throw ErrSeqArray(ERR_DIM, Path);
			DLen[0] = DLen[1] = DLen[2] = 0;
			GDS_Array_GetDim(Node, DLen, 2);
			if (DLen[0] != nSample)
				throw ErrSeqArray(ERR_DIM, Path);

			VariantStart = 0; VariantCount = DLen[1];
			CellCount = 1;
			SelPtr[0] = NeedTRUEs(1);
			if (DimCnt > 1)
			{
				CellCount *= DLen[1];
				SelPtr[1] = NeedTRUEs(DLen[1]);
			}
			break;

		default:
			throw ErrSeqArray("Internal Error in 'CVarApplyBySample::InitObject'.");
	}

	Reset();
}

void CVarApplyBySample::Reset()
{
	Position = 0;
	if (!MarginalSelect[0]) Next();
}

bool CVarApplyBySample::Next()
{
	Position ++;
	while ((Position<TotalNum_Sample) && !MarginalSelect[Position])
		Position ++;
	return (Position < TotalNum_Sample);
}

void CVarApplyBySample::ReadGenoData(int *Base)
{
	C_Int32 st[3] = { Position, VariantStart, 0 };
	C_Int32 cn[3] = { 1, VariantCount, DLen[2] };
	C_UInt8 *s = &GENO_BUFFER[0];

	GDS_Array_ReadDataEx(Node, st, cn, SelPtr, s, svUInt8);
	const int bit_mask = ~((-1) << NumOfBits);

	for (int i=0; i < Num_Variant; i++)
	{
		int *p;
		int missing = bit_mask;

		// the first 2 bits
		p = Base;
		for (int j=DLen[2]; j > 0; j--)
			*p++ = *s++;

		/// the left bits
		C_UInt8 shift = NumOfBits;
		for (int m=GenoCellCnt[i]; m > 1; m--)
		{
			p = Base;
			for (int j=DLen[2]; j > 0; j--)
			{
				*p |= int(*s) << shift;
				p ++; s ++;
			}
			shift += NumOfBits;
			missing = (missing << NumOfBits) | bit_mask;
		}

		for (int j=DLen[2]; j > 0; j--)
		{
			if (*Base == missing) *Base = NA_INTEGER;
			Base ++;
		}
	}
}

void CVarApplyBySample::ReadGenoData(C_UInt8 *Base)
{
	C_Int32 st[3] = { Position, VariantStart, 0 };
	C_Int32 cn[3] = { 1, VariantCount, DLen[2] };
	C_UInt8 *s = &GENO_BUFFER[0];

	GDS_Array_ReadDataEx(Node, st, cn, SelPtr, s, svUInt8);
	const int bit_mask = ~((-1) << NumOfBits);

	for (int i=0; i < Num_Variant; i++)
	{
		C_UInt8 *p;
		int missing = bit_mask;

		// the first 2 bits
		p = Base;
		for (int j=DLen[2]; j > 0; j--)
			*p++ = *s++;

		if (NumOfBits == 2)
		{
			if (GenoCellCnt[i] > 4)
				Rf_warning("RAW type may not be sufficient to store genotypes.");
		} else {
			if (GenoCellCnt[i] > 1)
				Rf_warning("RAW type may not be sufficient to store genotypes.");
		}

		/// the left bits
		C_UInt8 shift = NumOfBits;
		for (int m=GenoCellCnt[i]; m > 1; m--)
		{
			p = Base;
			for (int j=DLen[2]; j > 0; j--)
			{
				*p |= int(*s) << shift;
				p ++; s ++;
			}
			shift += NumOfBits;
			missing = (missing << NumOfBits) | bit_mask;
		}

		for (int j=DLen[2]; j > 0; j--)
		{
			if (*Base == missing) *Base = NA_RAW;
			Base ++;
		}
	}
}

void CVarApplyBySample::ReadData(SEXP Val)
{
	if (fVarType == ctGenotype)
	{
		if (UseRaw)
			ReadGenoData(RAW(Val));
		else
			ReadGenoData(INTEGER(Val));
	} else {
		C_Int32 st[3] = { Position, VariantStart, 0 };
		C_Int32 cn[3] = { 1, VariantCount, DLen[2] };

		if (COREARRAY_SV_INTEGER(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, cn, SelPtr, INTEGER(Val), svInt32);
		} else if (COREARRAY_SV_FLOAT(SVType))
		{
			GDS_Array_ReadDataEx(Node, st, cn, SelPtr, REAL(Val), svFloat64);
		} else if (COREARRAY_SV_STRING(SVType))
		{
			vector<string> buffer(CellCount);
			GDS_Array_ReadDataEx(Node, st, cn, SelPtr, &buffer[0], svStrUTF8);
			for (size_t i=0; i < buffer.size(); i++)
				SET_STRING_ELT(Val, i, Rf_mkChar(buffer[i].c_str()));
		}
	}
}

SEXP CVarApplyBySample::NeedRData(int &nProtected)
{
	map<size_t, SEXP>::iterator it = VarList.find(CellCount);
	if (it == VarList.end())
	{
		SEXP ans = R_NilValue, dim;
		if (COREARRAY_SV_INTEGER(SVType))
		{
			if (fVarType == ctGenotype)
			{
				if (UseRaw)
					PROTECT(ans = NEW_RAW(CellCount));
				else
					PROTECT(ans = NEW_INTEGER(CellCount));
				nProtected ++;
			} else {
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

		int *p;
		switch (fVarType)
		{
		case ctGenotype:
			p = INTEGER(dim = NEW_INTEGER(2));
			p[0] = DLen[2]; p[1] = Num_Variant;
			SET_DIM(ans, dim);
			break;

		case ctPhase:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = DLen[2]; p[1] = Num_Variant;
				SET_DIM(ans, dim);
			}
			break;

		case ctFormat:
			if (DimCnt > 2)  // DimCnt = 2 or 3 only
			{
				p = INTEGER(dim = NEW_INTEGER(2));
				p[0] = DLen[2]; p[1] = CellCount / DLen[2];
				SET_DIM(ans, dim);
			}
			break;

		default:
			break;
		}

		VarList.insert(pair<size_t, SEXP>(CellCount, ans));
		return ans;
	} else
		return it->second;
}

}


extern "C"
{
using namespace SeqArray;

// ===========================================================
// Apply functions over margins on a working space
// ===========================================================

COREARRAY_DLL_LOCAL extern const char *Txt_Apply_AsIs[];
COREARRAY_DLL_LOCAL extern const char *Txt_Apply_VarIdx[];


/// Apply functions over margins on a working space
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Sample(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP use_raw, SEXP rho)
{
	int use_raw_flag = Rf_asLogical(use_raw);
	if (use_raw_flag == NA_LOGICAL)
		Rf_error("'.useraw' must be TRUE or FALSE.");

	COREARRAY_TRY

		// the selection
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);

		// the number of calling PROTECT
		int nProtected = 0;

		// the number of selected variants
		int nSample = File.SampleSelNum();
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
			CVarApplyBySample::TVarType VarType;

			if (s == "sample.id")
			{
				VarType = CVarApplyBySample::ctBasic;
			} else if (s == "genotype")
			{
				VarType = CVarApplyBySample::ctGenotype;
				s.append("/data");
			} else if (s == "phase")
			{
				VarType = CVarApplyBySample::ctPhase;
				s.append("/data");
			} else if (strncmp(s.c_str(), "annotation/format/", 18) == 0)
			{
				VarType = CVarApplyBySample::ctFormat;
				s.append("/data");
			} else if (strncmp(s.c_str(), "sample.annotation/", 18) == 0)
			{
				VarType = CVarApplyBySample::ctSampleAnnot;
			} else {
				throw ErrSeqArray(
					"'%s' is not a standard variable name, and the standard format:\n"
					"    sample.id, genotype, annotation/format/VARIABLE_NAME\n"
					"    sample.annotation/VARIABLE_NAME",
					s.c_str());
			}

			NodeList[i].InitObject(VarType, s.c_str(), Root, File.VariantNum(),
				Sel.pVariant, File.SampleNum(), Sel.pSample,
				use_raw_flag != FALSE);
		}

		// ===============================================================
		// as.is
		int DatType = MatchText(CHAR(STRING_ELT(as_is, 0)), Txt_Apply_AsIs);
		if (DatType < 0)
			throw ErrSeqArray("'as.is' is not valid!");
		switch (DatType)
		{
		case 1:
			PROTECT(rv_ans = NEW_LIST(nSample)); nProtected ++; break;
		case 2:
			PROTECT(rv_ans = NEW_INTEGER(nSample)); nProtected ++; break;
		case 3:
			PROTECT(rv_ans = NEW_NUMERIC(nSample)); nProtected ++; break;
		case 4:
			PROTECT(rv_ans = NEW_CHARACTER(nSample)); nProtected ++; break;
		case 5:
			PROTECT(rv_ans = NEW_LOGICAL(nSample)); nProtected ++; break;
		case 6:
			PROTECT(rv_ans = NEW_RAW(nSample)); nProtected ++; break;
		default:
			rv_ans = R_NilValue;
		}

		// ===============================================================
		// rho
		if (!Rf_isEnvironment(rho))
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

		// ===============================================================
		// var.index
		int VarIdx = MatchText(CHAR(STRING_ELT(var_index, 0)), Txt_Apply_VarIdx);
		if (VarIdx < 0)
			throw ErrSeqArray("'var.index' is not valid!");

		SEXP R_fcall;
		SEXP R_Index = NULL;
		if (VarIdx > 0)
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
			case 1:  // relative
				INTEGER(R_Index)[0] = ans_index + 1; break;
			case 2:
				INTEGER(R_Index)[0] = NodeList.begin()->Position + 1; break;
			}

			if (NodeList.size() <= 1)
			{
				SEXP pm = NodeList[0].NeedRData(nProtected);
				if (pm != R_call_param)
				{
					if (VarIdx > 0)
						SETCADDR(R_fcall, pm);
					else
						SETCADR(R_fcall, pm);
					R_call_param = pm;
				}
				NodeList[0].ReadData(pm);
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
			SEXP val = Rf_eval(R_fcall, rho);
			switch (DatType)
			{
			case 1:
				SET_ELEMENT(rv_ans, ans_index, Rf_duplicate(val)); break;
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
				if (!it->Next())
					{ ifend = true; break; }
			}

		} while (!ifend);

		// finally
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

} // extern "C"
