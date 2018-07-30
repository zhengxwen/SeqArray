// ===========================================================
//
// ConvToGDS.cpp: format conversion
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
#include <vector>

using namespace std;
using namespace SeqArray;


extern "C"
{
// ======================================================================
// PLINK BED --> SeqArray GDS
// ======================================================================

/// to detect PLINK BED
COREARRAY_DLL_EXPORT SEXP SEQ_ConvBEDFlag(SEXP File, SEXP ReadBinFun, SEXP Rho)
{
	// 'readBin(File, raw(), 3)'
	SEXP R_Read_Call = PROTECT(
		LCONS(ReadBinFun, LCONS(File,
		LCONS(NEW_RAW(0), LCONS(ScalarInteger(3), R_NilValue)))));

	// call ...
	SEXP val = PROTECT(eval(R_Read_Call, Rho));
	unsigned char *prefix = RAW(val);

	if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
		error("Invalid prefix in the bed file.");

	UNPROTECT(2);
	return ScalarInteger((C_UInt8)prefix[2]);
}


/// to convert from PLINK BED to GDS
COREARRAY_DLL_EXPORT SEXP SEQ_ConvBED2GDS(SEXP GenoNode, SEXP Num, SEXP File,
	SEXP ReadBinFun, SEXP Rho)
{
	COREARRAY_TRY

		PdAbstractArray Mat = GDS_R_SEXP2Obj(GenoNode, FALSE);
		int n = Rf_asInteger(Num);
		int DLen[3];
		GDS_Array_GetDim(Mat, DLen, 3);

		int nGeno = DLen[1]*2;
		int nRe = DLen[1] % 4;
		int nRe4 = DLen[1] / 4;
		int nPack = (nRe > 0) ? (nRe4 + 1) : nRe4;

		// 'readBin(File, raw(), 3)'
		SEXP R_Read_Call = PROTECT(
			LCONS(ReadBinFun, LCONS(File,
			LCONS(NEW_RAW(0), LCONS(ScalarInteger(nPack), R_NilValue)))));

		vector<C_UInt8> dstgeno(nGeno);
		static const C_UInt8 cvt1[4] = { 0, 3, 1, 1 };
		static const C_UInt8 cvt2[4] = { 0, 3, 0, 1 };

		for (int i=0; i < n; i++)
		{
			// read genotypes
			SEXP val = eval(R_Read_Call, Rho);
			unsigned char *srcgeno = RAW(val);

			// unpacked
			C_UInt8 *p = &dstgeno[0];
			for (int k=0; k < nRe4; k++)
			{
				C_UInt8 g = srcgeno[k];
				p[0] = cvt1[g & 0x03]; p[1] = cvt2[g & 0x03];
				g >>= 2; p += 2;
				p[0] = cvt1[g & 0x03]; p[1] = cvt2[g & 0x03];
				g >>= 2; p += 2;
				p[0] = cvt1[g & 0x03]; p[1] = cvt2[g & 0x03];
				g >>= 2; p += 2;
				p[0] = cvt1[g & 0x03]; p[1] = cvt2[g & 0x03];
				g >>= 2; p += 2;
			}
			if (nRe > 0)
			{
				C_UInt8 g = srcgeno[nRe4];
				for (int k=0; k < nRe; k++)
				{
					p[0] = cvt1[g & 0x03]; p[1] = cvt2[g & 0x03];
					g >>= 2; p += 2;
				}
			}

			// append
			GDS_Array_AppendData(Mat, nGeno, &dstgeno[0], svUInt8);
		}

		UNPROTECT(1);

	COREARRAY_CATCH
}


// ======================================================================
// SNP GDS --> SeqArray GDS
// ======================================================================

static bool UseMajorAsRef = false;

COREARRAY_DLL_EXPORT SEXP FC_SNP2GDS_Ref(SEXP MajorRef)
{
	UseMajorAsRef = (Rf_asLogical(MajorRef) == TRUE);
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP FC_SNP2GDS(SEXP X)
{
	SEXP Geno = VECTOR_ELT(X, 0);
	size_t n  = Rf_length(Geno);

	SEXP Allele = VECTOR_ELT(X, 1);
	int sign_pos = -1;
	const char *base = CHAR(STRING_ELT(Allele, 0));
	{
		const char *p = base;
		for (; *p != 0; p++)
		{
			if (*p == '/')  // format A/B
				{ sign_pos = p - base; break; }
		}
	}

	// check allele
	bool rev_flag = false;
	if (UseMajorAsRef && sign_pos>=0)
	{
		int *s=INTEGER(Geno), sum=0, nvalid=0;
		for (size_t i=0; i < n; i++, s++)
		{
			if (0 <= *s && *s <= 2)
			{
				nvalid ++;
				sum += *s;
			}
		}
		rev_flag = (sum < nvalid);
	}

	SEXP rv_ans = PROTECT(NEW_LIST(2));
	SEXP Dest = NEW_INTEGER(2*n);
	SET_ELEMENT(rv_ans, 0, Dest);
	SET_ELEMENT(rv_ans, 1, Allele);

	int *s = INTEGER(Geno), *p = INTEGER(Dest);

	if (rev_flag)
	{
		string ss(strlen(base), 0);
		size_t nn = strlen(base+sign_pos+1);
		memcpy(&ss[0], base+sign_pos+1, nn);
		ss[nn] = ',';
		memcpy(&ss[nn+1], base, sign_pos);
		memcpy((void*)base, &ss[0], ss.size());
		for (; (n--) > 0; p+=2)
		{
			switch (*s++)
			{
				case 0: p[0] = p[1] = 0; break;
				case 1: p[0] = 0; p[1] = 1; break;
				case 2: p[0] = p[1] = 1; break;
				default: p[0] = p[1] = -1;
			}
		}
	} else {
		if (sign_pos >= 0)
			((char*)base)[sign_pos] = ',';
		for (; (n--) > 0; p+=2)
		{
			switch (*s++)
			{
				case 0: p[0] = p[1] = 1; break;
				case 1: p[0] = 0; p[1] = 1; break;
				case 2: p[0] = p[1] = 0; break;
				default: p[0] = p[1] = -1;
			}
		}
	}

	UNPROTECT(1);
	return rv_ans;
}


COREARRAY_DLL_EXPORT SEXP FC_Dosage2GDS(SEXP X)
{
	SEXP Geno = VECTOR_ELT(X, 0);
	size_t n  = Rf_length(Geno);

	SEXP Allele = VECTOR_ELT(X, 1);
	int sign_pos = -1;
	const char *base = CHAR(STRING_ELT(Allele, 0));
	{
		const char *p = base;
		for (; *p != 0; p++)
		{
			if (*p == '/')  // format A/B
				{ sign_pos = p - base; break; }
		}
	}

	// check allele
	bool rev_flag = false;
	if (UseMajorAsRef && sign_pos>=0)
	{
		double *s=REAL(Geno), sum=0;
		int nvalid=0;
		for (size_t i=0; i < n; i++, s++)
		{
			if (R_FINITE(*s) && 0 <= *s && *s <= 2)
			{
				nvalid ++;
				sum += *s;
			}
		}
		rev_flag = (sum < nvalid);
	}

	SEXP rv_ans = PROTECT(NEW_LIST(2));
	SEXP Dest = NEW_NUMERIC(n);
	SET_ELEMENT(rv_ans, 0, Dest);
	SET_ELEMENT(rv_ans, 1, Allele);

	double *s = REAL(Geno), *p = REAL(Dest);

	if (rev_flag)
	{
		// allele
		string ss(strlen(base), 0);
		size_t nn = strlen(base+sign_pos+1);
		memcpy(&ss[0], base+sign_pos+1, nn);
		ss[nn] = ',';
		memcpy(&ss[nn+1], base, sign_pos);
		memcpy((void*)base, &ss[0], ss.size());
		// 2 - dosage
		for (; n > 0; n--, p++, s++)
			*p = (R_FINITE(*s) && 0 <= *s && *s <= 2) ? (*s) : R_NaN;
	} else {
		if (sign_pos >= 0)
			((char*)base)[sign_pos] = ',';
		// copy dosage
		for (; n > 0; n--, p++, s++)
			*p = (R_FINITE(*s) && 0 <= *s && *s <= 2) ? (2 - *s) : R_NaN;
	}

	UNPROTECT(1);
	return rv_ans;
}


// ======================================================================
// SeqArray GDS --> SNP GDS
// ======================================================================

COREARRAY_DLL_EXPORT SEXP FC_GDS2SNP(SEXP geno)
{
	C_UInt8 *p = (C_UInt8*)RAW(geno);
	for (size_t n = XLENGTH(geno); n > 0; n--)
	{
		if (*p > 3) *p = 3;
		p ++;
	}
	return geno;
}


// ======================================================================
// SeqArray GDS --> Dosage GDS
// ======================================================================

static int FC_Num_Sample = 0;

COREARRAY_DLL_EXPORT SEXP FC_SetNumSamp(SEXP num)
{
	FC_Num_Sample = Rf_asInteger(num);
	return R_NilValue;
}

COREARRAY_DLL_EXPORT SEXP FC_GDS2Dosage(SEXP dosage)
{
	int n = LENGTH(dosage);
	if (n < FC_Num_Sample)
	{
		dosage = NEW_NUMERIC(FC_Num_Sample);
		double *dst = REAL(dosage);
		for (int i=0; i < FC_Num_Sample; i++)
			dst[i] = R_NaN;
	} else if (n > FC_Num_Sample)
	{
		double *src = REAL(dosage);
		dosage = NEW_NUMERIC(FC_Num_Sample);
		double *dst = REAL(dosage);
		memcpy(dst, src, sizeof(double)*FC_Num_Sample);
	}
	return dosage;
}

} // extern "C"
