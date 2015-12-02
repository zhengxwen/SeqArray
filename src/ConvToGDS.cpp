// ===========================================================
//
// ConvToGDS.cpp: format conversion
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

#include "Common.h"


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

} // extern "C"
