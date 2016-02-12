// ===========================================================
//
// Methods.cpp: the C/C++ codes for the SeqArray package
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

#include <Common.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


extern "C"
{
// ======================================================================

/// Calculate the missing rate per variant
COREARRAY_DLL_EXPORT SEXP FC_SNP2GDS(SEXP Geno)
{
	size_t n = Rf_length(Geno);
	SEXP Dest = NEW_INTEGER(2*n);
	int *s = INTEGER(Geno), *p = INTEGER(Dest);
	for (; (n--) > 0; p+=2)
	{
		switch (*s++)
		{
			case  0: p[0] = p[1] = 0; break;
			case  1: p[0] = 1; p[1] = 0; break;
			case  2: p[0] = p[1] = 1; break;
			default: p[0] = p[1] = -1;
		}
	}
	return Dest;
}


// ======================================================================

/// Calculate the missing rate per variant
COREARRAY_DLL_EXPORT SEXP FC_NumAllele(SEXP AlleleStr)
{
	return ScalarInteger(GetNumOfAllele(CHAR(STRING_ELT(AlleleStr, 0))));
}


// ======================================================================

/// Calculate the missing rate per variant
COREARRAY_DLL_EXPORT SEXP FC_Missing_PerVariant(SEXP Geno)
{
	size_t n = XLENGTH(Geno);
	size_t m = vec_int32_count(INTEGER(Geno), n, NA_INTEGER);
	return ScalarReal((n > 0) ? (double(m) / n) : R_NaN);
}

/// Calculate the missing rate per sample
COREARRAY_DLL_EXPORT SEXP FC_Missing_PerSample(SEXP Geno, SEXP sum)
{
	int *pdim = INTEGER(getAttrib(Geno, R_DimSymbol));
	int num_ploidy=pdim[0], num_sample=pdim[1];
	int *pG = INTEGER(Geno);
	int *pS = INTEGER(sum);

	for (int i=0; i < num_sample; i++)
	{
		for (int j=0; j < num_ploidy; j++)
		{
			if (*pG++ == NA_INTEGER)
				pS[i] ++;
		}
	}

	return R_NilValue;
}


// ======================================================================

/// Get a list of allele frequencies
COREARRAY_DLL_EXPORT SEXP FC_AF_List(SEXP List)
{
	int nAllele = GetNumOfAllele(CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));
	SEXP rv = NEW_NUMERIC(nAllele);
	double *pV = REAL(rv);
	memset((void*)pV, 0, sizeof(double)*nAllele);

	SEXP Geno = VECTOR_ELT(List, 0);
	int *p = INTEGER(Geno);
	int num = 0;

	for (size_t n=XLENGTH(Geno); n > 0; n--)
	{
		int g = *p ++;
		if (g != NA_INTEGER)
		{
			if ((0 <= g) && (g < nAllele))
			{
				num ++; pV[g] ++;
			} else
				warning("Invalid value in 'genotype/data'.");
		}
	}

	if (num > 0)
	{
		const double scale = 1.0 / num;
		for (; (nAllele--) > 0;) (*pV++) *= scale;
	} else {
		for (; (nAllele--) > 0;) (*pV++) = R_NaN;
	}

	return rv;
}


// ======================================================================

static ssize_t AlleleFreq_Index = 0;
static int *AlleleFreq_RefPtr = NULL;
static SEXP AlleleFreq_Allele = R_NilValue;

/// Set the reference allele with an index
COREARRAY_DLL_EXPORT SEXP FC_AF_SetIndex(SEXP RefIndex)
{
	if (XLENGTH(RefIndex) == 1)
	{
		AlleleFreq_Index = Rf_asInteger(RefIndex);
		AlleleFreq_RefPtr = NULL;
	} else {
		AlleleFreq_Index = 0;
		AlleleFreq_RefPtr = INTEGER(RefIndex);
	}
	return R_NilValue;
}

/// Get allele frequencies
COREARRAY_DLL_EXPORT SEXP FC_AF_Index(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	const int nAllele = GetNumOfAllele(CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	const size_t N = XLENGTH(Geno);
	size_t n = 0, m = 0;
	int A = (AlleleFreq_RefPtr==NULL) ?
		AlleleFreq_Index : AlleleFreq_RefPtr[AlleleFreq_Index++];

	if (A < nAllele)
	{
		vec_int32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
		n = N - n;
	}

	return ScalarReal((n > 0) ? (double(m) / n) : R_NaN);
}

/// Set the reference allele with string
COREARRAY_DLL_EXPORT SEXP FC_AF_SetAllele(SEXP RefAllele)
{
	AlleleFreq_Allele = RefAllele;
	AlleleFreq_Index = 0;
	return R_NilValue;
}

/// Get allele frequencies
COREARRAY_DLL_EXPORT SEXP FC_AF_Allele(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	int A = GetIndexOfAllele(
		CHAR(STRING_ELT(AlleleFreq_Allele, AlleleFreq_Index++)),
		CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	size_t n = 0, m = 0;
	if (A >= 0)
	{
		const size_t N = XLENGTH(Geno);
		vec_int32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
		n = N - n;
	}

	return ScalarReal((n > 0) ? (double(m) / n) : R_NaN);
}


// ======================================================================

/// Convert a Sequencing GDS file to a SNP GDS file in `seqGDS2SNP()`
COREARRAY_DLL_EXPORT SEXP FC_AlleleStr(SEXP allele)
{
	const R_xlen_t n = XLENGTH(allele);
	for (R_xlen_t i=0; i < n; i++)
	{
		char *s = (char*)CHAR(STRING_ELT(allele, i));
		while (*s)
		{
			if (*s == ',')
				{ *s = '/'; break; }
			s ++;
		}
	}
	return allele;
}

/// Convert a Sequencing GDS file to a SNP GDS file in `seqGDS2SNP()`
COREARRAY_DLL_EXPORT SEXP FC_AlleleStr2(SEXP allele)
{
	const R_xlen_t n = XLENGTH(allele);
	for (R_xlen_t i=0; i < n; i++)
	{
		char *s = (char*)CHAR(STRING_ELT(allele, i));
		while (*s)
		{
			if (*s == '/')
				{ *s = ','; break; }
			s ++;
		}
	}
	return allele;
}


// ======================================================================

/// Get a list of allele frequencies
COREARRAY_DLL_EXPORT SEXP FC_AlleleCount(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	int *p = INTEGER(Geno);
	int nAllele = GetNumOfAllele(CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	SEXP rv = NEW_INTEGER(nAllele);
	int *pV = INTEGER(rv);
	memset((void*)pV, 0, sizeof(int)*nAllele);

	for (size_t n=XLENGTH(Geno); n > 0; n--)
	{
		int g = *p ++;
		if (g != NA_INTEGER)
		{
			if ((0 <= g) && (g < nAllele))
			{
				pV[g] ++;
			} else
				warning("Invalid value in 'genotype/data'.");
		}
	}

	return rv;
}


// ======================================================================

static const char *pkg_digest = "digest";

#define PKG_LOAD(name)	\
	*(DL_FUNC*)(&name) = R_FindSymbol(#name, pkg_digest, NULL); \
	if (!name) \
		error("No function '%s' in the %s package", #name, pkg_digest);

typedef struct _md5_context
{
	C_UInt32 total[2];
	C_UInt32 state[4];
	C_UInt8  buffer[64];
	char tmp[1024];
} md5_context;


static void (*md5_starts)(md5_context*) = NULL;
static void (*md5_update)(md5_context*, void*, C_UInt32) = NULL;
static void (*md5_finish)(md5_context*, C_UInt8[16]) = NULL;
static md5_context md5_ctx;

/// Initialize digest method
COREARRAY_DLL_EXPORT SEXP FC_DigestInit(SEXP Algo)
{
	if (!md5_starts) PKG_LOAD(md5_starts);
	if (!md5_update) PKG_LOAD(md5_update);
	if (!md5_finish) PKG_LOAD(md5_finish);
	(*md5_starts)(&md5_ctx);
	return R_NilValue;
}

/// Finalize digest method
COREARRAY_DLL_EXPORT SEXP FC_DigestDone(SEXP Algo)
{
	C_UInt8 digest[16];
	(*md5_finish)(&md5_ctx, digest);

	int Len = sizeof(digest);
	char buffer[1024 + 1];
	char *p = buffer;
	C_UInt8 *Code = digest;
	for (; Len > 0; Code++, Len--)
	{
		C_UInt8 v1 = (*Code) & 0x0F;
		C_UInt8 v2 = (*Code) >> 4;
		*p++ = (v2 < 10) ? (v2 + '0') : (v2 - 10 + 'a');
		*p++ = (v1 < 10) ? (v1 + '0') : (v1 - 10 + 'a');
	}
	*p = 0;

	return mkString(buffer);
}

/// Applied digest function
COREARRAY_DLL_EXPORT SEXP FC_DigestScan(SEXP Data)
{
	if (Rf_isInteger(Data))
	{
		const size_t n = XLENGTH(Data);
		(*md5_update)(&md5_ctx, INTEGER(Data), n*sizeof(int));
	} else if (Rf_isLogical(Data))
	{
		const size_t n = XLENGTH(Data);
		(*md5_update)(&md5_ctx, LOGICAL(Data), n*sizeof(int));
	} else if (Rf_isReal(Data))
	{
		const size_t n = XLENGTH(Data);
		(*md5_update)(&md5_ctx, REAL(Data), n*sizeof(double));
	} else if (Rf_isString(Data))
	{
		const size_t n = XLENGTH(Data);
		for (size_t i=0; i < n; i++)
		{
			const char *s = CHAR(STRING_ELT(Data, i));
			(*md5_update)(&md5_ctx, (void*)s, strlen(s)+1);
		}
	} else
		error("Not support data type.");

	return R_NilValue;
}


} // extern "C"
