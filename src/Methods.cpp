// ===========================================================
//
// Methods.cpp: the C/C++ codes for the SeqArray package
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

#include <Common.h>


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
	int *p = INTEGER(Geno);
	size_t N = XLENGTH(Geno), m = 0;
	for (size_t n=N; n > 0; n--)
	{
		if (*p++ == NA_INTEGER)
			m ++;
	}
	return ScalarReal((N > 0) ? (double(m) / N) : R_NaN);
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
	SEXP Geno = VECTOR_ELT(List, 0);
	int *p = INTEGER(Geno);
	int nAllele = GetNumOfAllele(CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	SEXP rv = NEW_NUMERIC(nAllele);
	double *pV = REAL(rv);
	memset((void*)pV, 0, sizeof(double)*nAllele);
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
		for (int i=0; i < nAllele; i++) pV[i] /= num;
	} else {
		for (int i=0; i < nAllele; i++) pV[i] = R_NaN;
	}

	return rv;
}


// ======================================================================

static int AlleleFreq_Index = 0;
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
	int *p = INTEGER(Geno);
	int nAllele = GetNumOfAllele(CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));
	int n = 0, m = 0;

	if (AlleleFreq_RefPtr == NULL)
	{
		if (AlleleFreq_Index < nAllele)
		{
			for (size_t N=XLENGTH(Geno); N > 0; N--)
			{
				int g = *p ++;
				if (g != NA_INTEGER)
				{
					n ++;
					if (g == AlleleFreq_Index)
						m ++;
				}
			}
		} else {
			return ScalarReal(R_NaN);
		}
	} else {
		if (AlleleFreq_RefPtr[AlleleFreq_Index] < nAllele)
		{
			for (size_t N=XLENGTH(Geno); N > 0; N--)
			{
				int g = *p ++;
				if (g != NA_INTEGER)
				{
					n ++;
					if (g == AlleleFreq_RefPtr[AlleleFreq_Index])
						m ++;
				}
			}
			AlleleFreq_Index ++;
		} else {
			AlleleFreq_Index ++;
			return ScalarReal(R_NaN);
		}
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
	int *p = INTEGER(Geno);
	int idx = GetIndexOfAllele(
		CHAR(STRING_ELT(AlleleFreq_Allele, AlleleFreq_Index)),
		CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));
	AlleleFreq_Index ++;

	if (idx >= 0)
	{
		int n = 0, m = 0;
		for (size_t N=XLENGTH(Geno); N > 0; N--)
		{
			int g = *p ++;
			if (g != NA_INTEGER)
			{
				n ++;
				if (g == idx)
					m ++;
			}
		}
		return ScalarReal((n > 0) ? (double(m) / n) : R_NaN);
	} else
		return ScalarReal(R_NaN);
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

/// Get a matrix from the numerators and denominators
COREARRAY_DLL_EXPORT SEXP FC_IBD_Div(SEXP NumeratorDenominator, SEXP N)
{
	size_t n = Rf_asInteger(N);
	size_t size = n * (n + 1) / 2;
	if ((size_t)XLENGTH(NumeratorDenominator) != 2*size)
		error("Invalid 'numerator' and 'denominator'.");

	SEXP rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, n));
	double *base = REAL(rv_ans);
	double *pN = REAL(NumeratorDenominator);
	double *pD = REAL(NumeratorDenominator) + size;
	for (size_t i=0; i < n; i++)
	{
		for (size_t j=i; j < n; j++)
		{
			base[i*n + j] = base[j*n + i] = (*pN) / (*pD);
			pN ++; pD ++;
		}
	}
	UNPROTECT(1);

	return rv_ans;
}

#define MISSING    0x7F

/// Calculate average IBD over loci
COREARRAY_DLL_EXPORT SEXP FC_IBD_OneLocus(SEXP Geno, SEXP NumeratorDenominator,
	SEXP M_ij)
{
	int *pdim = INTEGER(getAttrib(Geno, R_DimSymbol));
	int num_ploidy=pdim[0], num_sample=pdim[1];
	if (num_ploidy != 2)
		error("Should be diploid.");

	C_Int8 *pM = (C_Int8*)RAW(M_ij);
	Rbyte *g_i = RAW(Geno);
	C_Int64 Sum = 0;
	int nSum = 0;

	for (int i=0; i < num_sample; i++, g_i+=2)
	{
		if ((g_i[0] != NA_RAW) && (g_i[1] != NA_RAW))
		{
			*pM ++ = (g_i[0] == g_i[1]) ? 4 : 0;  // 4 * M_i
			Rbyte *g_j = g_i + 2;

			for (int j=i+1; j < num_sample; j++, g_j+=2)
			{
				if ((g_j[0] != NA_RAW) && (g_j[1] != NA_RAW))
				{
					C_Int8 val = 0;  // 4 * M_{ij}
					if (g_i[0] == g_j[0]) val++;
					if (g_i[0] == g_j[1]) val++;
					if (g_i[1] == g_j[0]) val++;
					if (g_i[1] == g_j[1]) val++;
					*pM ++ = val;
					Sum += val; nSum ++;
				} else
					*pM ++ = MISSING;
			}
		} else {
			for (int j=i; j < num_sample; j++)
				*pM ++ = MISSING;
		}
	}

	if (nSum > 0)
	{
		size_t n = size_t(num_sample) * (num_sample + 1) / 2;
		double Mb = (double)Sum / nSum * 0.25;
		double OneMb = 1 - Mb;
		double *pN = REAL(NumeratorDenominator);
		double *pD = REAL(NumeratorDenominator) + n;
		C_Int8 *pM = (C_Int8*)RAW(M_ij);
		for (; n > 0; n--)
		{
			if (*pM != MISSING)
			{
				*pN += (*pM) * 0.25 - Mb;
				*pD += OneMb;
			}
			pN ++; pD ++; pM ++;
		}
	}

	return R_NilValue;
}


static vector<Rbyte> IBD_TwoLoci_GenoBuffer;
static int IBD_TwoLoci_Interval;
static int IBD_TwoLoci_Interval_Start;
static int IBD_TwoLoci_Interval_Index;

COREARRAY_DLL_EXPORT SEXP FC_IBD_TwoLoci_Init(SEXP interval, SEXP num_samp)
{
	IBD_TwoLoci_Interval = Rf_asInteger(interval);
	if (IBD_TwoLoci_Interval <= 0)
		error("Invalid 'interval'.");
	IBD_TwoLoci_Interval_Start = 0;
	IBD_TwoLoci_Interval_Index = 0;
	IBD_TwoLoci_GenoBuffer.resize(Rf_asInteger(num_samp)*2*IBD_TwoLoci_Interval);
	return R_NilValue;
}

/// Calculate average IBD over loci
COREARRAY_DLL_EXPORT SEXP FC_IBD_TwoLoci(SEXP Geno, SEXP NumeratorDenominator,
	SEXP M_ij)
{
	const int *pdim = INTEGER(getAttrib(Geno, R_DimSymbol));
	const int num_ploidy=pdim[0], num_sample=pdim[1];
	const size_t size = num_sample * 2;
	if (num_ploidy != 2)
		error("Should be diploid.");

	if (IBD_TwoLoci_Interval_Start < IBD_TwoLoci_Interval)
	{
		memcpy(&IBD_TwoLoci_GenoBuffer[size * IBD_TwoLoci_Interval_Start],
			RAW(Geno), size);
		IBD_TwoLoci_Interval_Start ++;
		return R_NilValue;
	}

	C_Int8 *pM = (C_Int8*)RAW(M_ij);
	Rbyte *g1_i = &IBD_TwoLoci_GenoBuffer[size * IBD_TwoLoci_Interval_Index];
	Rbyte *g2_i = RAW(Geno);
	C_Int64 Sum = 0;
	int nSum = 0;

	for (int i=0; i < num_sample; i++, g1_i+=2, g2_i+=2)
	{
		if ((g1_i[0] != NA_RAW) && (g1_i[1] != NA_RAW) &&
			(g2_i[0] != NA_RAW) && (g2_i[1] != NA_RAW))
		{
			const int h1_i = g1_i[0] | (int(g2_i[0]) << 16);
			const int h2_i = g1_i[1] | (int(g2_i[1]) << 16);
			*pM ++ = (h1_i == h2_i) ? 4 : 0;  // 4 * M_i
			Rbyte *g1_j = g1_i + 2;
			Rbyte *g2_j = g2_i + 2;

			for (int j=i+1; j < num_sample; j++, g1_j+=2, g2_j+=2)
			{
				if ((g1_j[0] != NA_RAW) && (g1_j[1] != NA_RAW) &&
					(g2_j[0] != NA_RAW) && (g2_j[1] != NA_RAW))
				{
					const int h1_j = g1_j[0] | (int(g2_j[0]) << 16);
					const int h2_j = g1_j[1] | (int(g2_j[1]) << 16);
					C_Int8 val = 0;  // 4 * M_{ij}
					if (h1_i == h1_j) val++;
					if (h1_i == h2_j) val++;
					if (h2_i == h1_j) val++;
					if (h2_i == h2_j) val++;
					*pM ++ = val;
					Sum += val; nSum ++;
				} else
					*pM ++ = MISSING;
			}
		} else {
			for (int j=i; j < num_sample; j++)
				*pM ++ = MISSING;
		}
	}

	if (nSum > 0)
	{
		size_t n = size_t(num_sample) * (num_sample + 1) / 2;
		double Mb = (double)Sum / nSum * 0.25;
		double OneMb = 1 - Mb;
		double *pN = REAL(NumeratorDenominator);
		double *pD = REAL(NumeratorDenominator) + n;
		C_Int8 *pM = (C_Int8*)RAW(M_ij);
		for (; n > 0; n--)
		{
			if (*pM != MISSING)
			{
				*pN += (*pM) * 0.25 - Mb;
				*pD += OneMb;
			}
			pN ++; pD ++; pM ++;
		}
	}

	// copy genotype to the buffer
	memcpy(&IBD_TwoLoci_GenoBuffer[size * IBD_TwoLoci_Interval_Index],
		RAW(Geno), size);
	IBD_TwoLoci_Interval_Index ++;
	if (IBD_TwoLoci_Interval_Index >= IBD_TwoLoci_Interval)
		IBD_TwoLoci_Interval_Index = 0;

	return R_NilValue;
}

} // extern "C"
