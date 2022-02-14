// ===========================================================
//
// Methods.cpp: the C/C++ codes for the SeqArray package
//
// Copyright (C) 2015-2022    Xiuwen Zheng
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

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "Index.h"

using namespace SeqArray;


// used in FC_SetPackedGenoSxV() & FC_SetPackedGenoVxS()

static inline unsigned char g2b(Rbyte g)
{
	return (0<=g && g<=2) ? g : 3;
}

static inline unsigned char g2b(int g)
{
	return (0<=g && g<=2) ? g : 3;
}

static inline unsigned char g2b(double g)
{
	return R_FINITE(g) ? g2b((int)round(g)) : 3;
}

template<typename TYPE>
	static void packed_geno_SxV(Rbyte *p, TYPE *s, size_t n)
{
	unsigned char b0, b1, b2, b3;
	for (; n >= 4; n-=4, s+=4)
	{
		b0 = g2b(s[0]); b1 = g2b(s[1]); b2 = g2b(s[2]); b3 = g2b(s[3]);
		*p++ = b0 | (b1<<2) | (b2<<4) | (b3<<6);
	}
	switch (n)
	{
		case 3:
			b0 = g2b(s[0]); b1 = g2b(s[1]); b2 = g2b(s[2]);
			*p++ = b0 | (b1<<2) | (b2<<4) | (0x03<<6);
			break;
		case 2:
			b0 = g2b(s[0]); b1 = g2b(s[1]);
			*p++ = b0 | (b1<<2) | (0x03<<4) | (0x03<<6);
			break;
		case 1:
			b0 = g2b(s[0]);
			*p++ = b0 | (0x03<<2) | (0x03<<4) | (0x03<<6);
	}
}


template<typename TYPE>
	static void packed_geno_VxS(Rbyte *p, TYPE *s, size_t n, size_t m,
		unsigned char bit_shift)
{
	// bit_shift: 0, 1, 2, 3
	const static int mask[4] = { ~0x3, ~(0x3<<2), ~(0x3<<4), ~(0x3<<6) };
	const unsigned char bit_mask = mask[bit_shift];
	bit_shift *= 2;
	for (size_t i=0; i < n; i++)
	{
		unsigned char b = g2b(s[i]);
		unsigned char g = *p;
		*p = (g & bit_mask) | (b << bit_shift);
		p += m;
	}
}


// ======================================================================

extern "C"
{

static const char *ERR_DS_TYPE = "Invalid type of dosage.";
static const char *ERR_DS_DIM  = "Invalid dimension of dosage.";

static void get_ds_n_m(SEXP DS, int &n, int &m)
{
	n = XLENGTH(DS); m = 1;
	SEXP dm = GET_DIM(DS);
	if (!isNull(dm))
	{
		if (XLENGTH(dm) != 2)
			throw ErrSeqArray("# of dimensions should be 2 for dosages.");
		m = INTEGER(dm)[1];
	}
}


// ======================================================================

/// Calculate the missing rate per variant
COREARRAY_DLL_EXPORT SEXP FC_Missing_PerVariant(SEXP Geno)
{
	size_t n = XLENGTH(Geno), n_miss;
	switch (TYPEOF(Geno))
	{
	case RAWSXP:
		n_miss = vec_i8_count((char*)RAW(Geno), n, NA_RAW); break;
	case INTSXP:
		n_miss = vec_i32_count(INTEGER(Geno), n, NA_INTEGER); break;
	case REALSXP:
		n_miss = vec_f64_num_notfinite(REAL(Geno), n); break;
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}
	return ScalarReal((n > 0) ? (double(n_miss) / n) : R_NaN);
}

/// Calculate the missing rate per sample
COREARRAY_DLL_EXPORT SEXP FC_Missing_PerSamp(SEXP Geno, SEXP sum)
{
	int *pdim = INTEGER(GET_DIM(Geno));
	int num_ploidy=pdim[0], num_sample=pdim[1];
	int *pS = INTEGER(sum);

	if (TYPEOF(Geno) == RAWSXP)
	{
		const Rbyte *pG = RAW(Geno), MISSING=NA_RAW;
		for (int i=0; i < num_sample; i++)
		{
			for (int j=0; j < num_ploidy; j++)
				if (*pG++ == MISSING) pS[i]++;
		}
	} else {
		const int *pG = INTEGER(Geno);
		for (int i=0; i < num_sample; i++)
		{
			for (int j=0; j < num_ploidy; j++)
				if (*pG++ == NA_INTEGER) pS[i]++;
		}
	}

	return R_NilValue;
}

/// Calculate the missing rate per sample for annotation/format/DS
COREARRAY_DLL_EXPORT SEXP FC_Missing_DS_PerSamp(SEXP DS, SEXP sum, SEXP tmp)
{
	int num_samp = XLENGTH(sum), n, m;
	get_ds_n_m(DS, n, m);
	if (m * num_samp != n) throw ErrSeqArray(ERR_DS_DIM);

	int *pT = INTEGER(tmp);
	memset(pT, 0, sizeof(int)*num_samp);
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		{
			const Rbyte *pG = RAW(DS), MISSING=NA_RAW;
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (*pG++ == MISSING) pT[i]++;
			break;
		}
	case INTSXP:
		{
			const int *pG = INTEGER(DS), MISSING=NA_INTEGER;
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (*pG++ == MISSING) pT[i]++;
			break;
		}
	case REALSXP:
		{
			const double *pG = REAL(DS);
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (!R_FINITE(*pG++)) pT[i]++;
			break;
		}
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	int *pS = INTEGER(sum);
	for (int i=0; i < num_samp; i++)
		pS[i] += (pT[i] > 0) ? 1 : 0;

	return R_NilValue;
}

/// Calculate the missing rate per sample and variant
COREARRAY_DLL_EXPORT SEXP FC_Missing_SampVariant(SEXP Geno, SEXP sum)
{
	int *pdim = INTEGER(GET_DIM(Geno));
	int num_ploidy=pdim[0], num_sample=pdim[1];
	int *pS = INTEGER(sum), n=0;

	if (TYPEOF(Geno) == RAWSXP)
	{
		const Rbyte *pG = RAW(Geno), MISSING=NA_RAW;
		for (int i=0; i < num_sample; i++)
		{
			for (int j=0; j < num_ploidy; j++)
				if (*pG++ == MISSING) { pS[i]++; n++; }
		}
	} else {
		const int *pG = INTEGER(Geno);
		for (int i=0; i < num_sample; i++)
		{
			for (int j=0; j < num_ploidy; j++)
				if (*pG++ == NA_INTEGER) { pS[i]++; n++; }
		}
	}

	return ScalarReal((double)n / (num_ploidy*num_sample));
}

/// Calculate the missing rate per sample and variant
COREARRAY_DLL_EXPORT SEXP FC_Missing_DS_SampVariant(SEXP DS, SEXP sum, SEXP tmp)
{
	int num_samp = XLENGTH(sum), n, m;
	get_ds_n_m(DS, n, m);
	if (m * num_samp != n) throw ErrSeqArray(ERR_DS_DIM);

	int *pT = INTEGER(tmp);
	memset(pT, 0, sizeof(int)*num_samp);
	int n_miss;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		{
			const Rbyte *pG = RAW(DS), MISSING=NA_RAW;
			n_miss = vec_i8_count((char*)pG, n, NA_RAW);
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (*pG++ == MISSING) pT[i]++;
			break;
		}
	case INTSXP:
		{
			const int *pG = INTEGER(DS), MISSING=NA_INTEGER;
			n_miss = vec_i32_count(pG, n, NA_INTEGER);
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (*pG++ == MISSING) pT[i]++;
			break;
		}
	case REALSXP:
		{
			const double *pG = REAL(DS);
			n_miss = vec_f64_num_notfinite(pG, n);
			for (int j=0; j < m; j++)
				for (int i=0; i < num_samp; i++)
					if (!R_FINITE(*pG++)) pT[i]++;
			break;
		}
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	int *pS = INTEGER(sum);
	for (int i=0; i < num_samp; i++)
		pS[i] += (pT[i] > 0) ? 1 : 0;

	return ScalarReal((n > 0) ? (double(n_miss) / n) : R_NaN);
}


// ======================================================================

static int  AFreq_Index  = 0;
static int* AFreq_RefPtr = NULL;
static SEXP AFreq_Allele = NULL;
static bool AFreq_Minor  = false;  //< true for return MAF/MAC
static int  AFreq_Ploidy = 0;      //< ploidy, e.g., 2 by default

/// Set the reference allele with an index
COREARRAY_DLL_EXPORT SEXP FC_AF_SetIndex(SEXP idx, SEXP minor, SEXP ploidy)
{
	if (XLENGTH(idx) == 1)
	{
		AFreq_Index = Rf_asInteger(idx);
		AFreq_RefPtr = NULL;
	} else {
		AFreq_Index = 0;
		AFreq_RefPtr = INTEGER(idx);
	}
	AFreq_Minor = (Rf_asLogical(minor)==TRUE);
	AFreq_Ploidy = Rf_asInteger(ploidy);
	return R_NilValue;
}

/// Set the reference allele with string
COREARRAY_DLL_EXPORT SEXP FC_AF_SetAllele(SEXP allele, SEXP minor, SEXP ploidy)
{
	AFreq_Allele = allele;
	AFreq_Index = 0;
	AFreq_Minor = (Rf_asLogical(minor)==TRUE);
	AFreq_Ploidy = Rf_asInteger(ploidy);
	return R_NilValue;
}

// ========

/// Get a list of allele frequencies
COREARRAY_DLL_EXPORT SEXP FC_AF_List(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	const size_t N = XLENGTH(Geno);

	int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));
	SEXP rv = NEW_NUMERIC(nAllele);
	double *pV = REAL(rv);

	size_t n1, n2, n3;
	switch (nAllele)
	{
	case 2:
		if (TYPEOF(Geno) == RAWSXP)
			vec_i8_count3((const char*)RAW(Geno), N, 0, 1, NA_RAW, &n1, &n2, &n3);
		else
			vec_i32_count3(INTEGER(Geno), N, 0, 1, NA_INTEGER, &n1, &n2, &n3);
		n3 = N - n3;
		if (n3 > 0)
		{
			pV[0] = (double)n1 / n3;
			pV[1] = (double)n2 / n3;
		} else
			pV[0] = pV[1] = R_NaN;
		break;

	case 1:
		if (TYPEOF(Geno) == RAWSXP)
			vec_i8_count2((const char*)RAW(Geno), N, 0, NA_RAW, &n1, &n2);
		else
			vec_i32_count2(INTEGER(Geno), N, 0, NA_INTEGER, &n1, &n2);
		n2 = N - n2;
		pV[0] = (n2 > 0) ? (double)n1 / n2 : R_NaN;
		break;

	default:
		int num = 0;
		memset((void*)pV, 0, sizeof(double)*nAllele);
		if (TYPEOF(Geno) == RAWSXP)
		{
			C_UInt8 *p = (C_UInt8*)RAW(Geno);
			for (size_t n=N; n > 0; n--)
			{
				C_UInt8 g = *p++;
				if (g != NA_RAW)
				{
					num ++;
					if (g < nAllele) pV[g] ++;
				}
			}
		} else {
			int *p = INTEGER(Geno);
			for (size_t n=N; n > 0; n--)
			{
				int g = *p++;
				if (g != NA_INTEGER)
				{
					num ++;
					if ((0 <= g) && (g < nAllele))
						pV[g] ++;
				}
			}
		}
		if (num > 0)
		{
			const double scale = 1.0 / num;
			for (; (nAllele--) > 0;) (*pV++) *= scale;
		} else {
			for (; (nAllele--) > 0;) (*pV++) = R_NaN;
		}
	}

	return rv;
}

/// Get reference allele frequency
COREARRAY_DLL_EXPORT SEXP FC_AF_Ref(SEXP Geno)
{
	const size_t N = XLENGTH(Geno);
	size_t n, m;
	if (TYPEOF(Geno) == RAWSXP)
		vec_i8_count2((const char*)RAW(Geno), N, 0, NA_RAW, &m, &n);
	else
		vec_i32_count2(INTEGER(Geno), N, 0, NA_INTEGER, &m, &n);
	n = N - n;
	if (n > 0)
	{
		double p = double(m) / n;
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}

#define GET_SUM_NUM(TYPE, TYPEPTR, START, NA_VAL, N)  \
	{ \
		const TYPE *p = TYPEPTR(DS) + START; \
		for (int i=0; i < N; i++, p++) \
			if (*p != NA_VAL) { sum += *p; num++; } \
		break; \
	}
#define GET_SUM_NUM_F(START, N)  \
	{ \
		const double *p = REAL(DS) + START; \
		for (int i=0; i < N; i++, p++) \
			if (R_FINITE(*p)) { sum += *p; num++; } \
		break; \
	}

/// Get reference allele frequency from dosage
COREARRAY_DLL_EXPORT SEXP FC_AF_DS_Ref(SEXP DS)
{
	int n, m, num=0;
	get_ds_n_m(DS, n, m);

	double sum=0;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		GET_SUM_NUM(Rbyte, RAW, 0, NA_RAW, n)
	case INTSXP:
		GET_SUM_NUM(int, INTEGER, 0, NA_INTEGER, n)
	case REALSXP:
		GET_SUM_NUM_F(0, n)
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	if (num > 0)
	{
		double p = 1 - sum * m / (num * AFreq_Ploidy);
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele frequency
COREARRAY_DLL_EXPORT SEXP FC_AF_Index(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	const int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));

	const size_t N = XLENGTH(Geno);
	size_t n = 0, m = 0;
	int A = (AFreq_RefPtr==NULL) ?
		AFreq_Index : AFreq_RefPtr[AFreq_Index++];

	if (A < nAllele)
	{
		if (TYPEOF(Geno) == RAWSXP)
			vec_i8_count2((const char*)RAW(Geno), N, A, NA_RAW, &m, &n);
		else
			vec_i32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
		n = N - n;
	}

	if (n > 0)
	{
		double p = double(m) / n;
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele frequency
COREARRAY_DLL_EXPORT SEXP FC_AF_DS_Index(SEXP List)
{
	SEXP DS = VECTOR_ELT(List, 0);
	const int A = (AFreq_RefPtr==NULL) ?
		AFreq_Index : AFreq_RefPtr[AFreq_Index++];
	if (A == 0) return FC_AF_DS_Ref(DS);

	const int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));
	if (A >= nAllele) return ScalarReal(R_NaN);

	int n, m;
	get_ds_n_m(DS, n, m);
	if (A > m) return ScalarReal(R_NaN);

	double sum=0;
	int num = 0, nrow = n/m;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		GET_SUM_NUM(Rbyte, RAW, nrow*(A-1), NA_RAW, nrow)
	case INTSXP:
		GET_SUM_NUM(int, INTEGER, nrow*(A-1), NA_INTEGER, nrow)
	case REALSXP:
		GET_SUM_NUM_F(nrow*(A-1), nrow)
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	if (num > 0)
	{
		double p = sum / (num * AFreq_Ploidy);
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele frequency
COREARRAY_DLL_EXPORT SEXP FC_AF_Allele(SEXP List)
{
	SEXP Ref = STRING_ELT(AFreq_Allele, AFreq_Index++);
	int A = -1;
	if (Ref != NA_STRING)
		A = GetIndexOfAllele(CHAR(Ref), CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	size_t n = 0, m = 0;
	if (A >= 0)
	{
		SEXP Geno = VECTOR_ELT(List, 0);
		const size_t N = XLENGTH(Geno);
		if (TYPEOF(Geno) == RAWSXP)
		{
			if (A < 255)
				vec_i8_count2((const char*)RAW(Geno), N, A, NA_RAW, &m, &n);
			else
				n = N;
		} else
			vec_i32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
		n = N - n;
	}

	if (n > 0)
	{
		double p = double(m) / n;
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele frequency
COREARRAY_DLL_EXPORT SEXP FC_AF_DS_Allele(SEXP List)
{
	SEXP Ref = STRING_ELT(AFreq_Allele, AFreq_Index++);
	int A = -1;
	if (Ref != NA_STRING)
		A = GetIndexOfAllele(CHAR(Ref), CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	SEXP DS = VECTOR_ELT(List, 0);
	if (A == 0) return FC_AF_DS_Ref(DS);

	int n, m;
	get_ds_n_m(DS, n, m);  // get total # and # of columns
	double sum = 0;
	int num = 0, nrow = n/m;
	if (A > 0 && A <= m)
	{
		switch (TYPEOF(DS))
		{
		case RAWSXP:
			GET_SUM_NUM(Rbyte, RAW, nrow*(A-1), NA_RAW, nrow)
		case INTSXP:
			GET_SUM_NUM(int, INTEGER, nrow*(A-1), NA_INTEGER, nrow)
		case REALSXP:
			GET_SUM_NUM_F(nrow*(A-1), nrow)
		default:
			throw ErrSeqArray(ERR_DS_TYPE);
		}
	}

	if (num > 0)
	{
		double p = sum / (num * AFreq_Ploidy);
		if (AFreq_Minor && p>0.5) p = 1 - p;
		return ScalarReal(p);
	} else
		return ScalarReal(R_NaN);
}


// ======================================================================

/// Get reference allele count
COREARRAY_DLL_EXPORT SEXP FC_AC_Ref(SEXP Geno)
{
	const size_t N = XLENGTH(Geno);
	size_t n, m;
	if (TYPEOF(Geno) == RAWSXP)
		vec_i8_count2((const char*)RAW(Geno), N, 0, NA_RAW, &m, &n);
	else
		vec_i32_count2(INTEGER(Geno), N, 0, NA_INTEGER, &m, &n);
	if (AFreq_Minor)
	{
		n = N - n - m;  // allele count for alternative
		if (n < m) m = n;
	}
	return ScalarInteger(m);
}

/// Get reference allele frequency from dosage
COREARRAY_DLL_EXPORT SEXP FC_AC_DS_Ref(SEXP DS)
{
	int n, m, num=0;
	get_ds_n_m(DS, n, m);

	double sum=0;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		GET_SUM_NUM(Rbyte, RAW, 0, NA_RAW, n)
	case INTSXP:
		GET_SUM_NUM(int, INTEGER, 0, NA_INTEGER, n)
	case REALSXP:
		GET_SUM_NUM_F(0, n)
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	if (num > 0)
	{
		double totac = double(num * AFreq_Ploidy) / m;
		double ac = totac - sum;
		if (AFreq_Minor && ac>0.5*totac) ac = totac - ac;
		return ScalarReal(ac);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele count
COREARRAY_DLL_EXPORT SEXP FC_AC_Index(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	const int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));

	const size_t N = XLENGTH(Geno);
	int A = (AFreq_RefPtr==NULL) ?
		AFreq_Index : AFreq_RefPtr[AFreq_Index++];

	int ans;
	if (A < nAllele)
	{
		size_t n, m;
		if (TYPEOF(Geno) == RAWSXP)
			vec_i8_count2((const char*)RAW(Geno), N, A, NA_RAW, &m, &n);
		else
			vec_i32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
		if (AFreq_Minor)
		{
			n = N - n - m;  // allele count for alternative
			if (n < m) m = n;
		}
		ans = m;
	} else
		ans = NA_INTEGER;
	
	return ScalarInteger(ans);
}

/// Get allele count
COREARRAY_DLL_EXPORT SEXP FC_AC_DS_Index(SEXP List)
{
	SEXP DS = VECTOR_ELT(List, 0);
	const int A = (AFreq_RefPtr==NULL) ?
		AFreq_Index : AFreq_RefPtr[AFreq_Index++];
	if (A == 0) return FC_AC_DS_Ref(DS);

	const int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));
	if (A >= nAllele) return ScalarReal(R_NaN);

	int n, m;
	get_ds_n_m(DS, n, m);
	if (A > m) return ScalarReal(R_NaN);

	double sum=0;
	int num = 0, nrow = n/m;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		GET_SUM_NUM(Rbyte, RAW, nrow*(A-1), NA_RAW, nrow)
	case INTSXP:
		GET_SUM_NUM(int, INTEGER, nrow*(A-1), NA_INTEGER, nrow)
	case REALSXP:
		GET_SUM_NUM_F(nrow*(A-1), nrow)
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	if (num > 0)
	{
		double sum2 = num * AFreq_Ploidy - sum;
		if (AFreq_Minor && sum>sum2) sum = sum2;
		return ScalarReal(sum);
	} else
		return ScalarReal(R_NaN);
}

/// Get allele count
COREARRAY_DLL_EXPORT SEXP FC_AC_Allele(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	int A = GetIndexOfAllele(
		CHAR(STRING_ELT(AFreq_Allele, AFreq_Index++)),
		CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	int ans = NA_INTEGER;
	if (A >= 0)
	{
		const size_t N = XLENGTH(Geno);
		size_t n, m;
		if (TYPEOF(Geno) == RAWSXP)
		{
			if (A < 255)
			{
				vec_i8_count2((const char*)RAW(Geno), N, A, NA_RAW, &m, &n);
				if (AFreq_Minor)
				{
					n = N - n - m;  // allele count for alternative
					if (n < m) m = n;
				}
				ans = m;
			}
		} else {
			vec_i32_count2(INTEGER(Geno), N, A, NA_INTEGER, &m, &n);
			if (AFreq_Minor)
			{
				n = N - n - m;  // allele count for alternative
				if (n < m) m = n;
			}
			ans = m;
		}
	}

	return ScalarInteger(ans);
}

/// Get allele count
COREARRAY_DLL_EXPORT SEXP FC_AC_DS_Allele(SEXP List)
{
	SEXP Ref = STRING_ELT(AFreq_Allele, AFreq_Index++);
	int A = -1;
	if (Ref != NA_STRING)
		A = GetIndexOfAllele(CHAR(Ref), CHAR(STRING_ELT(VECTOR_ELT(List, 1), 0)));

	SEXP DS = VECTOR_ELT(List, 0);
	if (A == 0) return FC_AC_DS_Ref(DS);

	int n, m;
	get_ds_n_m(DS, n, m);  // get total # and # of columns
	double sum = 0;
	int num = 0, nrow = n/m;
	if (A > 0 && A <= m)
	{
		switch (TYPEOF(DS))
		{
		case RAWSXP:
			GET_SUM_NUM(Rbyte, RAW, nrow*(A-1), NA_RAW, nrow)
		case INTSXP:
			GET_SUM_NUM(int, INTEGER, nrow*(A-1), NA_INTEGER, nrow)
		case REALSXP:
			GET_SUM_NUM_F(nrow*(A-1), nrow)
		default:
			throw ErrSeqArray(ERR_DS_TYPE);
		}
	}

	if (num > 0)
	{
		double sum2 = num * AFreq_Ploidy - sum;
		if (AFreq_Minor && sum>sum2) sum = sum2;
		return ScalarReal(sum);
	} else
		return ScalarReal(R_NaN);
}


// ======================================================================

/// Convert a SeqArray GDS file to a SNP GDS file in `seqGDS2SNP()`
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


// ======================================================================

/// Get a list of allele counts
COREARRAY_DLL_EXPORT SEXP FC_AlleleCount(SEXP List)
{
	SEXP Geno = VECTOR_ELT(List, 0);
	const size_t N = XLENGTH(Geno);

	int nAllele = Rf_asInteger(VECTOR_ELT(List, 1));
	SEXP rv = NEW_INTEGER(nAllele);
	int *pV = INTEGER(rv);

	size_t n1, n2;
	switch (nAllele)
	{
	case 2:
		if (TYPEOF(Geno) == RAWSXP)
			vec_i8_count2((const char*)RAW(Geno), N, 0, 1, &n1, &n2);
		else
			vec_i32_count2(INTEGER(Geno), N, 0, 1, &n1, &n2);
		pV[0] = n1; pV[1] = n2;
		break;

	case 1:
		if (TYPEOF(Geno) == RAWSXP)
			pV[0] = vec_i8_count((const char*)RAW(Geno), N, 0);
		else
			pV[0] = vec_i32_count(INTEGER(Geno), N, 0);
		break;

	default:
		memset((void*)pV, 0, sizeof(int)*nAllele);
		if (TYPEOF(Geno) == RAWSXP)
		{
			C_UInt8 *p = (C_UInt8*)RAW(Geno);
			for (size_t n=N; n > 0; n--)
			{
				C_UInt8 g = *p++;
				if (g < nAllele) pV[g] ++;
			}
		} else {
			int *p = INTEGER(Geno);
			for (size_t n=N; n > 0; n--)
			{
				int g = *p++;
				if ((0 <= g) && (g < nAllele)) pV[g] ++;
			}
		}
	}

	return rv;
}


// ======================================================================

// val[0] -- AF/MAF, val[1] -- AC/MAC, val[2] -- missing
static double *af_ac_miss_ptr = NULL;
// ploidy
static int af_ac_miss_ploidy = 0;
// whether return minor AF/AC
static bool af_ac_miss_minor = false;

/// Initialize maf/mac/missing internal variable
COREARRAY_DLL_EXPORT SEXP FC_AF_AC_MISS_Init(SEXP M_3n, SEXP ploidy, SEXP minor)
{
	af_ac_miss_ptr = REAL(M_3n);  // M_3n is a matrix with three rows
	af_ac_miss_ploidy = Rf_asInteger(ploidy);
	if (af_ac_miss_ploidy < 0) af_ac_miss_ploidy = 2;
	af_ac_miss_minor = (Rf_asLogical(minor) == TRUE);
	return R_NilValue;
}

/// Get MAF/MAC/missing rate from integer genotypes
COREARRAY_DLL_EXPORT SEXP FC_AF_AC_MISS_Geno(SEXP Geno)
{
	const size_t N = XLENGTH(Geno);
	size_t n0, nmiss;
	if (TYPEOF(Geno) == RAWSXP)
		vec_i8_count2((const char*)RAW(Geno), N, 0, 0xFF, &n0, &nmiss);
	else
		vec_i32_count2(INTEGER(Geno), N, 0, NA_INTEGER, &n0, &nmiss);
	size_t n = N - nmiss;
	// AF/MAF
	double af = R_NaN;
	if (n > 0)
	{
		af = (double)n0 / n;
		if (af_ac_miss_minor && (af > 0.5))
			af = 1 - af;
	}
	af_ac_miss_ptr[0] = af;
	// AC/MAC
	double ac = n0, ac2 = n - n0;
	if (af_ac_miss_minor && (ac > ac2)) ac = ac2;
	af_ac_miss_ptr[1] = ac;
	// missing rate
	af_ac_miss_ptr[2] = (double)nmiss / N;
	af_ac_miss_ptr += 3;
	// output nothing
	return R_NilValue;
}

/// Get MAF/MAC/missing rate from dosages
COREARRAY_DLL_EXPORT SEXP FC_AF_AC_MISS_DS(SEXP DS)
{
	int n, m, num=0;
	get_ds_n_m(DS, n, m);

	double sum=0;
	switch (TYPEOF(DS))
	{
	case RAWSXP:
		GET_SUM_NUM(Rbyte, RAW, 0, NA_RAW, n)
	case INTSXP:
		GET_SUM_NUM(int, INTEGER, 0, NA_INTEGER, n)
	case REALSXP:
		GET_SUM_NUM_F(0, n)
	default:
		throw ErrSeqArray(ERR_DS_TYPE);
	}

	double af=R_NaN, ac=R_NaN;
	if (num > 0)
	{
		af = sum * m / (num * af_ac_miss_ploidy);
		if (af_ac_miss_minor && (af > 0.5))
			af = 1 - af;
		double totac = double(num * af_ac_miss_ploidy) / m;
		ac = totac - sum;
		if (af_ac_miss_minor && (ac > 0.5*totac))
			ac = totac - ac;
	}

	// AF/MAF
	af_ac_miss_ptr[0] = af;
	// AC/MAC
	af_ac_miss_ptr[1] = ac;
	// missing rate
	af_ac_miss_ptr[2] = (double)(n - num) / n;
	af_ac_miss_ptr += 3;
	// output nothing
	return R_NilValue;
}


// ======================================================================
// ======================================================================

static const char *pkg_digest = "digest";
static int digest_data_type = -1;

#define PKG_LOAD(name)	\
	{ \
		DL_FUNC f =  R_FindSymbol(#name, pkg_digest, NULL); \
		if (!f) \
			error("No function '%s' in the %s package", #name, pkg_digest); \
		memcpy(&name, &f, sizeof(f)); \
	}

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
	digest_data_type = -1;
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
	if (digest_data_type < 0)
	{
		if (TYPEOF(Data) == RAWSXP)
			digest_data_type = 0;
		else if (TYPEOF(Data) == INTSXP)
			digest_data_type = !inherits(Data, "factor") ? 1 : 2;
		else if (Rf_isLogical(Data))
			digest_data_type = 3;
		else if (Rf_isReal(Data))
			digest_data_type = 4;
		else if (Rf_isString(Data))
			digest_data_type = 5;
		else if (!Rf_isNull(Data))
			error("Not support data type.");
	}
	switch (digest_data_type)
	{
		case 0:
			(*md5_update)(&md5_ctx, RAW(Data), XLENGTH(Data));
			break;
		case 1:
			(*md5_update)(&md5_ctx, INTEGER(Data), XLENGTH(Data)*sizeof(int));
			break;
		case 2:
			{
				int *p = INTEGER(Data);
				SEXP ls = getAttrib(Data, R_LevelsSymbol);
				int nls = LENGTH(ls);
				for (size_t n=XLENGTH(Data); n > 0; n--, p++)
				{
					const char *s;
					if ((0 < *p) && (*p <= nls))
						s = CHAR(STRING_ELT(ls, *p - 1));
					else
						s = "";
					(*md5_update)(&md5_ctx, (void*)s, strlen(s)+1);
				}
			}
			break;
		case 3:
			(*md5_update)(&md5_ctx, LOGICAL(Data), XLENGTH(Data)*sizeof(int));
			break;
		case 4:
			(*md5_update)(&md5_ctx, REAL(Data), XLENGTH(Data)*sizeof(double));
			break;
		case 5:
			{
				const size_t n = XLENGTH(Data);
				for (size_t i=0; i < n; i++)
				{
					const char *s = CHAR(STRING_ELT(Data, i));
					(*md5_update)(&md5_ctx, (void*)s, strlen(s)+1);
				}
			}
	}

	return R_NilValue;
}


// ======================================================================

static Rbyte *geno_raw_ptr = NULL;
static size_t geno_nrow = 0;
static size_t geno_ncol = 0;
static size_t geno_index = 0;

static const char *ERR_PACKED_GENO_N =
	"Internal error: store genotype in packed raw format!";
static const char *ERR_PACKED_GENO_TYPE =
	"Internal error: invalid data type of dosage!";

/// initialize the packed raw genotype matrix
COREARRAY_DLL_EXPORT SEXP FC_InitPackedGeno(SEXP geno)
{
	geno_raw_ptr = RAW(geno);
	geno_nrow = INTEGER(GET_DIM(geno))[0];
	geno_ncol = INTEGER(GET_DIM(geno))[1];
	geno_index = 0;
	return R_NilValue;
}

/// store dosage in a 2-bit packed matrix (sample by variant)
COREARRAY_DLL_EXPORT SEXP FC_SetPackedGenoSxV(SEXP dosage)
{
	size_t n = Rf_xlength(dosage);
	if (n > geno_nrow*4) Rf_error(ERR_PACKED_GENO_N);

	Rbyte *p = geno_raw_ptr + geno_nrow * geno_index;
	geno_index ++;
	switch (TYPEOF(dosage))
	{
	case RAWSXP:
		packed_geno_SxV<Rbyte>(p, RAW(dosage), n);
		break;
	case INTSXP:
		packed_geno_SxV<int>(p, INTEGER(dosage), n);
		break;
	case REALSXP:
		packed_geno_SxV<double>(p, REAL(dosage), n);
		break;
	default:
		Rf_error(ERR_PACKED_GENO_TYPE);
	}

	return R_NilValue;
}

/// store dosage in a 2-bit packed matrix (variant by sample)
COREARRAY_DLL_EXPORT SEXP FC_SetPackedGenoVxS(SEXP dosage)
{
	size_t n = Rf_xlength(dosage);
	if (n !=  geno_ncol) Rf_error(ERR_PACKED_GENO_N);

	Rbyte *p = geno_raw_ptr + (geno_index >> 2);
	unsigned char bit_shift = geno_index & 0x03;
	geno_index ++;
	switch (TYPEOF(dosage))
	{
	case RAWSXP:
		packed_geno_VxS<Rbyte>(p, RAW(dosage), n, geno_nrow, bit_shift);
		break;
	case INTSXP:
		packed_geno_VxS<int>(p, INTEGER(dosage), n, geno_nrow, bit_shift);
		break;
	case REALSXP:
		packed_geno_VxS<double>(p, REAL(dosage), n, geno_nrow, bit_shift);
		break;
	default:
		Rf_error(ERR_PACKED_GENO_TYPE);
	}

	return R_NilValue;
}

} // extern "C"
