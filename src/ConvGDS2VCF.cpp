// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// ConvGDS2VCF.cpp: the C++ code for the conversion from GDS to VCF
//
// Copyright (C) 2013 - 2014	Xiuwen Zheng [zhengx@u.washington.edu]
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
#include <vector>


// double quote the text if it is needed
static string QuoteText(const char *p)
{
	string rv;

	rv.clear();
	bool flag = false;
	for (; *p != 0; p ++)
	{
		switch (*p)
		{
			case ',': case ';':
				flag = true; rv.push_back(*p); break;
			case '\"':
				flag = true; rv.append("\\\""); break;
			case '\'':
				flag = true; rv.append("\\\'"); break;
			case ' ':
				flag = true; rv.push_back(' '); break;
			default:
				rv.push_back(*p);
		}
	}
	if (flag) // add double quote
	{
		rv.insert(0, "\"");
		rv.push_back('\"');
	}

	return rv;
}

// convert to a string
static const string TO_TEXT(SEXP X, int Start=0, int MaxCnt=-1,
	bool VarLength=false, bool NoBlank=true, int Step=1)
{
	char buffer[64];
	string ans;

	if (MaxCnt < 0)
		MaxCnt = (Rf_length(X) - Start) / Step;

	if (IS_INTEGER(X) || IS_LOGICAL(X))
	{
		int *Base = (IS_INTEGER(X) ? INTEGER(X) : LOGICAL(X)) + Start;
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
				if (Base[(MaxCnt-1)*Step] != NA_INTEGER) break;
		}
		for (int i=0; i < MaxCnt; i++, Base += Step)
		{
			if (i > 0) ans.push_back(',');
			if (*Base != NA_INTEGER)
			{
				snprintf(buffer, sizeof(buffer), "%d", *Base);
				ans.append(buffer);
			} else
				ans.push_back('.');
		}
	} else if (IS_NUMERIC(X))
	{
		double *Base = REAL(X) + Start;
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
				if (R_finite(Base[(MaxCnt-1)*Step])) break;
		}
		for (int i=0; i < MaxCnt; i++, Base += Step)
		{
			if (i > 0) ans.push_back(',');
			if (R_finite(*Base))
			{
				snprintf(buffer, sizeof(buffer), "%0.6g", *Base);
				ans.append(buffer);
			} else
				ans.push_back('.');
		}
	} else if (IS_CHARACTER(X) || Rf_isFactor(X))
	{
		if (Rf_isFactor(X))
			X = Rf_asCharacterFactor(X);
		if (VarLength || !NoBlank)
		{
			for (; MaxCnt > 0; MaxCnt --)
			{
				SEXP s = STRING_ELT(X, Start + (MaxCnt-1)*Step);
				if ((s != NA_STRING) && (CHAR(s)[0] != 0)) break;
			}
		}
		for (int i=0; i < MaxCnt; i ++, Start += Step)
		{
			if (i > 0) ans.push_back(',');
			if (STRING_ELT(X, Start) != NA_STRING)
				ans.append(QuoteText(CHAR(STRING_ELT(X, Start))));
			else
				ans.push_back('.');
		}
	}

	if (NoBlank)
	{
		if (ans.empty()) ans = ".";
	}

	return ans;
}


/// used in seq_OutVCF4
static vector<int> _VCF4_INFO_Number;    //< 
static vector<int> _VCF4_FORMAT_Number;  //< 



extern "C"
{
// ###########################################################
// Convert to VCF4: GDS -> VCF4
// ###########################################################

/// double quote text if needed
COREARRAY_DLL_EXPORT SEXP seq_Quote(SEXP text, SEXP dQuote)
{
	SEXP NewText, ans;
	PROTECT(NewText = AS_CHARACTER(text));
	PROTECT(ans = NEW_CHARACTER(Rf_length(NewText)));

	for (int i=0; i < Rf_length(NewText); i++)
	{
		string tmp = QuoteText(CHAR(STRING_ELT(NewText, i)));
		if (LOGICAL(dQuote)[0] == TRUE)
		{
			if ((tmp[0] != '\"') || (tmp[tmp.size()-1] != '\"'))
			{
				tmp.insert(0, "\"");
				tmp.push_back('\"');
			}
		}
		SET_STRING_ELT(ans, i, mkChar(tmp.c_str()));
	}

	UNPROTECT(2);
	return ans;
}


/// convert to VCF4
COREARRAY_DLL_EXPORT SEXP seq_InitOutVCF4(SEXP Info, SEXP Format)
{
	int *pInfo = INTEGER(Info);
	_VCF4_INFO_Number.assign(pInfo, pInfo + Rf_length(Info));
	int *pFmt = INTEGER(Format);
	_VCF4_FORMAT_Number.assign(pFmt, pFmt + Rf_length(Format));

	return R_NilValue;
}

/// convert to VCF4
COREARRAY_DLL_EXPORT SEXP seq_OutVCF4(SEXP X)
{
	const char *p, *s;
	string txt, tmp;
	int n;

	// variable list
	SEXP VarNames = getAttrib(X, R_NamesSymbol);

	// ************************************************************************
	// the first seven columns: chr, pos, id, allele (REF/ALT), qual, filter

	// CHROM
	txt.append(TO_TEXT(VECTOR_ELT(X, 0)));
	txt.push_back('\t');
	// POS
	txt.append(TO_TEXT(VECTOR_ELT(X, 1)));
	txt.push_back('\t');
	// ID
	txt.append(TO_TEXT(VECTOR_ELT(X, 2)));
	txt.push_back('\t');

	// allele -- REF/ALT
	s = p = CHAR(STRING_ELT(AS_CHARACTER(VECTOR_ELT(X, 3)), 0));
	n = 0;
	while ((*p != 0) && (*p != ','))
	{
		n ++; p ++;
	}

	// REF
	if (n > 0) txt.append(s, n); else txt.push_back('.');
	txt.push_back('\t');
	// ALT
	if (*p != 0)
	{
		p ++; txt.append((*p) ? p : ".");
	} else {
		txt.push_back('.');
	}
	txt.push_back('\t');

	// QUAL
	txt.append(TO_TEXT(VECTOR_ELT(X, 4)));
	txt.push_back('\t');
	// FILTER
	txt.append(TO_TEXT(VECTOR_ELT(X, 5)));
	txt.push_back('\t');


	// ************************************************************************
	// INFO

	bool NeedSeparator = false;
	n = 0;
	for (int i=0; i < (int)_VCF4_INFO_Number.size(); i++)
	{
		// name
		const char *nm = CHAR(STRING_ELT(VarNames, i + 8)) + 5;
		// SEXP
		SEXP D = VECTOR_ELT(X, i + 8);

		if (IS_LOGICAL(D))  // FLAG type
		{
			if (LOGICAL(D)[0] == TRUE)
			{
				if (NeedSeparator) txt.push_back(';');
				NeedSeparator = true;
				txt.append(nm);
				n ++;
			}
		} else {
			int L = _VCF4_INFO_Number[i];
			tmp = TO_TEXT(D, 0, (L < 0) ? -1 : L, (L < 0), false);
			if (!tmp.empty())
			{
				if (NeedSeparator) txt.push_back(';');
				NeedSeparator = true;
				txt.append(nm);
				txt.push_back('=');
				txt.append(tmp);
				n ++;
			}
		}
	}
	if (n <= 0) txt.push_back('.');	
	txt.push_back('\t');


	// ************************************************************************
	// FORMAT

	vector< pair<SEXP, int> > fmt_list;
	txt.append("GT");
	for (int i=0; i < (int)_VCF4_FORMAT_Number.size(); i ++)
	{
		const char *nm = CHAR(STRING_ELT(VarNames, i + 8 + _VCF4_INFO_Number.size()));
		SEXP D = VECTOR_ELT(X, i + 8 + _VCF4_INFO_Number.size());
		if (!isNull(D))
		{
			txt.push_back(':');
			txt.append(nm + 4);
			fmt_list.push_back(pair<SEXP, int>(D, i));
		}
	}
	txt.push_back('\t');


	// ************************************************************************
	// Genotypic data

	// genotype
	SEXP geno = VECTOR_ELT(X, 6);
	SEXP geno_dim = GET_DIM(geno);
	if (Rf_length(geno_dim) != 2)
		error("Invalid dimension of genotypes.");

	const int NumAllele = INTEGER(geno_dim)[0];
	const int NumSample = INTEGER(geno_dim)[1];

	// phase information
	SEXP phase = VECTOR_ELT(X, 7);
	if (Rf_length(phase) != (NumAllele-1)*NumSample)
		error("Invalid dimension of phasing information.");

	int *pSamp = INTEGER(geno);
	int *pAllele = INTEGER(phase);

	// for-loop of samples
	for (int i=0; i < NumSample; i ++)
	{
		// genotypes
		for (int j=0; j < NumAllele; j++, pSamp++)
		{
			if (j > 0)
			{
				txt.push_back(*pAllele ? '|' : '/');
				pAllele ++;
			}
			if (*pSamp != NA_INTEGER)
			{
				char buf[32];
				snprintf(buf, sizeof(buf), "%d", *pSamp);
				txt.append(buf);
			} else
				txt.push_back('.');
		}

		// annotation
		vector< pair<SEXP, int> >::iterator it;
		for (it=fmt_list.begin(); it != fmt_list.end(); it ++)
		{
			txt.push_back(':');
			int nTotal = Rf_length(it->first);
			int nColumn = nTotal / NumSample;
			if ((nTotal % NumSample) != 0)
				error("Internal Error: invalid dimension.");

			int L = _VCF4_FORMAT_Number[it->second];
			tmp = (L < 0) ? TO_TEXT(it->first, i, nColumn, true, true, NumSample) :
				TO_TEXT(it->first, i, L, false, true, NumSample);
			txt.append(tmp);
		}

		// add '\t'
		if (i < (NumSample-1)) txt.push_back('\t');
	}


	// append '\n'
	txt.push_back('\n');

	// return
	SEXP ans;
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(txt.c_str()));
	UNPROTECT(1);
	return ans;
}

} // extern "C"
