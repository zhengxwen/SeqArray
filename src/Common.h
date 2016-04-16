// ===========================================================
//
// Common.h: the C++ header file of SeqArray
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


#ifndef _HEADER_SEQ_COMMON_
#define _HEADER_SEQ_COMMON_

#include <R_GDS_CPP.h>
#include <dTrait.h>

#include <string>
#include <vector>
#include <list>
#include <map>

#include <cctype>
#include <cstring>

#include "Index.h"
#include "vectorization.h"

using namespace std;
using namespace CoreArray;


#define LongBool int

#define NA_RAW     0xFF


// ===========================================================
// The Initialized Object
// ===========================================================

/// the initial data
class COREARRAY_DLL_LOCAL TInitObject
{
public:
	/// constructor
	TInitObject();

	/// a vector of TRUE
	C_BOOL TRUE_ARRAY[1024];

	/// the buffer of genotypes
	vector<C_UInt8> GENO_BUFFER;
	/// allocator buffer according to size at least
	void Need_GenoBuffer(size_t size);
};

extern TInitObject Init;





// ===========================================================
// Library Functions
// ===========================================================

/// Get the number of TRUEs
#define GetNumOfTRUE(ptr, n)    vec_i8_cnt_nonzero((C_Int8*)(ptr), n)

/// Get the total count requiring the number of dimension is one
COREARRAY_DLL_LOCAL int GetGDSObjCount(PdAbstractArray Obj, const char *varname);

/// Get the number of alleles
COREARRAY_DLL_LOCAL int GetNumOfAllele(const char *allele_list);

/// Get the index in an allele list
COREARRAY_DLL_LOCAL int GetIndexOfAllele(const char *allele, const char *allele_list);

/// Get strings split by comma
COREARRAY_DLL_LOCAL void GetAlleles(const char *alleles, vector<string> &out);



// ===========================================================
// Private functions
// ===========================================================

/// get PdGDSObj from a SEXP object
inline static void GDS_PATH_PREFIX_CHECK(const char *path)
{
	for (; *path != 0; path++)
	{
		if ((*path == '~') || (*path == '@'))
		{
			throw SeqArray::ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.",
				*path);
		}
	}
}

/// check variable name
inline static void GDS_VARIABLE_NAME_CHECK(const char *p)
{
	for (; *p != 0; p++)
	{
		if ((*p == '~') || (*p == '@') || (*p == '/'))
		{
			throw SeqArray::ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.", *p);
		}
	}
}

/// get PdGDSObj from a SEXP object
inline static string GDS_PATH_PREFIX(const string &path, char prefix)
{
	string s = path;
	for (int i=s.size()-1; i >= 0; i--)
	{
		if (s[i] == '/')
		{
			if (((int)s.size() > i+1) && (s[i+1] == '~'))
				s[i+1] = prefix;
			else
				s.insert(i+1, &prefix, 1);
			return s;
		}
	}

	if ((s.size() > 0) && (s[0] == '~'))
		s[0] = prefix;
	else
		s.insert(s.begin(), prefix);

	return s;
}

/// get PdGDSObj from a SEXP object
inline static string GDS_UP_PATH(const char *path)
{
	const char *p = path + strlen(path) - 1;
	while ((p!=path) && (*p != '/')) p --;
	return string(path, p);
}

#endif /* _HEADER_SEQ_COMMON_ */
