// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// Common.h: the header file of SeqArray C/C++ codes
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

#include <dType.h>
#include <dString.h>
#include <CoreGDSLink.h>

#include <string>
#include <vector>
#include <list>
#include <map>

#include <ctype.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>

// R_XLEN_T_MAX is defined, R >= v3.0
#ifndef R_XLEN_T_MAX
#  define R_xlen_t    R_len_t
#  define XLENGTH     Rf_length
#endif


using namespace std;
using namespace CoreArray;
using namespace GDSInterface;


#define LongBool int


#ifdef COREARRAY_GNUG
#  ifdef COREARRAY_WINDOWS
#    define DLLEXPORT __attribute__((dllexport))
#  else
#    define DLLEXPORT
#  endif
#else
#  define DLLEXPORT __declspec(dllexport)
#endif



// ###########################################################
// defined macro
// ###########################################################

#define CORETRY			try {
#define CORETRY_CALL	bool has_error = false; CORETRY

#define CORECATCH(cmd)	} \
	catch (exception &E) { \
		gds_LastError() = E.what(); \
		cmd; \
	} \
	catch (const char *E) { \
		gds_LastError() = E; \
		cmd; \
	} \
	catch (...) { \
		gds_LastError() = "unknown error!"; \
		cmd; \
	}
#define CORECATCH_CALL	CORECATCH(has_error = true); \
	if (has_error) error(gds_LastError().c_str());




// ###########################################################
// The initialized object
// ###########################################################

/// the initial data
class TInitObject
{
public:
	struct TSelection
	{
		vector<CBOOL> Sample;
		vector<CBOOL> Variant;
	};

	typedef list<TSelection> TSelList;

	TInitObject();

	TSelection &Selection(SEXP gds);
	void Check_TrueArray(int Cnt);

	vector<CBOOL> TRUE_ARRAY;
	vector<UInt8> GENO_BUFFER;
	map<int, TSelList> _Map;
};

extern TInitObject Init;



// ###########################################################
// define exception
// ###########################################################

class ErrSeqArray: public ErrCoreArray
{
public:
	ErrSeqArray(): ErrCoreArray()
		{ }
	ErrSeqArray(const char *fmt, ...): ErrCoreArray()
		{ _COREARRAY_ERRMACRO_(fmt); }
	ErrSeqArray(const std::string &msg): ErrCoreArray()
		{ fMessage = msg; }
};



// ###########################################################
// private functions
// ###########################################################

/// get the list element named str, or return NULL
COREARRAY_INLINE static SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (R_len_t i = 0; i < length(list); i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

/// check CoreArray function
COREARRAY_INLINE static void CHECK(bool retval)
{
	if (!retval)
		throw ErrSeqArray(gds_LastError());
}
/// check CoreArray function
COREARRAY_INLINE static void* CHECK(void* retval)
{
	if (retval == NULL)
		throw ErrSeqArray(gds_LastError());
	return retval;
}
/// check CoreArray function
COREARRAY_INLINE static Int64 CHECK(Int64 retval)
{
	if (retval < 0)
		throw ErrSeqArray(gds_LastError());
	return retval;
}
/// check CoreArray function
COREARRAY_INLINE static size_t CHECK_SIZE(size_t retval)
{
	if (retval == ((size_t)-1))
		throw ErrSeqArray(gds_LastError());
	return retval;
}

/// check CoreArray function
COREARRAY_INLINE static const char *SKIP(const char *p)
{
	while (isspace(*p)) p ++;
	return p;
}

/// check CoreArray function
COREARRAY_INLINE static string SHORT_TEXT(const char *p, int MaxNum=16)
{
	if ((int)strlen(p) <= MaxNum)
		return string(p);
	else
		return string(p, MaxNum) + "...";
}


/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static PdGDSObj GDS_OBJECT(SEXP obj)
{
	PdGDSObj N;
	memcpy(&N, INTEGER(obj), sizeof(N));
	return N;
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static void GDS_PATH_PREFIX_CHECK(const char *path)
{
	for (; *path != 0; path++)
	{
		if ((*path == '~') || (*path == '@'))
		{
			throw ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.",
				*path);
		}
	}
}

/// check variable name
COREARRAY_INLINE static void GDS_VARIABLE_NAME_CHECK(const char *p)
{
	for (; *p != 0; p++)
	{
		if ((*p == '~') || (*p == '@') || (*p == '/'))
		{
			throw ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.", *p);
		}
	}
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static string GDS_PATH_PREFIX(const string &path, char prefix)
{
	string s = path;
	for (int i=s.size()-1; i >= 0; i--)
	{
		if (s[i] == '/')
		{
			s.insert(i+1, &prefix, 1);
			return s;
		}
	}
	s.insert(s.begin(), prefix);
	return s;
}

/// get PdGDSObj from a SEXP object
COREARRAY_INLINE static string GDS_UP_PATH(const char *path)
{
	const char *p = path + strlen(path) - 1;
	while ((p!=path) && (*p != '/')) p --;
	return string(path, p);
}

/// convert _SEXP to SEXP
COREARRAY_INLINE static SEXP _(_SEXP_ v)
{
	union {
		_SEXP_ f;
		SEXP t;
	} u;
	u.f = v;
	return u.t;
}
