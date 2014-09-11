// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// dString.h: Character type and operation
//
// Copyright (C) 2007 - 2014	Xiuwen Zheng
//
// This file is part of CoreArray.
//
// CoreArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// CoreArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     dString.h
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2007 - 2014
 *	\brief    Character type and operation
 *	\details
**/


#ifndef _H_COREARRAY_STRING_
#define _H_COREARRAY_STRING_

#include <CoreDEF.h>
#include <dType.h>


namespace CoreArray
{
	// String

	/// UTF-8 character
	typedef char UTF8;

#if (WCHAR_MAX == UINT16_MAX) || (WCHAR_MAX == INT16_MAX)
#  define COREARRAY_SIZEOF_WCHAR 2
	/// UTF-16 character
	typedef wchar_t UTF16;
	/// UTF-32 character
	typedef Int32 UTF32;

#elif (WCHAR_MAX == UINT32_MAX) || (WCHAR_MAX == INT32_MAX)
#  define COREARRAY_SIZEOF_WCHAR 4
	/// UTF-16 character
	typedef Int16 UTF16;
	/// UTF-32 character
	typedef wchar_t UTF32;

#else
#  error "Unable to determine sizeof(wchar_t)."
#endif


	/// The string type
	typedef std::string                 RawString;

	/// UTF-8 string
	// class UTF8String: public std::basic_string<UTF8>
	//{
	//public:
	//	UTF8String()
	//}



	/// UTF-8 string
	typedef std::basic_string<UTF8>     UTF8String;
	/// UTF-16 string
	typedef std::basic_string<UTF16>    UTF16String;
	/// UTF-32 string
	typedef std::basic_string<UTF32>    UTF32String;


	// String Traits

	template<> struct TdTraits<UTF8String>
	{
		typedef UTF8String TType;
		typedef UTF8 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 8u;
		static const bool isClass = true;
		static const TSVType SVType = svStrUTF8;

		static const char * TraitName() { return "Str8"; }
	};

	template<> struct TdTraits<UTF8*>
	{
		typedef UTF8String TType;
		typedef UTF8 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 8u;
		static const bool isClass = false;
		static const TSVType SVType = svStrUTF8;

		static const char * TraitName() { return "Str8"; }
	};

	template<> struct TdTraits<UTF16String>
	{
		typedef UTF16String TType;
		typedef UTF16 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 16u;
		static const bool isClass = true;
		static const TSVType SVType = svStrUTF16;

		static const char * TraitName() { return "Str16"; }
	};

	template<> struct TdTraits<UTF16*>
	{
		typedef UTF16String TType;
		typedef UTF16 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 16u;
		static const bool isClass = false;
		static const TSVType SVType = svStrUTF16;

		static const char * TraitName() { return "Str16"; }
	};

	template<> struct TdTraits<UTF32String>
	{
		typedef UTF32String TType;
		typedef UTF32 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 32u;
		static const bool isClass = true;
		static const TSVType SVType = svCustomStr;

		static const char * TraitName() { return "Str32"; }
	};

	template<> struct TdTraits<UTF32*>
	{
		typedef UTF32String TType;
		typedef UTF32 ElmType;
		static const int trVal = COREARRAY_TR_STRING;
		static const unsigned BitOf = 32u;
		static const bool isClass = false;
		static const TSVType SVType = svCustomStr;

		static const char * TraitName() { return "Str32"; }
	};
}

#endif /* _H_COREARRAY_STRING_ */
