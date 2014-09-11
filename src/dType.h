// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// dType.h: Template classes for elementary types
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
 *	\file     dType.h
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2007 - 2014
 *	\brief    Template classes for elementary types
 *	\details
**/


#ifndef _H_COREARRAY_TYPE_
#define _H_COREARRAY_TYPE_

#include <CoreDEF.h>

#ifndef __STDC_LIMIT_MACROS
#  define __STDC_LIMIT_MACROS
#endif

#ifdef COREARRAY_UNIX
#  include <unistd.h>
#endif

#ifdef COREARRAY_MSC
#  include <msvc/stdint.h>
#else
#  include <stdint.h>
#endif

#include <cfloat>
#include <limits>
#include <string>


namespace CoreArray
{
	// ******************************************************************** //
	// ******************************************************************** //

	// define little-endian atomic type

	/// little-endian int8_t
	#define COREARRAY_LE_INT8_TRAIT_ID            0
	/// little-endian uint8_t
	#define COREARRAY_LE_UINT8_TRAIT_ID           1
	/// little-endian int16_t
	#define COREARRAY_LE_INT16_TRAIT_ID           2
	/// little-endian uint16_t
	#define COREARRAY_LE_UINT16_TRAIT_ID          3
	/// little-endian int32_t
	#define COREARRAY_LE_INT32_TRAIT_ID           4
	/// little-endian uint32_t
	#define COREARRAY_LE_UINT32_TRAIT_ID          5
	/// little-endian int64_t
	#define COREARRAY_LE_INT64_TRAIT_ID           6
	/// little-endian uint64_t
	#define COREARRAY_LE_UINT64_TRAIT_ID          7

	/// IEEE 32-bit floating point
	#define COREARRAY_IEEE_LE_FLOAT32_TRAIT_ID    10
	/// IEEE 64-bit floating point
	#define COREARRAY_IEEE_LE_FLOAT64_TRAIT_ID    11


	// define big-endian atomic type

	/// big-endian int8_t
	#define COREARRAY_BE_INT8_TRAIT_ID            20
	/// big-endian uint8_t
	#define COREARRAY_BE_UINT8_TRAIT_ID           21
	/// big-endian int16_t
	#define COREARRAY_BE_INT16_TRAIT_ID           22
	/// big-endian uint16_t
	#define COREARRAY_BE_UINT16_TRAIT_ID          23
	/// big-endian int32_t
	#define COREARRAY_BE_INT32_TRAIT_ID           24
	/// big-endian uint32_t
	#define COREARRAY_BE_UINT32_TRAIT_ID          25
	/// big-endian int64_t
	#define COREARRAY_BE_INT64_TRAIT_ID           26
	/// big-endian uint64_t
	#define COREARRAY_BE_UINT64_TRAIT_ID          27

	/// IEEE 32-bit floating point
	#define COREARRAY_IEEE_BE_FLOAT32_TRAIT_ID    30
	/// IEEE 64-bit floating point
	#define COREARRAY_IEEE_BE_FLOAT64_TRAIT_ID    31



	// define native atomic type

	#if defined(COREARRAY_LITTLE_ENDIAN)

		/// native int8_t
		#define COREARRAY_NATIVE_INT8_TRAIT_ID      COREARRAY_LE_INT8_TRAIT_ID
		/// native uint8_t
		#define COREARRAY_NATIVE_UINT8_TRAIT_ID     COREARRAY_LE_UINT8_TRAIT_ID
		/// native int16_t
		#define COREARRAY_NATIVE_INT16_TRAIT_ID     COREARRAY_LE_INT16_TRAIT_ID
		/// native uint16_t
		#define COREARRAY_NATIVE_UINT16_TRAIT_ID    COREARRAY_LE_UINT16_TRAIT_ID
		/// native int32_t
		#define COREARRAY_NATIVE_INT32_TRAIT_ID     COREARRAY_LE_INT32_TRAIT_ID
		/// native uint32_t
		#define COREARRAY_NATIVE_UINT32_TRAIT_ID    COREARRAY_LE_UINT32_TRAIT_ID
		/// native int64_t
		#define COREARRAY_NATIVE_INT64_TRAIT_ID     COREARRAY_LE_INT64_TRAIT_ID
		/// native uint64_t
		#define COREARRAY_NATIVE_UINT64_TRAIT_ID    COREARRAY_LE_UINT64_TRAIT_ID

		/// native 32-bit floating point
		#define COREARRAY_IEEE_NATIVE_FLOAT32_TRAIT_ID    COREARRAY_IEEE_LE_FLOAT32_TRAIT_ID
		/// native 64-bit floating point
		#define COREARRAY_IEEE_NATIVE_FLOAT64_TRAIT_ID    COREARRAY_IEEE_LE_FLOAT64_TRAIT_ID

	#elif defined(COREARRAY_BIG_ENDIAN)

		/// native int8_t
		#define COREARRAY_NATIVE_INT8_TRAIT_ID      COREARRAY_BE_INT8_TRAIT_ID
		/// native uint8_t
		#define COREARRAY_NATIVE_UINT8_TRAIT_ID     COREARRAY_BE_UINT8_TRAIT_ID
		/// native int16_t
		#define COREARRAY_NATIVE_INT16_TRAIT_ID     COREARRAY_BE_INT16_TRAIT_ID
		/// native uint16_t
		#define COREARRAY_NATIVE_UINT16_TRAIT_ID    COREARRAY_BE_UINT16_TRAIT_ID
		/// native int32_t
		#define COREARRAY_NATIVE_INT32_TRAIT_ID     COREARRAY_BE_INT32_TRAIT_ID
		/// native uint32_t
		#define COREARRAY_NATIVE_UINT32_TRAIT_ID    COREARRAY_BE_UINT32_TRAIT_ID
		/// native int64_t
		#define COREARRAY_NATIVE_INT64_TRAIT_ID     COREARRAY_BE_INT64_TRAIT_ID
		/// native uint64_t
		#define COREARRAY_NATIVE_UINT64_TRAIT_ID    COREARRAY_BE_UINT64_TRAIT_ID

		/// native 32-bit floating point
		#define COREARRAY_IEEE_NATIVE_FLOAT32_TRAIT_ID    COREARRAY_IEEE_BE_FLOAT32_TRAIT_ID
		/// native 64-bit floating point
		#define COREARRAY_IEEE_NATIVE_FLOAT64_TRAIT_ID    COREARRAY_IEEE_BE_FLOAT64_TRAIT_ID

	#else
	#  error "Unknown endianness"
    #endif




	/// Memory data type id
	enum TSVType {
		svCustom = 0,   ///< Unknown or customized type
		svCustomInt,    ///< Customized signed integer
		svCustomUInt,   ///< Customized unsigned integer
		svCustomFloat,  ///< Customized float number
		svCustomStr,    ///< Customized string type
		svInt8,         ///< Signed integer of 8 bits
		svUInt8,        ///< Unsigned integer of 8 bits
		svInt16,        ///< Signed integer of 16 bits
		svUInt16,       ///< Unsigned integer of 16 bits
		svInt32,        ///< Signed integer of 32 bits
		svUInt32,       ///< Unsigned integer of 32 bits
		svInt64,        ///< Signed integer of 64 bits
		svUInt64,       ///< Unsigned integer of 64 bits
		svFloat32,      ///< Float number of single precision (32 bits)
		svFloat64,      ///< Float number of double precision (64 bits)
		svStrUTF8,      ///< UTF-8 string
		svStrUTF16      ///< UTF-16 string
	};


	/// Whether x (TSVType) is an integer or not
	#define COREARRAY_SV_INTEGER(x) \
		((svInt8<=(x) && (x)<=svUInt64) || (x)==svCustomInt || (x)==svCustomUInt)

	/// Whether x (TSVType) is a signed integer or not
	#define COREARRAY_SV_SINT(x) \
		((x)==svInt8 || (x)==svInt16 || (x)==svInt32 || (x)==svInt64 || (x)==svCustomInt)

	/// Whether x (TSVType) is an unsigned integer or not
	#define COREARRAY_SV_UINT(x) \
		((x)==svUInt8 || (x)==svUInt16 || (x)==svUInt32 || (x)==svUInt64 || (x)==svCustomUInt)

	/// Whether x (TSVType) is a float number or not
	#define COREARRAY_SV_FLOAT(x) \
		((x)==svFloat32 || (x)==svFloat64 || (x)==svCustomFloat)

	/// Whether x (TSVType) is a string or not
	#define COREARRAY_SV_STRING(x) \
		((x)==svStrUTF8 || (x)==svStrUTF16 || (x)==svCustomStr)


	#define COREARRAY_TR_UNKNOWN                  -1
	#define COREARRAY_TR_CUSTOM                    0

	#define COREARRAY_TR_INTEGER                   1
	#define COREARRAY_TR_BIT_INTEGER               2

	#define COREARRAY_TR_FLOAT                     3

	#define COREARRAY_TR_STRING                    4
	#define COREARRAY_TR_FIXED_LENGTH_STRING       5
	#define COREARRAY_TR_VARIABLE_LENGTH_STRING    6





	// ******************************************************************** //
	// ******************************************************************** //

	// Integers

	/// Signed integer of 8 bits
	typedef int8_t      Int8;
	typedef Int8*       PInt8;		

	/// Unsigned integer of 8 bits
	typedef uint8_t     UInt8;
	typedef UInt8*      PUInt8;

	/// Signed integer of 16 bits
	typedef int16_t     Int16;
	typedef Int16*      PInt16;

	/// Unsigned integer of 16 bits
	typedef uint16_t    UInt16;
	typedef UInt16*     PUInt16;

	/// Signed integer of 32 bits
	typedef int32_t     Int32;
	typedef Int32*      PInt32;

	/// Unsigned integer of 32 bits
	typedef uint32_t    UInt32;
	typedef UInt32*     PUInt32;

	/// Signed integer of 64 bits
	typedef int64_t     Int64;
	typedef Int64*      PInt64;

	/// Unsigned integer of 64 bits
	typedef uint64_t    UInt64;
	typedef UInt64*     PUInt64;


	/// CoreArray Boolean
	typedef int8_t      CBOOL;
	typedef CBOOL*      PCBOOL;


	#if defined(COREARRAY_MSC) && !defined(ssize_t)
	typedef ptrdiff_t	ssize_t;
	#endif


	// ******************************************************************** //
	// ******************************************************************** //

	// Integer Traits

	template<typename T> struct TdTraits
	{
    	typedef T TType;
		typedef T ElmType;
		static const int trVal = COREARRAY_TR_UNKNOWN;
		static const unsigned BitOf = sizeof(T)*8u;
		static const bool isClass = false;
		static const TSVType SVType = svCustom;
	};

	template<> struct TdTraits<Int8>
	{
		typedef Int8 TType;
		typedef Int8 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 8u;
		static const bool isClass = false;
		static const TSVType SVType = svInt8;

		static const char * TraitName() { return "Int8"; }
		static const char * StreamName() { return "dInt8"; }

		COREARRAY_INLINE static short Min() { return INT8_MIN; }
		COREARRAY_INLINE static short Max() { return INT8_MAX; }
	};

	template<> struct TdTraits<UInt8>
	{
		typedef UInt8 TType;
		typedef UInt8 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 8u;
		static const bool isClass = false;
		static const TSVType SVType = svUInt8;
		enum {
			isNumeric = true
		};
		static const char * TraitName() { return "UInt8"; }
		static const char * StreamName() { return "dUInt8"; }

		COREARRAY_INLINE static unsigned short Min() { return 0; }
		COREARRAY_INLINE static unsigned short Max() { return UINT8_MAX; }
	};

	template<> struct TdTraits<Int16>
	{
		typedef Int16 TType;
		typedef Int16 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 16u;
		static const bool isClass = false;
		static const TSVType SVType = svInt16;

		static const char * TraitName() { return "Int16"; }
		static const char * StreamName() { return "dInt16"; }

		COREARRAY_INLINE static Int16 Min() { return INT16_MIN; }
		COREARRAY_INLINE static Int16 Max() { return INT16_MAX; }
	};

	template<> struct TdTraits<UInt16>
	{
		typedef UInt16 TType;
		typedef UInt16 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 16u;
		static const bool isClass = false;
		static const TSVType SVType = svUInt16;

		static const char * TraitName() { return "UInt16"; }
		static const char * StreamName() { return "dUInt16"; }

		COREARRAY_INLINE static UInt16 Min() { return 0; }
		COREARRAY_INLINE static UInt16 Max() { return UINT16_MAX; }
	};

	template<> struct TdTraits<Int32>
	{
		typedef Int32 TType;
		typedef Int32 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 32u;
		static const bool isClass = false;
		static const TSVType SVType = svInt32;

		static const char * TraitName() { return "Int32"; }
		static const char * StreamName() { return "dInt32"; }

		COREARRAY_INLINE static Int32 Min() { return INT32_MIN; }
		COREARRAY_INLINE static Int32 Max() { return INT32_MAX; }
	};

	template<> struct TdTraits<UInt32>
	{
		typedef UInt32 TType;
		typedef UInt32 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 32u;
		static const bool isClass = false;
		static const TSVType SVType = svUInt32;

		static const char * TraitName() { return "UInt32"; }
		static const char * StreamName() { return "dUInt32"; }

		COREARRAY_INLINE static UInt32 Min() { return 0; }
		COREARRAY_INLINE static UInt32 Max() { return UINT32_MAX; }
	};

	template<> struct TdTraits<Int64>
	{
		typedef Int64 TType;
		typedef Int64 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 64u;
		static const bool isClass = false;
		static const TSVType SVType = svInt64;

		static const char * TraitName() { return "Int64"; }
		static const char * StreamName() { return "dInt64"; }

		COREARRAY_INLINE static Int64 Min() { return std::numeric_limits<Int64>::min(); }
		COREARRAY_INLINE static Int64 Max() { return std::numeric_limits<Int64>::max(); }
	};

	template<> struct TdTraits<UInt64>
	{
		typedef UInt64 TType;
		typedef UInt64 ElmType;
		static const int trVal = COREARRAY_TR_INTEGER;
		static const unsigned BitOf = 64u;
		static const bool isClass = false;
		static const TSVType SVType = svUInt64;

		static const char * TraitName() { return "UInt64"; }
		static const char * StreamName() { return "dUInt64"; }

		COREARRAY_INLINE static UInt64 Min() { return 0; }
		COREARRAY_INLINE static UInt64 Max() { return std::numeric_limits<UInt64>::max(); }
	};



	// Float

	/// Float number of single precision (32 bits)
	typedef float Float32;
	/// Float number of double precision (64 bits)
	typedef double Float64;
	/// Float number of long precision
	typedef long double LongFloat;


  	// Float Traits

	template<> struct TdTraits<Float32>
	{
		typedef Float32 TType;
		typedef Float32 ElmType;
		static const int trVal = COREARRAY_TR_FLOAT;
		static const unsigned BitOf = 32u;
		static const bool isClass = false;
		static const TSVType SVType = svFloat32;

		static const char * TraitName() { return "Float32"; }
		static const char * StreamName() { return "dFloat32"; }

		COREARRAY_INLINE static Float32 Min() { return FLT_MIN; }
		COREARRAY_INLINE static Float32 Max() { return FLT_MAX; }
		COREARRAY_INLINE static Float32 Epsilon() { return FLT_EPSILON; }
		COREARRAY_INLINE static int Digits() { return FLT_MANT_DIG; }
	};

	template<> struct TdTraits<Float64>
	{
		typedef Float64 TType;
		typedef Float64 ElmType;
		static const int trVal = COREARRAY_TR_FLOAT;
		static const unsigned BitOf = 64u;
		static const bool isClass = false;
		static const TSVType SVType = svFloat64;

		static const char * TraitName() { return "Float64"; }
		static const char * StreamName() { return "dFloat64"; }

		COREARRAY_INLINE static Float64 Min() { return DBL_MIN; }
		COREARRAY_INLINE static Float64 Max() { return DBL_MAX; }
		COREARRAY_INLINE static Float64 Epsilon() { return DBL_EPSILON; }
		COREARRAY_INLINE static int Digits() { return DBL_MANT_DIG; }
	};

	template<> struct TdTraits<LongFloat>
	{
		typedef LongFloat TType;
		typedef LongFloat ElmType;
		static const int trVal = COREARRAY_TR_FLOAT;
		static const unsigned BitOf = sizeof(LongFloat)*8u;
		static const bool isClass = false;
		static const TSVType SVType = svCustomFloat;

		static const char * StreamName()
		{
		#if defined(COREARRAY_HAVE_FLOAT128)
			return "dFloat128";
		#elif defined(COREARRAY_LONGFLOAT_IS_DOUBLE)
			return "dFloat64";
		#else
        	return "dFloat80";
		#endif
		}
		static const char * TraitName() { return StreamName()+1; }

		COREARRAY_INLINE static LongFloat Min() { return LDBL_MIN; }
		COREARRAY_INLINE static LongFloat Max() { return LDBL_MAX; }
		COREARRAY_INLINE static LongFloat Epsilon() { return LDBL_EPSILON; }
		COREARRAY_INLINE static int Digits() { return LDBL_MANT_DIG; }
	};



	/// Customized integer type
	/** \tparam TYPE  any data type, e.g integer or float number
	 *  \tparam SIZE  to specify the structure size, can be != sizeof(TYPE)
	**/
	template<typename TYPE, ssize_t SIZE> struct TdNumber
	{
	public:
		/// The size of this type
		static const ssize_t size = SIZE;

		TdNumber() {}
		TdNumber(TYPE val) { fVal = val; }

		TdNumber<TYPE, SIZE> & operator+= (TYPE val)
			{ fVal += val; return *this; }
		TdNumber<TYPE, SIZE> & operator-= (TYPE val)
			{ fVal -= val; return *this; }
		TdNumber<TYPE, SIZE> & operator++ ()
			{ fVal++; return *this; }
		TdNumber<TYPE, SIZE> & operator-- ()
			{ fVal--; return *this; }
		/// Assignment
		TdNumber<TYPE, SIZE> & operator= (TYPE val)
			{ fVal = val; return *this; }

		bool operator== (const TdNumber<TYPE, SIZE> &val) const
			{ return fVal == val.fVal; }
		bool operator!= (const TdNumber<TYPE, SIZE> &val) const
			{ return fVal != val.fVal; }

		operator TYPE() const { return fVal; }
		TYPE &get() { return fVal; }
		const TYPE &get() const { return fVal; }

		/// Return minimum value of the type
		static const TYPE min() { return TdTraits<TYPE>::Min(); }
		/// Return maximum value of the type
		static const TYPE max() { return TdTraits<TYPE>::Max(); }

	private:
		TYPE fVal;
    };
}

#endif /* _H_COREARRAY_TYPE_ */
