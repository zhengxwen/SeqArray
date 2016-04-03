// ===========================================================
//
// vectorization.h: compiler optimization with vectorization
//
// Copyright (C) 2016    Xiuwen Zheng
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
// You should have received a copy of the GNU General Public
// License along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     vectorization.h
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2016
 *	\brief    compiler optimization with vectorization
 *	\details
**/


#ifndef _HEADER_COREARRAY_VECTORIZATION_
#define _HEADER_COREARRAY_VECTORIZATION_

#include "CoreDEF.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#if (defined(__SSE__) && defined(__SSE2__))

#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2

#   if defined(__SSE3__)  // SSE3
#       include <pmmintrin.h>
#   endif

#   if defined(__SSSE3__)  // SSSE3
#       include <tmmintrin.h>
#   endif

#   if defined(__SSE4_2__) || defined(__POPCNT__)
#       define COREARRAY_HARDWARE_POPCNT
#       include <nmmintrin.h>  // SSE4_2, for POPCNT
#   endif

#   if defined(__AVX__) || defined(__AVX2__)
#       include <immintrin.h>  // AVX, AVX2
#   endif

#endif



#ifdef __cplusplus

namespace Vectorization
{
	/// an aligned pointer
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR
	{
	private:
		void *alloc_ptr, *base_ptr;

	public:
		ALIGN_PTR()
			{ alloc_ptr = base_ptr = NULL; }
		ALIGN_PTR(size_t n, size_t align)
			{
				alloc_ptr = base_ptr = NULL;
				reset(n, align);
			}
		~ALIGN_PTR()
			{
				if (alloc_ptr) free(alloc_ptr);
				alloc_ptr = base_ptr = NULL;
			}

		void reset(size_t n, size_t align)
		{
			if (n > 0)
			{
				if (align < 1) align = 1;
				alloc_ptr = realloc(alloc_ptr, n + align - 1);
				if (!alloc_ptr)
					throw "Insufficient memory.";
				size_t r = ((size_t)alloc_ptr) % align;
				base_ptr = r ? (void*)((char*)alloc_ptr + align - r) : alloc_ptr;
			} else {
				if (alloc_ptr) free(alloc_ptr);
				alloc_ptr = base_ptr = NULL;
			}
		}

		inline void *get() { return base_ptr; }
	};

	/// an aligned pointer with 16-byte alignment
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR_SSE: public ALIGN_PTR
	{
	public:
		ALIGN_PTR_SSE(): ALIGN_PTR()
			{ }
		ALIGN_PTR_SSE(size_t n): ALIGN_PTR(n, 16u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 16u); }
	};

	/// an aligned pointer with 32-byte alignment
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR_AVX: public ALIGN_PTR
	{
	public:
		ALIGN_PTR_AVX(): ALIGN_PTR()
			{ }
		ALIGN_PTR_AVX(size_t n): ALIGN_PTR(n, 32u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 32u); }
	};


	// auto pointer

#if defined(COREARRAY_SIMD_AVX)

	typedef ALIGN_PTR_AVX AUTO_PTR;

#elif defined(COREARRAY_SIMD_SSE)

	typedef ALIGN_PTR_SSE AUTO_PTR;

#else
	struct COREARRAY_DLL_DEFAULT AUTO_PTR: public ALIGN_PTR
	{
	public:
		AUTO_PTR(): ALIGN_PTR()
			{ }
		AUTO_PTR(size_t n): ALIGN_PTR(n, 1u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 1u); }
	};
#endif
}

#endif



#ifdef __cplusplus
extern "C" {
#endif


#ifdef COREARRAY_HARDWARE_POPCNT

#   define POPCNT_U32(x)    _mm_popcnt_u32((uint32_t)(x))
#   ifdef COREARRAY_REGISTER_BIT32
#       define POPCNT_U64(x)    \\
#           _mm_popcnt_u32((uint32_t)(x)) + _mm_popcnt_u32((uint64_t)(x) >> 32)
#   else
#       define POPCNT_U64(x)    _mm_popcnt_u64((uint64_t)(x))
#   endif

#else

inline static int POPCNT_U32(uint32_t x)
{
	x = x - ((x >> 1) & 0x55555555);
	x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
	return (((x + (x >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}

inline static int POPCNT_U64(uint64_t x)
{
	x -= ((x >> 1) & 0x5555555555555555LLU);
	x = (x & 0x3333333333333333LLU) + ((x >> 2) & 0x3333333333333333LLU);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
	return (x * 0x0101010101010101LLU) >> 56;
}

#endif



// ===========================================================

#ifdef __SSE__

#   ifdef __SSE3__
#       define MM_LOADU_128(p)    _mm_lddqu_si128((__m128i const*)(p))
#   else
#       define MM_LOADU_128(p)    _mm_loadu_si128((__m128i const*)(p))
#   endif

#   ifdef __SSE2__
#       define MM_BLEND_128(a, b, mask)  \
		    _mm_or_si128(_mm_and_si128(mask, a), _mm_andnot_si128(mask, b))
#   endif

#endif


#ifdef __AVX__

#   define MM_LOADU_256(p)    _mm256_loadu_si256((__m256i const *)(p))

#   ifdef __AVX2__
#       define MM_BLEND_256(a, b, mask)  \
		    _mm256_or_si256(_mm256_and_si256(mask, a), _mm256_andnot_si256(mask, b))
#   endif

#endif



/// get the number of non-zero
COREARRAY_DLL_DEFAULT size_t vec_i8_cnt_nonzero(const int8_t *p, size_t n);


// ===========================================================
// functions for int8
// ===========================================================

/// return the pointer to the non-zero character starting from p
COREARRAY_DLL_DEFAULT const char *vec_i8_ptr_nonzero(const char *p, size_t n);


/// count how many 'val' in 'p'
// COREARRAY_DLL_DEFAULT size_t vec_i8_count(int8_t *p, size_t n, int8_t val);


/// replace 'val' in the array of 'p' by 'substitute'
COREARRAY_DLL_DEFAULT void vec_i8_replace(int8_t *p, size_t n, int8_t val,
	int8_t substitute);

///
COREARRAY_DLL_DEFAULT void vec_i8_cnt_dosage2(const int8_t *p,
	int8_t *out, size_t n, int8_t val, int8_t missing,
	int8_t missing_substitute);



// ===========================================================
// functions for int32
// ===========================================================

/// count how many val in p, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT size_t vec_i32_count(const int32_t *p, size_t n, int32_t val);

/// count how many val1 and val2 in p, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_count2(const int32_t *p, size_t n,
	int32_t val1, int32_t val2, size_t *out_n1, size_t *out_n2);


COREARRAY_DLL_DEFAULT void vec_int32_set(int32_t *p, size_t n, int32_t val);

/// replace 'val' in the array of 'p' by 'substitute', assuming 'p' is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_replace(int32_t *p, size_t n, int32_t val,
	int32_t substitute);

/// assuming 'out' is 4-byte aligned, output (p[0]==val) + (p[1]==val) or missing_substitute
COREARRAY_DLL_DEFAULT void vec_i32_cnt_dosage2(const int32_t *p,
	int32_t *out, size_t n, int32_t val, int32_t missing,
	int32_t missing_substitute);




// ===========================================================
// functions for float64
// ===========================================================



// ===========================================================
// functions for char
// ===========================================================

COREARRAY_DLL_DEFAULT const char *vec_char_find_CRLF(const char *p, size_t n);



#ifdef __cplusplus
}
#endif

#endif /* _HEADER_COREARRAY_VECTORIZATION_ */
