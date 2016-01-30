// ===========================================================
//
// vectorization.h: compiler optimization with vectorization
//
// Copyright (C) 2016    Xiuwen Zheng
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
#include <stdint.h>
#include <string.h>

#if (defined(COREARRAY_SIMD_SSE) && defined(COREARRAY_SIMD_SSE2))

#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2

#   if defined(__SSE4_2__) || defined(__POPCNT__)
#       define COREARRAY_HARDWARE_POPCNT
#       include <nmmintrin.h>  // SSE4_2, for POPCNT
#   endif

#   if defined(COREARRAY_SIMD_AVX)
#       include <immintrin.h>  // AVX, AVX2
#   endif

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



/// get the number of non-zero
COREARRAY_DLL_DEFAULT size_t vec_byte_count(const uint8_t *p, size_t n);


COREARRAY_DLL_DEFAULT void vec_int32_set(int32_t *p, size_t n, int32_t val);


#ifdef __cplusplus
}
#endif

#endif /* _HEADER_COREARRAY_VECTORIZATION_ */
