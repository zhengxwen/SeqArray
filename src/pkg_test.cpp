// ===========================================================
//
// pkg_test.cpp: package testing with C/C++ codes
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
// You should have received a copy of the GNU General Public License
// along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.

#include "vectorization.h"
#include <R.h>
#include <Rdefines.h>

extern "C"
{

SEXP test_array_popcnt32(SEXP val)
{
	int n = length(val);
	int *s = INTEGER(val);
	SEXP rv_ans = NEW_INTEGER(n);
	int *p = INTEGER(rv_ans);
	for (int i=0; i < n; i++)
		p[i] = POPCNT_U32(s[i]);
	return rv_ans;
}


SEXP test_array_popcnt64(SEXP v1, SEXP v2)
{
	int n = length(v1);
	if (n != length(v2))
		error("error in 'test_popcnt64'.");
	int *s1 = INTEGER(v1), *s2 = INTEGER(v2);
	SEXP rv_ans = NEW_INTEGER(n);
	int *p = INTEGER(rv_ans);
	for (int i=0; i < n; i++)
	{
		uint64_t v = ((uint64_t)s1[i]) << 32 | (uint64_t)s2[i];
		p[i] = POPCNT_U64(v);
	}
	return rv_ans;
}


SEXP test_byte_count(SEXP val, SEXP start)
{
	int st = Rf_asInteger(start) - 1;
	int8_t *p = (int8_t *)RAW(val);
	int n = length(val);
	return ScalarInteger(vec_i8_cnt_nonzero(p + st, n - st));
}


SEXP test_int8_replace(SEXP val, SEXP start, SEXP find, SEXP substitute)
{
	int st = Rf_asInteger(start) - 1;
	int v1 = Rf_asInteger(find);
	int v2 = Rf_asInteger(substitute);
	int n = length(val);

	SEXP rv_ans = duplicate(val);
	int8_t *p = (int8_t *)RAW(rv_ans);
	vec_i8_replace(p + st, n - st, v1, v2);

	return rv_ans;
}


SEXP test_int32_count(SEXP val, SEXP start, SEXP find)
{
	int st = Rf_asInteger(start) - 1;
	int fd = Rf_asInteger(find);
	int *p = INTEGER(val);
	int n = length(val);
	return ScalarInteger(vec_i32_count(p + st, n - st, fd));
}


SEXP test_int32_count2(SEXP val, SEXP start, SEXP find1, SEXP find2)
{
	int st = Rf_asInteger(start) - 1;
	int fd1 = Rf_asInteger(find1);
	int fd2 = Rf_asInteger(find2);
	int *p = INTEGER(val);
	int n = length(val);

	size_t n1, n2;
	vec_i32_count2(p + st, n - st, fd1, fd2, &n1, &n2);
	SEXP rv_ans = NEW_INTEGER(2);
	INTEGER(rv_ans)[0] = n1;
	INTEGER(rv_ans)[1] = n2;

	return rv_ans;
}


SEXP test_int32_replace(SEXP val, SEXP start, SEXP find, SEXP substitute)
{
	int st = Rf_asInteger(start) - 1;
	int v1 = Rf_asInteger(find);
	int v2 = Rf_asInteger(substitute);
	int n = length(val);

	SEXP rv_ans = duplicate(val);
	int *p = INTEGER(rv_ans);
	vec_i32_replace(p + st, n - st, v1, v2);

	return rv_ans;
}

}
