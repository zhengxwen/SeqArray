// ===========================================================
//
// ReadByUnit.cpp: Read data variant by units of selected variants
//
// Copyright (C) 2019-2024    Xiuwen Zheng
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

#include "Index.h"


extern "C"
{
using namespace SeqArray;

/// 
COREARRAY_DLL_EXPORT SEXP SEQ_Unit_SlidingWindows(SEXP Pos, SEXP Idx,
	SEXP WinSize, SEXP WinShift, SEXP WinStart, SEXP DupFlag, SEXP Tmp)
{
	int n = Rf_length(Pos);
	int *pos = INTEGER(Pos), *idx = INTEGER(Idx);
	int winsize  = Rf_asInteger(WinSize);
	int winshift = Rf_asInteger(WinShift);
	int winstart = Rf_asInteger(WinStart);
	int duprmflag = Rf_asLogical(DupFlag);
	if (duprmflag == NA_LOGICAL)
		Rf_error("'dup.rm' must be TRUE or FALSE.");
	int *tmp = INTEGER(Tmp);

	// get the total number of windows
	int i=0, num=0, w_start=winstart;
	int old_i=i, old_i2=i;
	while (i < n)
	{
		while (i<n && pos[i]<w_start) i++;
		int wend = w_start + winsize;
		int i2 = i;
		while (i<n && pos[i]<wend) i++;
		if (duprmflag)
		{
			if (i > i2)
			{
				if (i!=old_i || i2!=old_i2)
				{
					old_i = i; old_i2 = i2;
					num++;
				}
			}
		} else {
			num++;
		}
		w_start += winshift;
		if (winshift < winsize) i = i2;
	}

	// save the results
	SEXP rv_ans = PROTECT(NEW_LIST(2));
	SEXP rv_st  = PROTECT(NEW_INTEGER(num));
	SEXP rv_lst = PROTECT(NEW_LIST(num));
	SET_VECTOR_ELT(rv_ans, 0, rv_st);
	SET_VECTOR_ELT(rv_ans, 1, rv_lst);

	// get list
	i=0; num=0; w_start=winstart;
	old_i = old_i2 = i;
	while (i < n)
	{
		while (i<n && pos[i]<w_start) i++;
		int wend = w_start + winsize;
		int i2 = i;
		while (i<n && pos[i]<wend)
		{
			tmp[i - i2] = idx[i];
			i++;
		}
		if (duprmflag)
		{
			if (i > i2)
			{
				if (i!=old_i || i2!=old_i2)
				{
					old_i = i; old_i2 = i2;
					INTEGER(rv_st)[num] = w_start;
					SEXP v = NEW_INTEGER(i - i2);
					memcpy(INTEGER(v), tmp, sizeof(int)*(i - i2));
					SET_ELEMENT(rv_lst, num, v);
					num++;
				}
			}
		} else {
			INTEGER(rv_st)[num] = w_start;
			SET_ELEMENT(rv_lst, num, NEW_INTEGER(0));
			num++;
		}
		w_start += winshift;
		if (winshift < winsize) i = i2;
	}

	UNPROTECT(3);
	return rv_ans;
}


/// Apply functions over margins on a working space
/*
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Unit(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho);
*/
}
