// ===========================================================
//
// ReadByUnit.cpp: Read data variant by units of selected variants
//
// Copyright (C) 2019    Xiuwen Zheng
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


namespace SeqArray
{

}


extern "C"
{
using namespace SeqArray;

/// 
COREARRAY_DLL_EXPORT SEXP SEQ_Unit_SlidingWindows(SEXP Pos, SEXP Idx,
	SEXP WinSize, SEXP WinShift, SEXP WinStart, SEXP Tmp)
{
	int n = Rf_length(Pos);
	int *pos = INTEGER(Pos), *idx = INTEGER(Idx);
	int winsize  = Rf_asInteger(WinSize);
	int winshift = Rf_asInteger(WinShift);
	int winstart = Rf_asInteger(WinStart);
	int *tmp = INTEGER(Tmp);
	// get # of windows
	int i=0, num=0, w_start=winstart;
	while (i < n)
	{
		while (i<n && pos[i]<w_start) i++;
		int wend = w_start + winsize;
		int i2 = i;
		while (i<n && pos[i]<wend) i++;
		if (i > i2) num++;
		w_start += winshift;
		if (winshift < winsize) i = i2;
	}
	// get list
	SEXP rv = PROTECT(NEW_LIST(num));
	i=0; num=0; w_start=winstart;
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
		if (i > i2)
		{
			SEXP v = NEW_INTEGER(i - i2);
			memcpy(INTEGER(v), tmp, sizeof(int)*(i - i2));
			SET_ELEMENT(rv, num, v);
			num++;
		}
		w_start += winshift;
		if (winshift < winsize) i = i2;
	}
	UNPROTECT(1);
	return rv;
}


/// Apply functions over margins on a working space
/*
COREARRAY_DLL_EXPORT SEXP SEQ_Apply_Unit(SEXP gdsfile, SEXP var_name,
	SEXP FUN, SEXP as_is, SEXP var_index, SEXP param, SEXP rho);
*/
}
