// ===========================================================
//
// Index.h: Indexing Objects
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


#ifndef _HEADER_SEQ_INDEX_
#define _HEADER_SEQ_INDEX_

#include <R_GDS_CPP.h>

#include <string>
#include <vector>
#include <map>

#include <cctype>
#include <cstring>
#include "vectorization.h"

using namespace std;
using namespace CoreArray;


namespace SeqArray
{

// ===========================================================
// Chromosome indexing
// ===========================================================

/// Chromosome indexing object
class COREARRAY_DLL_LOCAL CChromIndex
{
public:
	/// range object
	struct TRange
	{
		int Start;   ///< the starting position
		int Length;  ///< the length
	};

	typedef vector<TRange> TRangeList;

	/// constructor
	CChromIndex() { }

	/// represent chromosome codes as a RLE object in Map
	void AddChrom(PdGDSFolder Root);

	/// the total length of a TRangeList object
	size_t RangeTotalLength(const TRangeList &RngList);

	/// 
	void ReadData(PdAbstractArray Obj, const TRangeList &RngList,
		void *OutBuf, C_SVType SVType);

	/// 
	map<string, TRangeList> Map;
};




// ===========================================================
// Define Exception
// ===========================================================

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

}

#endif /* _HEADER_SEQ_INDEX_ */
