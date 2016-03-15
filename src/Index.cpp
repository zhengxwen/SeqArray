// ===========================================================
//
// Index.cpp: Indexing Objects
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

#include "Index.h"

using namespace SeqArray;


// ===========================================================
// Chromosome Indexing
// ===========================================================

void CChromIndex::AddChrom(PdGDSFolder Root)
{
	PdAbstractArray varVariant = GDS_Node_Path(Root, "variant.id", TRUE);
	C_Int32 NumVariant = GDS_Array_GetTotalCount(varVariant);

	PdAbstractArray varChrom = GDS_Node_Path(Root, "chromosome", TRUE);
	C_Int32 NumChrom = GDS_Array_GetTotalCount(varChrom);

	if ((GDS_Array_DimCnt(varChrom) != 1) || (NumVariant != NumChrom))
		throw ErrSeqArray("Invalid dimension of 'chromosome'.");
	if (NumChrom <= 0) return;

	C_Int32 idx=0, len=1;
	string last;
	GDS_Array_ReadData(varChrom, &idx, &len, &last, svStrUTF8);
	idx ++;

	TRange rng;
	rng.Start = 0;
	rng.Length = 1;

	Map.clear();
	string txt[1024];

	while (idx < NumChrom)
	{
		len = NumChrom - idx;
		if (len > 1024) len = 1024;
		GDS_Array_ReadData(varChrom, &idx, &len, &txt, svStrUTF8);
		for (int i=0; i < len; i++)
		{
			if (txt[i] == last)
			{
				rng.Length ++;
			} else {
				Map[last].push_back(rng);
				last = string(txt[i].begin(), txt[i].end());
				rng.Start = idx + i;
				rng.Length = 1;
			}
		}
		idx += len;
	}

	Map[last].push_back(rng);
}

size_t CChromIndex::RangeTotalLength(const TRangeList &RngList)
{
	size_t ans = 0;
	vector<TRange>::const_iterator it;
	for (it=RngList.begin(); it != RngList.end(); it ++)
		ans += it->Length;
	return ans;
}



// ===========================================================
// Genomic Range Set
// ===========================================================

bool CRangeSet::less_range::operator()(const TRange &lhs, const TRange &rhs) const
{
	// -1 for two possible adjacent regions
	return (lhs.End < rhs.Start-1);
}

void CRangeSet::Clear()
{
	_RangeSet.clear();
}

void CRangeSet::AddRange(int start, int end)
{
	if (end < start) end = start;
	TRange rng;
	rng.Start = start; rng.End = end;

	do {
		set<TRange, less_range>::iterator it = _RangeSet.find(rng);
		if (it != _RangeSet.end())
		{
			if ((rng.Start < it->Start) || (rng.End > it->End))
			{
				if (rng.Start > it->Start) rng.Start = it->Start;
				if (rng.End < it->End) rng.End = it->End;
				_RangeSet.erase(it);
			} else
				break;
		} else {
			_RangeSet.insert(rng);
			break;
		}
	} while (1);
}

bool CRangeSet::IsIncluded(int point)
{
	TRange rng;
	rng.Start = rng.End = point;
	set<TRange, less_range>::iterator it = _RangeSet.find(rng);
	return it != _RangeSet.end();
}
