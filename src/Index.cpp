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

using namespace std;


namespace SeqArray
{
// ===========================================================
// Chromosome Indexing
// ===========================================================

CChromIndex::CChromIndex() { }

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

void CChromIndex::Clear()
{
	Map.clear();
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



// ===========================================================
// SeqArray GDS file information
// ===========================================================

static const char *ERR_DIM = "Invalid dimension of '%s'.";
static const char *ERR_FILE_ROOT = "CFileInfo::FileRoot should be initialized.";

CFileInfo::CFileInfo(PdGDSFolder root)
{
	_Root = NULL;
	_SampleNum = _VariantNum = 0;
	ResetRoot(root);
}

CFileInfo::~CFileInfo()
{
	_Root = NULL;
	_SampleNum = _VariantNum = 0;
}

void CFileInfo::ResetRoot(PdGDSFolder root)
{
	if (_Root != root)
	{
		// initialize
		_Root = root;
		SelList.clear();
		_Chrom.Clear();
		_Position.clear();
		// sample.id
		PdAbstractArray N = GDS_Node_Path(root, "sample.id", TRUE);
		C_Int64 n = GDS_Array_GetTotalCount(N);
		if ((n < 0) || (n > 2147483647))
			throw ErrSeqArray(ERR_DIM, "sample.id");
		_SampleNum = n;
		// variant.id
		N = GDS_Node_Path(root, "variant.id", TRUE);
		n = GDS_Array_GetTotalCount(N);
		if ((n < 0) || (n > 2147483647))
			throw ErrSeqArray(ERR_DIM, "variant.id");
		_VariantNum = n;
	}
}

TSelection &CFileInfo::Selection()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	if (SelList.empty())
		SelList.push_back(TSelection());

	TSelection &s = SelList.back();
	if (s.Sample.empty())
		s.Sample.resize(_SampleNum, TRUE);
	if (s.Variant.empty())
		s.Variant.resize(_VariantNum, TRUE);

	return s;
}

CChromIndex &CFileInfo::Chromosome()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	if (_Chrom.Map.empty())
		_Chrom.AddChrom(_Root);
	return _Chrom;
}

vector<C_Int32> &CFileInfo::Position()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	if (_Position.empty())
	{
		PdAbstractArray N = GetObj("position", TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != _VariantNum))
			throw ErrSeqArray(ERR_DIM, "position");
		// read
		_Position.resize(_VariantNum, 0);
		GDS_Array_ReadData(N, NULL, NULL, &_Position[0], svInt32);
	}
	return _Position;
}

PdAbstractArray CFileInfo::GetObj(const char *name, C_BOOL MustExist)
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	return GDS_Node_Path(_Root, name, MustExist);
}

int CFileInfo::SampleSelNum()
{
	TSelection &sel = Selection();
	return vec_i8_cnt_nonzero((C_Int8*)&sel.Sample[0], _SampleNum);
}

int CFileInfo::VariantSelNum()
{
	TSelection &sel = Selection();
	return vec_i8_cnt_nonzero((C_Int8*)&sel.Variant[0], _VariantNum);
}


// ===========================================================

/// File info list
std::map<int, CFileInfo> COREARRAY_DLL_LOCAL GDSFile_ID_Info;

/// get the associated CFileInfo
COREARRAY_DLL_LOCAL CFileInfo &GetFileInfo(SEXP gdsfile)
{
	SEXP ID = RGetListElement(gdsfile, "id");
	if (Rf_isNull(ID))
		throw ErrSeqArray("Invalid gds object.");

	int id = Rf_asInteger(ID);
	PdGDSFolder root = GDS_R_SEXP2FileRoot(gdsfile);

	map<int, CFileInfo>::iterator p = GDSFile_ID_Info.find(id);
	if (p == GDSFile_ID_Info.end())
	{
		GDSFile_ID_Info[id].ResetRoot(root);
		p = GDSFile_ID_Info.find(id);
	} else {
		if (p->second.Root() != root)
			p->second.ResetRoot(root);
	}

	return p->second;
}



// ===========================================================
// GDS Variable Type
// ===========================================================

static C_BOOL TRUEs[64] = {
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1
};

C_BOOL *CVarApply::NeedTRUEs(size_t size)
{
	if (size <= sizeof(TRUEs))
	{
		return TRUEs;
	} else if (size > _TRUE.size())
	{
		_TRUE.resize(size, TRUE);
	}
	return &_TRUE[0];
}






// ===========================================================
// Define Functions
// ===========================================================

// the buffer of TRUEs
static vector<C_BOOL> TrueBuffer;

COREARRAY_DLL_LOCAL size_t RLength(SEXP val)
{
	return (!Rf_isNull(val)) ? XLENGTH(val) : 0;
}

COREARRAY_DLL_LOCAL SEXP RGetListElement(SEXP list, const char *name)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	size_t n = RLength(names);
	for (size_t i = 0; i < n; i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

COREARRAY_DLL_LOCAL C_BOOL *NeedArrayTRUEs(size_t len)
{
	if (len <= sizeof(TRUEs))
		return TRUEs;
	else if (len > TrueBuffer.size())
		TrueBuffer.resize(len, TRUE);
	return &TrueBuffer[0];
}


static char pretty_num_buffer[32];

/// Get pretty text for an integer with comma
COREARRAY_DLL_LOCAL const char *PrettyInt(int val)
{
	char *p = pretty_num_buffer + sizeof(pretty_num_buffer);
	*(--p) = 0;

	bool sign = (val < 0);
	if (sign) val = -val;

	int digit = 0;
	do {
		*(--p) = (val % 10) + '0';
		val /= 10;
		if (((++digit) >= 3) && (val > 0))
		{
			*(--p) = ',';
			digit = 0;
		}
	} while (val > 0);

	if (sign) *(--p) = '-';
	return p;
}

/// Text matching, return -1 when no maching
COREARRAY_DLL_LOCAL int MatchText(const char *txt, const char *list[])
{
	for (int i=0; *list; list++, i++)
	{
		if (strcmp(txt, *list) == 0)
			return i;
	}
	return -1;
}

}
