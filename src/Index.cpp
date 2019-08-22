// ===========================================================
//
// Index.cpp: Indexing Objects
//
// Copyright (C) 2016-2019    Xiuwen Zheng
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
#include <stdio.h>
#include <string.h>

using namespace std;


#ifdef NO_R_v3_3
extern "C"
{
extern Rconnection getConnection(int n);
COREARRAY_DLL_LOCAL Rconnection My_R_GetConnection(SEXP x)
	{ return getConnection(Rf_asInteger(x)); }
}
#endif


namespace SeqArray
{

static const char *ERR_INDEX_VALUE =
	"'%s' should not contain negative values or NA (replaced by zero).";


// ===========================================================
// Indexing object
// ===========================================================

CIndex::CIndex()
{
	TotalLength = 0;
	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
}

void CIndex::Init(PdContainer Obj, const char *varname)
{
	Values.clear();
	Lengths.clear();
	int Buffer[65536];
	C_Int64 n = GDS_Array_GetTotalCount(Obj);
	if (n > INT_MAX)
		throw ErrSeqArray("Invalid dimension in CIndex.");

	CdIterator it;
	GDS_Iter_GetStart(Obj, &it);
	TotalLength = n;
	int last = -1;
	C_UInt32 repeat = 0;
	bool if_neg_val = false;

	while (n > 0)
	{
		ssize_t m = (n <= 65536) ? n : 65536;
		GDS_Iter_RData(&it, Buffer, m, svInt32);
		n -= m;
		for (int *p = Buffer; m > 0; m--)
		{
			int v = *p++;
			if (v < 0) { v = 0; if_neg_val = true; }
			if (v == last)
			{
				repeat ++;
			} else {
				if (repeat > 0)
				{
					Values.push_back(last);
					Lengths.push_back(repeat);					
				}
				last = v; repeat = 1;
			}
		}
	}
	if (repeat > 0)
	{
		Values.push_back(last);
		Lengths.push_back(repeat);					
	}

	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
	if (if_neg_val && varname)
		warning(ERR_INDEX_VALUE, varname);
}

void CIndex::InitOne(int num)
{
	Values.clear();
	Values.push_back(1);
	Lengths.clear();
	Lengths.push_back(num);
	TotalLength = num;
	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
}

void CIndex::GetInfo(size_t pos, C_Int64 &Sum, int &Value)
{
	if (pos >= TotalLength)
		throw ErrSeqArray("Invalid position in CIndex.");
	if (pos < Position)
	{
		Position = 0;
		AccSum = 0;
		AccIndex = AccOffset = 0;
	}
	for (; Position < pos; )
	{
		size_t L = Lengths[AccIndex];
		size_t n = L - AccOffset;
		if ((Position + n) <= pos)
		{
			AccSum += (Values[AccIndex] * n);
			AccIndex ++; AccOffset = 0;
		} else {
			n = pos - Position;
			AccSum += (Values[AccIndex] * n);
			AccOffset += n;
		}
		Position += n;
	}
	Sum = AccSum;
	Value = Values[AccIndex];
}

SEXP CIndex::GetLen_Sel(const C_BOOL sel[])
{
	size_t n;
	const C_BOOL *p = (C_BOOL *)vec_i8_cnt_nonzero_ptr((const int8_t *)sel,
		TotalLength, &n);
	SEXP ans = NEW_INTEGER(n);
	if (n > 0)
	{
		int *pV = &Values[0];
		C_UInt32 *pL = &Lengths[0];
		size_t L = *pL;
		// skip non-selection
		for (size_t m=p-sel; m > 0; )
		{
			if (L == 0)
			{
				L = *(++pL); pV ++;
				continue;  // in case, L = 0
			}
			if (L <= m)
			{
				m -= L; L = 0;
			} else {
				L -= m; m = 0;
			}
		}
		// get lengths
		int *pAns = INTEGER(ans);
		while (n > 0)
		{
			if (L == 0)
			{
				L = *(++pL); pV ++;
				continue;  // in case, L = 0
			}
			L--;
			if (*p++)
			{
				*pAns++ = *pV;
				n --;
			}
		}
	}
	return ans;
}

SEXP CIndex::GetLen_Sel(const C_BOOL sel[], int &out_var_start,
	int &out_var_count, vector<C_BOOL> &out_var_sel)
{
	size_t n;
	const C_BOOL *p = (C_BOOL *)vec_i8_cnt_nonzero_ptr((const int8_t *)sel,
		TotalLength, &n);
	SEXP ans = NEW_INTEGER(n);
	out_var_start = 0;
	out_var_count = 0;

	if (n > 0)
	{
		int *pV = &Values[0];
		C_UInt32 *pL = &Lengths[0];
		size_t L = *pL;
		// skip non-selection
		for (size_t m=p-sel; m > 0; )
		{
			if (L == 0)
			{
				L = *(++pL); pV ++;
				continue;  // in case, L = 0
			}
			if (L <= m)
			{
				m -= L; out_var_start += L * (*pV); L = 0;
			} else {
				L -= m; out_var_start += m * (*pV); m = 0;
			}
		}
		sel = p;
		// get the total length
		int *pVV = pV;
		C_UInt32 *pLL = pL;
		size_t LL = L;
		int *pAns = INTEGER(ans);
		for (size_t m=n; m > 0; )
		{
			if (L == 0)
			{
				L = *(++pL); pV ++;
				continue;  // in case, L = 0
			}
			L--;
			out_var_count += (*pV);
			if (*p++)
			{
				*pAns++ = *pV;
				m --;
			}
		}
		// set bool selection
		out_var_sel.resize(out_var_count, TRUE);
		C_BOOL *pB = &out_var_sel[0];
		p = sel; pV = pVV; pL = pLL; L = LL;
		while (n > 0)
		{
			if (L == 0)
			{
				L = *(++pL); pV ++;
				continue;  // in case, L = 0
			}
			L--;
			if (*p++)
			{
				pB += *pV; n --;
			} else {
				for (size_t m=*pV; m > 0; m--) *pB++ = FALSE;
			}
		}
	} else {
		out_var_sel.clear();
	}

	return ans;
}



// ===========================================================

CGenoIndex::CGenoIndex()
{
	TotalLength = 0;
	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
}

void CGenoIndex::Init(PdContainer Obj, const char *varname)
{
	Values.clear();
	Lengths.clear();
	C_UInt16 Buffer[65536];
	C_Int64 n = GDS_Array_GetTotalCount(Obj);
	if (n > INT_MAX)
		throw ErrSeqArray("Invalid dimension in CIndex.");

	CdIterator it;
	GDS_Iter_GetStart(Obj, &it);
	TotalLength = n;
	C_UInt16 last = 0xFFFF;
	C_UInt32 repeat = 0;
	bool if_neg_val = false;

	while (n > 0)
	{
		ssize_t m = (n <= 65536) ? n : 65536;
		GDS_Iter_RData(&it, Buffer, m, svUInt16);
		n -= m;
		for (C_UInt16 *p = Buffer; m > 0; m--)
		{
			C_UInt16 v = *p++;
			if (v < 0) { v = 0; if_neg_val = true; }
			if (v == last)
			{
				repeat ++;
			} else {
				if (repeat > 0)
				{
					Values.push_back(last);
					Lengths.push_back(repeat);					
				}
				last = v; repeat = 1;
			}
		}
	}

	if (repeat > 0)
	{
		Values.push_back(last);
		Lengths.push_back(repeat);					
	}

	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
	if (if_neg_val && varname)
		warning(ERR_INDEX_VALUE, varname);
}

void CGenoIndex::GetInfo(size_t pos, C_Int64 &Sum, C_UInt8 &Value)
{
	if (pos >= TotalLength)
		throw ErrSeqArray("Invalid position in CIndex.");
	if (pos < Position)
	{
		Position = 0;
		AccSum = 0;
		AccIndex = AccOffset = 0;
	}
	for (; Position < pos; )
	{
		size_t L = Lengths[AccIndex];
		size_t n = L - AccOffset;
		if ((Position + n) <= pos)
		{
			AccSum += (Values[AccIndex] * n);
			AccIndex ++; AccOffset = 0;
		} else {
			n = pos - Position;
			AccSum += (Values[AccIndex] * n);
			AccOffset += n;
		}
		Position += n;
	}
	Sum = AccSum;
	Value = Values[AccIndex] & 0x0F;
}



// ===========================================================
// Chromosome Indexing
// ===========================================================

CChromIndex::CChromIndex() { }

void CChromIndex::AddChrom(PdGDSFolder Root)
{
	static const char *warn_msg =
		"Warning: @chrom_rle_val and @chrom_rle_len are not correct, "
		"please call 'seqOptimize(..., target='chromosome') to update "
		"the chromosome indexing.\n";

	PdAbstractArray varVariant = GDS_Node_Path(Root, "variant.id", TRUE);
	C_Int32 NumVariant = GDS_Array_GetTotalCount(varVariant);

	PdAbstractArray varChrom = GDS_Node_Path(Root, "chromosome", TRUE);
	C_Int32 NumChrom = GDS_Array_GetTotalCount(varChrom);

	if ((GDS_Array_DimCnt(varChrom) != 1) || (NumVariant != NumChrom))
		throw ErrSeqArray("Invalid dimension of 'chromosome'.");
	if (NumChrom <= 0) return;

	// initialize Map and RleChr
	Map.clear();
	RleChr.Clear();

	// check whether RLE representation of chromosome
	PdAbstractArray varChrVal = GDS_Node_Path(Root, "@chrom_rle_val", FALSE);
	PdAbstractArray varChrLen = GDS_Node_Path(Root, "@chrom_rle_len", FALSE);
	if (varChrVal && varChrLen)
	{
		int n1 = GDS_Array_GetTotalCount(varChrVal);
		int n2 = GDS_Array_GetTotalCount(varChrLen);
		if (GDS_Array_DimCnt(varChrVal)==1 && GDS_Array_DimCnt(varChrLen)==1 &&
			n1 == n2)
		{
			string val[n1];
			C_Int32 len[n1], idx=0;
			GDS_Array_ReadData(varChrVal, &idx, &n1, &val[0], svStrUTF8);
			GDS_Array_ReadData(varChrLen, &idx, &n1, &len[0], svInt32);

			int ntot = 0;
			for (int i=0; i < n1; i++) ntot += len[i];
			if (ntot == NumVariant)
			{
				TRange rng;
				rng.Start = 0;
				for (int i=0; i < n1; i++)
				{
					rng.Length = len[i];
					Map[val[i]].push_back(rng);
					rng.Start += rng.Length;
					RleChr.Add(val[i], len[i]);
				}
				RleChr.Init();
				return;
			} else
				REprintf(warn_msg);

		} else
			REprintf(warn_msg);
	} else if (varChrVal || varChrLen)
	{
		REprintf(warn_msg);
	}

	C_Int32 idx=0, len=1;
	string last;
	GDS_Array_ReadData(varChrom, &idx, &len, &last, svStrUTF8);
	idx ++;

	TRange rng;
	rng.Start = 0; rng.Length = 1;

	const C_Int32 NMAX = 65536;
	string txt[NMAX];

	while (idx < NumChrom)
	{
		len = NumChrom - idx;
		if (len > NMAX) len = NMAX;
		GDS_Array_ReadData(varChrom, &idx, &len, &txt, svStrUTF8);
		for (int i=0; i < len; i++)
		{
			if (txt[i] == last)
			{
				rng.Length ++;
			} else {
				Map[last].push_back(rng);
				RleChr.Add(last, rng.Length);
				last = string(txt[i].begin(), txt[i].end());
				rng.Start = idx + i;
				rng.Length = 1;
			}
		}
		idx += len;
	}

	Map[last].push_back(rng);
	RleChr.Add(last, rng.Length);
	RleChr.Init();
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
	set<TRange, less_range>::iterator p = _RangeSet.find(rng);
	return (p != _RangeSet.end()) && (p->Start <= point) && (point <= p->End);
}

void CRangeSet::GetRanges(int Start[], int End[])
{
	set<TRange, less_range>::const_iterator it = _RangeSet.begin();
	for (size_t n=_RangeSet.size(); n > 0; n--, it++)
	{
		*Start++ = it->Start;
		*End++   = it->End;
	}
}



// ===========================================================
// SeqArray GDS file information
// ===========================================================

TSelection::TSelection(CFileInfo &File, bool init)
{
	Link = NULL;
	if (File.Ploidy() <= 0)
		throw ErrSeqArray("Unable to determine ploidy.");
	numPloidy = File.Ploidy();
	numSamp = File.SampleNum(); pSample = new C_BOOL[numSamp];
	if (init) memset(pSample, TRUE, numSamp);
	numVar = File.VariantNum(); pVariant = new C_BOOL[numVar];
	if (init) memset(pVariant, TRUE, numVar);
	pFlagGenoSel = NULL;
	varTrueNum = -1; varStart = varEnd = 0;
}

TSelection::~TSelection()
{
	if (pSample)
		{ delete[] pSample; pSample = NULL; }
	if (pVariant)
		{ delete[] pVariant; pVariant = NULL; }
	ClearStructSample();
	Link = NULL;
}

TSelection::TSampStruct *TSelection::GetStructSample()
{
	// the block size considered in the block reading
	static ptrdiff_t block_size = 512;

	if (!pFlagGenoSel)
	{
		const size_t SIZE = numSamp * numPloidy;
		pFlagGenoSel = new C_BOOL[SIZE];  // set the output
		C_BOOL *p = pFlagGenoSel, *s = pSample;
		memset(p, TRUE, SIZE);
		for (size_t n=numSamp; n > 0; n--)
		{
			if (*s++ == FALSE)
			{
				for (size_t m=numPloidy; m > 0; m--) *p++ = FALSE;
			} else {
				p += numPloidy;
			}
		}
	}

	if (pSampList.empty())
	{
		TSampStruct ss, last;
		last.length = last.offset = 0; last.sel = NULL;
		C_BOOL *pSt = pSample, *pEnd = pSample + numSamp;
		// find the first TRUE
		while (pSt<pEnd && !*pSt) pSt++;
		// for-loop all TRUE blocks
		for (; pSt < pEnd; )
		{
			// find the last TRUE of the TRUE block
			C_BOOL *pE = pSt;
			while (pE<pEnd && *pE) pE++;
			// find the first TRUE of the next TRUE block
			C_BOOL *pSt2 = pE;
			while (pSt2<pEnd && !*pSt2) pSt2++;
			//
			if (pE-pSt >= block_size)
			{
				// TRUE block: long enough
				if (last.length > 0)
				{
					pSampList.push_back(last);
					last.length = 0;
				}
				ss.length = (pE - pSt) * numPloidy;
				ss.offset = (pSt - pSample) * numPloidy;
				ss.sel = NULL;
				pSampList.push_back(ss);
			} else if (pSt2-pE >= block_size)
			{
				// TRUE block: short, but far away from the next TRUE block
				if (last.length > 0)
				{
					// merge with the last block
					last.length = (pE - pSample) * numPloidy - last.offset;
					if (!last.sel)
						last.sel = pFlagGenoSel + last.offset;
					pSampList.push_back(last);
					last.length = 0;
				} else {
					ss.length = (pE - pSt) * numPloidy;
					ss.offset = (pSt - pSample) * numPloidy;
					ss.sel = NULL;
					pSampList.push_back(ss);
				}
			} else {
				// sparse
				if (last.length > 0)
				{
					// merge with the last block
					last.length = (pE - pSample) * numPloidy - last.offset;
					if (!last.sel)
						last.sel = pFlagGenoSel + last.offset;
				} else {
					last.length = (pE - pSt) * numPloidy;
					last.offset = (pSt - pSample) * numPloidy;
					last.sel = NULL;
				}
			}
			// at last
			pSt = pSt2;
		}

		// the ending one
		if (last.length > 0)
			pSampList.push_back(last);
		ss.length = ss.offset = 0; ss.sel = NULL;
		pSampList.push_back(ss);
	}

	{
		// check
		ssize_t num = 0;
		TSampStruct *p = &pSampList[0];
		for (; p->length > 0; p++)
		{
			if (p->sel)
				num += GetNumOfTRUE(p->sel, p->length);
			else
				num += p->length;
		}
		ssize_t num_samp = GetNumOfTRUE(pSample, numSamp);
		if (num_samp*int(numPloidy) != num)
			throw ErrSeqArray("Internal error when preparing structure for selected samples, please email to zhengxwen@gmail.com.");
	}

	return &pSampList[0];
}

void TSelection::ClearStructSample()
{
	if (pFlagGenoSel)
	{
		delete[] pFlagGenoSel;
		pFlagGenoSel = NULL;
	}
	pSampList.clear();
}

void TSelection::GetStructVariant()
{
	if (varTrueNum < 0)
	{
		C_BOOL *end = pVariant + numVar;
		C_BOOL *p = VEC_BOOL_FIND_TRUE(pVariant, end);
		varStart = p - pVariant;
		C_BOOL *last = end - 1;
		size_t num = 0;
		for (; p < end; p++)
			if (*p) { num++; last = p; }
		varTrueNum = num;
		varEnd = last + 1 - pVariant;
	}
}

void TSelection::ClearSelectVariant()
{
	if (varTrueNum < 0)
	{
		memset(pVariant, 0, numVar);
	} else {
		memset(pVariant+varStart, 0, varEnd-varStart);
	}
	varTrueNum = varStart = varEnd = 0;
}

void TSelection::ClearStructVariant()
{
	varTrueNum = -1;
	varStart = varEnd = 0;
}


//

static const char *ERR_DIM = "Invalid dimension of '%s'.";
static const char *ERR_FILE_ROOT = "CFileInfo::FileRoot should be initialized.";

CFileInfo::CFileInfo(PdGDSFolder root)
{
	_Root = NULL;
	_SelList = NULL;
	_SampleNum = _VariantNum = 0;
	ResetRoot(root);
}

CFileInfo::~CFileInfo()
{
	_Root = NULL;
	_SampleNum = _VariantNum = 0;
	clear_selection();
}

void CFileInfo::clear_selection()
{
	for (TSelection *p=_SelList; p != NULL; )
	{
		TSelection *n = p;
		p = p->Link;
		delete n;
	}
	_SelList = NULL;
}

void CFileInfo::ResetRoot(PdGDSFolder root)
{
	if (_Root != root)
	{
		// initialize
		_Root = root;
		_Chrom.Clear();
		_Position.clear();
		clear_selection();

		// sample.id
		PdAbstractArray Node = GDS_Node_Path(root, "sample.id", TRUE);
		C_Int64 n = GDS_Array_GetTotalCount(Node);
		if ((n < 0) || (n > 2147483647))
			throw ErrSeqArray(ERR_DIM, "sample.id");
		_SampleNum = n;

		// variant.id
		Node = GDS_Node_Path(root, "variant.id", TRUE);
		n = GDS_Array_GetTotalCount(Node);
		if ((n < 0) || (n > 2147483647))
			throw ErrSeqArray(ERR_DIM, "variant.id");
		_VariantNum = n;

		// genotypes
		_Ploidy = -1;
		Node = GDS_Node_Path(root, "genotype/data", FALSE);
		if (Node != NULL)
		{
			if (GDS_Array_DimCnt(Node) == 3)
			{
				C_Int32 DLen[3];
				GDS_Array_GetDim(Node, DLen, 3);
				_Ploidy = DLen[2];
			}
		} else
			_Ploidy = 2;

		// initialize selection
		_SelList = new TSelection(*this, true);
	}
}

TSelection &CFileInfo::Selection()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	return *_SelList;
}

TSelection &CFileInfo::Push_Selection(bool init_samp, bool init_var)
{
	TSelection *n = new TSelection(*this, false);
	n->Link = _SelList;
	if (init_samp)
		memcpy(n->pSample, _SelList->pSample, _SampleNum);
	if (init_var)
		memcpy(n->pVariant, _SelList->pVariant, _VariantNum);
	_SelList = n;
	return *n;
}

void CFileInfo::Pop_Selection()
{
	if (_SelList==NULL || _SelList->Link==NULL)
		throw ErrSeqArray("No filter can be pop up.");
	TSelection *n = _SelList;
	_SelList = n->Link;
	delete n;
}

CChromIndex &CFileInfo::Chromosome()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	if (_Chrom.Empty())
		_Chrom.AddChrom(_Root);
	return _Chrom;
}

void CFileInfo::ResetChromosome()
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	_Chrom.Clear();
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
		_Position.resize(_VariantNum);
		GDS_Array_ReadData(N, NULL, NULL, &_Position[0], svInt32);
	}
	return _Position;
}

CGenoIndex &CFileInfo::GenoIndex()
{
	if (_GenoIndex.Empty())
	{
		static const char *geno_var = "genotype/@data";
		PdAbstractArray I = GetObj(geno_var, TRUE);
		_GenoIndex.Init(I, geno_var);
	}
	return _GenoIndex;
}

CIndex &CFileInfo::VarIndex(const string &varname)
{
	CIndex &I = _VarIndex[varname];
	if (I.Empty())
	{
		PdAbstractArray N = GDS_Node_Path(_Root, varname.c_str(), FALSE);
		if (N == NULL)
			I.InitOne(_VariantNum);
		else
			I.Init(N, varname.c_str());
	}
	return I;
}

PdAbstractArray CFileInfo::GetObj(const char *name, C_BOOL MustExist)
{
	if (!_Root)
		throw ErrSeqArray(ERR_FILE_ROOT);
	return GDS_Node_Path(_Root, name, MustExist);
}

int CFileInfo::SampleSelNum()
{
	TSelection &s = Selection();
	return vec_i8_cnt_nonzero((C_Int8*)s.pSample, _SampleNum);
}

int CFileInfo::VariantSelNum()
{
	TSelection &s = Selection();
	s.GetStructVariant();
	return s.varTrueNum;
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

static C_BOOL ArrayTRUEs[64] = {
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1
};

CVarApply::CVarApply()
{
	fVarType = ctNone;
	MarginalStart = MarginalEnd = 0;
	MarginalSelect = NULL;
	Node = NULL;
	Position = 0;
}

CVarApply::~CVarApply()
{ }

void CVarApply::Reset()
{
	Position = MarginalStart;
	if ((MarginalEnd - MarginalStart) > 0)
		if (!MarginalSelect[MarginalStart]) Next();
}

bool CVarApply::Next()
{
	C_BOOL *p = MarginalSelect;
	Position = VEC_BOOL_FIND_TRUE(p+Position+1, p+MarginalEnd) - p;
	return (Position < MarginalEnd);
}

C_BOOL *CVarApply::NeedTRUEs(size_t size)
{
	if (size <= sizeof(ArrayTRUEs))
	{
		return ArrayTRUEs;
	} else if (size > _TRUE.size())
	{
		_TRUE.resize(size, TRUE);
	}
	return &_TRUE[0];
}


CApply_Variant::CApply_Variant(): CVarApply()
{ }

CApply_Variant::CApply_Variant(CFileInfo &File): CVarApply()
{
	InitMarginal(File);
}

void CApply_Variant::InitMarginal(CFileInfo &File)
{
	TSelection &Sel = File.Selection();
	Sel.GetStructVariant();
	MarginalStart = Sel.varStart;
	MarginalEnd = Sel.varEnd;
	MarginalSelect = Sel.pVariant;
}


CVarApplyList::~CVarApplyList()
{
	for (iterator p = begin(); p != end(); p++)
	{
		CVarApply *v = (*p);
		*p = NULL;
		delete v;
	}
}

bool CVarApplyList::CallNext()
{
	bool has_next = true;
	for (iterator p = begin(); p != end(); p++)
	{
		if (!(*p)->Next())
			has_next = false;
	}
	return has_next;
}



// ===========================================================
// Progress object
// ===========================================================

static const int PROGRESS_BAR_CHAR_NUM = 50;
static const int PROGRESS_LINE_NUM = 100000;

static const double S_MIN  =  60;
static const double S_HOUR =  60 * S_MIN;
static const double S_DAY  =  24 * S_HOUR;
static const double S_YEAR = 365 * S_DAY;

static const char *time_str(double s)
{
	if (R_FINITE(s))
	{
		static char buffer[64];
		if (s < S_MIN)
			sprintf(buffer, "%.0fs", s);
		else if (s < S_HOUR)
			sprintf(buffer, "%.1fm", s/S_MIN);
		else if (s < S_DAY)
			sprintf(buffer, "%.1fh", s/S_HOUR);
		else if (s < S_YEAR)
			sprintf(buffer, "%.1fd", s/S_DAY);
		else
			sprintf(buffer, "%.1f years", s/S_YEAR);
		return buffer;
	} else
		return "---";
}


CProgress::CProgress(C_Int64 start, C_Int64 count, SEXP conn, bool newline)
{
	TotalCount = count;
	Counter = (start >= 0) ? start : 0;
	double percent;
	File = NULL;
	if (conn && !Rf_isNull(conn))
		File = R_GetConnection(conn);
	NewLine = newline;

	if (count > 0)
	{
		int n = 100;
		if (n > count) n = count;
		if (n < 1) n = 1;
		_start = _step = (double)count / n;
		_hit = (C_Int64)(_start);
		if (Counter > count) Counter = count;
		percent = (double)Counter / count;
	} else {
		_start = _step = 0;
		_hit = PROGRESS_LINE_NUM;
		percent = 0;
	}

	time_t s; time(&s);
	_start_time = s;
	_timer.reserve(128);
	_timer.push_back(pair<double, time_t>(percent, s));

	ShowProgress();
}

CProgress::~CProgress()
{ }

void CProgress::Forward(C_Int64 Inc)
{
	Counter += Inc;
	if (TotalCount > 0 && Counter > TotalCount)
		Counter = TotalCount;
	if (Counter >= _hit)
	{
		if (TotalCount > 0)
		{
			while (Counter >= _hit)
			{
				_start += _step;
				_hit = (C_Int64)(_start);
			}
			if (_hit > TotalCount) _hit = TotalCount;
		} else {
			while (Counter >= _hit) _hit += PROGRESS_LINE_NUM;
		}
		ShowProgress();
	}
}

void CProgress::ShowProgress()
{
	if (File)
	{
		if (TotalCount > 0)
		{
			char bar[PROGRESS_BAR_CHAR_NUM + 1];
			double p = (double)Counter / TotalCount;
			int n = (int)round(p * PROGRESS_BAR_CHAR_NUM);
			memset(bar, '.', sizeof(bar));
			memset(bar, '=', n);
			if ((Counter > 0) && (n < PROGRESS_BAR_CHAR_NUM))
				bar[n] = '>';
			bar[PROGRESS_BAR_CHAR_NUM] = 0;

			// ETC: estimated time to complete
			n = (int)_timer.size() - 20;  // 20% as a sliding window size
			if (n < 0) n = 0;
			time_t now; time(&now);
			_timer.push_back(pair<double, time_t>(p, now));

			// in seconds
			double s = difftime(now, _timer[n].second);
			double diff = p - _timer[n].first;
			if (diff > 0)
				s = s / diff * (1 - p);
			else
				s = R_NaN;
			p *= 100;

			// show
			if (NewLine)
			{
				ConnPutText(File, "[%s] %2.0f%%, ETC: %s", bar, p, time_str(s));
				if (R_Process_Count && R_Process_Index && (*R_Process_Count>1))
					ConnPutText(File, " (process %d)", *R_Process_Index);
				ConnPutText(File, "\n");
			} else {
				ConnPutText(File, "\r[%s] %2.0f%%, ETC: %s", bar, p, time_str(s));
				if (R_Process_Count && R_Process_Index && (*R_Process_Count>1))
					ConnPutText(File, " (process %d)", *R_Process_Index);
				ConnPutText(File, "    ");
				if (Counter >= TotalCount) ConnPutText(File, "\n");
			}
		} else {
			int n = Counter / PROGRESS_LINE_NUM;
			n = (n / 10) + (n % 10 ? 1 : 0);
			string s(n, '.');
			if (NewLine)
			{
				if (Counter > 0)
					ConnPutText(File, "[:%s (%dk lines)]", s.c_str(), Counter/1000);
				else
					ConnPutText(File, "[: (0 line)]");
				if (R_Process_Count && R_Process_Index && (*R_Process_Count>1))
					ConnPutText(File, " (process %d)", *R_Process_Index);
				ConnPutText(File, "\n");
			} else {
				if (Counter > 0)
					ConnPutText(File, "\r[:%s (%dk lines)]", s.c_str(), Counter/1000);
				else
					ConnPutText(File, "\r[: (0 line)]");
				if (R_Process_Count && R_Process_Index && (*R_Process_Count>1))
					ConnPutText(File, " (process %d)", *R_Process_Index);
			}
		}
		(*File->fflush)(File);
	}
}


CProgressStdOut::CProgressStdOut(C_Int64 count, bool verbose):
	CProgress(0, count, NULL, false)
{
	if (count < 0)
		throw ErrSeqArray("%s, 'count' should be greater than zero.", __func__);
	_last_time = _timer.back().second;
	Verbose = verbose;
	ShowProgress();
}

void CProgressStdOut::ShowProgress()
{
	if (Verbose && (TotalCount > 0))
	{
		char bar[PROGRESS_BAR_CHAR_NUM + 1];
		double p = (double)Counter / TotalCount;
		int n = (int)round(p * PROGRESS_BAR_CHAR_NUM);
		memset(bar, '.', sizeof(bar));
		memset(bar, '=', n);
		if ((Counter > 0) && (n < PROGRESS_BAR_CHAR_NUM))
			bar[n] = '>';
		bar[PROGRESS_BAR_CHAR_NUM] = 0;

		// ETC: estimated time to complete
		n = (int)_timer.size() - 20;  // 20% as a sliding window size
		if (n < 0) n = 0;
		time_t now; time(&now);
		_timer.push_back(pair<double, time_t>(p, now));

		// in seconds
		double interval = difftime(now, _last_time);
		double s = difftime(now, _timer[n].second);
		double diff = p - _timer[n].first;
		if (diff > 0)
			s = s / diff * (1 - p);
		else
			s = R_NaN;
		p *= 100;

		// show
		_last_time = now;
		if (Counter >= TotalCount)
		{
			char buffer[512];
			s = difftime(_last_time, _start_time);
			int n = sprintf(buffer, "\r[%s] 100%%, completed, %s", bar, time_str(s));
			if (R_Process_Count && R_Process_Index && (*R_Process_Count>1))
				sprintf(buffer+n, " (process %d)", *R_Process_Index);
			Rprintf("%s\n", buffer);
		} else if ((interval >= 5) || (Counter <= 0))
		{
			char buffer[512];
			_last_time = now;
			int n = sprintf(buffer, "\r[%s] %2.0f%%, ETC: %s", bar, p, time_str(s));
			if ((Counter>0) && R_Process_Count && R_Process_Index && (*R_Process_Count>1))
				n += sprintf(buffer+n, " (process %d)", *R_Process_Index);
			strcpy(buffer+n, "    ");
			Rprintf("%s", buffer);
		}
	}
}



// ===========================================================
// Pre-defined R objects
// ===========================================================

SEXP R_Geno_Dim2_Name = R_NilValue;
SEXP R_Geno_Dim3_Name = R_NilValue;
SEXP R_Dosage_Name    = R_NilValue;
SEXP R_Data_Name      = R_NilValue;
SEXP R_Data_Dim2_Name = R_NilValue;

int* R_Process_Count = NULL;
int* R_Process_Index = NULL;



// ===========================================================
// Define Functions
// ===========================================================

// the buffer of ArrayTRUEs
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


/// Allocate R object given by SVType
COREARRAY_DLL_LOCAL SEXP RObject_GDS(PdAbstractArray Node, size_t n,
	int &nProtected, bool bit1_is_logical)
{
	SEXP ans = R_NilValue;
	C_SVType SVType = GDS_Array_GetSVType(Node);

	if (COREARRAY_SV_INTEGER(SVType))
	{
		char classname[128];
		GDS_Node_GetClassName(Node, classname, sizeof(classname));
		if (strcmp(classname, "dBit1") == 0)
		{
			if (bit1_is_logical)
				PROTECT(ans = NEW_LOGICAL(n));
			else
				PROTECT(ans = NEW_INTEGER(n));
		} else if (GDS_R_Is_Logical(Node))
		{
			PROTECT(ans = NEW_LOGICAL(n));
		} else {
			PROTECT(ans = NEW_INTEGER(n));
			nProtected += GDS_R_Set_IfFactor(Node, ans);
		}
		nProtected ++;
	} else if (COREARRAY_SV_FLOAT(SVType))
	{
		PROTECT(ans = NEW_NUMERIC(n));
		nProtected ++;
	} else if (COREARRAY_SV_STRING(SVType))
	{
		PROTECT(ans = NEW_CHARACTER(n));
		nProtected ++;
	}

	return ans;
}


/// Append data to a GDS node
COREARRAY_DLL_LOCAL void RAppendGDS(PdAbstractArray Node, SEXP Val)
{
	switch (TYPEOF(Val))
	{
	case LGLSXP:  // logical
		GDS_Array_AppendData(Node, XLENGTH(Val), LOGICAL(Val), svInt32);
		break;
	case INTSXP:  // integer
		GDS_Array_AppendData(Node, XLENGTH(Val), INTEGER(Val), svInt32);
		break;
	case REALSXP:  // numeric
		GDS_Array_AppendData(Node, XLENGTH(Val), REAL(Val), svFloat64);
		break;
	case RAWSXP:  // RAW
		GDS_Array_AppendData(Node, XLENGTH(Val), RAW(Val), svInt8);
		break;
	case STRSXP:  // character
		{
			R_xlen_t n = XLENGTH(Val);
			vector<string> buf(n);
			for (R_xlen_t i=0; i < n; i++)
			{
				SEXP s = STRING_ELT(Val, i);
				if (s != NA_STRING) buf[i] = translateCharUTF8(s);
			}
			GDS_Array_AppendData(Node, n, &buf[0], svStrUTF8);
		}
		break;
	default:
		throw ErrSeqArray("the user-defined function should return a vector.");
	}
}


COREARRAY_DLL_LOCAL C_BOOL *NeedArrayTRUEs(size_t len)
{
	if (len <= sizeof(ArrayTRUEs))
		return ArrayTRUEs;
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


/// Get the number of alleles
COREARRAY_DLL_LOCAL int GetNumOfAllele(const char *allele_list)
{
	int n = 0;
	while (*allele_list)
	{
		if (*allele_list != ',')
		{
			n ++;
			while ((*allele_list != ',') && (*allele_list != 0))
				allele_list ++;
			if (*allele_list == ',')
			{
				allele_list ++;
				if (*allele_list == 0)
				{
					n ++;
					break;
				}
			}
		}
	}
	return n;
}


/// Get the index in an allele list
COREARRAY_DLL_LOCAL int GetIndexOfAllele(const char *allele, const char *allele_list)
{
	const size_t len = strlen(allele);
	const char *st = allele_list;
	int idx = 0;
	while (*allele_list)
	{
		while ((*allele_list != ',') && (*allele_list != 0))
			allele_list ++;
		size_t n = allele_list - st;
		if ((len==n) && (strncmp(allele, st, n)==0))
			return idx;
		if (*allele_list == ',')
		{
			idx ++;
			allele_list ++;
			st = allele_list;
		}
	}
	return -1;
}


/// Get strings split by comma
COREARRAY_DLL_LOCAL void GetAlleles(const char *alleles, vector<string> &out)
{
	out.clear();
	const char *p, *s;
	p = s = alleles;
	do {
		if ((*p == 0) || (*p == ','))
		{
			out.push_back(string(s, p));
			if (*p == ',') p ++;
			s = p;
			if (*p == 0) break;
		}
		p ++;
	} while (1);
}


/// get PdGDSObj from a SEXP object
COREARRAY_DLL_LOCAL void GDS_PATH_PREFIX_CHECK(const char *path)
{
	for (; *path != 0; path++)
	{
		if ((*path == '~') || (*path == '@'))
		{
			throw SeqArray::ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.",
				*path);
		}
	}
}


COREARRAY_DLL_LOCAL void GDS_VARIABLE_NAME_CHECK(const char *p)
{
	for (; *p != 0; p++)
	{
		if ((*p == '~') || (*p == '@') || (*p == '/'))
		{
			throw ErrSeqArray(
				"the variable name contains an invalid prefix '%c'.", *p);
		}
	}
}


/// get PdGDSObj from a SEXP object
COREARRAY_DLL_LOCAL string GDS_PATH_PREFIX(const string &path, char prefix)
{
	string s = path;
	for (int i=s.size()-1; i >= 0; i--)
	{
		if (s[i] == '/')
		{
			if (((int)s.size() > i+1) && (s[i+1] == '~'))
				s[i+1] = prefix;
			else
				s.insert(i+1, &prefix, 1);
			return s;
		}
	}

	if ((s.size() > 0) && (s[0] == '~'))
		s[0] = prefix;
	else
		s.insert(s.begin(), prefix);

	return s;
}


/// output to a connection
COREARRAY_DLL_LOCAL void ConnPutText(Rconnection file, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*file->vfprintf)(file, fmt, args);
	va_end(args);
}

}
