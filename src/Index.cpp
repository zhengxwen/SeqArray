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

void CIndex::Init(PdContainer Obj)
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

	while (n > 0)
	{
		ssize_t m = (n <= 65536) ? n : 65536;
		GDS_Iter_RData(&it, Buffer, m, svInt32);
		n -= m;
		for (int *p = Buffer; m > 0; m--)
		{
			int v = *p++;
			if (v < 0) v = 0;
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


// ===========================================================

CGenoIndex::CGenoIndex()
{
	TotalLength = 0;
	Position = 0;
	AccSum = 0;
	AccIndex = AccOffset = 0;
}

void CGenoIndex::Init(PdContainer Obj)
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

	while (n > 0)
	{
		ssize_t m = (n <= 65536) ? n : 65536;
		GDS_Iter_RData(&it, Buffer, m, svUInt16);
		n -= m;
		for (C_UInt16 *p = Buffer; m > 0; m--)
		{
			C_UInt16 v = *p++;
			if (v < 0) v = 0;
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
	PosToChr.Clear();

	const C_Int32 NMAX = 4096;
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
				PosToChr.Add(last, rng.Length);
				last = string(txt[i].begin(), txt[i].end());
				rng.Start = idx + i;
				rng.Length = 1;
			}
		}
		idx += len;
	}

	Map[last].push_back(rng);
	PosToChr.Add(last, rng.Length);
	PosToChr.Init();
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
		}
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
	if (_Chrom.Empty())
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
		_Position.resize(_VariantNum);
		GDS_Array_ReadData(N, NULL, NULL, &_Position[0], svInt32);
	}
	return _Position;
}

CGenoIndex &CFileInfo::GenoIndex()
{
	if (_GenoIndex.Empty())
	{
		PdAbstractArray I = GetObj("genotype/@data", TRUE);
		_GenoIndex.Init(I);
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
			I.Init(N);
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

static C_BOOL ArrayTRUEs[64] = {
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1
};

CVarApply::CVarApply()
{
	fVarType = ctNone;
	MarginalSize = 0;
	MarginalSelect = NULL;
	Node = NULL;
	Position = 0;
}

CVarApply::~CVarApply()
{ }

void CVarApply::Reset()
{
	Position = 0;
	if (MarginalSize > 0)
		if (!MarginalSelect[0]) Next();
}

bool CVarApply::Next()
{
	C_BOOL *p = MarginalSelect + Position;
	while (Position < MarginalSize)
	{
		Position ++;
		if (*(++p)) break;
	}
	return (Position < MarginalSize);
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
	MarginalSize = File.VariantNum();
	MarginalSelect = File.Selection().pVariant();
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

static int Progress_ShowNum = 50;
static int Progress_Line_Num = 100000;

inline static void put_text(Rconnection conn, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*conn->vfprintf)(conn, fmt, args);
	va_end(args);
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
		_hit = Progress_Line_Num;
		percent = 0;
	}

	time_t s; time(&s);
	_timer.reserve(128);
	_timer.push_back(pair<double, time_t>(percent, s));

	ShowProgress();
}

CProgress::~CProgress()
{ }

void CProgress::Forward()
{
	Counter ++;
	if (Counter >= _hit)
	{
		if (TotalCount > 0)
		{
			_start += _step;
			_hit = (C_Int64)(_start);
		} else {
			_hit += Progress_Line_Num;
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
			char ss[Progress_ShowNum + 1];
			double percent = (double)Counter / TotalCount;
			int n = (int)round(percent * Progress_ShowNum);
			memset(ss, '.', sizeof(ss));
			memset(ss, '=', n);
			if (n < Progress_ShowNum) ss[n] = '>';
			ss[Progress_ShowNum] = 0;

			// ETC: estimated time to complete
			n = (int)_timer.size() - 20;  // 20% as a sliding window size
			if (n < 0) n = 0;
			time_t now; time(&now);
			_timer.push_back(pair<double, time_t>(percent, now));

			double sec = difftime(now, _timer[n].second);
			double diff = percent - _timer[n].first;
			if (diff > 0)
				sec = sec / diff * (1 - percent);
			else
				sec = 999.9 * 60 * 60;
			percent *= 100;

			// show
			if (NewLine)
			{
				if (sec < 60)
				{
					put_text(File, "[%s] %2.0f%%, ETC: %.0fs\n", ss,
						percent, sec);
				} else if (sec < 3600)
				{
					put_text(File, "[%s] %2.0f%%, ETC: %.1fm\n", ss,
						percent, sec/60);
				} else {
					put_text(File, "[%s] %2.0f%%, ETC: %.1fh\n", ss,
						percent, sec/(60*60));
				}
			} else {
				if (sec < 60)
				{
					put_text(File, "\r[%s] %2.0f%%, ETC: %.0fs  ", ss,
						percent, sec);
				} else if (sec < 3600)
				{
					put_text(File, "\r[%s] %2.0f%%, ETC: %.1fm  ", ss,
						percent, sec/60);
				} else {
					put_text(File, "\r[%s] %2.0f%%, ETC: %.1fh  ", ss,
						percent, sec/(60*60));
				}
				if (Counter >= TotalCount)
					put_text(File, "\n");
			}
		} else {
			int n = Counter / Progress_Line_Num;
			string s(n, '.');
			if (NewLine)
			{
				if (Counter > 0)
					put_text(File, "[:%s (%dk lines)]\n", s.c_str(), Counter/1000);
				else
					put_text(File, "[: (0 line)]\n");
			} else {
				if (Counter > 0)
					put_text(File, "\r[:%s (%dk lines)]", s.c_str(), Counter/1000);
				else
					put_text(File, "\r[: (0 line)]");
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
		char ss[Progress_ShowNum + 1];
		double percent = (double)Counter / TotalCount;
		int n = (int)round(percent * Progress_ShowNum);
		memset(ss, '.', sizeof(ss));
		memset(ss, '=', n);
		if (n < Progress_ShowNum) ss[n] = '>';
		ss[Progress_ShowNum] = 0;

		// ETC: estimated time to complete
		n = (int)_timer.size() - 20;  // 20% as a sliding window size
		if (n < 0) n = 0;
		time_t now; time(&now);
		_timer.push_back(pair<double, time_t>(percent, now));

		double int_sec = difftime(now, _last_time);
		double sec = difftime(now, _timer[n].second);
		double diff = percent - _timer[n].first;
		if (diff > 0)
			sec = sec / diff * (1 - percent);
		else
			sec = 999.9 * 60 * 60;
		percent *= 100;

		// show
		if (Counter >= TotalCount)
		{
			Rprintf("\r[%s] 100%%, completed  \n", ss);
		} else if ((int_sec >= 5) || (Counter <= 0))
		{
			_last_time = now;
			if (sec < 60)
			{
				Rprintf("\r[%s] %2.0f%%, ETC: %.0fs  ", ss, percent, sec);
			} else if (sec < 3600)
			{
				Rprintf("\r[%s] %2.0f%%, ETC: %.1fm  ", ss, percent, sec/60);
			} else {
				if (sec >= 999 * 60 * 60)
					Rprintf("\r[%s] %2.0f%%, ETC: NA    ", ss, percent);
				else
					Rprintf("\r[%s] %2.0f%%, ETC: %.1fh  ", ss, percent, sec/(60*60));
			}
		}
	}
}



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

}
