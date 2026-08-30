// ===========================================================
//
// CompressIO.cpp: reading and writing the compressed files
//
// Copyright (C) 2026    Xiuwen Zheng
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

#include "CompressIO.h"
#include "Index.h"

#include <string.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <zlib.h>

#include <Rdefines.h>


namespace SeqArray
{

using namespace std;

// ===========================================================
// The BGZF block header
// ===========================================================

/// the fixed part of a BGZF block header, the last two bytes are BSIZE
static const unsigned char BGZF_HEAD[18] = {
	0x1F, 0x8B, 0x08, 0x04, 0, 0, 0, 0, 0, 0xFF, 0x06, 0x00,
	'B', 'C', 0x02, 0x00, 0, 0
};

/// the 28-byte EOF marker, i.e., an empty BGZF block
static const unsigned char BGZF_EOF[28] = {
	0x1F, 0x8B, 0x08, 0x04, 0, 0, 0, 0, 0, 0xFF, 0x06, 0x00,
	'B', 'C', 0x02, 0x00, 0x1B, 0x00, 0x03, 0x00,
	0, 0, 0, 0, 0, 0, 0, 0
};

/// whether 'head' is a valid BGZF block header
static inline bool bgzf_check_head(const unsigned char *h)
{
	return h[0]==0x1F && h[1]==0x8B && h[2]==0x08 && (h[3] & 0x04) &&
		h[12]=='B' && h[13]=='C' && h[14]==0x02 && h[15]==0x00;
}

bool BGZF_IsValid(const char *fn)
{
	FILE *f = fopen(fn, "rb");
	if (!f) return false;
	// the header of the first block
	unsigned char h[28];
	bool rv = (fread(h, 1, 18, f) == 18) && bgzf_check_head(h);
	// the EOF marker at the end of file: a file consisting of BGZF blocks
	//   always ends with it, so a file mixing BGZF and other gzip members
	//   (e.g., 'cat a.bgz b.gz') is not treated as a BGZF file
	if (rv)
	{
	#ifdef _WIN32
		rv = (_fseeki64(f, -(int64_t)sizeof(BGZF_EOF), SEEK_END) == 0);
	#else
		rv = (fseeko(f, -(off_t)sizeof(BGZF_EOF), SEEK_END) == 0);
	#endif
		if (rv)
		{
			rv = (fread(h, 1, sizeof(BGZF_EOF), f) == sizeof(BGZF_EOF)) &&
				(memcmp(h, BGZF_EOF, sizeof(BGZF_EOF)) == 0);
		}
	}
	fclose(f);
	return rv;
}


// ===========================================================
// CBgzfReader
// ===========================================================

CBgzfReader::CBgzfReader()
{
	fFile = NULL; fAddr = fNext = 0;
}

CBgzfReader::~CBgzfReader()
{
	Close();
}

void CBgzfReader::Open(const char *fn, int64_t addr)
{
	Close();
	fFile = fopen(fn, "rb");
	if (!fFile)
		throw ErrSeqArray("Cannot open '%s'.", fn);
	fBuffer.resize(BGZF_MAX_BLOCK);
	fAddr = fNext = addr;
}

void CBgzfReader::Close()
{
	if (fFile)
	{
		fclose(fFile); fFile = NULL;
		fBuffer.clear();
		vector<unsigned char>().swap(fBuffer);
	}
	fAddr = fNext = 0;
}

size_t CBgzfReader::ReadHead(int64_t addr, unsigned char head[18])
{
#ifdef _WIN32
	if (_fseeki64(fFile, addr, SEEK_SET) != 0)
#else
	if (fseeko(fFile, addr, SEEK_SET) != 0)
#endif
		throw ErrSeqArray("Failed to seek the BGZF file.");
	size_t n = fread(head, 1, 18, fFile);
	if (n == 0) return 0;  // the end of file
	if (n < 18 || !bgzf_check_head(head))
	{
		throw ErrSeqArray("Invalid BGZF block at the offset %lld.",
			(long long int)addr);
	}
	// BSIZE: the total size of the block minus one
	return (size_t)(head[16] | (head[17] << 8)) + 1;
}

size_t CBgzfReader::ReadBlock(void *dst)
{
	unsigned char head[18];
	while (true)
	{
		fAddr = fNext;
		size_t bsize = ReadHead(fAddr, head);
		if (bsize == 0) return 0;  // the end of file
		fNext = fAddr + bsize;
		const int64_t addr = fAddr;
		// read the rest of the block
		size_t n = bsize - 18;
		if (n < 8 || fread(&fBuffer[0], 1, n, fFile) != n)
		{
			throw ErrSeqArray("Truncated BGZF block at the offset %lld.",
				(long long int)fAddr);
		}
		// the last 8 bytes are CRC32 and ISIZE
		z_stream zs;
		memset(&zs, 0, sizeof(zs));
		zs.next_in = &fBuffer[0];         zs.avail_in = n - 8;
		zs.next_out = (unsigned char*)dst;  zs.avail_out = BGZF_MAX_BLOCK;
		if (inflateInit2(&zs, -15) != Z_OK)
			throw ErrSeqArray("Failed to initialize the zlib inflation.");
		int r = inflate(&zs, Z_FINISH);
		size_t rv = zs.total_out;
		inflateEnd(&zs);
		if (r != Z_STREAM_END)
		{
			throw ErrSeqArray("Failed to decompress the BGZF block at the "
				"offset %lld.", (long long int)fAddr);
		}
		if (rv > 0) return rv;
		// an empty block: if it is the last one (i.e., the EOF marker), stop
		//   here, so that the virtual offset is the same as htslib's
		//   bgzf_tell(); otherwise skip it and read the next block
		unsigned char h[18];
		if (ReadHead(fNext, h) == 0)
		{
			fAddr = fNext = addr;
			return 0;
		}
	}
}


// ===========================================================
// CBgzfWriter
// ===========================================================

CBgzfWriter::CBgzfWriter()
{
	fFile = NULL; fLevel = Z_DEFAULT_COMPRESSION; fSize = 0;
}

CBgzfWriter::~CBgzfWriter()
{
	if (fFile) { fclose(fFile); fFile = NULL; }
}

void CBgzfWriter::Create(const char *fn, int level)
{
	Close();
	fFile = fopen(fn, "wb");
	if (!fFile)
		throw ErrSeqArray("Cannot create '%s'.", fn);
	fLevel = level;
	fBuffer.resize(BGZF_BLOCK_SIZE);
	fCBuffer.resize(BGZF_MAX_BLOCK);
	fSize = 0;
}

void CBgzfWriter::Flush()
{
	if (!fFile || fSize<=0) return;

	// compress the raw deflate stream
	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	zs.next_in = &fBuffer[0];
	zs.avail_in = fSize;
	zs.next_out = &fCBuffer[0] + 18;
	zs.avail_out = fCBuffer.size() - 18 - 8;
	if (deflateInit2(&zs, fLevel, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw ErrSeqArray("Failed to initialize the zlib deflation.");
	int r = deflate(&zs, Z_FINISH);
	size_t clen = zs.total_out;
	deflateEnd(&zs);
	if (r != Z_STREAM_END)
		throw ErrSeqArray("Failed to compress a BGZF block.");

	// the block header, with BSIZE = the total size of the block minus one
	const size_t bsize = clen + 18 + 8;
	memcpy(&fCBuffer[0], BGZF_HEAD, 18);
	fCBuffer[16] = (unsigned char)((bsize-1) & 0xFF);
	fCBuffer[17] = (unsigned char)(((bsize-1) >> 8) & 0xFF);
	// the trailer: CRC32 and ISIZE
	unsigned int crc = crc32(crc32(0L, NULL, 0), &fBuffer[0], fSize);
	unsigned char *p = &fCBuffer[0] + 18 + clen;
	p[0] = crc & 0xFF;  p[1] = (crc >> 8) & 0xFF;
	p[2] = (crc >> 16) & 0xFF;  p[3] = (crc >> 24) & 0xFF;
	p[4] = fSize & 0xFF;  p[5] = (fSize >> 8) & 0xFF;
	p[6] = (fSize >> 16) & 0xFF;  p[7] = (fSize >> 24) & 0xFF;

	if (fwrite(&fCBuffer[0], 1, bsize, fFile) != bsize)
		throw ErrSeqArray("Failed to write a BGZF block.");
	fSize = 0;
}

void CBgzfWriter::Write(const void *buf, size_t size)
{
	const unsigned char *s = (const unsigned char *)buf;
	while (size > 0)
	{
		size_t n = fBuffer.size() - fSize;
		if (n > size) n = size;
		memcpy(&fBuffer[fSize], s, n);
		fSize += n; s += n; size -= n;
		if (fSize >= fBuffer.size()) Flush();
	}
}

void CBgzfWriter::Close()
{
	if (fFile)
	{
		Flush();
		// the EOF marker
		fwrite(BGZF_EOF, 1, sizeof(BGZF_EOF), fFile);
		fclose(fFile); fFile = NULL;
		fBuffer.clear();  vector<unsigned char>().swap(fBuffer);
		fCBuffer.clear(); vector<unsigned char>().swap(fCBuffer);
		fSize = 0;
	}
}


// ===========================================================
// Splitting a BGZF file at the block boundaries
// ===========================================================

void BGZF_Split(const char *fn, size_t num, vector<TBgzfPart> &parts)
{
	parts.clear();
	if (num < 1) num = 1;
	FILE *f = fopen(fn, "rb");
	if (!f) throw ErrSeqArray("Cannot open '%s'.", fn);

	// the file size
#ifdef _WIN32
	_fseeki64(f, 0, SEEK_END);
	const int64_t fsize = _ftelli64(f);
#else
	fseeko(f, 0, SEEK_END);
	const int64_t fsize = ftello(f);
#endif

	// scan the block headers, and record the block boundaries closest to
	//   fsize/num, 2*fsize/num, ...; only 18 bytes are read for each block
	size_t k = 1;
	int64_t addr = 0, prev = 0;
	while (k < num)
	{
		const int64_t target = (int64_t)((double)fsize * k / num);
		// move to the first block at or after 'target'
		while (addr < target)
		{
		#ifdef _WIN32
			if (_fseeki64(f, addr, SEEK_SET) != 0) break;
		#else
			if (fseeko(f, addr, SEEK_SET) != 0) break;
		#endif
			unsigned char h[18];
			if (fread(h, 1, 18, f) != 18) { addr = fsize; break; }
			if (!bgzf_check_head(h))
			{
				fclose(f);
				throw ErrSeqArray("Invalid BGZF block at the offset %lld.",
					(long long int)addr);
			}
			prev = addr;
			addr += (int64_t)(h[16] | (h[17] << 8)) + 1;
		}
		if (addr >= fsize) break;
		TBgzfPart p;
		p.open_addr = prev; p.first_addr = addr; p.end_addr = 0;
		parts.push_back(p);
		k ++;
	}
	fclose(f);

	// the first part starts from the beginning of the file
	TBgzfPart p0;
	p0.open_addr = p0.first_addr = 0; p0.end_addr = fsize;
	parts.insert(parts.begin(), p0);
	// the end of each part is the beginning of the next one
	for (size_t i=0; i+1 < parts.size(); i++)
		parts[i].end_addr = parts[i+1].first_addr;
	parts[parts.size()-1].end_addr = fsize;
}


// ===========================================================
// Reading a BGZF file line by line
// ===========================================================

/// Reading a BGZF file line by line, tracking the virtual offsets
class CBgzfLineReader
{
public:
	CBgzfLineReader(const char *fn)
	{
		fReader.Open(fn);
		fBuffer.resize(BGZF_MAX_BLOCK);
		fPos = fEnd = 0; fEOF = false;
	}

	/// the virtual offset of the current position; the next block is not
	///   read here, so that the offsets are the same as those of htslib
	inline uint64_t VOffset() const
	{
		return ((uint64_t)fReader.Addr() << 16) | (uint64_t)fPos;
	}

	/// read a line without the end-of-line character(s), return false at EOF
	bool ReadLine(string &line)
	{
		line.clear();
		if ((fPos >= fEnd) && !Fill()) return false;
		while (true)
		{
			size_t i = fPos;
			while ((i < fEnd) && (fBuffer[i] != '\n')) i ++;
			line.append((const char*)&fBuffer[fPos], i - fPos);
			if (i < fEnd)
			{
				fPos = i + 1;  // skip '\n'
				break;
			}
			fPos = i;
			if (!Fill()) break;
		}
		// remove the trailing '\r', if any
		if (!line.empty() && (line[line.size()-1] == '\r'))
			line.resize(line.size()-1);
		return true;
	}

private:
	CBgzfReader fReader;
	std::vector<unsigned char> fBuffer;
	size_t fPos, fEnd;
	bool fEOF;

	bool Fill()
	{
		if (fEOF) return false;
		size_t n = fReader.ReadBlock(&fBuffer[0]);
		fPos = 0; fEnd = n;
		if (n <= 0) { fEOF = true; return false; }
		return true;
	}
};


// ===========================================================
// The CSI index
// ===========================================================

// the CSI parameters used by tabix, see htslib: n_lvls =
//   (TBX_MAX_SHIFT - min_shift + 2) / 3, where TBX_MAX_SHIFT = 37
static const int CSI_MIN_SHIFT = 14;
static const int CSI_N_LVLS = 8;
static const uint32_t CSI_NULL_BIN = 0xFFFFFFFFu;
/// the pseudo-bin storing the range and the number of the records of a
///   reference sequence, see htslib META_BIN()
static const uint32_t CSI_META_BIN =
	(uint32_t)((((int64_t)1 << (3*(CSI_N_LVLS+1))) - 1) / 7) + 1;

/// a chunk of a BGZF file, in the virtual offsets [u, v)
struct TChunk { uint64_t u, v; };

/// the index of a reference sequence
struct TRefIndex
{
	map< uint32_t, vector<TChunk> > bins;  ///< bin --> the chunks
};

/// the bin of the region [beg, end), see htslib hts_reg2bin()
static inline uint32_t csi_reg2bin(int64_t beg, int64_t end)
{
	int l, s = CSI_MIN_SHIFT;
	int t = ((1 << (CSI_N_LVLS*3)) - 1) / 7;
	for (--end, l=CSI_N_LVLS; l > 0; --l, s += 3, t -= 1 << (l*3))
		if ((beg>>s) == (end>>s)) return t + (uint32_t)(beg>>s);
	return 0;
}

/// add the chunk [u, v) to 'bin'
static inline void csi_add_chunk(TRefIndex &ref, uint32_t bin, uint64_t u,
	uint64_t v)
{
	TChunk c; c.u = u; c.v = v;
	ref.bins[bin].push_back(c);
}

/// get the interval [beg, end) and the chromosome of a VCF line, return false
///   if the line should be skipped; see htslib get_intv() for TBX_VCF:
///   beg is POS-1, end is beg plus the length of REF, and END= in the INFO
///   field overrides end
static bool csi_vcf_intv(const string &line, string &chr, int64_t &beg,
	int64_t &end)
{
	const size_t n = line.size();
	chr.clear(); beg = 0; end = 1;
	size_t st = 0;
	for (int id=1; (id <= 8) && (st <= n); id++)
	{
		size_t ed = line.find('\t', st);
		if (ed == string::npos) ed = n;
		switch (id)
		{
		case 1:  // CHROM
			chr.assign(line, st, ed-st);
			break;
		case 2:  // POS
			beg = strtoll(line.c_str()+st, NULL, 10) - 1;
			end = beg + 1;
			break;
		case 4:  // REF
			end = beg + (int64_t)(ed - st);
			break;
		case 8:  // INFO, look for "END="
			{
				string info(line, st, ed-st);
				size_t k = string::npos;
				if (info.compare(0, 4, "END=") == 0)
				{
					k = 0;
				} else {
					k = info.find(";END=");
					if (k != string::npos) k ++;
				}
				if (k != string::npos)
					end = strtoll(info.c_str()+k+4, NULL, 0);
			}
			break;
		}
		if (ed >= n) break;
		st = ed + 1;
	}
	return !chr.empty();
}

/// write a little-endian 32-bit integer
static inline void csi_w32(CBgzfWriter &wr, uint32_t v)
{
	unsigned char b[4];
	b[0]=v&0xFF; b[1]=(v>>8)&0xFF; b[2]=(v>>16)&0xFF; b[3]=(v>>24)&0xFF;
	wr.Write(b, 4);
}

/// write a little-endian 64-bit integer
static inline void csi_w64(CBgzfWriter &wr, uint64_t v)
{
	unsigned char b[8];
	for (int i=0; i < 8; i++) b[i] = (v >> (i*8)) & 0xFF;
	wr.Write(b, 8);
}

int64_t BGZF_BuildCSI_VCF(const char *fn, const char *fnidx)
{
	CBgzfLineReader rd(fn);

	vector<string> names;          // the reference sequence names
	map<string, int> name2id;      // name --> tid
	vector<TRefIndex> refs;        // the index of each reference sequence

	// the state of the index building, see htslib hts_idx_push()
	int last_tid = -1, save_tid = -1;
	uint32_t last_bin = CSI_NULL_BIN, save_bin = CSI_NULL_BIN;
	uint64_t save_off = 0, off_beg = 0;
	int64_t num_var = 0, n_mapped = 0;

	string line, chr;
	uint64_t off = rd.VOffset();
	while (rd.ReadLine(line))
	{
		const uint64_t cur_off = off;
		off = rd.VOffset();  // the beginning of the next line
		if (line.empty() || line[0]=='#') continue;

		int64_t beg, end;
		if (!csi_vcf_intv(line, chr, beg, end)) continue;
		if (beg < 0) beg = 0;
		if (end <= 0) end = 1;
		num_var ++;

		// the id of the reference sequence
		int tid;
		map<string,int>::iterator it = name2id.find(chr);
		if (it != name2id.end())
		{
			tid = it->second;
		} else {
			tid = (int)names.size();
			names.push_back(chr);
			name2id[chr] = tid;
			refs.resize(names.size());
		}

		if (last_tid != tid)  // the change of a chromosome
		{
			if (last_tid >= 0)
			{
				// the pseudo-bin of the previous chromosome
				csi_add_chunk(refs[last_tid], CSI_META_BIN, off_beg, cur_off);
				csi_add_chunk(refs[last_tid], CSI_META_BIN, n_mapped, 0);
			}
			off_beg = cur_off; n_mapped = 0;
			last_tid = tid; last_bin = CSI_NULL_BIN;
		}
		n_mapped ++;

		const uint32_t bin = csi_reg2bin(beg, end);
		if (last_bin != bin)
		{
			if (save_bin != CSI_NULL_BIN)
				csi_add_chunk(refs[save_tid], save_bin, save_off, cur_off);
			save_off = cur_off;
			save_bin = last_bin = bin;
			save_tid = tid;
		}
	}
	// the virtual offset at the end of file, the same as htslib's bgzf_tell()
	//   after the last (failed) reading
	off = rd.VOffset();

	// finish, see htslib hts_idx_finish()
	if (save_tid >= 0)
	{
		csi_add_chunk(refs[save_tid], save_bin, save_off, off);
		csi_add_chunk(refs[last_tid], CSI_META_BIN, off_beg, off);
		csi_add_chunk(refs[last_tid], CSI_META_BIN, n_mapped, 0);
	}
	for (size_t i=0; i < refs.size(); i++)
	{
		TRefIndex &r = refs[i];
		// merge the adjacent chunks starting from the same BGZF block
		for (map< uint32_t, vector<TChunk> >::iterator p = r.bins.begin();
			p != r.bins.end(); p++)
		{
			if (p->first == CSI_META_BIN) continue;  // not a real bin
			vector<TChunk> &v = p->second;
			size_t m = 0;
			for (size_t l=1; l < v.size(); l++)
			{
				if ((v[m].v >> 16) >= (v[l].u >> 16))
				{
					if (v[m].v < v[l].v) v[m].v = v[l].v;
				} else
					v[++m] = v[l];
			}
			v.resize(m + 1);
		}
	}

	// output
	CBgzfWriter wr;
	wr.Create(fnidx);
	wr.Write("CSI\1", 4);
	csi_w32(wr, CSI_MIN_SHIFT);
	csi_w32(wr, CSI_N_LVLS);
	// the auxiliary data: the tabix configuration and the sequence names
	size_t l_nm = 0;
	for (size_t i=0; i < names.size(); i++) l_nm += names[i].size() + 1;
	csi_w32(wr, (uint32_t)(28 + l_nm));
	csi_w32(wr, 2);    // format: TBX_VCF
	csi_w32(wr, 1);    // col_seq
	csi_w32(wr, 2);    // col_beg
	csi_w32(wr, 0);    // col_end
	csi_w32(wr, '#');  // meta
	csi_w32(wr, 0);    // skip
	csi_w32(wr, (uint32_t)l_nm);
	for (size_t i=0; i < names.size(); i++)
		wr.Write(names[i].c_str(), names[i].size()+1);
	// the index of each reference sequence
	csi_w32(wr, (uint32_t)refs.size());
	for (size_t i=0; i < refs.size(); i++)
	{
		TRefIndex &r = refs[i];
		csi_w32(wr, (uint32_t)r.bins.size());
		for (map< uint32_t, vector<TChunk> >::iterator p = r.bins.begin();
			p != r.bins.end(); p++)
		{
			csi_w32(wr, p->first);
			// 'loffset': the virtual offset of the first record in the bin
			csi_w64(wr, (p->first != CSI_META_BIN) && !p->second.empty() ?
				p->second[0].u : 0);
			csi_w32(wr, (uint32_t)p->second.size());
			for (size_t k=0; k < p->second.size(); k++)
			{
				csi_w64(wr, p->second[k].u);
				csi_w64(wr, p->second[k].v);
			}
		}
	}
	csi_w64(wr, 0);  // n_no_coor
	wr.Close();

	return num_var;
}

}


// ===========================================================
// The functions called by R
// ===========================================================

extern "C"
{
using namespace SeqArray;
using std::string;

// note: 'xprivate' is 'private' of Rconnection, renamed in Index.h for C++

static void bgzf_con_close(Rconnection con)
{
	CBgzfWriter *p = (CBgzfWriter *)con->xprivate;
	con->xprivate = NULL;
	con->isopen = FALSE;
	if (p)
	{
		// no C++ exception should be thrown when R closes a connection
		const char *err = NULL;
		try {
			p->Close();
		} catch (std::exception &e) {
			err = e.what();
		} catch (...) {
			err = "Failed to close the BGZF file.";
		}
		delete p;
		if (err) Rf_warning("%s", err);
	}
}

static size_t bgzf_con_write(const void *ptr, size_t size, size_t nitems,
	Rconnection con)
{
	CBgzfWriter *p = (CBgzfWriter *)con->xprivate;
	if (!p) Rf_error("The BGZF file has been closed.");
	try {
		p->Write(ptr, size*nitems);
	} catch (std::exception &e) {
		Rf_error("%s", e.what());
	}
	return nitems;
}

/// create a BGZF file, and return a connection object for writing
COREARRAY_DLL_EXPORT SEXP SEQ_bgzip_create(SEXP filename)
{
	const char *fn = CHAR(STRING_ELT(filename, 0));
	Rconnection con;
	SEXP r_con = R_new_custom_connection(fn, "wb", "bgzip_file", &con);

	CBgzfWriter *p = new CBgzfWriter;
	try {
		p->Create(R_ExpandFileName(fn));
	} catch (std::exception &e) {
		delete p;
		Rf_error("%s", e.what());
	}

	con->xprivate = p;
	con->isopen = TRUE;
	con->canwrite = TRUE;
	con->canread = FALSE;
	con->text = FALSE;
	con->close = &bgzf_con_close;
	con->write = &bgzf_con_write;

	return r_con;
}

/// whether 'filename' is a BGZF file
COREARRAY_DLL_EXPORT SEXP SEQ_bgzip_is(SEXP filename)
{
	return Rf_ScalarLogical(
		BGZF_IsValid(R_ExpandFileName(CHAR(STRING_ELT(filename, 0)))));
}

/// split a BGZF file into 'num' parts, return a matrix of 3 columns
COREARRAY_DLL_EXPORT SEXP SEQ_bgzip_split(SEXP filename, SEXP num)
{
	const string fn(R_ExpandFileName(CHAR(STRING_ELT(filename, 0))));
	COREARRAY_TRY
		vector<TBgzfPart> parts;
		BGZF_Split(fn.c_str(), (size_t)Rf_asInteger(num), parts);
		const size_t n = parts.size();
		rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, 3));
		for (size_t i=0; i < n; i++)
		{
			REAL(rv_ans)[i] = (double)parts[i].open_addr;
			REAL(rv_ans)[i + n] = (double)parts[i].first_addr;
			REAL(rv_ans)[i + 2*n] = (double)parts[i].end_addr;
		}
		UNPROTECT(1);
	COREARRAY_CATCH
}

/// build a CSI index for a BGZF-compressed VCF file
COREARRAY_DLL_EXPORT SEXP SEQ_bgzip_index(SEXP filename, SEXP idxfilename)
{
	COREARRAY_TRY
		// note: R_ExpandFileName() may return a pointer to a static buffer
		const string fn(R_ExpandFileName(CHAR(STRING_ELT(filename, 0))));
		const string fnidx(R_ExpandFileName(CHAR(STRING_ELT(idxfilename, 0))));
		int64_t n = BGZF_BuildCSI_VCF(fn.c_str(), fnidx.c_str());
		rv_ans = Rf_ScalarReal((double)n);
	COREARRAY_CATCH
}

}
