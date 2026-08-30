// ===========================================================
//
// CompressIO.h: reading and writing the compressed files
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

#ifndef _HEADER_SEQ_COMPRESS_IO_
#define _HEADER_SEQ_COMPRESS_IO_

#include <stdio.h>
#include <stdint.h>
#include <vector>


namespace SeqArray
{

// ===========================================================
// BGZF: the blocked gzip format (i.e., created by bgzip)
//
//     A BGZF file is a series of independent gzip members ('blocks'), each of
// which has an extra field 'BC' storing the size of the block, and stores at
// most 64KiB of the uncompressed data. So a block boundary can be found
// without decompressing anything, and the decompression can start from any
// block. A 'virtual offset' is a pair (the file offset of the block, the
// offset within the uncompressed block), which is used for seeking.
// ===========================================================

/// the maximum size of the uncompressed data in a BGZF block
const size_t BGZF_MAX_BLOCK = 65536;
/// the maximum size of the uncompressed data written to a BGZF block
const size_t BGZF_BLOCK_SIZE = 0xFF00;

/// whether 'fn' is a BGZF file
bool BGZF_IsValid(const char *fn);


/// Reading a BGZF file, one block at a time
class CBgzfReader
{
public:
	CBgzfReader();
	~CBgzfReader();

	/// open a BGZF file, starting from the block at the file offset 'addr'
	void Open(const char *fn, int64_t addr=0);
	/// close the file, if any
	void Close();
	/// decompress the next block to 'dst' (at least BGZF_MAX_BLOCK bytes),
	///   return the number of the uncompressed bytes (0 for the end of file)
	size_t ReadBlock(void *dst);

	/// the file offset of the current block
	inline int64_t Addr() const { return fAddr; }
	inline bool IsOpen() const { return fFile != NULL; }

private:
	FILE *fFile;      ///< the file handler
	int64_t fAddr;    ///< the file offset of the current block
	int64_t fNext;    ///< the file offset of the next block
	std::vector<unsigned char> fBuffer;  ///< the compressed data of a block

	/// read the header of the block at 'addr', return the total size of the
	///   block, or 0 for the end of file
	size_t ReadHead(int64_t addr, unsigned char head[18]);
};


/// Writing a BGZF file
class CBgzfWriter
{
public:
	CBgzfWriter();
	~CBgzfWriter();

	/// create a BGZF file; 'level' is the zlib compression level
	void Create(const char *fn, int level=6);
	/// write the remaining data, append the EOF marker and close the file
	void Close();
	/// write 'size' bytes
	void Write(const void *buf, size_t size);
	/// compress and write the buffered data, if any
	void Flush();

	inline bool IsOpen() const { return fFile != NULL; }

private:
	FILE *fFile;      ///< the file handler
	int fLevel;       ///< the zlib compression level
	std::vector<unsigned char> fBuffer;  ///< the uncompressed data
	std::vector<unsigned char> fCBuffer; ///< the compressed block
	size_t fSize;     ///< the number of bytes stored in fBuffer
};


// ===========================================================
// The CSI index (coordinate-sorted index) for a BGZF-compressed VCF file
// ===========================================================

/// build a CSI index for the BGZF-compressed VCF file 'fn', and save it to
///   'fnidx'; return the number of the indexed variants
int64_t BGZF_BuildCSI_VCF(const char *fn, const char *fnidx);


/// Split a BGZF file into 'num' parts at the block boundaries, by scanning
///   the block headers only (i.e., nothing is decompressed). For each part,
///   'open_addr' is the file offset of the block where the reading starts
///   (i.e., the block before 'first_addr', so that the first complete line
///   can be found), and the part covers the blocks in
///   [first_addr, end_addr).
struct TBgzfPart { int64_t open_addr, first_addr, end_addr; };
void BGZF_Split(const char *fn, size_t num, std::vector<TBgzfPart> &parts);

}

#endif /* _HEADER_SEQ_COMPRESS_IO_ */
