// ===========================================================
//
// samtools_ext: the functions interacting with Rsamtools
//
// Copyright (C) 2017-2024    Xiuwen Zheng
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

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Connections.h>
#include <R_ext/Utils.h>


// ======================================================================
// ======================================================================

static const char *pkg_samtools = "Rsamtools";

#define BGZF void
#define PKG_LOAD(name)	{ \
		DL_FUNC f = R_FindSymbol(#name, pkg_samtools, NULL); \
		if (!f) \
			Rf_error("No function '%s' in the %s package", #name, pkg_samtools); \
		memcpy(&name, &f, sizeof(f)); \
	}

static BGZF* (*bgzf_open)(const char* path, const char *mode) = NULL;
static int (*bgzf_close)(BGZF *fp) = NULL;
static ssize_t (*bgzf_write)(BGZF *fp, const void *data, ssize_t length) = NULL;


static void init_bgzf(void)
{
	PKG_LOAD(bgzf_open);
	PKG_LOAD(bgzf_close);
	PKG_LOAD(bgzf_write);
}


static void bzfile_close(Rconnection con)
{
	if (con->private)
	{
		(*bgzf_close)((BGZF*)con->private);
		con->private = NULL;
	}
	con->isopen = FALSE;
}

static size_t bzfile_write(const void *ptr, size_t size, size_t nitems,
	Rconnection con)
{
	BGZF *fp = (BGZF *)con->private;
	/* uses 'unsigned' for len */
	if ((double) size * (double) nitems > UINT_MAX)
		Rf_error("too large a block specified");
	return (*bgzf_write)(fp, ptr, (unsigned int)(size*nitems)) / size;
}


/// Create a bgzip connection object
SEXP SEQ_bgzip_create(SEXP filename)
{
	init_bgzf();

	const char *fn = CHAR(STRING_ELT(filename, 0));
	Rconnection con;
	SEXP r_con = R_new_custom_connection(fn, "wb", "bgzip_file", &con);
	BGZF *bz = (*bgzf_open)(R_ExpandFileName(fn), "wb");
	if (!bz)
		Rf_error("Cannot open '%s'.", fn);

	con->private = bz;
	con->isopen = TRUE;
	con->canwrite = TRUE;
	con->canread = FALSE;
	con->text = FALSE;
	con->close = &bzfile_close;
	con->write = &bzfile_write;

	return r_con;
}
