#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Reading and writing the compressed files
#



#######################################################################
# BGZF: the blocked gzip format (i.e., created by bgzip)
#

# whether 'fn' is a BGZF file
.bgzf_is <- function(fn) .Call(SEQ_bgzip_is, fn)

# create a BGZF file, and return a connection object for writing
.bgzf_create <- function(fn) .Call(SEQ_bgzip_create, fn)

# build a .csi index for a BGZF-compressed VCF file, return the file name of
# the index; 'idxfn' is <fn>.csi if NULL
.bgzf_index <- function(fn, idxfn=NULL)
{
    stopifnot(is.character(fn), length(fn)==1L)
    if (is.null(idxfn)) idxfn <- paste0(fn, ".csi")
    if (!.bgzf_is(fn))
        stop("'", fn, "' should be a BGZF file (i.e., created by bgzip).")
    .Call(SEQ_bgzip_index, fn, idxfn)
    invisible(idxfn)
}
