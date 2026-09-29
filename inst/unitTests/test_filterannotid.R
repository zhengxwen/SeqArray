###########################################################################
#
# Unit tests for seqSetFilterAnnotID
#

library(SeqArray)
library(RUnit)


# compare with the expected selection (the variants with an ID in 'id',
#   ignoring NA, "" and ".") and the expected ret.idx (the same as
#   match(id, seqGetData(f, "annotation/id")) after the filter is set)
.check_annotid <- function(f, all_id, id, msg)
{
    sel <- all_id %in% id[!is.na(id) & id!="" & id!="."]

    v <- seqSetFilterAnnotID(f, id, verbose=FALSE)
    checkTrue(is.null(v), paste0(msg, ": ret.idx=FALSE returns NULL"))
    checkIdentical(seqGetFilter(f)$variant.sel, sel,
        paste0(msg, ": selected variants"))

    v <- seqSetFilterAnnotID(f, id, ret.idx=TRUE, verbose=FALSE)
    checkIdentical(seqGetFilter(f)$variant.sel, sel,
        paste0(msg, ": selected variants with ret.idx=TRUE"))
    checkIdentical(v, match(id, seqGetData(f, "annotation/id")),
        paste0(msg, ": ret.idx"))

    invisible()
}


test_filterannotid_basic <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    all_id <- seqGetData(f, "annotation/id")
    ids <- unique(all_id)

    .check_annotid(f, all_id, ids[1:10], "basic: in order")
    .check_annotid(f, all_id, rev(ids[1:20]), "basic: reverse order")
    set.seed(1000)
    .check_annotid(f, all_id, sample(ids, 50L), "basic: random order")
    .check_annotid(f, all_id, ids, "basic: all IDs")
    .check_annotid(f, all_id, c(ids[3], ids[3], "rs_none", ids[1], ids[10]),
        "basic: duplicated and nonexistent IDs")
    .check_annotid(f, all_id, character(0), "basic: empty query")

    # the IDs are unique in the example file
    v <- seqSetFilterAnnotID(f, c(ids[2], NA, "", ".", ids[2], ids[1]),
        ret.idx=TRUE, verbose=FALSE)
    checkIdentical(v, c(2L, NA, NA, NA, 2L, 1L), "basic: ret.idx")

    invisible()
}


test_filterannotid_missing_dup <- function()
{
    # a copy of the example file with duplicated and missing IDs
    fn <- tempfile(fileext=".gds")
    on.exit(unlink(fn))
    file.copy(seqExampleFileName("gds"), fn)
    f <- seqOpen(fn, readonly=FALSE)
    nv <- length(seqGetData(f, "variant.id"))
    id <- paste0("rs", seq_len(nv) %% 300L)
    id[c(5L, 50L, 500L)] <- ""
    id[c(7L, 70L, 700L)] <- "."
    seqAddValue(f, "annotation/id", id, replace=TRUE, verbose=FALSE)
    seqClose(f)
    f <- seqOpen(fn)
    on.exit(seqClose(f), add=TRUE, after=FALSE)

    all_id <- seqGetData(f, "annotation/id")
    checkIdentical(all_id, id, "missing_dup: IDs in the GDS file")

    # NA, "" and "." do not select the variants without an ID
    v <- seqSetFilterAnnotID(f, c(NA, "", "."), ret.idx=TRUE, verbose=FALSE)
    checkTrue(!any(seqGetFilter(f)$variant.sel),
        "missing_dup: missing IDs select no variant")
    checkIdentical(v, rep(NA_integer_, 3L),
        "missing_dup: ret.idx of missing IDs")

    # ret.idx is the index of the first variant with the ID
    v <- seqSetFilterAnnotID(f, c("rs2", "rs1", "rs2"), ret.idx=TRUE,
        verbose=FALSE)
    checkIdentical(seqGetFilter(f)$variant.sel, all_id %in% c("rs1", "rs2"),
        "missing_dup: duplicated IDs in the GDS file")
    checkIdentical(v, c(2L, 1L, 2L), "missing_dup: ret.idx")

    # random queries
    ids <- c(unique(all_id), NA, "rs_none")
    set.seed(1000)
    for (i in 1:50)
    {
        q <- sample(ids, sample(0:60, 1L), replace=TRUE)
        .check_annotid(f, all_id, q, paste("missing_dup: random query", i))
    }

    invisible()
}


test_filterannotid_ret_idx_na <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    ids <- seqGetData(f, "annotation/id")
    seqSetFilterAnnotID(f, ids[1:3], verbose=FALSE)

    # ret.idx=NA is treated as ret.idx=FALSE
    v <- seqSetFilterAnnotID(f, ids[4:6], ret.idx=NA, verbose=FALSE)
    checkTrue(is.null(v), "ret_idx_na: returns NULL")
    checkIdentical(seqGetFilter(f)$variant.sel, ids %in% ids[4:6],
        "ret_idx_na: selected variants")

    invisible()
}
