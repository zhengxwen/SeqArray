###########################################################################
#
# Unit tests for seqSetFilterPos
#

library(SeqArray)
library(RUnit)


test_filterpos_basic <- function()
{
    # open the GDS file
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    # get all positions and chromosomes
    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")
    vid <- seqGetData(f, "variant.id")

    # test 1: single chromosome, all positions
    chr1 <- chr[1L]
    idx <- which(chr == chr1)
    seqSetFilterPos(f, chr1, pos[idx], verbose=FALSE)
    v <- seqGetData(f, "variant.id")
    checkTrue(all(v %in% vid[idx]), "basic: single chr, all positions")

    # test 2: subset of positions
    sub_idx <- idx[1:5]
    seqSetFilterPos(f, chr1, pos[sub_idx], verbose=FALSE)
    v <- seqGetData(f, "variant.id")
    checkTrue(length(v) >= 1L, "basic: subset positions returns results")
    checkTrue(all(seqGetData(f, "position") %in% pos[sub_idx]),
        "basic: returned positions are in query")

    invisible()
}


test_filterpos_multi_chr <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")

    # pick positions from multiple chromosomes
    chr_lst <- unique(chr)
    if (length(chr_lst) >= 2L)
    {
        # take first 3 from each of first 2 chromosomes
        sel <- c(which(chr == chr_lst[1L])[1:3],
            which(chr == chr_lst[2L])[1:3])
        sel <- sel[!is.na(sel)]

        seqSetFilterPos(f, chr[sel], pos[sel], verbose=FALSE)
        v_pos <- seqGetData(f, "position")
        v_chr <- seqGetData(f, "chromosome")
        checkTrue(length(v_pos) >= 1L, "multi_chr: returns results")
    }

    invisible()
}


test_filterpos_ret_idx <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")

    # single chromosome
    chr1 <- unique(chr)[1L]
    idx <- which(chr == chr1)
    sub_idx <- idx[1:5]

    ri <- seqSetFilterPos(f, chr1, pos[sub_idx], ret.idx=TRUE, verbose=FALSE)
    checkEquals(length(ri), length(sub_idx),
        "ret.idx: length matches input length")
    # non-NA entries should be valid indices into current filter
    vi <- seqGetData(f, "$variant_index")
    checkTrue(all(ri[!is.na(ri)] %in% seq_along(vi)),
        "ret.idx: indices are valid")

    invisible()
}


test_filterpos_ref_alt <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")
    allele <- seqGetData(f, "allele")

    # parse ref and alt from allele strings
    sp <- strsplit(allele, ",")
    ref_all <- vapply(sp, `[`, character(1L), 1L)
    alt_all <- vapply(sp, `[`, character(1L), 2L)

    # test with correct ref/alt - should find matches
    chr1 <- unique(chr)[1L]
    idx <- which(chr == chr1)[1:5]
    idx <- idx[!is.na(idx)]

    seqSetFilterPos(f, chr1, pos[idx], ref=ref_all[idx], alt=alt_all[idx],
        verbose=FALSE)
    v_pos <- seqGetData(f, "position")
    checkTrue(length(v_pos) >= 1L, "ref_alt: correct alleles find matches")

    # test with wrong ref - should find fewer or no matches
    wrong_ref <- rep("ZZZ", length(idx))
    seqSetFilterPos(f, chr1, pos[idx], ref=wrong_ref, alt=alt_all[idx],
        verbose=FALSE)
    v_pos2 <- seqGetData(f, "position")
    checkTrue(length(v_pos2) <= length(v_pos),
        "ref_alt: wrong ref finds fewer matches")

    # test with NA ref (matches any ref)
    na_ref <- rep(NA_character_, length(idx))
    seqSetFilterPos(f, chr1, pos[idx], ref=na_ref, alt=alt_all[idx],
        verbose=FALSE)
    v_pos3 <- seqGetData(f, "position")
    checkTrue(length(v_pos3) >= length(v_pos),
        "ref_alt: NA ref matches at least as many")

    invisible()
}


test_filterpos_multi_pos_false <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")

    chr1 <- unique(chr)[1L]
    idx <- which(chr == chr1)

    # multi.pos=TRUE (default) may include duplicates
    seqSetFilterPos(f, chr1, pos[idx], multi.pos=TRUE, verbose=FALSE)
    n_true <- length(seqGetData(f, "variant.id"))

    # multi.pos=FALSE should give at most as many
    seqSetFilterPos(f, chr1, pos[idx], multi.pos=FALSE, verbose=FALSE)
    n_false <- length(seqGetData(f, "variant.id"))

    checkTrue(n_false <= n_true,
        "multi.pos=FALSE gives <= multi.pos=TRUE variants")

    invisible()
}


test_filterpos_intersect <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")
    vid_all <- seqGetData(f, "variant.id")

    chr1 <- unique(chr)[1L]
    idx <- which(chr == chr1)

    # set a partial filter first
    seqSetFilter(f, variant.sel=idx[1:min(10, length(idx))], verbose=FALSE)
    vid_before <- seqGetData(f, "variant.id")

    # intersect=TRUE should give subset of current filter
    seqSetFilterPos(f, chr1, pos[idx], intersect=TRUE, verbose=FALSE)
    vid_after <- seqGetData(f, "variant.id")
    checkTrue(all(vid_after %in% vid_before),
        "intersect: result is subset of previous filter")

    invisible()
}


test_filterpos_no_match <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")

    # query a position that doesn't exist
    chr1 <- unique(chr)[1L]
    seqSetFilterPos(f, chr1, 999999999L, verbose=FALSE)
    v <- seqGetData(f, "variant.id")
    checkEquals(length(v), 0L,
        "no_match: returns empty when no positions match")

    invisible()
}


test_filterpos_duplicated_input <- function()
{
    f <- seqOpen(seqExampleFileName("gds"))
    on.exit(seqClose(f))

    chr <- seqGetData(f, "chromosome")
    pos <- seqGetData(f, "position")

    # duplicate the same position in the query
    chr1 <- unique(chr)[1L]
    idx <- which(chr == chr1)[1:3]
    idx <- idx[!is.na(idx)]

    dup_pos <- rep(pos[idx], each=2L)
    seqSetFilterPos(f, chr1, dup_pos, verbose=FALSE)
    v_pos <- seqGetData(f, "position")
    checkTrue(length(v_pos) >= 1L,
        "duplicated_input: handles duplicate query positions")

    # with ret.idx
    ri <- seqSetFilterPos(f, chr1, dup_pos, ret.idx=TRUE, verbose=FALSE)
    checkEquals(length(ri), length(dup_pos),
        "duplicated_input: ret.idx length matches duplicated input")

    invisible()
}
