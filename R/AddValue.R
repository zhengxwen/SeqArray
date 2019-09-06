#######################################################################
#
# Package Name: SeqArray
#
# Description: Add data to an existing SeqArray GDS file
#


#######################################################################
# Add or modify values in a GDS file
#
seqAddValue <- function(gdsfile, varnm, val, desp="", replace=FALSE, compress="LZMA_RA",
    packed=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(varnm), length(varnm)==1L)
    stopifnot(is.character(desp), length(desp)==1L, !is.na(desp))
    stopifnot(is.logical(replace), length(replace)==1L)
    stopifnot(is.character(compress))
    stopifnot(is.logical(packed), length(packed)==1L)

    # dm[1] -- ploidy, dm[2] -- # of total samples, dm[3] -- # of total variants
    dm <- .dim(gdsfile)
    nsamp <- dm[2L]
    nvar  <- dm[3L]

    if (varnm == "sample.id")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nsamp)
        if (anyDuplicated(val))
            stop("'val' should be unique.")
        n <- add.gdsn(gdsfile, "sample.id", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)

    } else if (varnm == "variant.id")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nvar)
        if (anyDuplicated(val))
            stop("'val' should be unique.")
        n <- add.gdsn(gdsfile, "variant.id", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)

    } else if (varnm == "position")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.numeric(val), length(val)==nvar)
        n <- add.gdsn(gdsfile, "position", as.integer(val), compress=compress,
            closezip=TRUE, replace=TRUE)
        .DigestCode(n, TRUE, FALSE)

    } else if (varnm == "chromosome")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nvar)
        n <- add.gdsn(gdsfile, "chromosome", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        .optim_chrom(gdsfile)  # RLE-coded chromosome
        .Call(SEQ_ResetChrom, gdsfile)

    } else if (varnm == "allele")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val), length(val)==nvar)
        n <- add.gdsn(gdsfile, "allele", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)

    } else if (varnm == "sample.annotation")
    {
        stopifnot(is.data.frame(val), NROW(val)==nsamp)
        n <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
        if (!is.null(n)) stopifnot(replace)
        n <- add.gdsn(gdsfile, "sample.annotation", val, compress=compress,
            closezip=TRUE, replace=TRUE)
        delete.attr.gdsn(n, "R.class")
        for (nm in ls.gdsn(n))
            .DigestCode(index.gdsn(n, nm), TRUE, FALSE)

    } else if (substr(varnm, 1L, 18L) == "sample.annotation/")
    {
        stopifnot(is.vector(val), length(val)==nsamp)
        n <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
        if (is.null(n)) n <- addfolder.gdsn(gdsfile, "sample.annotation")
        varnm <- substring(varnm, 19L)
        n <- add.gdsn(n, varnm, val, compress=compress, closezip=TRUE, replace=replace)
        .DigestCode(n, TRUE, FALSE)

    } else {
        stop("Invalid `varnm`.")
    }

    invisible()
}
