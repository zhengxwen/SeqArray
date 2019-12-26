#######################################################################
#
# Package Name: SeqArray
#
# Description: Add data to an existing SeqArray GDS file
#


#######################################################################
# Add or modify values in a GDS file
#
seqAddValue <- function(gdsfile, varnm, val, desp=character(), replace=FALSE,
    compress="LZMA_RA", packed=TRUE, verbose=TRUE, verbose.attr=TRUE)
{
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(varnm), length(varnm)==1L)
    stopifnot(is.character(desp))
    stopifnot(is.logical(replace), length(replace)==1L)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(packed), length(packed)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.attr), length(verbose.attr)==1L)

    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1L)
        if (verbose)
            cat("Open '", gdsfile, "' ...\n", sep="")
        gdsfile <- seqOpen(gdsfile, readonly=FALSE)
        on.exit(seqClose(gdsfile))
    }

    # dm[1] -- ploidy, dm[2] -- # of total samples, dm[3] -- # of total variants
    dm <- .dim(gdsfile)
    nsamp <- dm[2L]; nvar  <- dm[3L]

    if (varnm == "sample.id")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nsamp)
        if (anyDuplicated(val)) stop("'val' should be unique.")
        if (length(desp)) warning("'desp' is not used.")
        n <- add.gdsn(gdsfile, "sample.id", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        if (verbose) print(n, attribute=verbose.attr)

    } else if (varnm == "variant.id")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nvar)
        if (anyDuplicated(val)) stop("'val' should be unique.")
        if (length(desp)) warning("'desp' is not used.")
        n <- add.gdsn(gdsfile, "variant.id", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        if (verbose) print(n, attribute=verbose.attr)

    } else if (varnm == "position")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.numeric(val), length(val)==nvar)
        if (length(desp)) warning("'desp' is not used.")
        n <- add.gdsn(gdsfile, "position", as.integer(val), compress=compress,
            closezip=TRUE, replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        if (verbose) print(n, attribute=verbose.attr)

    } else if (varnm == "chromosome")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val) | is.numeric(val), length(val)==nvar)
        if (length(desp)) warning("'desp' is not used.")
        n <- add.gdsn(gdsfile, "chromosome", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        .optim_chrom(gdsfile)  # RLE-coded chromosome
        .Call(SEQ_ResetChrom, gdsfile)
        if (verbose) print(n, attribute=verbose.attr)

    } else if (varnm == "allele")
    {
        stopifnot(replace)
        stopifnot(is.vector(val), is.character(val), length(val)==nvar)
        if (length(desp)) warning("'desp' is not used.")
        n <- add.gdsn(gdsfile, "allele", val, compress=compress, closezip=TRUE,
            replace=TRUE)
        .DigestCode(n, TRUE, FALSE)
        if (verbose) print(n, attribute=verbose.attr)

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
        if (verbose) print(n, attribute=verbose.attr)

    } else if (substr(varnm, 1L, 16L) == "annotation/info/")
    {
        if (nchar(varnm) <= 16L) stop("Invalid 'varnm'.")
        n <- index.gdsn(gdsfile, varnm, silent=TRUE)
        if (!is.null(n)) stopifnot(replace)
        node <- index.gdsn(gdsfile, dirname(varnm))
        nm <- basename(varnm)
        if (is.data.frame(val))
        {
            stopifnot(nrow(val) == nvar)
            if (length(desp)) stopifnot(ncol(val) == length(desp))
            n <- addfolder.gdsn(node, nm, replace=TRUE)
            if (verbose) print(n, attribute=verbose.attr)
            for (i in seq_len(ncol(val)))
            {
                seqAddValue(gdsfile, paste0(varnm, "/", names(val)[i]), val[[i]],
                    desp=if (length(desp)) desp[i] else desp,
                    replace=replace, compress=compress, packed=packed,
                    verbose=verbose, verbose.attr=verbose.attr)
            }
        } else if (is.null(val))
        {
            n <- addfolder.gdsn(node, nm, replace=TRUE)
            if (length(desp))
                put.attr.gdsn(n, "Description", desp[1L])
            if (verbose) print(n, attribute=verbose.attr)
        } else if ((is.vector(val) || is.factor(val) || is.matrix(val)) && !is.list(val))
        {
            isvec <- is.vector(val) || is.factor(val)
            if (isvec)
                stopifnot(length(val) == nvar)
            else if (is.matrix(val))
                stopifnot(ncol(val) == nvar)
            if (packed && isvec && any(is.na(val)))
            {
                x <- !is.na(val)
                n <- add.gdsn(node, nm, val[x], compress=compress, closezip=TRUE,
                    replace=TRUE)
                nidx <- add.gdsn(node, paste0("@", nm), x, storage="bit1",
                    compress=compress, closezip=TRUE, replace=TRUE, visible=FALSE)
                num <- "."
            } else {
                n <- add.gdsn(node, nm, val, compress=compress, closezip=TRUE,
                    replace=TRUE)
                nidx <- NULL
                num <- as.character(ifelse(isvec, 1L, nrow(val)))
            }
            put.attr.gdsn(n, "Number", num)
            put.attr.gdsn(n, "Type", .vcf_type(n))
            if (!length(desp)) desp <- ""
            put.attr.gdsn(n, "Description", desp[1L])
            .DigestCode(n, TRUE, FALSE)
            if (!is.null(nidx)) .DigestCode(nidx, TRUE, FALSE)
            if (verbose) print(n, attribute=verbose.attr)
        } else {
            stop("Invalid type of 'val': ", typeof(val))
        }

    } else {
        stop("Invalid `varnm`.")
    }

    invisible()
}
