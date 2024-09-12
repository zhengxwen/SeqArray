#######################################################################
#
# Package Name: SeqArray
#
# Description: Add data to an existing SeqArray GDS file
#


# get matched GDS type for indexing
.maxlen_bit_type <- function(nmax)
{
    if (nmax < 2L)
        "bit1"
    else if (nmax < 4L)
        "bit2"
    else if (nmax < 16L)
        "bit4"
    else if (nmax < 256L)
        "uint8"
    else if (nmax < 65536L)
        "uint16"
    else
        "int32"
}

# modify sample.id
.r_sample_id <- function(gdsfile, val, replace, nsamp, desp, compress, verbose,
    verbose.attr)
{
    stopifnot(replace)
    stopifnot(is.vector(val), is.character(val) | is.numeric(val))
    if (anyDuplicated(val))
        stop("'val' should be unique.")
    if (nsamp>0L && length(val)!=nsamp)
        stop("length(val) should be ", nsamp, ", but received ", length(val))
    if (length(desp))
        warning("'desp' is not used.")
    n <- add.gdsn(gdsfile, "sample.id", val, compress=compress, closezip=TRUE,
        replace=TRUE)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify variant.id
.r_variant_id <- function(gdsfile, val, replace, nvar, desp, compress,
    verbose, verbose.attr)
{
    stopifnot(replace)
    stopifnot(is.vector(val), is.character(val) | is.numeric(val))
    if (anyDuplicated(val))
        stop("'val' should be unique.")
    # if (nvar>0L && length(val)!=nvar)
    #     stop("length(val) should be ", nvar, ", but received ", length(val))
    if (length(desp))
        warning("'desp' is not used.")
    n <- add.gdsn(gdsfile, "variant.id", val, compress=compress, closezip=TRUE,
        replace=TRUE)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify position
.r_position <- function(gdsfile, val, replace, nvar, desp, compress, verbose,
    verbose.attr)
{
    stopifnot(replace)
    stopifnot(is.vector(val), is.numeric(val))
    if (length(val) != nvar)
        stop("length(val) should be ", nvar, ", but received ", length(val))
    if (length(desp))
        warning("'desp' is not used.")
    n <- add.gdsn(gdsfile, "position", as.integer(val), compress=compress,
        closezip=TRUE, replace=TRUE)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify chromosome
.r_chrom <- function(gdsfile, val, replace, nvar, desp, compress, verbose,
    verbose.attr)
{
    stopifnot(replace)
    stopifnot(is.vector(val), is.character(val) | is.numeric(val))
    if (length(val) != nvar)
        stop("length(val) should be ", nvar, ", but received ", length(val))
    if (length(desp))
        warning("'desp' is not used.")
    n <- add.gdsn(gdsfile, "chromosome", val, compress=compress, closezip=TRUE,
        replace=TRUE)
    .DigestCode(n, TRUE, FALSE)
    .optim_chrom(gdsfile)  # RLE-coded chromosome
    .Call(SEQ_ResetChrom, gdsfile)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify allele
.r_allele <- function(gdsfile, val, replace, nvar, desp, compress, verbose,
    verbose.attr)
{
    stopifnot(replace)
    stopifnot(is.vector(val), is.character(val))
    if (length(val) != nvar)
        stop("length(val) should be ", nvar, ", but received ", length(val))
    if (length(desp))
        warning("'desp' is not used.")
    n <- add.gdsn(gdsfile, "allele", val, compress=compress, closezip=TRUE,
        replace=TRUE)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify sample.annotation
.r_samp_annot <- function(gdsfile, val, replace, nsamp, compress, verbose)
{
    stopifnot(is.data.frame(val), NROW(val)==nsamp)
    n <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
    if (!is.null(n))
        stopifnot(replace)
    n <- add.gdsn(gdsfile, "sample.annotation", val, compress=compress,
        closezip=TRUE, replace=TRUE)
    delete.attr.gdsn(n, "R.class")
    for (nm in ls.gdsn(n))
        .DigestCode(index.gdsn(n, nm), TRUE, FALSE)
    TRUE
}

# modify sample.annotation/VARNAME
.r_samp_annot_sub <- function(gdsfile, varnm, val, replace, nsamp, desp,
    compress, verbose, verbose.attr)
{
    stopifnot(is.vector(val), length(val)==nsamp)
    n <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
    if (is.null(n))
        n <- addfolder.gdsn(gdsfile, "sample.annotation")
    varnm <- substring(varnm, 19L)
    if (varnm == "")
        stop("Invalid input variable name in 'sample.annotation/'.")
    n <- add.gdsn(n, varnm, val, compress=compress, closezip=TRUE,
        replace=replace)
    if (length(desp)>0L && !is.na(desp[1L]))
        put.attr.gdsn(n, "Description", desp[1L])
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=verbose.attr)
    TRUE
}

# modify annotation/id
.r_annot_id <- function(gdsfile, val, replace, nvar, compress, verbose)
{
    stopifnot(is.vector(val), length(val)==nvar)
    if (!is.character(val))
    {
        warning("It is suggested to use character for 'annotation/id'.",
            immediate.=TRUE)
    }
    if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
        stopifnot(replace)
    n <- add.gdsn(index.gdsn(gdsfile, "annotation"), "id", val,
        compress=compress, closezip=TRUE, replace=replace)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=TRUE)
    TRUE
}

# modify annotation/qual
.r_annot_qual <- function(gdsfile, val, replace, nvar, compress, verbose)
{
    stopifnot(is.vector(val), length(val)==nvar)
    if (!is.numeric(val))
    {
        warning("It is suggested to use numeric for 'annotation/qual'.",
            immediate.=TRUE)
    }
    if (!is.null(index.gdsn(gdsfile, "annotation/qual", silent=TRUE)))
        stopifnot(replace)
    n <- add.gdsn(index.gdsn(gdsfile, "annotation"), "qual", val,
        compress=compress, closezip=TRUE, replace=replace)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=TRUE)
    TRUE
}

# modify annotation/filter
.r_annot_flt <- function(gdsfile, val, replace, nvar, compress, verbose)
{
    stopifnot(is.vector(val) || is.factor(val), length(val)==nvar)
    if (!is.null(index.gdsn(gdsfile, "annotation/filter", silent=TRUE)))
        stopifnot(replace)
    n <- add.gdsn(index.gdsn(gdsfile, "annotation"), "filter", val,
        compress=compress, closezip=TRUE, replace=replace)
    .DigestCode(n, TRUE, FALSE)
    if (verbose)
        print(n, attribute=TRUE)
    TRUE
}

# add annotation/info
.r_annot_info <- function(gdsfile)
{
    n <- index.gdsn(gdsfile, "annotation/info", silent=TRUE)
    if (is.null(n))
        addfolder.gdsn(index.gdsn(gdsfile, "annotation"), "info")
    TRUE
}

# add annotation/format
.r_annot_fmt <- function(gdsfile)
{
    n <- index.gdsn(gdsfile, "annotation/format", silent=TRUE)
    if (is.null(n))
        addfolder.gdsn(index.gdsn(gdsfile, "annotation"), "format")
    TRUE
}

# add annotation/info/VARNAME
.r_annot_info_sub <- function(gdsfile, varnm, val, replace, nvar, desp,
    compress, packed, packed.idx, verbose, verbose.attr)
{
    n <- index.gdsn(gdsfile, varnm, silent=TRUE)
    if (!is.null(n))
        stopifnot(replace)
    node <- index.gdsn(gdsfile, dirname(varnm))
    nm <- basename(varnm)
    if (is.data.frame(val))
    {
        # store a data.frame
        stopifnot(nrow(val) == nvar)
        if (length(desp))
            stopifnot(ncol(val) == length(desp))
        n <- addfolder.gdsn(node, nm, replace=TRUE)
        if (verbose)
            print(n, attribute=verbose.attr)
        for (i in seq_len(ncol(val)))
        {
            seqAddValue(gdsfile, paste0(varnm, "/", names(val)[i]), val[[i]],
                desp=if (length(desp)) desp[i] else desp,
                replace=replace, compress=compress, packed=packed,
                verbose=verbose, verbose.attr=verbose.attr)
        }
        n <- nidx <- NULL
    } else if (is.null(val))
    {
        # it is a folder
        n <- addfolder.gdsn(node, nm, replace=TRUE)
        if (length(desp))
            put.attr.gdsn(n, "Description", desp[1L])
        if (verbose)
            print(n, attribute=verbose.attr)
        n <- nidx <- NULL
    } else if ((is.vector(val) || is.factor(val) || is.matrix(val)) && !is.list(val))
    {
        # a vector or matrix
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
                compress=compress, closezip=TRUE, replace=TRUE,
                visible=FALSE)
            num <- "."
        } else {
            n <- add.gdsn(node, nm, val, compress=compress, closezip=TRUE,
                replace=TRUE)
            nidx <- index.gdsn(node, paste0("@", nm), silent=TRUE)
            if (!is.null(nidx)) delete.gdsn(nidx)
            nidx <- NULL
            num <- as.character(ifelse(isvec, 1L, nrow(val)))
        }
    } else if (inherits(val, "SeqVarDataList"))
    {
        # it is an SeqVarDataList object list(length, data)
        ns <- val$length
        ns[is.na(ns) | ns<0L] <- 0L
        if (length(ns) != nvar)
            stop("length(val$length) should be ", nvar)
        if (!isTRUE(sum(ns)==length(val$data)))
            stop("Invalid 'SeqVarDataList' input.")
        n <- add.gdsn(node, nm, val$data, compress=compress, closezip=TRUE,
            replace=TRUE)
        st <- if (packed.idx) .maxlen_bit_type(max(ns)) else "int"
        nidx <- add.gdsn(node, paste0("@", nm), ns, storage=st,
            compress=compress, closezip=TRUE, replace=TRUE, visible=FALSE)
    } else if (is.list(val))
    {
        # store lists
        stopifnot(length(val) == nvar)
        val <- lapply(val, function(x) unlist(x, use.names=FALSE))
        ns <- lengths(val)
        val <- unlist(val, use.names=FALSE)
        if (is.null(val)) val <- logical()
        n <- add.gdsn(node, nm, val, compress=compress, closezip=TRUE,
            replace=TRUE)
        st <- if (packed.idx) .maxlen_bit_type(max(ns)) else "int"
        nidx <- add.gdsn(node, paste0("@", nm), ns, storage=st,
            compress=compress, closezip=TRUE, replace=TRUE, visible=FALSE)
    } else
        stop("Invalid type of 'val': ", typeof(val))

    if (!is.null(n))
    {
        put.attr.gdsn(n, "Number", ".")
        put.attr.gdsn(n, "Type", .vcf_type(n))
        if (!length(desp)) desp <- ""
        put.attr.gdsn(n, "Description", desp[1L])
        .DigestCode(n, TRUE, FALSE)
        if (verbose)
            print(n, attribute=verbose.attr)
    }
    if (!is.null(nidx))
        .DigestCode(nidx, TRUE, FALSE)

    TRUE
}



#######################################################################
# Add or modify values in a GDS file
#
seqAddValue <- function(gdsfile, varnm, val, desp=character(), replace=FALSE,
    compress="LZMA_RA", packed=TRUE, packed.idx=TRUE, verbose=TRUE,
    verbose.attr=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(varnm), length(varnm)==1L)
    stopifnot(is.character(desp))
    stopifnot(is.logical(replace), length(replace)==1L)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(packed), length(packed)==1L)
    stopifnot(is.logical(packed.idx), length(packed.idx)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.attr), length(verbose.attr)==1L)

    # open the file if needed
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
    nsamp <- objdesp.gdsn(index.gdsn(gdsfile, "sample.id"))$dim
    nvar  <- objdesp.gdsn(index.gdsn(gdsfile, "variant.id"))$dim

    ret <- switch(varnm,
        "sample.id"  = .r_sample_id(gdsfile, val, replace, nsamp, desp,
                compress, verbose, verbose.attr),
        "variant.id" = .r_variant_id(gdsfile, val, replace, nvar, desp,
                compress, verbose, verbose.attr),
        "position"   = .r_position(gdsfile, val, replace, nvar, desp,
                compress, verbose, verbose.attr),
        "chromosome" = .r_chrom(gdsfile, val, replace, nvar, desp,
                compress, verbose, verbose.attr),
        "allele"     = .r_allele(gdsfile, val, replace, nvar, desp,
                compress, verbose, verbose.attr),
        "sample.annotation" = .r_samp_annot(gdsfile, val, replace, nsamp,
                compress, verbose),
        "annotation/id"   = .r_annot_id(gdsfile, val, replace, nvar,
                compress, verbose),
        "annotation/qual" = .r_annot_qual(gdsfile, val, replace, nvar,
                compress, verbose),
        "annotation/filter" = .r_annot_flt(gdsfile, val, replace, nvar,
                compress, verbose),
        "annotation/info"   = .r_annot_info(gdsfile),
        "annotation/format" = .r_annot_fmt(gdsfile),
        FALSE
    )
    if (ret) return(invisible())

    # others
    if (substr(varnm, 1L, 18L) == "sample.annotation/")
    {
        .r_samp_annot_sub(gdsfile, varnm, val, replace, nsamp, desp, compress,
            verbose, verbose.attr)
    } else if (substr(varnm, 1L, 16L) == "annotation/info/")
    {
        if (nchar(varnm) <= 16L)
            stop("Invalid 'varnm'.")
        .r_annot_info_sub(gdsfile, varnm, val, replace, nvar, desp, compress,
            packed, packed.idx, verbose, verbose.attr)
    } else
        stop("Invalid `varnm`.")

    invisible()
}
