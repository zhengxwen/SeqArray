#######################################################################
#
# Package Name: SeqArray
#
# Description: Methods
#


#######################################################################
# Open a SeqArray GDS file
#
seqOpen <- function(gds.fn, readonly=TRUE, allow.duplicate=FALSE)
{
    # check
    stopifnot(is.character(gds.fn), length(gds.fn)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)

    # open the file
    ans <- openfn.gds(gds.fn, readonly=readonly, allow.fork=TRUE,
        allow.duplicate=allow.duplicate)

    # FileFormat
    at <- get.attr.gdsn(ans$root)
    if (!is.null(at$FileFormat))
    {
        # it does not throw any warning or error if FileFormat does not exist,
        # but it is encouraged to add this attribute
        if (identical(at$FileFormat, "SNP_ARRAY"))
        {
            closefn.gds(ans)
            stop(sprintf(
                "'%s' is a SNP GDS file, please use SNPRelate::snpgdsOpen().",
                gds.fn), "\n",
                "Or use SeqArray::seqSNP2GDS() for converting SNP GDS to SeqArray GDS.")
        }
        if (!identical(at$FileFormat, "SEQ_ARRAY"))
        {
            closefn.gds(ans)
            stop(sprintf("'%s' is not a SeqArray GDS file (%s).",
                gds.fn, "'FileFormat' should be 'SEQ_ARRAY'"))
        }
    }

    # FileVersion
    version <- at$FileVersion
    if (!is.null(version))
    {
        # it does not throw any warning or error if FileVersion does not exist,
        # but it is encouraged to add this attribute
        if (!identical(version, "v1.0"))
        {
            closefn.gds(ans)
            stop(sprintf(
                "Invalid FileVersion '%s' (should be v1.0).",
                as.character(version)))
        }
    }

    # check header
    n <- index.gdsn(ans, "description", silent=TRUE)
    if (is.null(n))
    {
        closefn.gds(ans)
        stop(sprintf("'%s' is not a SeqArray GDS file.", gds.fn))
    }

    .Call(SEQ_File_Init, ans)
    new("SeqVarGDSClass", ans)
}



#######################################################################
# Close a SeqArray GDS file
#
setMethod("seqClose", signature(object="gds.class"),
    function(object)
    {
        closefn.gds(object)
    }
)

setMethod("seqClose", signature(object="SeqVarGDSClass"),
    function(object)
    {
        .Call(SEQ_File_Done, object)
        closefn.gds(object)
        invisible()
    }
)



#######################################################################
# Set a working space with selected samples and variants
#
setMethod("seqSetFilter", signature(object="SeqVarGDSClass", variant.sel="ANY"),
    function(object, variant.sel, sample.sel=NULL, variant.id=NULL,
        sample.id=NULL, action=c("set", "intersect", "push", "push+set",
        "push+intersect", "pop"), ret.idx=FALSE, warn=TRUE, verbose=TRUE)
    {
        # check
        if (missing(variant.sel)) variant.sel <- NULL
        action <- match.arg(action)
        stopifnot(is.logical(ret.idx), length(ret.idx)==1L)
        stopifnot(is.logical(warn), length(warn)==1L)
        stopifnot(is.logical(verbose), length(verbose)==1L)

        # action behavior
        setflag <- FALSE
        switch(action,
            "set" = NULL,
            "intersect" = { setflag <- TRUE },
            "push" = {
                if (!all(is.null(sample.id), is.null(variant.id),
                        is.null(sample.sel), is.null(variant.sel)))
                {
                    stop("The arguments 'sample.id', 'variant.id', ",
                        "'sample.sel' and 'variant.sel' should be NULL.")
                }
                .Call(SEQ_FilterPushLast, object)
                return(invisible())
            },
            "push+set" = {
                .Call(SEQ_FilterPushLast, object)
            },
            "push+intersect" = {
                .Call(SEQ_FilterPushLast, object)
                setflag <- TRUE
            },
            "pop" = {
                if (!all(is.null(sample.id), is.null(variant.id),
                        is.null(sample.sel), is.null(variant.sel)))
                {
                    stop("The arguments 'sample.id', 'variant.id', ",
                        "'sample.sel' and 'variant.sel' should be NULL.")
                }
                .Call(SEQ_FilterPop, object)
                return(invisible())
            }
        )

        # set sample filter
        ii_samp <- NULL
        if (!is.null(sample.id))
        {
            stopifnot(is.vector(sample.id))
            stopifnot(is.numeric(sample.id) | is.character(sample.id))
            .Call(SEQ_SetSpaceSample, object, sample.id, setflag, verbose)
            if (ret.idx)
                ii_samp <- match(sample.id, seqGetData(object, "sample.id"))
        } else if (!is.null(sample.sel))
        {
            stopifnot(is.vector(sample.sel))
            stopifnot(is.logical(sample.sel) | is.raw(sample.sel) |
                is.numeric(sample.sel))
            .Call(SEQ_SetSpaceSample2, object, sample.sel, setflag, warn,
                verbose)
            if (ret.idx)
            {
                if (is.numeric(sample.sel))
                {
                    ii_samp <- order(sample.sel)
                    if (anyNA(sample.sel))
                        ii_samp[is.na(sample.sel)] <- NA_integer_
                } else {
                    ii_samp <- seq_len(SeqArray:::.seldim(object)[2L])
                }
            }
        }

        # set variant filter
        ii_var <- NULL
        if (!is.null(variant.id))
        {
            stopifnot(is.vector(variant.id))
            stopifnot(is.numeric(variant.id) | is.character(variant.id))
            .Call(SEQ_SetSpaceVariant, object, variant.id, setflag, verbose)
            if (ret.idx)
                ii_var <- match(variant.id, seqGetData(object, "variant.id"))
        } else if (!is.null(variant.sel))
        {
            stopifnot(is.vector(variant.sel))
            stopifnot(is.logical(variant.sel) | is.raw(variant.sel) |
                is.numeric(variant.sel))
            .Call(SEQ_SetSpaceVariant2, object, variant.sel, setflag, warn,
                verbose)
            if (ret.idx)
            {
                if (is.numeric(variant.sel))
                {
                    ii_var <- order(variant.sel)
                    if (anyNA(variant.sel))
                        ii_var[is.na(variant.sel)] <- NA_integer_
                } else {
                    ii_var <- seq_len(SeqArray:::.seldim(object)[3L])
                }
            }
        } else {
            if (is.null(sample.id) & is.null(sample.sel))
            {
                .Call(SEQ_SetSpaceSample, object, NULL, setflag, verbose)
                .Call(SEQ_SetSpaceVariant, object, NULL, setflag, verbose)
            }
        }

        # output
        if (ret.idx)
            list(sample_idx=ii_samp, variant_idx=ii_var)
        else
            invisible()
    }
)

setMethod("seqSetFilter", signature(object="SeqVarGDSClass",
    variant.sel="GRanges"),
    function(object, variant.sel, rm.txt="chr", intersect=FALSE, verbose=TRUE)
    {
        z <- seqnames(variant.sel)
        levels(z) <- sub(rm.txt, "", levels(z))

        seqSetFilterChrom(object,
            include = as.character(z),
            from.bp = BiocGenerics::start(variant.sel),
            to.bp   = BiocGenerics::end(variant.sel),
            intersect = intersect,
            verbose = verbose)
        invisible()
    }
)

setMethod("seqSetFilter", signature(object="SeqVarGDSClass",
    variant.sel="GRangesList"),
    function(object, variant.sel, rm.txt="chr", intersect=FALSE, verbose=TRUE)
    {
        seqSetFilter(object, unlist(variant.sel), rm.txt, intersect, verbose)
        invisible()
    }
)

setMethod("seqSetFilter", signature(object="SeqVarGDSClass",
    variant.sel="IRanges"),
    function(object, variant.sel, chr, intersect=FALSE, verbose=TRUE)
    {
        stopifnot(is.vector(chr))
        if (length(chr) > 1L)
            stopifnot(length(chr) == length(variant.sel))
        else
            chr <- rep(chr, length(variant.sel))

        seqSetFilterChrom(object,
            include = chr,
            from.bp = BiocGenerics::start(variant.sel),
            to.bp   = BiocGenerics::end(variant.sel),
            intersect = intersect,
            verbose = verbose)
        invisible()
    }
)



#######################################################################
# Reset a working space without selected samples and variants
#
seqFilterPush <- function(object)
{
    stopifnot(inherits(object, "SeqVarGDSClass"))
    seqSetFilter(object, action="push", verbose=FALSE)
}

seqFilterPop <- function(object)
{
    stopifnot(inherits(object, "SeqVarGDSClass"))
    seqSetFilter(object, action="pop", verbose=FALSE)
}


#######################################################################
# Reset a working space without selected samples and variants
#
seqResetFilter <- function(object, sample=TRUE, variant=TRUE, verbose=TRUE)
{
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.logical(sample), length(sample)==1L)
    stopifnot(is.logical(variant), length(variant)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (sample)
        .Call(SEQ_SetSpaceSample, object, NULL, FALSE, verbose)
    if (variant)
        .Call(SEQ_SetSpaceVariant, object, NULL, FALSE, verbose)

    invisible()
}



#######################################################################
# Set a filter according to specified chromosomes
#
seqSetFilterChrom <- function(object, include=NULL, is.num=NA,
    from.bp=NULL, to.bp=NULL, intersect=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.null(include) | is.numeric(include) | is.character(include))
    stopifnot(is.logical(is.num), length(is.num)==1L)
    stopifnot(is.null(from.bp) | is.numeric(from.bp))
    stopifnot(is.null(to.bp) | is.numeric(to.bp))
    stopifnot(is.logical(intersect), length(intersect)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # call C function
    .Call(SEQ_SetSpaceChrom, object, include, is.num, from.bp, to.bp, intersect,
        verbose)

    invisible()
}



#######################################################################
# Set a filter according to specified chromosomes, positions, ref & alt alleles
#
seqSetFilterPos <- function(object, chr, pos, ref=NULL, alt=NULL,
    intersect=FALSE, multi.pos=TRUE, ret.idx=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(chr)==1L || length(chr)==length(pos))
    has_ref_alt <- !is.null(ref) && !is.null(alt)
    if (!is.null(ref))
    {
        stopifnot(is.vector(ref), length(ref)==length(pos), is.character(ref))
        if (is.null(alt)) stop("'alt' is missing.")
    }
    if (!is.null(alt))
    {
        stopifnot(is.vector(alt), length(alt)==length(pos), is.character(alt))
        if (is.null(ref)) stop("'ref' is missing.")
    }
    stopifnot(is.logical(intersect), length(intersect)==1L)
    stopifnot(is.logical(multi.pos), length(multi.pos)==1L)
    stopifnot(is.logical(ret.idx), length(ret.idx)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # user-defined set
    d0 <- data.frame(
        chr=if (length(chr)==1L) rep(chr, length(pos)) else chr,
        pos=pos, stringsAsFactors=FALSE)
    if (ret.idx) d0$i0 <- seq_along(pos)

    # set filter on chromosome first
    # need reset the fitler on variants?
    seqSetFilterChrom(object, chr, intersect=intersect, verbose=FALSE)

    # gds variant set
    d1 <- data.frame(
        chr=seqGetData(object, "chromosome"),
        pos=seqGetData(object, "position"),
        i1=seqGetData(object, "$variant_index"), stringsAsFactors=FALSE)

    if (has_ref_alt)
    {
        # set filter on chromosomes and positions first
        d <- merge(d0, d1, sort=FALSE)
        seqSetFilter(object, variant.sel=d$i1, warn=FALSE, verbose=FALSE)
        d1 <- data.frame(
            chr=seqGetData(object, "chromosome"),
            pos=seqGetData(object, "position"),
            i1=seqGetData(object, "$variant_index"), stringsAsFactors=FALSE)
        # then on alleles
        d0$allele <- paste0(ref, ",", alt)
        d1$allele <- seqGetData(object, "allele")
    }

    # match
    d <- merge(d0, d1, all.x=TRUE, sort=FALSE)

    # output
    if (!isFALSE(multi.pos))
    {
        # multi.pos = TRUE
        i1 <- d$i1
        seqSetFilter(object, variant.sel=i1, warn=FALSE, verbose=verbose)
        if (ret.idx)
        {
            if (nrow(d) <= nrow(d0))
            {
                # no duplicated d$i0
                i1 <- i1[order(d$i0)]
            } else {
                # find the smallest d$i1 (the first variant in GDS)
                j <- order(d$i0, d$i1)
                i1 <- d$i1[j][!duplicated(d$i0[j])]
            }
            match(i1, seqGetData(object, "$variant_index"))
        } else {
            invisible()
        }
    } else {
        # multi.pos = FALSE
        if (ret.idx)
        {
            j <- order(d$i0, d$i1)
            j <- j[!duplicated(d$i0[j])]
            i <- d$i1[j]
            seqSetFilter(object, variant.sel=i, warn=FALSE, verbose=verbose)
            match(i, seqGetData(object, "$variant_index"))
        } else {
            i <- d$i1
            j <- order(i)
            i <- i[j[!duplicated(d$i0[j])]]
            seqSetFilter(object, variant.sel=i, warn=FALSE, verbose=verbose)
            invisible()
        }
    }
}



#######################################################################
# Set a filter according to specified conditions of MAF, MAC and missing rates
#
seqSetFilterCond <- function(gdsfile, maf=NaN, mac=1L, missing.rate=NaN,
    parallel=seqGetParallel(), .progress=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.numeric(maf), length(maf) %in% 1:2)
    stopifnot(is.numeric(mac), length(mac) %in% 1:2)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.logical(.progress), length(.progress)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    verbose <- .progress || verbose

    if (!all(is.na(maf), is.na(mac), is.na(missing.rate)))
    {
        # get MAF/MAC/missing rate
        v <- .Get_MAF_MAC_Missing(gdsfile, parallel, verbose)
        # selection
        sel <- rep(TRUE, length(v$maf))
        # check mac[1] <= ... < mac[2]
        if (!is.na(mac[1L]))
            sel <- sel & (mac[1L] <= v$mac)
        if (!is.na(mac[2L]))
            sel <- sel & (v$mac < mac[2L])
        # check maf[1] <= ... < maf[2]
        if (any(!is.na(maf)))
        {
            if (!is.na(maf[1L]))
                sel <- sel & (maf[1L] <= v$maf)
            if (!is.na(maf[2L]))
                sel <- sel & (v$maf < maf[2L])
        }
        # check ... <= missing.rate
        if (!is.na(missing.rate))
            sel <- sel & (v$miss <= missing.rate)
        # set filter
        seqSetFilter(gdsfile, variant.sel=sel, action="intersect",
            verbose=verbose)
    } else if (verbose)
        cat("No action in the filter of MAF, MAC and missing rate.\n")

    invisible()
}



#######################################################################
# Set a filter with RS ID (stored in annotation/id)
#
seqSetFilterAnnotID <- function(object, id, ret.idx=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(object, "SeqVarGDSClass"))
    stopifnot(is.character(id))
    stopifnot(is.logical(ret.idx), length(ret.idx)==1L)
    # call C function
    .Call(SEQ_SetSpaceAnnotID, object, id, verbose)
    # output
    if (ret.idx)
        match(id, seqGetData(object, "annotation/id"))
    else
        invisible()
}



#######################################################################
# To get a working space
#
seqGetFilter <- function(gdsfile, .useraw=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    # call C function
    .Call(SEQ_GetSpace, gdsfile, .useraw)
}



#######################################################################
# Get data from a working space with selected samples and variants
#
seqGetData <- function(gdsfile, var.name, .useraw=FALSE, .padNA=TRUE,
    .tolist=FALSE, .envir=NULL)
{
    # check
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }
    .Call(SEQ_GetData, gdsfile, var.name, .useraw, .padNA, .tolist, .envir)
}

print.SeqVarDataList <- function(x, ...) str(x)

seqNewVarData <- function(len, data)
{
    # check
    stopifnot(is.vector(len), is.numeric(len) | is.raw(len))
    stopifnot(is.vector(data) | is.matrix(data) | is.factor(data))
    len <- as.integer(len)
    if (any(len < 0L) || anyNA(len))
        stop("'len' should be a non-negative vector.'")
    if (!is.matrix(data))
    {
        if (sum(len) != length(data))
            stop("Invalid length of data.")
    } else {
        if (sum(len) != ncol(data))
            stop("Invalid column number of data.")
    }
    # output
    rv <- list(length=len, data=data)
    class(rv) <- "SeqVarDataList"
    rv
}

seqListVarData <- function(obj)
{
    # check
    stopifnot(inherits(obj, "SeqVarDataList"))
    # call
    .Call(SEQ_ListVarData, obj$length, obj$data)
}



#######################################################################
# Apply functions over margins on a working space with selected samples and variants
#
seqApply <- function(gdsfile, var.name, FUN,
    margin=c("by.variant", "by.sample"),
    as.is=c("none", "list", "integer", "double", "character", "logical", "raw"),
    var.index=c("none", "relative", "absolute"), parallel=FALSE,
    .useraw=FALSE, .progress=FALSE, .list_dup=TRUE, ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name)>0L)

    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    var.index <- match.arg(var.index)
    njobs <- .NumParallel(parallel)
    param <- list(useraw=.useraw, progress=.progress, list_dup=.list_dup)

    if (inherits(as.is, "connection") | inherits(as.is, "gdsn.class"))
    {
        if (njobs > 1L)
        {
            stop("the parallel function is not available ",
                "when 'as.is' is 'connection' or 'gdsn.class'.")
        }
    } else {
        as.is <- match.arg(as.is)
    }

    dm <- .seldim(gdsfile)
    if (margin == "by.variant")
    {
        if ((njobs <= 1L) || (dm[3L] <= 0L))
        {
            # C call, by.variant
            rv <- .Call(SEQ_Apply_Variant, gdsfile, var.name, FUN, as.is,
                var.index, param, new.env())
        } else {
            rv <- seqParallel(parallel, gdsfile,
                FUN=function(gdsfile, .vn, .FUN, .as.is, .varidx, .param, ...)
                {
                    .Call(SEQ_Apply_Variant, gdsfile, .vn, .FUN, .as.is,
                        .varidx, .param, new.env())
                }, split=margin, .vn=var.name, .FUN=FUN, .as.is=as.is,
                .varidx=var.index, .param=param, ...)
        }
    } else {
        if ((njobs <= 1L) || (dm[2L] <= 0L))
        {
            # C call, by.sample
            rv <- .Call(SEQ_Apply_Sample, gdsfile, var.name, FUN, as.is,
                var.index, .useraw, new.env())
        } else {
            rv <- seqParallel(parallel, gdsfile,
                FUN=function(gdsfile, .vn, .FUN, .as.is, .varidx, .param, ...)
                {
                    .Call(SEQ_Apply_Sample, gdsfile, .vn, .FUN, .as.is,
                        .varidx, .param, new.env())
                }, split=margin, .vn=var.name, .FUN=FUN, .as.is=as.is,
                .varidx=var.index, .param=param, ...)
        }
    }

    if (!is.character(as.is) | identical(as.is, "none"))
        return(invisible())
    rv
}



#######################################################################
# Apply functions over margins with chunks
#
seqBlockApply <- function(gdsfile, var.name, FUN, margin=c("by.variant"),
    as.is=c("none", "list", "unlist"),
    var.index=c("none", "relative", "absolute"), bsize=1024L, parallel=FALSE,
    .useraw=FALSE, .padNA=TRUE, .tolist=FALSE, .progress=FALSE, ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name)>0L)
    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    var.index <- match.arg(var.index)
    stopifnot(is.numeric(bsize), length(bsize)==1L)
    njobs <- .NumParallel(parallel)
    param <- list(bsize=bsize, useraw=.useraw, padNA=.padNA, tolist=.tolist,
        progress=.progress)

    if (!inherits(as.is, "connection") & !inherits(as.is, "gdsn.class"))
    {
        as.is <- match.arg(as.is)
    }

    dm <- .seldim(gdsfile)
    if (margin == "by.variant")
    {
        if ((njobs <= 1L) || (dm[3L] <= 0L))
        {
            # C call, by.variant
            rv <- .Call(SEQ_BApply_Variant, gdsfile, var.name, FUN, as.is,
                var.index, param, new.env())
        } else {
            rv <- seqParallel(parallel, gdsfile,
                FUN=function(gdsfile, .vn, .FUN, .as.is, .varidx, .param, ...)
                {
                    .Call(SEQ_BApply_Variant, gdsfile, .vn, .FUN, .as.is,
                        .varidx, .param, new.env())
                }, split=margin, .vn=var.name, .FUN=FUN, .as.is=as.is,
                .varidx=var.index, .param=param, ...)
        }
    }

    if (!is.character(as.is) | identical(as.is, "none"))
        return(invisible())
    else if (identical(as.is, "unlist"))
        rv <- unlist(rv, recursive = FALSE)
    rv
}



#######################################################################
# Data Analysis
#######################################################################

#######################################################################
# The number of alleles per site
#
seqNumAllele <- function(gdsfile)
{
    seqGetData(gdsfile, "$num_allele")
}



#######################################################################
# Missing rate
#
seqMissing <- function(gdsfile, per.variant=TRUE, parallel=seqGetParallel(),
    verbose=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(per.variant), length(per.variant)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check genotypes
    nm <- "genotype"
    gv <- .has_geno(gdsfile)
    if (!gv) nm <- .has_dosage(gdsfile)

    # calculate
    if (is.na(per.variant))
    {
        dm <- .seldim(gdsfile)
        if (gv)
        {
            sv <- seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, num, pg)
                {
                    ssum <- integer(num)
                    v <- seqApply(f, "genotype", as.is="double",
                        FUN=.cfunction2("FC_Missing_SampVariant"), y=ssum,
                        .progress=pg & (process_index==1L))
                    list(v, ssum)
                }, .combine=function(v1, v2) {
                    list(c(v1[[1L]], v2[[1L]]), v1[[2L]]+v2[[2L]])
                }, num=dm[2L], pg=verbose)
            sv[[2L]] <- sv[[2L]] / (dm[1L] * dm[3L])
        } else {
            sv <- seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, num, pg, nm)
                {
                    ssum <- integer(num)
                    tmpsum <- integer(num)
                    v <- seqApply(f, nm, as.is="double",
                        FUN=.cfunction3("FC_Missing_DS_SampVariant"),
                        y=ssum, z=tmpsum, .progress=pg & (process_index==1L))
                    list(v, ssum)
                }, .combine=function(v1, v2) {
                    list(c(v1[[1L]], v2[[1L]]), v1[[2L]]+v2[[2L]])
                }, num=dm[2L], pg=verbose, nm=nm)
            sv[[2L]] <- sv[[2L]] / dm[3L]
        }
        names(sv) <- c("variant", "sample")
        sv

    } else if (per.variant)
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg, nm)
            {
               seqApply(f, nm, as.is="double",
                   FUN=.cfunction("FC_Missing_PerVariant"), .useraw=NA,
                   .progress=pg & (process_index==1L))
            }, pg=verbose, nm=nm)

    } else {
        dm <- .seldim(gdsfile)
        if (gv)
        {
            sv <- seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, num, pg)
                {
                    ssum <- integer(num)
                    seqApply(f, "genotype", as.is="none",
                        FUN=.cfunction2("FC_Missing_PerSamp"), y=ssum,
                        .useraw=NA, .progress=pg & (process_index==1L))
                    ssum
                }, .combine="+", num=dm[2L], pg=verbose)
            sv / (dm[1L] * dm[3L])
        } else {
            sv <- seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, num, pg, nm)
                {
                    ssum <- integer(num)
                    tmpsum <- integer(num)
                    seqApply(f, nm, as.is="none",
                        FUN=.cfunction3("FC_Missing_DS_PerSamp"),
                        y=ssum, z=tmpsum, .useraw=NA,
                        .progress=pg & (process_index==1L))
                    ssum
                }, .combine="+", num=dm[2L], pg=verbose, nm=nm)
            sv / dm[3L]
        }
    }
}



#######################################################################
# Allele frequency
#
seqAlleleFreq <- function(gdsfile, ref.allele=0L, minor=FALSE,
    parallel=seqGetParallel(), verbose=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | is.numeric(ref.allele) |
        is.character(ref.allele))
    stopifnot(is.logical(minor), length(minor)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check genotypes
    gv <- .has_geno(gdsfile)
    if (gv)
    {
        nm <- "genotype"
        ploidy <- .dim(gdsfile)[1L]
    } else {
        nm <- .has_dosage(gdsfile)
        ploidy <- getOption("seqarray.ploidy", 2L)
        err <- "No 'genotype/data', try seqAlleleFreq(, ref.allele=0)"
    }
    if (is.na(ploidy)) ploidy <- 2L

    # calculate
    if (is.null(ref.allele))
    {
        if (!gv) stop(err)
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
                seqApply(f, c("genotype", "$num_allele"), as.is="list",
                    FUN = .cfunction("FC_AF_List"), .list_dup=FALSE,
                    .useraw=NA, .progress=pg & process_index==1L)
            }, pg=verbose)
    } else if (is.numeric(ref.allele))
    {
        if (length(ref.allele) == 1L)
        {
            if (ref.allele == 0L)
            {
                seqParallel(parallel, gdsfile, split="by.variant",
                    FUN = function(f, pg, nm, mi, pl, cn)
                    {
                        .cfunction3("FC_AF_SetIndex")(0L, mi, pl)
                        seqApply(f, nm, as.is="double", FUN=.cfunction(cn),
                            .useraw=NA, .progress=pg & process_index==1L)
                    }, pg=verbose, nm=nm, mi=minor, pl=ploidy,
                        cn=ifelse(gv, "FC_AF_Ref", "FC_AF_DS_Ref"))
            } else {
                seqParallel(parallel, gdsfile, split="by.variant",
                    FUN = function(f, ref, pg, nm, mi, pl, cn)
                    {
                        .cfunction3("FC_AF_SetIndex")(ref, mi, pl)
                        seqApply(f, c(nm, "$num_allele"), as.is="double",
                            FUN=.cfunction(cn), .useraw=NA,
                            .progress=pg & process_index==1L)
                    }, ref=ref.allele, pg=verbose, nm=nm, mi=minor, pl=ploidy,
                        cn=ifelse(gv, "FC_AF_Index", "FC_AF_DS_Index"))
            }
        } else {
            dm <- .seldim(gdsfile)
            if (length(ref.allele) != dm[3L])
                stop("'length(ref.allele)' should be 1 or the number of selected variants.")

            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(f, selflag, ref, pg, nm, mi, pl, cn)
                {
                    s <- ref[selflag]
                    .cfunction3("FC_AF_SetIndex")(s, mi, pl)
                    seqApply(f, c(nm, "$num_allele"), as.is="double",
                        FUN=.cfunction(cn), .useraw=NA,
                        .progress=pg & process_index==1L)
                }, ref=ref.allele, pg=verbose, nm=nm, mi=minor, pl=ploidy,
                    cn=ifelse(gv, "FC_AF_Index", "FC_AF_DS_Index"))
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (length(ref.allele) != dm[3L])
            stop("'length(ref.allele)' should be the number of selected variants.")

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(f, selflag, ref, pg, nm, mi, pl, cn)
            {
                s <- ref[selflag]
                .cfunction3("FC_AF_SetAllele")(s, mi, pl)
                seqApply(f, c(nm, "allele"), as.is="double", FUN=.cfunction(cn),
                    .useraw=NA, .progress=pg & process_index==1L)
            }, ref=ref.allele, pg=verbose, nm=nm, mi=minor, pl=ploidy,
                cn=ifelse(gv, "FC_AF_Allele", "FC_AF_DS_Allele"))
    } else
        stop("Invalid 'ref.allele'.")
}



#######################################################################
# Allele counts
#
seqAlleleCount <- function(gdsfile, ref.allele=0L, minor=FALSE,
    parallel=seqGetParallel(), verbose=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(minor), length(minor)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check genotypes
    gv <- .has_geno(gdsfile)
    if (gv)
    {
        nm <- "genotype"
        ploidy <- .dim(gdsfile)[1L]
        tp <- "integer"
    } else {
        nm <- .has_dosage(gdsfile)
        ploidy <- getOption("seqarray.ploidy", 2L)
        tp <- "double"
        err <- "No 'genotype/data', try seqAlleleFreq(, ref.allele=0)"
    }
    if (is.na(ploidy)) ploidy <- 2L

    # calculate
    if (is.null(ref.allele))
    {
        if (!gv) stop(err)
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
                seqApply(f, c("genotype", "$num_allele"), margin="by.variant",
                    as.is="list", FUN = .cfunction("FC_AlleleCount"),
                    .useraw=NA, .list_dup=FALSE,
                    .progress=pg & process_index==1L)
            }, pg=verbose)
    } else if (is.numeric(ref.allele))
    {
        if (length(ref.allele) == 1L)
        {
            if (ref.allele == 0L)
            {
                seqParallel(parallel, gdsfile, split="by.variant",
                    FUN = function(f, pg, tp, nm, mi, pl, cn)
                    {
                        .cfunction3("FC_AF_SetIndex")(0L, mi, pl)
                        seqApply(f, nm, as.is=tp, FUN=.cfunction(cn),
                            .useraw=NA, .progress=pg & process_index==1L)
                    }, pg=verbose, tp=tp, nm=nm, mi=minor, pl=ploidy,
                        cn=ifelse(gv, "FC_AC_Ref", "FC_AC_DS_Ref"))
            } else {
                seqParallel(parallel, gdsfile, split="by.variant",
                    FUN = function(f, ref, pg, tp, nm, mi, pl, cn)
                    {
                        .cfunction3("FC_AF_SetIndex")(ref, mi, pl)
                        seqApply(f, c(nm, "$num_allele"), as.is=tp,
                            FUN=.cfunction(cn), .useraw=NA,
                            .progress=pg & process_index==1L)
                    }, ref=ref.allele, pg=verbose, tp=tp, nm=nm, mi=minor, pl=ploidy,
                        cn=ifelse(gv, "FC_AC_Index", "FC_AC_DS_Index"))
            }
        } else {
            dm <- .seldim(gdsfile)
            if (length(ref.allele) != dm[3L])
                stop("'length(ref.allele)' should be 1 or the number of selected variants.")

            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(f, selflag, ref, pg, tp, nm, mi, pl, cn)
                {
                    s <- ref[selflag]
                    .cfunction3("FC_AF_SetIndex")(s, mi, pl)
                    seqApply(f, c(nm, "$num_allele"), as.is=tp,
                        FUN=.cfunction(cn), .useraw=NA,
                        .progress=pg & process_index==1L)
                }, ref=ref.allele, pg=verbose, tp=tp, nm=nm, mi=minor, pl=ploidy,
                    cn=ifelse(gv, "FC_AC_Index", "FC_AC_DS_Index"))
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (length(ref.allele) != dm[3L])
            stop("'length(ref.allele)' should be the number of selected variants.")

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(f, selflag, ref, pg, tp, nm, mi, pl, cn)
            {
                s <- ref[selflag]
                .cfunction3("FC_AF_SetAllele")(s, mi, pl)
                seqApply(f, c(nm, "allele"), as.is=tp, FUN=.cfunction(cn),
                    .useraw=NA, .progress=pg & process_index==1L)
            }, ref=ref.allele, pg=verbose, tp=tp, nm=nm, mi=minor, pl=ploidy,
                cn=ifelse(gv, "FC_AC_Allele", "FC_AC_DS_Allele"))
    } else
        stop("Invalid 'ref.allele'.")
}



#######################################################################
# Get AF/MAF, AC/MAC and missing rate for variants
#
# [deprecated]
.Get_MAF_MAC_Missing <- function(gdsfile, parallel, verbose)
{
    v <- seqGetAF_AC_Missing(gdsfile, minor=TRUE, parallel=parallel,
        verbose=verbose)
    list(maf=v$af, mac=v$ac, miss=v$miss)
}

seqGetAF_AC_Missing <- function(gdsfile, minor=FALSE, parallel=seqGetParallel(),
    verbose=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(minor), length(minor)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    .NumParallel(parallel)

    # check genotypes
    gv <- .has_geno(gdsfile)
    if (gv)
    {
        nm <- "genotype"
        ploidy <- .dim(gdsfile)[1L]
    } else {
        nm <- .has_dosage(gdsfile)
        ploidy <- getOption("seqarray.ploidy", 2L)
    }
    if (is.na(ploidy)) ploidy <- 2L

    # calculate
    m3s <- seqParallel(parallel, gdsfile, split="by.variant", .combine="list",
        FUN = function(f, pg, nm, pl, minor, cn)
        {
            m3 <- matrix(0, nrow=3L, ncol=.seldim(f)[3L])
            .cfunction3("FC_AF_AC_MISS_Init")(m3, pl, minor)
            seqApply(f, nm, as.is="none", FUN=.cfunction(cn),
                .useraw=NA, .progress=pg & (process_index==1L))
            m3
        }, pg=verbose, nm=nm, pl=ploidy, minor=minor,
            cn=ifelse(gv, "FC_AF_AC_MISS_Geno", "FC_AF_AC_MISS_DS"))

    # merge
    if (is.list(m3s)) m3s <- do.call(cbind, m3s)
    # output
    data.frame(af=m3s[1L,], ac=m3s[2L,], miss=m3s[3L,])
}



#######################################################################
# get 2-bit packed genotypes in a raw matrix
#

# deprecated
.seqGet2bGeno <- function(gdsfile, verbose=TRUE)
{
    seqGet2bGeno(gdsfile, samp_by_var=TRUE, verbose=verbose)
}

seqGet2bGeno <- function(gdsfile, samp_by_var=TRUE, verbose=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(samp_by_var), length(samp_by_var)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # get gds node
    nd <- index.gdsn(gdsfile, "genotype/data", silent=TRUE)
    if (!is.null(nd))
        varnm <- "$dosage_alt"
    else if (!is.null(index.gdsn(gdsfile, "annotation/format/DS", silent=TRUE)))
        varnm <- "annotation/format/DS"
    else
        stop("No 'genotype' or 'annotation/format/DS' is available.")

    # get num of samples and variants
    dm <- .seldim(gdsfile)
    nsamp <- dm[2L]
    nvar  <- dm[3L]
    if (isTRUE(samp_by_var))
    {
        geno <- matrix(as.raw(0xFF), nrow=ceiling(nsamp/4), ncol=nvar)
        cfunc <- .cfunction("FC_SetPackedGenoSxV")
    } else {
        geno <- matrix(as.raw(0xFF), nrow=ceiling(nvar/4), ncol=nsamp)
        cfunc <- .cfunction("FC_SetPackedGenoVxS")
    }
    if (length(geno) <= 0) return(geno)

    # initialize
    .cfunction("FC_InitPackedGeno")(geno)
    # fill
    seqApply(gdsfile, varnm, FUN=cfunc, as.is="none", .useraw=NA,
        .progress=verbose)

    # output
    geno
}
