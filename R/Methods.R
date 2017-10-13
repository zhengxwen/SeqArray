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
        if (!identical(at$FileFormat, "SEQ_ARRAY"))
        {
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
        "push+intersect", "pop"), verbose=TRUE)
    {
        if (missing(variant.sel)) variant.sel <- NULL

        # check
        action <- match.arg(action)
        stopifnot(is.logical(verbose))

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
                .Call(SEQ_FilterPushEmpty, object)
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

        if (!is.null(sample.id))
        {
            stopifnot(is.vector(sample.id))
            stopifnot(is.numeric(sample.id) | is.character(sample.id))
            .Call(SEQ_SetSpaceSample, object, sample.id, setflag, verbose)
        } else if (!is.null(sample.sel))
        {
            stopifnot(is.vector(sample.sel))
            stopifnot(is.logical(sample.sel) | is.raw(sample.sel) | is.numeric(sample.sel))
            .Call(SEQ_SetSpaceSample2, object, sample.sel, setflag, verbose)
        }

        if (!is.null(variant.id))
        {
            stopifnot(is.vector(variant.id))
            stopifnot(is.numeric(variant.id) | is.character(variant.id))
            .Call(SEQ_SetSpaceVariant, object, variant.id, setflag, verbose)
        } else if (!is.null(variant.sel))
        {
            stopifnot(is.vector(variant.sel))
            stopifnot(is.logical(variant.sel) | is.raw(variant.sel) | is.numeric(variant.sel))
            .Call(SEQ_SetSpaceVariant2, object, variant.sel, setflag, verbose)
        } else {
            if (is.null(sample.id) & is.null(sample.sel))
            {
                .Call(SEQ_SetSpaceSample, object, NULL, setflag, verbose)
                .Call(SEQ_SetSpaceVariant, object, NULL, setflag, verbose)
            }
        }

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
    .Call(SEQ_SetChrom, object, include, is.num, from.bp, to.bp, intersect,
        verbose)

    invisible()
}



#######################################################################
# Set a filter according to specified conditions
#
seqSetFilterCond <- function(gdsfile, maf=NaN, mac=1L, missing.rate=NaN,
    parallel=seqGetParallel(), .progress=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.numeric(maf), length(maf) %in% 1L:2L)
    stopifnot(is.numeric(mac), length(mac) %in% 1L:2L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.logical(.progress), length(.progress)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (!all(c(is.na(maf), is.na(mac), is.na(missing.rate))))
    {
        # calculation
        n <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
                seqApply(f, "genotype", margin="by.variant", as.is="list",
                    FUN = .cfunction("FC_AlleleCount2"),
                    .useraw=NA, .progress=pg & (process_index==1L))
            }, pg=.progress)
        N <- .seldim(gdsfile)
        N <- N[1L] * N[2L]    # the total number of alleles per site
        n0 <- sapply(n, `[`, i=1L)
        nm <- sapply(n, `[`, i=2L)
        nn <- N - nm    # the number of non-missing alleles per site
        n0 <- pmin(n0, nn-n0)
        # selection
        sel <- rep(TRUE, length(n0))
        if (!is.na(mac[1L]))
            sel <- sel & (mac[1L] <= n0)
        if (!is.na(mac[2L]))
            sel <- sel & (n0 < mac[2L])
        if (any(!is.na(maf)))
        {
            p <- n0 / nn
            if (!is.na(maf[1L]))
                sel <- sel & (maf[1L] <= p)
            if (!is.na(maf[2L]))
                sel <- sel & (p < maf[2L])
        }
        if (!is.na(missing.rate))
            sel <- sel & (nm/N <= missing.rate)
        # set filter
        seqSetFilter(gdsfile, variant.sel=sel, action="intersect",
            verbose=verbose)
    }

    invisible()
}



#######################################################################
# To set a working space with selected variants
#
# seqSetFilterVariant <- function(gdsfile, autosome.only=TRUE,
#     remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, verbose=TRUE)
# {
#    # check
#    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
#    stopifnot(is.logical(verbose))
# }



#######################################################################
# To get a working space
#
seqGetFilter <- function(gdsfile, .useraw=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    .Call(SEQ_GetSpace, gdsfile, .useraw)
}



#######################################################################
# Get data from a working space with selected samples and variants
#
seqGetData <- function(gdsfile, var.name, .useraw=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name)==1L)

    .Call(SEQ_GetData, gdsfile, var.name, .useraw)
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

    if (margin == "by.variant")
    {
        if (njobs <= 1L)
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
        if (njobs <= 1L)
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
    .useraw=FALSE, .progress=FALSE, ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name)>0L)

    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    var.index <- match.arg(var.index)
    stopifnot(is.numeric(bsize), length(bsize)==1L)
    njobs <- .NumParallel(parallel)
    param <- list(bsize=bsize, useraw=.useraw, progress=.progress)

    if (!inherits(as.is, "connection") & !inherits(as.is, "gdsn.class"))
    {
        as.is <- match.arg(as.is)
    }

    if (margin == "by.variant")
    {
        if (njobs <= 1L)
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
seqMissing <- function(gdsfile, per.variant=TRUE, .progress=FALSE,
    parallel=seqGetParallel())
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(per.variant))
    stopifnot(is.logical(.progress), length(.progress)==1L)

    if (per.variant)
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
               seqApply(f, "genotype", margin="by.variant",
                   as.is="double", FUN=.cfunction("FC_Missing_PerVariant"),
                   .useraw=NA, .progress=pg & (process_index==1L))
            }, pg=.progress)
    } else {
        dm <- .seldim(gdsfile)
        sum <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, num, pg)
            {
                tmpsum <- integer(num)
                seqApply(f, "genotype", margin="by.variant",
                    as.is="none", FUN=.cfunction2("FC_Missing_PerSample"),
                    y=tmpsum, .progress=pg & (process_index==1L))
                tmpsum
            }, .combine="+", num=dm[2L], pg=.progress)
        sum / (dm[1L] * dm[3L])
    }
}



#######################################################################
# Allele frequency
#
seqAlleleFreq <- function(gdsfile, ref.allele=0L, .progress=FALSE,
    parallel=seqGetParallel())
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | is.numeric(ref.allele) |
        is.character(ref.allele))
    stopifnot(is.logical(.progress), length(.progress)==1L)

    if (is.null(ref.allele))
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
                seqApply(f, c("genotype", "$num_allele"), as.is="list",
                    FUN = .cfunction("FC_AF_List"), .list_dup=FALSE,
                    .useraw=NA, .progress=pg & (process_index==1L))
            }, pg=.progress)
    } else if (is.numeric(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (!(length(ref.allele) %in% c(1L, dm[3L])))
            stop("'length(ref.allele)' should be 1 or the number of selected variants.")

        if (length(ref.allele) == 1L)
        {
            seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, ref, pg)
                {
                    .cfunction("FC_AF_SetIndex")(ref)
                    seqApply(f, c("genotype", "$num_allele"),
                        as.is="double", FUN = .cfunction("FC_AF_Index"),
                        .useraw=NA, .progress=pg & (process_index==1L))
                }, ref=ref.allele, pg=.progress)
        } else {
            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(f, selflag, ref, pg)
                {
                    s <- ref[selflag]
                    .cfunction("FC_AF_SetIndex")(s)
                    seqApply(f, c("genotype", "$num_allele"),
                        as.is="double", FUN = .cfunction("FC_AF_Index"),
                        .useraw=NA, .progress=pg & (process_index==1L))
                }, ref=ref.allele, pg=.progress)
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (length(ref.allele) != dm[3L])
            stop("'length(ref.allele)' should be the number of selected variants.")

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(f, selflag, ref, pg)
            {
                s <- ref[selflag]
                .cfunction("FC_AF_SetAllele")(s)
                seqApply(f, c("genotype", "allele"), margin="by.variant",
                    as.is="double", FUN = .cfunction("FC_AF_Allele"),
                    .useraw=NA, .progress=pg & (process_index==1L))
            }, ref=ref.allele, pg=.progress)
    } else {
        stop("Invalid 'ref.allele'.")
    }
}



#######################################################################
# Allele counts
#
seqAlleleCount <- function(gdsfile, ref.allele=0L, .progress=FALSE,
    parallel=seqGetParallel())
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(.progress), length(.progress)==1L)

    if (is.null(ref.allele))
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, pg)
            {
                seqApply(f, c("genotype", "$num_allele"), margin="by.variant",
                    as.is="list", FUN = .cfunction("FC_AlleleCount"),
                    .useraw=NA, .progress=pg & (process_index==1L))
            }, pg=.progress)
    } else if (is.numeric(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (!(length(ref.allele) %in% c(1L, dm[3L])))
            stop("'length(ref.allele)' should be 1 or the number of selected variants.")

        if (length(ref.allele) == 1L)
        {
            seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, ref, pg)
                {
                    .cfunction("FC_AF_SetIndex")(ref)
                    seqApply(f, c("genotype", "$num_allele"),
                        as.is="integer", FUN = .cfunction("FC_AC_Index"),
                        .useraw=NA, .progress=pg & (process_index==1L))
                }, ref=ref.allele, pg=.progress)
        } else {
            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(f, selflag, ref, pg)
                {
                    s <- ref[selflag]
                    .cfunction("FC_AF_SetIndex")(s)
                    seqApply(f, c("genotype", "$num_allele"),
                        as.is="integer", FUN = .cfunction("FC_AC_Index"),
                        .useraw=NA, .progress=pg & (process_index==1L))
                }, ref=ref.allele, pg=.progress)
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        if (length(ref.allele) != dm[3L])
            stop("'length(ref.allele)' should be the number of selected variants.")

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(f, selflag, ref, pg)
            {
                s <- ref[selflag]
                .cfunction("FC_AF_SetAllele")(s)
                seqApply(f, c("genotype", "allele"), margin="by.variant",
                    as.is="double", FUN = .cfunction("FC_AC_Allele"),
                    .useraw=NA, .progress=pg & (process_index==1L))
            }, ref=ref.allele, pg=.progress)
    } else {
        stop("Invalid 'ref.allele'.")
    }
}
