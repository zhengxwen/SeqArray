#######################################################################
#
# Package Name: SeqArray
#
# The interface to the class "SeqVarGDSClass"
#


#######################################################################
# Open a sequencing-variant GDS file
#
seqOpen <- function(gds.fn, readonly=TRUE)
{
    # check
    stopifnot(is.character(gds.fn) & (length(gds.fn)==1))
    stopifnot(is.logical(readonly))

    # open the file
    ans <- openfn.gds(gds.fn, readonly=readonly, allow.fork=TRUE)

    # FileFormat
    at <- get.attr.gdsn(ans$root)
    if ("FileFormat" %in% names(at))
    {
        # it does not throw any warning or error if FileFormat does not exist,
        # but it is encouraged to add this attribute
        if (!identical(at$FileFormat, "SEQ_ARRAY"))
        {
            stop(sprintf("'%s' is not a sequencing-variant GDS file (%s).",
                gds.fn, "'FileFormat' should be 'SEQ_ARRAY'"))
        }
    }

    # check header
    n <- index.gdsn(ans, "description", silent=TRUE)
    if (is.null(n))
    {
        closefn.gds(ans)
        stop(sprintf("'%s' is not a sequencing-variant GDS file.", gds.fn))
    }

    at <- get.attr.gdsn(n)
    if (!("sequence.variant.format" %in% names(at)))
    {
        closefn.gds(ans)
        stop(sprintf("'%s' is not a sequencing-variant GDS file.", gds.fn))
    }

    .Call(SEQ_File_Init, ans)

    new("SeqVarGDSClass", ans)
}



#######################################################################
# Close a sequencing-variant GDS file
#
setMethod("seqClose", "SeqVarGDSClass", function(object)
    {
        .Call(SEQ_File_Done, object)
        closefn.gds(object)
        invisible()
    }
)



#######################################################################
# To set a working space with selected samples and variants
#
seqSetFilter <- function(gdsfile, sample.id=NULL, variant.id=NULL,
    samp.sel=NULL, variant.sel=NULL,
    action=c("set", "intersect", "push", "push+set", "push+intersect", "pop"),
    verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(verbose))

    action <- match.arg(action)
    intersect.flag <- FALSE
    switch(action,
        "set" = NULL,
        "intersect" = { intersect.flag <- TRUE },
        "push" = {
            if (!all(is.null(sample.id), is.null(variant.id),
                    is.null(samp.sel), is.null(variant.sel)))
            {
                stop("The arguments 'sample.id', 'variant.id', ",
                    "'samp.sel' and 'variant.sel' should be NULL.")
            }
            .Call(SEQ_FilterPushLast, gdsfile)
        },
        "push+set" = {
            .Call(SEQ_FilterPushEmpty, gdsfile)
        },
        "push+intersect" = {
            .Call(SEQ_FilterPushEmpty, gdsfile)
            intersect.flag <- TRUE
        },
        "pop" = {
            if (!all(is.null(sample.id), is.null(variant.id),
                    is.null(samp.sel), is.null(variant.sel)))
            {
                stop("The arguments 'sample.id', 'variant.id', ",
                    "'samp.sel' and 'variant.sel' should be NULL.")
            }
            .Call(SEQ_FilterPop, gdsfile)
            return(invisible())
        }
    )

    if (!is.null(sample.id))
    {
        stopifnot(is.vector(sample.id))
        stopifnot(is.numeric(sample.id) | is.character(sample.id))
        .Call(SEQ_SetSpaceSample, gdsfile, sample.id, intersect.flag, verbose)
    } else if (!is.null(samp.sel))
    {
        stopifnot(is.vector(samp.sel) & is.logical(samp.sel))
        .Call(SEQ_SetSpaceSample, gdsfile, samp.sel, intersect.flag, verbose)
    }

    if (!is.null(variant.id))
    {
        stopifnot(is.vector(variant.id))
        stopifnot(is.numeric(variant.id) | is.character(variant.id))
        .Call(SEQ_SetSpaceVariant, gdsfile, variant.id, intersect.flag, verbose)
    } else if (!is.null(variant.sel))
    {
        stopifnot(is.vector(variant.sel) & is.logical(variant.sel))
        .Call(SEQ_SetSpaceVariant, gdsfile, variant.sel, intersect.flag, verbose)
    } else {
        if (is.null(sample.id) & is.null(samp.sel))
        {
            .Call(SEQ_SetSpaceSample, gdsfile, NULL, intersect.flag, verbose)
            .Call(SEQ_SetSpaceVariant, gdsfile, NULL, intersect.flag, verbose)
        }
    }

    invisible()
}



#######################################################################
# Set a filter according to specified chromosomes
#
seqSetFilterChrom <- function(gdsfile, include=NULL, is.num=NA,
    from.bp=NaN, to.bp=NaN)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(include) | is.numeric(include) | is.character(include))
    stopifnot(is.logical(is.num))

    stopifnot(is.numeric(from.bp) & is.vector(from.bp))
    stopifnot(length(from.bp) == 1L)
    stopifnot(is.numeric(to.bp) & is.vector(to.bp))
    stopifnot(length(to.bp) == 1L)

    # call C function
    .Call(SEQ_SetChrom, gdsfile, include, is.num)

    if (is.finite(from.bp) | is.finite(to.bp))
    {
        pos <- seqGetData(gdsfile, "position")
        if (is.finite(from.bp))
        {
            flag <- (pos >= from.bp)
            if (is.finite(to.bp))
                flag <- flag & (pos <= to.bp)
        } else {
            flag <- (pos <= to.bp)
        }
        seqSetFilter(gdsfile, action="intersect", variant.sel=flag, verbose=FALSE)
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
seqGetFilter <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    .Call(SEQ_GetSpace, gdsfile)
}



#######################################################################
# Get data from a working space with selected samples and variants
#
seqGetData <- function(gdsfile, var.name)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name)==1))

    .Call(SEQ_GetData, gdsfile, var.name)
}



#######################################################################
# Apply functions over margins on a working space with selected samples and variants
#
seqApply <- function(gdsfile, var.name, FUN,
    margin = c("by.variant", "by.sample"), as.is = c("none", "list",
    "integer", "double", "character", "logical", "raw"),
    var.index = c("none", "relative", "absolute"), ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name) > 0))

    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    as.is <- match.arg(as.is)
    var.index <- match.arg(var.index)
    var.index <- match(var.index, c("none", "relative", "absolute"))

    if (margin == "by.variant")
    {
        # C call
        rv <- .Call(SEQ_Apply_Variant, gdsfile, var.name, FUN, as.is,
            var.index, new.env())
        if (as.is == "none") return(invisible())
    } else if (margin == "by.sample")
    {
        # C call
        rv <- .Call(SEQ_Apply_Sample, gdsfile, var.name, FUN, as.is,
            var.index, new.env())
        if (as.is == "none") return(invisible())
    }
    rv
}



#######################################################################
# Apply functions via a sliding window over variants
#
seqSlidingWindow <- function(gdsfile, var.name, win.size, shift=1, FUN,
    as.is = c("list", "integer", "double", "character", "none"),
    var.index = c("none", "relative", "absolute"), ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name) > 0))

    stopifnot(is.numeric(win.size) & (length(win.size)==1))
    win.size <- as.integer(win.size)
    stopifnot(is.finite(win.size))
    if (win.size < 1)
        stop("`win.size' should be greater than 0.")

    stopifnot(is.numeric(shift) & (length(shift)==1))
    shift <- as.integer(shift)
    stopifnot(is.finite(shift))
    if (shift < 1)
        stop("`shift' should be greater than 0.")

    as.is <- match.arg(as.is)
    var.index <- match.arg(var.index)
    var.index <- match(var.index, c("none", "relative", "absolute"))

    FUN <- match.fun(FUN)

    # C call
    rv <- .Call(SEQ_SlidingWindow, gdsfile, var.name, win.size, shift,
        FUN, as.is, var.index, new.env())

    if (as.is == "none") return(invisible())
    rv
}



#######################################################################
# Summarize the GDS file
#
seqSummary <- function(gdsfile, varname=NULL,
    check=c("check", "full.check", "none"), verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(varname) | is.character(varname))
    if (is.character(varname))
        stopifnot(length(varname) == 1)
    check <- match.arg(check)

    # initialize a GDS object
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1)
        gds <- seqOpen(gdsfile)
        on.exit(seqClose(gds))   
    } else {
        gds <- gdsfile
    }


    # whether specify the variable name or not
    if (is.null(varname))
    {
        ########################################################
        check.dim <- function(node, dim.len)
        {
            dm <- objdesp.gdsn(node)$dim
            if (check != "")
            {
                if (length(dm) != dim.len)
                {
                    print(node)
                    stop("Error in dimension!")
                }
            }
            dm
        }


        ########################################################
        # initialize
        ans <- list()
        ans$filename <- gds$filename
        if (verbose)
            cat("File: ", ans$filename, "\n", sep="")


        ########################################################
        # check header
        n <- index.gdsn(gds, "description")
        tmp <- get.attr.gdsn(n)
        if ("sequence.variant.format" %in% names(tmp))
        {
            ans$format.version <- tmp$sequence.variant.format
            if (verbose)
                cat("Format Version: ", ans$format.version, "\n", sep="")
        }

        ########################################################
        # the number of samples
        n <- index.gdsn(gds, "sample.id")
        dm <- check.dim(n, 1)
        n.samp <- dm[1]
        ans$num.of.sample <- n.samp

        ########################################################
        # the number of variants
        n <- index.gdsn(gds, "variant.id")
        dm <- check.dim(n, 1)
        n.variant <- dm[1]
        ans$num.of.variant <- n.variant

        if (verbose)
        {
            cat("The number of samples: ", n.samp, "\n", sep="")
            cat("The number of variants: ", n.variant, "\n", sep="")
        }


        ########################################################
        # position
        n <- index.gdsn(gds, "position")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'position' variable.")

        ########################################################
        # chromosome
        n <- index.gdsn(gds, "chromosome")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'chromosome' variable.")
        if (verbose)
        {
            chr <- seqGetData(gds, "chromosome")
            tab <- table(chr, exclude=NULL)
            names(dimnames(tab)) <- "The chromosomes:"
            print(tab)
            rm(chr)
        }

        ########################################################
        # allele
        n <- index.gdsn(gds, "allele")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'allele' variable.")
        if ((check != "") | verbose)
            nallele <- .Call(SEQ_NumOfAllele, n)
        if (verbose)
        {
            tab <- table(nallele)
            names(dimnames(tab)) <- "The number of alleles per site:"
            print(tab)
        }


        ########################################################
        # genotype
        n <- index.gdsn(gds, "genotype/data")
        dm <- check.dim(n, 3)
        if (dm[2] != n.samp)
        {
            print(n)
            stop("Invalid number of samples.")
        }


        ########################################################
        # phase
        n <- index.gdsn(gds, "phase/data")
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) > 2)
            dm <- dm[-c(length(dm)-1, length(dm))]
        if (dm[1] != n.samp)
        {
            print(n)
            stop("Invalid number of samples.")
        }
        if (dm[2] != n.variant)
        {
            print(n)
            stop("Invalid number of variants.")
        }


        ########################################################
        # annotation/id
        n <- index.gdsn(gds, "annotation/id", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }

        ########################################################
        # annotation/qual
        n <- index.gdsn(gds, "annotation/qual", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }

        ########################################################
        # annotation/filter
        n <- index.gdsn(gds, "annotation/filter", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }


        ########################################################
        # annotation/info
        n <- index.gdsn(gds, "annotation/info", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            vars <- unlist(strsplit(vars, "@", fixed=TRUE))
            vars <- unique(vars[vars != ""])
            dat <- NULL
            if (verbose)
                cat("Annotation, INFO variable(s):\n")
            for (nm in vars)
            {
                n <- index.gdsn(gds, paste("annotation/info/", nm, sep=""))
                a <- get.attr.gdsn(n)
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ", a$Description,
                        "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name = nm,
                    number = a$Number, type = a$Type, description = a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$info <- dat
        }


        ########################################################
        # annotation/format
        n <- index.gdsn(gds, "annotation/format", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- NULL
            if (verbose)
                cat("Annotation, FORMAT variable(s):\n")
            for (nm in vars)
            {
                n <- index.gdsn(gds, paste("annotation/format/", nm, sep=""))
                a <- get.attr.gdsn(n)
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ", a$Description,
                        "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name = nm,
                    number = a$Number, type = a$Type, description = a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$format <- dat
        }


        ########################################################
        # sample.annotation/format
        n <- index.gdsn(gds, "sample.annotation", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- data.frame(var.name=vars, stringsAsFactors=FALSE)
            if (verbose)
                cat("Annotation, sample variable(s):", vars, "\n")
            ans$sample.annot <- dat
        }

        # output
        return(invisible(ans))

    } else {

        # get a description of variable
        .Call(SEQ_Summary, gds, varname)
    }
}



#######################################################################
# Data Analysis
#######################################################################

#######################################################################
# The number of alleles per site
#
seqNumAllele <- function(gdsfile, parallel=getOption("gds.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    seqParallel(parallel, gdsfile, split="by.variant",
        FUN = function(gdsfile)
        {
            seqApply(gdsfile, "allele", margin="by.variant", as.is="integer",
                FUN = .cfunction("FC_NumAllele"))
        })
}



#######################################################################
# Missing rate
#
seqMissing <- function(gdsfile, per.variant=TRUE,
    parallel=getOption("gds.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(per.variant))

    if (per.variant)
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(gdsfile)
            {
               seqApply(gdsfile, "genotype", margin="by.variant",
                   as.is="double", FUN=.cfunction("FC_Missing_PerVariant"))
            })
    } else {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants

        sum <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(gdsfile, num)
            {
                tmpsum <- integer(num)
                seqApply(gdsfile, "genotype", margin="by.variant",
                   as.is="none", FUN=.cfunction2("FC_Missing_PerSample"),
                   y=tmpsum)
                tmpsum
            }, .combine="+", num=dm[1])
        sum / (2L * dm[2])
    }
}



#######################################################################
# Allele frequency
#
seqAlleleFreq <- function(gdsfile, ref.allele=0L,
    parallel=getOption("gds.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | is.numeric(ref.allele) | is.character(ref.allele))

    if (is.null(ref.allele))
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(gdsfile)
            {
                seqApply(gdsfile, c("genotype", "allele"), margin="by.variant",
                    as.is="list", FUN = .cfunction("FC_AF_List"))
            })
    } else if (is.numeric(ref.allele))
    {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
        if (!(length(ref.allele) %in% c(1L, dm[2L])))
            stop("'length(ref.allele)' should be 1 or the number of selected variants.")

        if (length(ref.allele) == 1L)
        {
            seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(gdsfile, ref)
                {
                    .cfunction("FC_AF_SetIndex")(ref)
                    seqApply(gdsfile, c("genotype", "allele"), margin="by.variant",
                        as.is="double", FUN = .cfunction("FC_AF_Index"))
                }, ref=ref.allele)
        } else {
            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(gdsfile, selflag, ref)
                {
                    .cfunction("FC_AF_SetIndex")(ref[selflag])
                    seqApply(gdsfile, c("genotype", "allele"), margin="by.variant",
                        as.is="double", FUN = .cfunction("FC_AF_Index"))
                }, ref=ref.allele)
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
        if (length(ref.allele) != dm[2L])
            stop("'length(ref.allele)' should be the number of selected variants.")

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(gdsfile, selflag, ref)
            {
                .cfunction("FC_AF_SetAllele")(ref[selflag])
                seqApply(gdsfile, c("genotype", "allele"), margin="by.variant",
                    as.is="double", FUN = .cfunction("FC_AF_Allele"))
            }, ref=ref.allele)
    }
}