#######################################################################
#
# Package Name: SeqArray
#
# The interface to the class "SeqVarGDSClass"
#


#######################################################################
# Open a SeqArray GDS file
#
seqOpen <- function(gds.fn, readonly=TRUE)
{
    # check
    stopifnot(is.character(gds.fn), length(gds.fn)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)

    # open the file
    ans <- openfn.gds(gds.fn, readonly=readonly, allow.fork=TRUE)

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
setMethod("seqSetFilter", signature(object="SeqVarGDSClass"),
    function(object, sample.id=NULL, variant.id=NULL, samp.sel=NULL,
        variant.sel=NULL, action=c("set", "intersect", "push", "push+set",
        "push+intersect", "pop"), verbose=TRUE)
    {
        # check
        action <- match.arg(action)
        stopifnot(is.logical(verbose))

        setflag <- FALSE
        switch(action,
            "set" = NULL,
            "intersect" = { setflag <- TRUE },
            "push" = {
                if (!all(is.null(sample.id), is.null(variant.id),
                        is.null(samp.sel), is.null(variant.sel)))
                {
                    stop("The arguments 'sample.id', 'variant.id', ",
                        "'samp.sel' and 'variant.sel' should be NULL.")
                }
                .Call(SEQ_FilterPushLast, object)
            },
            "push+set" = {
                .Call(SEQ_FilterPushEmpty, object)
            },
            "push+intersect" = {
                .Call(SEQ_FilterPushEmpty, object)
                setflag <- TRUE
            },
            "pop" = {
                if (!all(is.null(sample.id), is.null(variant.id),
                        is.null(samp.sel), is.null(variant.sel)))
                {
                    stop("The arguments 'sample.id', 'variant.id', ",
                        "'samp.sel' and 'variant.sel' should be NULL.")
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
        } else if (!is.null(samp.sel))
        {
            stopifnot(is.vector(samp.sel) & is.logical(samp.sel))
            .Call(SEQ_SetSpaceSample, object, samp.sel, setflag, verbose)
        }

        if (!is.null(variant.id))
        {
            stopifnot(is.vector(variant.id))
            stopifnot(is.numeric(variant.id) | is.character(variant.id))
            .Call(SEQ_SetSpaceVariant, object, variant.id, setflag, verbose)
        } else if (!is.null(variant.sel))
        {
            stopifnot(is.vector(variant.sel) & is.logical(variant.sel))
            .Call(SEQ_SetSpaceVariant, object, variant.sel, setflag, verbose)
        } else {
            if (is.null(sample.id) & is.null(samp.sel))
            {
                .Call(SEQ_SetSpaceSample, object, NULL, setflag, verbose)
                .Call(SEQ_SetSpaceVariant, object, NULL, setflag, verbose)
            }
        }

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
    stopifnot(is.character(var.name), length(var.name)==1L)

    .Call(SEQ_GetData, gdsfile, var.name)
}



#######################################################################
# Apply functions over margins on a working space with selected samples and variants
#
seqApply <- function(gdsfile, var.name, FUN,
    margin = c("by.variant", "by.sample"), as.is = c("none", "list",
    "integer", "double", "character", "logical", "raw"),
    var.index = c("none", "relative", "absolute"),
    .useraw=FALSE, .list_dup=TRUE, ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name)>0L)

    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    as.is <- match.arg(as.is)
    var.index <- match.arg(var.index)
    var.index <- match(var.index, c("none", "relative", "absolute"))

    if (margin == "by.variant")
    {
        # C call
        rv <- .Call(SEQ_Apply_Variant, gdsfile, var.name, FUN, as.is,
            var.index, .useraw, .list_dup, new.env())
        if (as.is == "none") return(invisible())
    } else if (margin == "by.sample")
    {
        # C call
        rv <- .Call(SEQ_Apply_Sample, gdsfile, var.name, FUN, as.is,
            var.index, .useraw, new.env())
        if (as.is == "none") return(invisible())
    }
    rv
}



#######################################################################
# Apply functions via a sliding window over variants
#
.seqSlidingWindow <- function(gdsfile, var.name, win.size, shift=1, FUN,
    as.is = c("list", "integer", "double", "character", "none"),
    var.index = c("none", "relative", "absolute"), ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name), length(var.name) > 0L)

    stopifnot(is.numeric(win.size), length(win.size)==1L)
    win.size <- as.integer(win.size)
    stopifnot(is.finite(win.size))
    if (win.size < 1)
        stop("`win.size' should be greater than 0.")

    stopifnot(is.numeric(shift), length(shift)==1L)
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
        stopifnot(length(varname) == 1L)
    check <- match.arg(check)

    # initialize a GDS object
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1L)
        gds <- seqOpen(gdsfile)
        on.exit(seqClose(gds))   
    } else {
        gds <- gdsfile
    }


    # whether specify the variable name or not
    if (is.null(varname))
    {
        ########
        check_dim <- function(node, dimlen)
        {
            dm <- objdesp.gdsn(node)$dim
            if (check != "none")
            {
                if (length(dm) != dimlen)
                {
                    print(node)
                    stop("Error in dimension!")
                }
            }
            dm
        }

        STOP <- function(...)
        {
            if (check != "none") stop(...)
        }
        STOP2 <- function(node, ...)
        {
            if (check != "none")
                { print(node); stop(...) }
        }


        ########
        # initialize
        ans <- list()
        ans$filename <- gds$filename
        if (verbose)
            cat("File: ", ans$filename, "\n", sep="")


        ########
        # check header
        tmp <- get.attr.gdsn(gds$root)
        if ("FileVersion" %in% names(tmp))
        {
            ans$format.version <- tmp$FileVersion
            if (verbose)
                cat("Format Version: ", ans$format.version, "\n", sep="")
        }

        ########
        # the number of samples
        n <- index.gdsn(gds, "sample.id")
        dm <- check_dim(n, 1L)
        n.samp <- dm[1L]
        ans$num.of.sample <- n.samp

        ########
        # the number of variants
        n <- index.gdsn(gds, "variant.id")
        dm <- check_dim(n, 1L)
        n.variant <- dm[1L]
        ans$num.of.variant <- n.variant

        if (verbose)
        {
            cat("The number of samples: ", n.samp, "\n", sep="")
            cat("The number of variants: ", n.variant, "\n", sep="")
        }


        ########
        # position
        n <- index.gdsn(gds, "position")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'position' variable.")

        ########
        # chromosome
        n <- index.gdsn(gds, "chromosome")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'chromosome' variable.")
        if (verbose & (check != "none"))
        {
            chr <- seqGetData(gds, "chromosome")
            tab <- table(chr, exclude=NULL)
            names(dimnames(tab)) <- "The chromosomes:"
            print(tab)
            rm(chr)
        }

        ########
        # allele
        n <- index.gdsn(gds, "allele")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'allele' variable.")
        if (check != "none")
        {
            nallele <- .Call(SEQ_NumOfAllele, n)
            if (verbose)
            {
                tab <- table(nallele)
                names(dimnames(tab)) <- "The number of alleles per site:"
                print(tab)
            }
        }


        ########
        # genotype
        n <- index.gdsn(gds, "genotype/data")
        dm <- check_dim(n, 3L)
        if (dm[2L] != n.samp)
            STOP2(n, "Invalid dimension of sample in 'genotype/data'.")
        if (dm[3L] < n.variant)
            STOP2(n, "Invalid dimension of variant in 'genotype/data'.")

        n <- index.gdsn(gds, "genotype/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 3L)
            if (dm[2L] < n.variant)
                STOP(n, "Invalid dimension of variant in 'genotype/~data'.")
            if (dm[3L] != n.samp)
                STOP(n, "Invalid dimension of sample in 'genotype/~data'.")
        }


        ########
        # phase
        n <- index.gdsn(gds, "phase/data")
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) > 2L)
            dm <- dm[-c(length(dm)-1L, length(dm))]
        if (dm[1L] != n.samp)
            STOP2(n, "Invalid dimension of sample in 'phase/data'.")
        if (dm[2L] != n.variant)
            STOP2(n, "Invalid dimension of variant in 'phase/data'.")

        n <- index.gdsn(gds, "phase/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- objdesp.gdsn(n)$dim
            if (length(dm) > 2L)
                dm <- dm[-c(length(dm)-1L, length(dm))]
            if (dm[1L] != n.variant)
                STOP2(n, "Invalid dimension of variant in 'phase/~data'.")
            if (dm[2L] != n.samp)
                STOP2(n, "Invalid dimension of sample in 'phase/~data'.")
        }


        ########
        # annotation/id
        n <- index.gdsn(gds, "annotation/id", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }

        ########
        # annotation/qual
        n <- index.gdsn(gds, "annotation/qual", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }

        ########
        # annotation/filter
        n <- index.gdsn(gds, "annotation/filter", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }


        ########
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
                n1 <- index.gdsn(n, nm)
                d <- objdesp.gdsn(n1)
                a <- get.attr.gdsn(n1)
                n2 <- index.gdsn(n, paste("@", nm, sep=""), silent=TRUE)
                if (is.null(n2))
                {
                    if (d$dim[length(d$dim)] != n.variant)
                        STOP2(n1, "Invalid dimension of variant.")
                } else {
                    dm <- check_dim(n2, 1L)
                    if (dm != n.variant)
                        STOP2(n1, "Invalid dimension of variant.")
                }
                if (is.null(a$Number))
                {
                    if (is.null(n2))
                        a$Number <- prod(d$dim[-length(d$dim)])
                    else
                        a$Number <- "A"
                }
                if (is.null(a$Type))
                {
                    a$Type <- .vcf_type(n1)
                }
                if (is.null(a$Description))
                {
                    a$Description <- "."
                }
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ",
                        a$Description, "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name=nm,
                    number=a$Number, type=a$Type, description=a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$info <- dat
        }


        ########
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


        ########
        # sample.annotation
        n <- index.gdsn(gds, "sample.annotation", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- data.frame(var.name=vars, type=rep(".", length(vars)),
                stringsAsFactors=FALSE)
            for (i in seq_len(length(vars)))
                dat$type[i] <- .vcf_type(index.gdsn(n, vars[i]))
            if (verbose)
            {
                cat("Annotation, sample variable(s):\n")
                cat(paste("\t", vars, sep=""), sep="\n")
            }
            ans$sample.annot <- dat
        }

        # output
        invisible(ans)

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
seqNumAllele <- function(gdsfile, parallel=getOption("seqarray.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    seqParallel(parallel, gdsfile, split="by.variant",
        FUN = function(f)
        {
            seqApply(f, "allele", margin="by.variant", as.is="integer",
                FUN = .cfunction("FC_NumAllele"))
        })
}



#######################################################################
# Missing rate
#
seqMissing <- function(gdsfile, per.variant=TRUE,
    parallel=getOption("seqarray.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(per.variant))

    if (per.variant)
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f)
            {
               seqApply(f, "genotype", margin="by.variant",
                   as.is="double", FUN=.cfunction("FC_Missing_PerVariant"))
            })
    } else {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants

        sum <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f, num)
            {
                tmpsum <- integer(num)
                seqApply(f, "genotype", margin="by.variant",
                    as.is="none", FUN=.cfunction2("FC_Missing_PerSample"),
                    y=tmpsum)
                tmpsum
            }, .combine="+", num=dm[1L])
        sum / (2L * dm[2L])
    }
}



#######################################################################
# Allele frequency
#
seqAlleleFreq <- function(gdsfile, ref.allele=0L,
    parallel=getOption("seqarray.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(ref.allele) | is.numeric(ref.allele) |
        is.character(ref.allele))

    if (is.null(ref.allele))
    {
        seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(f)
            {
                seqApply(f, c("genotype", "allele"), margin="by.variant",
                    as.is="list", FUN = .cfunction("FC_AF_List"),
                    .list_dup=FALSE)
            })
    } else if (is.numeric(ref.allele))
    {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
        if (!(length(ref.allele) %in% c(1L, dm[2L])))
        {
            stop("'length(ref.allele)' should be 1 or the number of selected variants.")
        }

        if (length(ref.allele) == 1L)
        {
            seqParallel(parallel, gdsfile, split="by.variant",
                FUN = function(f, ref)
                {
                    .cfunction("FC_AF_SetIndex")(ref)
                    seqApply(f, c("genotype", "allele"), margin="by.variant",
                        as.is="double", FUN = .cfunction("FC_AF_Index"))
                }, ref=ref.allele)
        } else {
            ref.allele <- as.integer(ref.allele)
            seqParallel(parallel, gdsfile, split="by.variant",
                .selection.flag=TRUE,
                FUN = function(f, selflag, ref)
                {
                    s <- ref[selflag]
                    .cfunction("FC_AF_SetIndex")(s)
                    seqApply(f, c("genotype", "allele"), margin="by.variant",
                        as.is="double", FUN = .cfunction("FC_AF_Index"))
                }, ref=ref.allele)
        }
    } else if (is.character(ref.allele))
    {
        dm <- .seldim(gdsfile)
        # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
        if (length(ref.allele) != dm[2L])
        {
            stop("'length(ref.allele)' should be the number of selected variants.")
        }

        seqParallel(parallel, gdsfile, split="by.variant",
            .selection.flag=TRUE,
            FUN = function(f, selflag, ref)
            {
                s <- ref[selflag]
                .cfunction("FC_AF_SetAllele")(s)
                seqApply(f, c("genotype", "allele"), margin="by.variant",
                    as.is="double", FUN = .cfunction("FC_AF_Allele"))
            }, ref=ref.allele)
    }
}



#######################################################################
# Allele Counts
#
seqAlleleCount <- function(gdsfile,
    parallel=getOption("seqarray.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    seqParallel(parallel, gdsfile, split="by.variant",
        FUN = function(f)
        {
            seqApply(f, c("genotype", "allele"), margin="by.variant",
                as.is="list", FUN = .cfunction("FC_AlleleCount"))
        })
}



#######################################################################
# IBD
#
ssIBD <- function(gdsfile, method=c("OneLocus", "TwoLoci"), interval=1L,
    parallel=getOption("seqarray.parallel", FALSE))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    method <- match.arg(method)

    if (method == "OneLocus")
    {
        v <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(gdsfile)
            {
                n <- .seldim(gdsfile)[1L]
                size <- n * (n + 1L) / 2L
                # m[,1] -- numerator, m[,2] -- denominator
                m <- matrix(0.0, nrow=size, ncol=2L)
                seqApply(gdsfile, "genotype", margin = "by.variant",
                    as.is = "none", FUN = .cfunction3("FC_IBD_OneLocus"),
                    y = m, z = raw(size), .useraw=TRUE)
                m
            }, .combine="+")
    } else if (method == "TwoLoci")
    {
        v <- seqParallel(parallel, gdsfile, split="by.variant",
            FUN = function(gdsfile, interval)
            {
                n <- .seldim(gdsfile)[1L]
                size <- n * (n + 1L) / 2L
                # m[,1] -- numerator, m[,2] -- denominator
                m <- matrix(0.0, nrow=size, ncol=2L)
                .cfunction2("FC_IBD_TwoLoci_Init")(interval, n)
                # apply
                seqApply(gdsfile, "genotype", margin = "by.variant",
                    as.is = "none", FUN = .cfunction3("FC_IBD_TwoLoci"),
                    y = m, z = raw(size), .useraw=TRUE)
                m
            }, .combine="+", interval=interval)
    }

    .cfunction2("FC_IBD_Div")(v, .seldim(gdsfile)[1L])
}
