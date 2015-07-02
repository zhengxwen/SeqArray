#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Big Data Management of Sequencing-based Genetic Variants
#



#######################################################################
# Get the file name of an example
#
seqExampleFileName <- function(type=c("gds", "vcf", "KG_Phase1"))
{
    type <- match.arg(type)
    switch(type,
        gds = system.file("extdata", "CEU_Exon.gds", package="SeqArray"),
        vcf = system.file("extdata", "CEU_Exon.vcf.gz", package="SeqArray"),
        KG_Phase1 =
            system.file("extdata", "1KG_phase1_release_v3_chr22.gds",
            package="SeqArray")
    )
}



#######################################################################
# Merge multiple GDS files
#
seqMerge <- function(gds.fn, out.fn, compress.option = seqCompress.Option(),
    verbose = TRUE)
{
    # check
    stopifnot(is.character(gds.fn) & (length(gds.fn)>0))
    stopifnot(is.character(out.fn) & (length(out.fn)==1))
    stopifnot(is.logical(verbose) & (length(verbose)==1))


    # define functions
    compress <- function(var.name)
    {
        if (var.name %in% names(compress.option))
            compress.option[[var.name]]
        else
            ""
    }


    ##################################################
    # open GDS files

    opfile <- vector("list", length(gds.fn))
    on.exit({
        for (i in 1:length(opfile))
        {
            if (!is.null(opfile[[i]]))
                closefn.gds(opfile[[i]])
        }
    })

    # for - loop
    samp.list <- list()
    variant.num <- 0
    for (i in 1:length(gds.fn))
    {
        if (verbose)
            cat(sprintf("Open (%02d): %s\n", i, gds.fn[i]))
        opfile[[i]] <- seqOpen(gds.fn[i])

        samp.list[[i]] <- read.gdsn(index.gdsn(opfile[[i]], "sample.id"))
        desp <- objdesp.gdsn(index.gdsn(opfile[[i]], "variant.id"))
        variant.num <- variant.num + prod(desp$dim)
    }

    # merge all sample id
    samp.id <- unique(unlist(samp.list))
    if (verbose)
    {
        cat(sprintf("Input: %d samples, %d variants.\n",
            length(samp.id), variant.num))
    }


    ##################################################
    # create a GDS file

    gfile <- createfn.gds(out.fn)
    on.exit({ closefn.gds(gfile) }, add=TRUE)
    if (verbose)
        cat("Output: ", out.fn, "\n", sep="")

    # add header
    n <- add.gdsn(gfile, name="description", storage="folder")
    put.attr.gdsn(n, "sequence.variant.format", "v1.0")

    # add sample id
    add.gdsn(gfile, "sample.id", samp.id, compress=compress("sample.id"), closezip=TRUE)

    # add variant.id
    # TODO: 
    node <- add.gdsn(gfile, "variant.id", storage="int32", valdim=c(0),
        compress=compress("variant.id"))
    for (i in 1:length(opfile))
        assign.gdsn(node, index.gdsn(opfile[[i]], "variant.id"), TRUE)
    readmode.gdsn(node)

    # add position
    # TODO: need to check whether position can be stored in 'int32'
    node <- add.gdsn(gfile, "position", storage="int32", valdim=c(0),
        compress=compress("position"))
    for (i in 1:length(opfile))
        assign.gdsn(node, index.gdsn(opfile[[i]], "position"), TRUE)
    readmode.gdsn(node)

    # add chromosome
    node <- add.gdsn(gfile, "chromosome", storage="string", valdim=c(0),
        compress=compress("chromosome"))
    for (i in 1:length(opfile))
        assign.gdsn(node, index.gdsn(opfile[[i]], "chromosome"), TRUE)
    readmode.gdsn(node)

    # add allele
    node <- add.gdsn(gfile, "allele", storage="string", valdim=c(0),
        compress=compress("allele"))
    for (i in 1:length(opfile))
        assign.gdsn(node, index.gdsn(opfile[[i]], "allele"), TRUE)
    readmode.gdsn(node)

    # sync file
    sync.gds(gfile)


    # add a folder for genotypes
    att <- get.attr.gdsn(index.gdsn(opfile[[1]], "genotype"))
    varGeno <- add.gdsn(gfile, name="genotype", storage="folder")
    for (i in 1:length(att))
        put.attr.gdsn(varGeno, names(att)[i], att[[i]])
    readmode.gdsn(node)

    # add length data to the folder of genotype
    node <- add.gdsn(varGeno, "length", storage="int32", valdim=c(0),
        compress=compress("genotype"))
    for (i in 1:length(opfile))
        assign.gdsn(node, index.gdsn(opfile[[i]], "genotype/length"), TRUE)
    readmode.gdsn(node)

    # add genotypic data to the folder of genotype
    desp <- objdesp.gdsn(index.gdsn(opfile[[1]], "genotype/data"))
    desp$dim[length(desp$dim)] <- 0
    node <- add.gdsn(varGeno, "data", storage="bit2", valdim=desp$dim,
        compress=compress("genotype"))
    for (i in 1:length(opfile))
    {
        assign.gdsn(node, index.gdsn(opfile[[i]], "genotype/data"), TRUE)
        if (verbose)
            cat(sprintf("\tGenotype (%02d) done.\n", i))
    }
    readmode.gdsn(node)


    # add phase folder
    if (!is.null(index.gdsn(opfile[[1]], "phase", FALSE)))
    {
        varPhase <- add.gdsn(gfile, name="phase", storage="folder")

        # add data
        desp <- objdesp.gdsn(index.gdsn(opfile[[1]], "phase/data"))
        desp$dim[length(desp$dim)] <- 0
        node <- add.gdsn(varPhase, "data", storage="bit1", valdim=desp$dim,
            compress=compress("phase"))
        for (i in 1:length(opfile))
        {
            assign.gdsn(node, index.gdsn(opfile[[i]], "phase/data"), TRUE)
            if (verbose)
                cat(sprintf("\tPhase (%02d) done.\n", i))
        }
        readmode.gdsn(node)
    }

    # return
    invisible()
}



#######################################################################
# Read the header of a VCF file
#
seqCompress.Option <- function(default="ZIP_RA.MAX", ...)
{
    rv <- list(...)
    n <- names(rv)

    # mandatory

    if (!("description" %in% n))
        rv$description <- default

    if (!("sample.id" %in% n))
        rv$sample.id <- default
    if (!("variant.id" %in% n))
        rv$variant.id <- default
    if (!("position" %in% n))
        rv$position <- default
    if (!("chromosome" %in% n))
        rv$chromosome <- default
    if (!("allele" %in% n))
        rv$allele <- default

    if (!("genotype" %in% n))
        rv$genotype <- default
    if (!("genotype.extra" %in% n))
        rv$genotype.extra <- default

    if (!("phase" %in% n))
        rv$phase <- default
    if (!("phase.extra" %in% n))
        rv$phase.extra <- default

    # options

    if (!("id" %in% n))
        rv$id <- default
    if (!("qual" %in% n))
        rv$qual <- default
    if (!("filter" %in% n))
        rv$filter <- default
    if (!("info" %in% n))
        rv$info <- default
    if (!("format" %in% n))
        rv$format <- default

    if (!("annotation" %in% n))
        rv$annotation <- default

    class(rv) <- "SeqGDSCompressFlagClass"
    return(rv)
}



#######################################################################
# Apply functions in parallel
#
seqParallel <- function(cl, gdsfile, FUN, split=c("by.variant", "by.sample", "none"),
    .combine="unlist", .selection.flag=FALSE, ...)
{
    # check
    stopifnot(is.null(cl) | is.logical(cl) | is.numeric(cl) | inherits(cl, "cluster"))
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.function(FUN))
    split <- match.arg(split)
    stopifnot(is.character(.combine) | is.function(.combine))
    stopifnot(is.logical(.selection.flag))

    if (is.character(.combine))
    {
        stopifnot(length(.combine) == 1L)
        if (!(.combine %in% c("unlist", "list", "none")))
            .combine <- match.fun(.combine)
    }

    if (is.null(cl) | identical(cl, FALSE) | identical(cl, 1L) | identical(cl, 1))
    {
        #################################################################
        # a single process

        if (.selection.flag)
        {
            dm <- .seldim(gdsfile)
            # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
            if (split == "by.variant")
                ans <- FUN(gdsfile, rep(TRUE, dm[2L]), ...)
            else if (split == "by.sample")
                ans <- FUN(gdsfile, rep(TRUE, dm[1L]), ...)
            else
                ans <- FUN(gdsfile, NULL, ...)
        } else
            ans <- FUN(gdsfile, ...)

    } else if (identical(cl, TRUE) | is.numeric(cl))
    {
        # library
        if (!requireNamespace("parallel"))
            stop("the 'parallel' package should be installed.")

        if (identical(cl, TRUE))
        {
            cl <- parallel::detectCores() - 1L
            if (cl <= 1L) cl <- 2L
        }
        stopifnot(length(cl) == 1L)
        if (cl <= 0L)
            cl <- getOption("mc.cores", 2L)
        if (.Platform$OS.type == "windows")
            cl <- 1L
        if (cl <= 0L)
        	stop("Invalid number of cores!")

        if (split %in% c("by.variant", "by.sample"))
        {
            dm <- .seldim(gdsfile)
            # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
            if (split == "by.variant")
            {
                if (dm[2L] <= 0) stop("No variants selected.")
                if (cl > dm[2L]) cl <- dm[2L]
            } else {
                if (dm[1L] <= 0) stop("No samples selected.")
                if (cl > dm[1L]) cl <- dm[1L]
            }
        }

        ans <- parallel::mclapply(seq_len(cl),
            mc.preschedule=FALSE, mc.cores=cl, mc.cleanup=TRUE,
            FUN = function(i, .fun)
            {
                sel <- .Call(SEQ_SplitSelection, gdsfile, split, i, cl,
                    .selection.flag)
                # call the user-defined function
                if (.selection.flag)
                    FUN(gdsfile, sel, ...)
                else
                    FUN(gdsfile, ...)
            }, .fun = FUN)

        if (is.list(ans))
        {
            if (identical(.combine, "unlist"))
            {
                ans <- unlist(ans, recursive=FALSE)
            } else if (is.function(.combine))
            {
                rv <- ans[[1]]
                for (i in seq_len(length(ans)-1L) + 1L)
                    rv <- .combine(rv, ans[[i]])
                ans <- rv
            }
        }

    } else if (inherits(cl, "cluster"))
    {
        #################################################################
        # multiple processes

        # library
        if (!requireNamespace("parallel"))
            stop("the 'parallel' package should be installed.")

        if (split %in% c("by.variant", "by.sample"))
        {
            dm <- .seldim(gdsfile)
            # dm[1] -- Num of selected samples, dm[2] -- Num of selected variants
            if (split == "by.variant")
            {
                if (dm[2L] <= 0) stop("No variants selected.")
                if (length(cl) > dm[2L]) cl <- cl[seq_len(dm[2L])]
            } else {
                if (dm[1L] <= 0) stop("No samples selected.")
                if (length(cl) > dm[1L]) cl <- cl[seq_len(dm[1L])]
            }
        }

        ans <- .DynamicClusterCall(cl, length(cl), .fun =
            function(.idx, .n_process, .gds.fn, .selection, FUN, .split, .selection.flag, ...)
        {
            # load the package
            library("SeqArray")

            # open the file
            gfile <- seqOpen(.gds.fn)
            on.exit({ closefn.gds(gfile) })

            # set filter
            seqSetFilter(gfile, samp.sel=.selection$sample.sel,
                variant.sel=.selection$variant.sel)

            sel <- .Call(SEQ_SplitSelection, gfile, .split, .idx, .n_process,
                .selection.flag)
            # call the user-defined function
            if (.selection.flag)
                FUN(gdsfile, sel, ...)
            else
                FUN(gdsfile, ...)

        }, .combinefun = .combine, .stopcluster=FALSE,
            .n_process = length(cl), .gds.fn = gdsfile$filename,
            .selection = seqGetFilter(gdsfile), FUN = FUN,
            .split = split, .selection.flag=.selection.flag, ...
        )

        if (is.list(ans) & identical(.combine, "unlist"))
            ans <- unlist(ans, recursive=FALSE)
    } else {
        stop("Invalid 'cl' in seqParallel.")
    }

    # output
    if (identical(.combine, "none"))
        invisible()
    else
        ans
}



#######################################################################
# Modify the SeqArray data structure
#######################################################################

#######################################################################
# Delete data variables
#

seqDelete <- function(gdsfile, info.varname=character(),
    format.varname=character(), verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    if (gdsfile$readonly)
        stop("The GDS file is read-only.")

    stopifnot(is.character(info.varname))
    stopifnot(is.character(format.varname))
    stopifnot(is.logical(verbose))

    if (verbose) cat("Deleting INFO variable(s):")
    for (nm in info.varname)
    {
        n <- index.gdsn(gdsfile, paste0("annotation/info/", nm))
        delete.gdsn(n, force=TRUE)
        n <- index.gdsn(gdsfile, paste0("annotation/info/@", nm), silent=TRUE)
        if (!is.null(n))
            delete.gdsn(n, force=TRUE)
        if (verbose) cat("", nm)
    }
    if (verbose) cat("\n")

    if (verbose) cat("Deleting FORMAT variable(s):")
    for (nm in format.varname)
    {
        n <- index.gdsn(gdsfile, paste0("annotation/format/", nm))
        delete.gdsn(n, force=TRUE)
        if (verbose) cat("", nm)
    }
    if (verbose) cat("\n")

    # return
    invisible()
}



#######################################################################
# Transpose data variable(s)
#

.Transpose <- function(gdsfile, src.fn, prefix, compress=NULL)
{
    dst.fn <- .var_path(src.fn, prefix)
    if (is.null(index.gdsn(gdsfile, dst.fn, silent=TRUE)))
    {
        node <- index.gdsn(gdsfile, src.fn)
        desp <- objdesp.gdsn(node)
        dm <- desp$dim
        if (length(dm) > 1L)
        {
            # dimension
            dm <- c(dm[-(length(dm)-1L)], 0L)
            # folder
            nm <- unlist(strsplit(src.fn, "/"))
            if (length(nm) <= 1)
                folder <- gdsfile$root
            else
                folder <- index.gdsn(gdsfile, index=nm[-length(nm)])
            # compress
            if (is.null(compress))
                compress <- desp$compress

            pm <- list(node = folder,
                name = paste(prefix, nm[length(nm)], sep=""),
                val = NULL, storage = desp$storage,
                valdim = dm, compress = compress)
            if (!is.null(desp$param))
                pm <- c(pm, desp$param)

            newnode <- do.call(add.gdsn, pm)
            moveto.gdsn(newnode, node, relpos="after")

            # write data
            apply.gdsn(node, margin=length(dm)-1L, as.is="gdsnode",
                FUN=`c`, target.node=newnode, .useraw=TRUE)

            readmode.gdsn(newnode)
        }
    }
    invisible()
}


seqTranspose <- function(gdsfile, var.name, compress=NULL, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & is.vector(var.name))
    stopifnot(length(var.name) == 1L)

    node <- index.gdsn(gdsfile, var.name)
    desp <- objdesp.gdsn(node)
    dm <- desp$dim
    if (length(dm) > 1L)
    {
        # dimension
        dm <- c(dm[-(length(dm)-1L)], 0L)
        # folder
        index <- unlist(strsplit(var.name, "/"))
        if (length(index) <= 1L)
            folder <- gdsfile$root
        else
            folder <- index.gdsn(gdsfile, index=index[-length(index)])
        # compress
        if (is.null(compress))
            compress <- desp$compress

        name <- paste("~", index[length(index)], sep="")
        newnode <- add.gdsn(folder, name, val=NULL, storage=desp$storage,
            valdim=dm, compress=compress)

        # write data
        apply.gdsn(node, margin=length(dm)-1L, as.is="none", FUN=function(g) {
            append.gdsn(newnode, g)
        }, .useraw=TRUE)

        readmode.gdsn(newnode)
    } else
        warning("It is a vector.")

    invisible()
}



#######################################################################
# Optimize data by transposing
#

seqOptimize <- function(gdsfn, target=c("by.sample"),
    format.var=TRUE, cleanup=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfn) & is.vector(gdsfn))
    target <- match.arg(target)
    stopifnot(is.logical(format.var) || is.character(format.var))
    stopifnot(is.logical(cleanup))
    stopifnot(is.logical(verbose))

    gdsfile <- seqOpen(gdsfn, FALSE)
    on.exit({ seqClose(gdsfile) })

    if (target == "by.sample")
    {
        # genotype
        if (verbose) cat("Working on 'genotype' ...\n")
        .Transpose(gdsfile, "genotype/data", "~")

        # phase
        if (verbose) cat("Working on 'phase' ...\n")
        .Transpose(gdsfile, "phase/data", "~")

        # annotation - format
        if (identical(format.var, TRUE) || is.character(format.var))
        {
            n <- index.gdsn(gdsfile, "annotation/format", silent=TRUE)
            if (!is.null(n))
            {
                nm <- ls.gdsn(n)
                if (identical(format.var, TRUE))
                    format.var <- nm
                for (i in nm)
                {
                    if (i %in% format.var)
                    {
                        if (verbose)
                        {
                            cat("Working on 'annotation/format/", i,
                                "' ...\n", sep="")
                        }
                        .Transpose(gdsfile,
                            paste("annotation/format", i, "data", sep="/"), "~")
                    }
                }
            }
        }
    }

    if (cleanup)
    {
        on.exit()
        seqClose(gdsfile)
        cleanup.gds(gdsfn, verbose=verbose)
    }

    invisible()
}
