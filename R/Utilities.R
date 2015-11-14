#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Big Data Management of Sequencing-based Genetic Variants
#


#######################################################################
.Last.lib <- function(libpath)
{
    cl <- getOption("seqarray.parallel", FALSE)
    if (inherits(cl, "cluster"))
    {
        if (requireNamespace("parallel", quietly=TRUE))
            parallel::stopCluster(cl)
        options(seqarray.parallel=NULL)
    }
}



#######################################################################
# Get the file name of an example
#
seqExampleFileName <- function(type=c("gds", "vcf", "KG_Phase1"))
{
    type <- match.arg(type)
    switch(type,
        gds =
            system.file("extdata", "CEU_Exon.gds", package="SeqArray"),
        vcf =
            system.file("extdata", "CEU_Exon.vcf.gz", package="SeqArray"),
        KG_Phase1 =
            system.file("extdata", "1KG_phase1_release_v3_chr22.gds",
                package="SeqArray")
    )
}



#######################################################################
# Setup the parallel parameters in SeqArray
#
seqParallelSetup <- function(cluster=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.null(cluster) | is.logical(cluster) |
        is.numeric(cluster) | inherits(cluster, "cluster"))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.null(cluster) || identical(cluster, FALSE))
    {
        opt <- getOption("seqarray.parallel", NULL)
        if (inherits(opt, "cluster"))
        {
            .LoadParallelPackage()
            parallel::stopCluster(opt)
        }
        if (verbose)
            cat("Stop the computing cluster.\n")
        options(seqarray.parallel=cluster)
        return(invisible())
    }

    # Windows platform does not support forking, we have to setup a cluster
    if (.Platform$OS.type == "windows")
    {
        setup <- function(num.cores)
        {
            cl <- parallel::makeCluster(num.cores)
            parallel::clusterCall(cl, function() { library(SeqArray); NULL })
            cl
        }

        if (is.logical(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster)
            {
                .LoadParallelPackage()
                cl <- parallel::detectCores() - 1L
                if (cl <= 1L) cl <- 2L
                cluster <- setup(cl)
                if (verbose)
                {
                    cat("Enable the computing cluster with", cl,
                        "R processes.\n")
                }
            } else {
                if (verbose) cat("No computing cluster.\n")
            }
        } else if (is.numeric(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster > 1L)
            {
                .LoadParallelPackage()
                cl <- cluster
                cluster <- setup(cluster)
                if (verbose)
                {
                    cat("Enable the computing cluster with", cl,
                        "R processes.\n")
                }
            }
        }
    } else {
        # unix forking technique
        if (identical(cluster, TRUE))
        {
            .LoadParallelPackage()
            n <- parallel::detectCores() - 1L
            if (n <= 1L) n <- 2L
            if (verbose)
            {
                cat("Enable the computing cluster with", n,
                    "forked R processes.\n")
            }
        } else if (is.numeric(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster > 1L)
            {
                .LoadParallelPackage()
                if (verbose)
                {
                    cat("Enable the computing cluster with", cluster,
                        "forked R processes.\n")
                }
            }
        }
    }

    options(seqarray.parallel=cluster)
    invisible()
}



#######################################################################
# Export to a GDS file
#
seqExport <- function(gdsfile, out.fn, info.var=NULL, fmt.var=NULL,
    samp.var=NULL, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))
    stopifnot(is.null(samp.var) | is.character(samp.var))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    #######################################################################

    cp <- function(folder, sel, name, name2=NULL, show=verbose)
    {
        if (is.null(name2)) name2 <- name
        if (show)
            cat("Exporting \"", name2, "\" ...\n", sep="")
        src <- index.gdsn(gdsfile, name2)
        dst <- add.gdsn(folder, name, storage=src)
        put.attr.gdsn(dst, val=src)
        if (!is.null(sel))
        {
            dm <- objdesp.gdsn(src)$dim
            ss <- vector("list", length(dm))
            ss[[length(dm)]] <- sel
            for (i in seq_len(length(dm)-1L))
                ss[[i]] <- rep(TRUE, dm[i])
            sel <- ss
        }
        assign.gdsn(dst, src, append=FALSE, seldim=sel)
    }

    cp2 <- function(folder, samp.sel, var.sel, name, show=verbose)
    {
        if (show)
            cat("Exporting \"", name, "\" ...\n", sep="")

        dat1 <- index.gdsn(gdsfile, paste(name, "data", sep="/"))
        dat2 <- index.gdsn(gdsfile, paste(name, "~data", sep="/"), silent=TRUE)
        idx <- index.gdsn(gdsfile, paste(name, "@data", sep="/"))
        num <- read.gdsn(idx)
        stopifnot(length(num) == length(var.sel))

        dst1 <- add.gdsn(folder, "data", storage=dat1)
        put.attr.gdsn(dst1, val=dat1)
        if (!is.null(dat2))
        {
            dst2 <- add.gdsn(folder, "~data", storage=dat2)
            put.attr.gdsn(dst2, val=dat2)
        }
        dstidx <- add.gdsn(folder, "@data", storage=idx)
        put.attr.gdsn(dstidx, val=idx)

        dm <- objdesp.gdsn(dat1)$dim
        sel <- vector("list", length(dm))
        for (i in seq_len(length(dm)-2L))
            sel[[i]] <- rep(TRUE, dm[i])
        flag <- rep(var.sel, num)

        sel[[length(dm)]] <- flag
        sel[[length(dm)-1L]] <- samp.sel
        assign.gdsn(dst1, dat1, append=FALSE, seldim=sel)

        if (!is.null(dat2))
        {
            sel[[length(dm)]] <- samp.sel
            sel[[length(dm)-1L]] <- flag
            assign.gdsn(dst2, dat2, append=FALSE, seldim=sel)
        }

        assign.gdsn(dstidx, idx, append=FALSE, seldim=var.sel)
    }

    cp.info <- function(folder, sel, name, name2, show=verbose)
    {
        if (show)
            cat("Exporting \"", name2, "\" ...\n", sep="")
        src <- index.gdsn(gdsfile, name2)
        idx <- index.gdsn(gdsfile, .var_path(name2, "@"), silent=TRUE)

        dst <- add.gdsn(folder, name, storage=src)
        put.attr.gdsn(dst, val=src)
        if (!is.null(idx))
        {
            dstidx <- add.gdsn(folder, paste("@", name, sep=""), storage=idx)
            put.attr.gdsn(dstidx, val=idx)
        }

        dm <- objdesp.gdsn(src)$dim
        ss <- vector("list", length(dm))
        ss[[length(dm)]] <- sel
        for (i in seq_len(length(dm)-1L))
            ss[[i]] <- rep(TRUE, dm[i])
        assign.gdsn(dst, src, append=FALSE, seldim=ss)

        if (!is.null(idx))
            assign.gdsn(dstidx, idx, append=FALSE, seldim=sel)
    }

    #######################################################################

    # create the GDS file
    outfile <- createfn.gds(out.fn)
    on.exit({ closefn.gds(outfile) })

    # copy folders and attributes
    put.attr.gdsn(outfile$root, val=gdsfile$root)
    copyto.gdsn(outfile, index.gdsn(gdsfile, "description"))

    # the selection
    S <- seqGetFilter(gdsfile)

    ## sample.id, etc
    cp(outfile, S$sample.sel, "sample.id")
    cp(outfile, S$variant.sel, "variant.id")
    cp(outfile, S$variant.sel, "position")
    cp(outfile, S$variant.sel, "chromosome")
    cp(outfile, S$variant.sel, "allele")

    ## genotype
    node <- addfolder.gdsn(outfile, "genotype")
    put.attr.gdsn(node, val=index.gdsn(gdsfile, "genotype"))
    cp2(node, S$sample.sel, S$variant.sel, "genotype")

    if (prod(objdesp.gdsn(index.gdsn(gdsfile, "genotype/extra.index"))$dim) <= 0)
    {
        copyto.gdsn(node, index.gdsn(gdsfile, "genotype/extra.index"))
        copyto.gdsn(node, index.gdsn(gdsfile, "genotype/extra"))
    } else  # TODO
        stop("Not implemented in \"genotype/extra.index\".")

    ## annotation
    node <- addfolder.gdsn(outfile, "annotation")
    put.attr.gdsn(node, val=index.gdsn(gdsfile, "annotation"))
    node.info <- NULL
    node.fmt <- NULL
    lst <- ls.gdsn(index.gdsn(gdsfile, "annotation"))
    for (nm in lst)
    {
        if (nm == "info")
        {
            if (is.null(node.info))
            {
                node.info <- addfolder.gdsn(node, "info")
                put.attr.gdsn(node.info,
                    val=index.gdsn(gdsfile, "annotation/info"))
            }
            lst.info <- ls.gdsn(index.gdsn(gdsfile, "annotation/info"))
            for (nm2 in lst.info)
            {
                cp.info(node.info, S$variant.sel, nm2,
                    paste("annotation", "info", nm2, sep="/"))
            }
        } else if (nm == "format")
        {
            if (is.null(node.fmt))
            {
                node.fmt <- addfolder.gdsn(node, "format")
                put.attr.gdsn(node.fmt,
                    val=index.gdsn(gdsfile, "annotation/format"))
            }
            lst.fmt <- ls.gdsn(index.gdsn(gdsfile, "annotation/format"))
            for (nm2 in lst.fmt)
            {
                s <- paste("annotation", "format", nm2, sep="/")
                n2 <- addfolder.gdsn(node.fmt, nm2)
                put.attr.gdsn(n2, val=index.gdsn(gdsfile, s))
                cp2(n2, S$sample.sel, S$variant.sel, s)
            }
        } else {
            cp(node, S$variant.sel, nm, paste("annotation", nm, sep="/"))
        }
    }

    ## sample.annotation
    node <- addfolder.gdsn(outfile, "sample.annotation")
    n2 <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
    if (!is.null(n2))
    {
        put.attr.gdsn(node, val=n2)
        lst <- ls.gdsn(n2)
        for (nm in lst)
            cp(node, S$sample.sel, nm, paste("sample.annotation", nm, sep="/"))
    }

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Merge multiple GDS files
#
seqMerge <- function(gds.fn, out.fn, storage.option=seqStorage.Option(),
    info.var=NULL, fmt.var=NULL, samp.var=NULL, optimize=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gds.fn))
    if (length(gds.fn) <= 1L)
        stop("'gds.fn' should have more than one files.")
    stopifnot(is.character(out.fn), length(out.fn)==1L)

    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))
    stopifnot(is.null(samp.var) | is.character(samp.var))

    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat(sprintf("Merging %d GDS files:\n", length(gds.fn)))
    }

    # open all GDS files
    flist <- vector("list", length(gds.fn))
    on.exit({ for (f in flist) seqClose(f) })
    for (i in seq_along(gds.fn))
        flist[[i]] <- seqOpen(gds.fn[i])

    # samples
    samp.id <- samp2.id <- seqGetData(flist[[1L]], "sample.id")
    for (f in flist[-1L])
    {
        s <- seqGetData(f, "sample.id")
        samp.id <- unique(c(samp.id, s))
        samp2.id <- intersect(samp2.id, s)
    }

    if (verbose)
    {
        cat(sprintf("\t%d sample(s) in total, %d sample(s) in common\n",
            length(samp.id), length(samp2.id)))
    }

    # variants
    variant <- function(f) {
        paste(seqGetData(f, "chromosome"), seqGetData(f, "position"), sep="-")
    }
    variant.id <- variant2.id <- variant(flist[[1L]])
    for (f in flist[-1L])
    {
        s <- variant(f)
        variant.id <- unique(c(variant.id, s))
        variant2.id <- intersect(variant2.id, s)
    }

    if (verbose)
    {
        cat(sprintf("\t%d variant(s) in total, %d variant(s) in common\n",
            length(variant.id), length(variant2.id)))
    }

    # common samples
    if (length(samp2.id) > 0L)
    {
        if (length(variant2.id) > 0L)
            stop("There are overlapping on both samples and variants.")
    } else {
        stop("not implemented yet!")
    }


    ## create a GDS file
    gfile <- createfn.gds(out.fn)
    on.exit({ if (!is.null(gfile)) closefn.gds(gfile) }, add=TRUE)

    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(gfile, "description")
    .MergeNodeAttr(n, flist, "description")

    v <- vector("list", length(gds.fn))
    for (i in seq_along(flist))
        v[[i]] <- seqSummary(flist[[i]], check="none", verbose=FALSE)$reference
    reference <- as.character(unique(unlist(v)))
    .AddVar(storage.option, n, "reference", reference, closezip=TRUE,
        visible=FALSE)
    
    alt <- NULL
    for (i in seq_along(flist))
    {
        z <- index.gdsn(flist[[i]], "description/vcf.alt", silent=TRUE)
        if (!is.null(z))
            alt <- rbind(alt, read.gdsn(z))
    }
    if (!is.null(alt))
    {
        .AddVar(storage.option, n, "vcf.alt", .UniqueDataFrame(alt),
            closezip=TRUE, visible=FALSE)
    }
        
    contig <- NULL
    for (i in seq_along(flist))
    {
        z <- index.gdsn(flist[[i]], "description/vcf.contig", silent=TRUE)
        if (!is.null(z))
            contig <- rbind(contig, read.gdsn(z))
    }
    if (!is.null(contig))
    {
        .AddVar(storage.option, n, "vcf.contig", .UniqueDataFrame(contig),
            closezip=TRUE, visible=FALSE)
    }

    header <- NULL
    for (i in seq_along(flist))
    {
        z <- index.gdsn(flist[[i]], "description/vcf.header", silent=TRUE)
        if (!is.null(z))
            header <- rbind(header, read.gdsn(z))
    }
    if (!is.null(header))
    {
        .AddVar(storage.option, n, "vcf.header", .UniqueDataFrame(header),
            closezip=TRUE, visible=FALSE)
    }
        

    ## add sample.id
    if (verbose) cat("sample.id\n")
    .AddVar(storage.option, gfile, "sample.id", samp.id, closezip=TRUE)

    ## add variant.id
    if (verbose) cat("variant.id\n")
    .AddVar(storage.option, gfile, "variant.id", seq_along(variant.id),
        storage="int32", closezip=TRUE)

    if (length(samp2.id) > 0L)
    {
        nSamp <- length(samp.id)

        ## add position, chromsome, allele
        # TODO: need to check whether position can be stored in 'int32'
        if (verbose) cat("position\n")
        n <- .AddVar(storage.option, gfile, "position", storage="int32")
        .append_gds(n, flist, "position")

        if (verbose) cat("chromosome\n")
        n <- .AddVar(storage.option, gfile, "chromosome", storage="string")
        .append_gds(n, flist, "chromosome")

        if (verbose) cat("allele\n")
        n <- .AddVar(storage.option, gfile, "allele", storage="string")
        .append_gds(n, flist, "allele")


        ## add a folder for genotypes
        if (verbose) cat("genotype, phase [")
        varGeno <- addfolder.gdsn(gfile, "genotype")
        .MergeNodeAttr(varGeno, flist, "genotype")

        ploidy <- seqSummary(flist[[1L]], check="none", verbose=FALSE)$ploidy
        for (i in seq_along(gds.fn))
        {
            if (!identical(ploidy, seqSummary(flist[[1L]], check="none",
                    verbose=FALSE)$ploidy))
                stop("ploidy should be the same!")
        }
        n1 <- .AddVar(storage.option, varGeno, "data", storage="bit2",
            valdim=c(ploidy, nSamp, 0L))

        ## add phase folder
        varPhase <- addfolder.gdsn(gfile, "phase")
        n2 <- .AddVar(storage.option, varPhase, "data", storage="bit1",
            valdim=c(header$ploidy-1L, nSamp, 0L))

        for (i in seq_along(flist))
        {
            if (verbose)
                cat(ifelse(i > 1L, ", ", ""), i, sep="")
            sid <- seqGetData(flist[[i]], "sample.id")
            n3 <- index.gdsn(flist[[i]], "genotype/data")
            n4 <- index.gdsn(flist[[i]], "phase/data")
            if (identical(sid, samp.id))
            {
                append.gdsn(n1, n3)
                append.gdsn(n2, n4)
            } else {
                k <- match(samp.id, sid)
                s <- vector("list", 3L)
                s[2L] <- k
                assign.gdsn(n1, n3, seldim=s, append=TRUE,
                    .value=NA_integer_, .substitute=3L)
                s <- vector("list", length(objdesp.gdsn(n4)$dim))
                s[length(s)-1L] <- k
                assign.gdsn(n2, n4, seldim=s, append=TRUE,
                    .value=NA_integer_, .substitute=0)
            }
        }

        readmode.gdsn(n1)
        readmode.gdsn(n2)

        n <- .AddVar(storage.option, varGeno, "@data", storage="uint8",
            visible=FALSE)
        .append_gds(n, flist, "genotype/@data")

        # TODO
        n <- .AddVar(storage.option, varGeno, "extra.index", storage="int32",
            valdim=c(3L,0L))
        put.attr.gdsn(n, "R.colnames",
            c("sample.index", "variant.index", "length"))
        .AddVar(storage.option, varGeno, "extra", storage="int16", closezip=TRUE)

        n <- .AddVar(storage.option, varPhase, "extra.index",
            storage="int32", valdim=c(3L,0L))
        put.attr.gdsn(n, "R.colnames",
            c("sample.index", "variant.index", "length"))
        .AddVar(storage.option, varPhase, "extra", storage="bit1", closezip=TRUE)

        # sync file
        sync.gds(gfile)
        if (verbose) cat("]\n")


        ## add annotation folder
        varAnnot <- addfolder.gdsn(gfile, "annotation")

        # add id
        n <- .AddVar(storage.option, varAnnot, "id", storage="string")
        .append_gds(n, flist, "annotation/id")
        # add qual
        n <- .AddVar(storage.option, varAnnot, "qual", storage="float")
        .append_gds(n, flist, "annotation/qual")
        # add filter
        n <- .AddVar(storage.option, varAnnot, "filter", storage="int32")
        .append_gds(n, flist, "annotation/filter")


        ## VCF INFO
        if (verbose) cat("annotation/info\n")
        varInfo <- addfolder.gdsn(varAnnot, "info")
        varnm <- NULL
        for (i in seq_along(flist))
        {
            n <- index.gdsn(flist[[i]], "annotation/info", silent=TRUE)
            if (!is.null(n))
                varnm <- unique(c(varnm, ls.gdsn(n)))
        }
        if (!is.null(info.var))
        {
            s <- setdiff(info.var, varnm)
            if (length(s) > 0L)
            {
                warning("No INFO variable(s): ", paste(s, collapse=", "),
                    immediate.=TRUE)
            }
            varnm <- intersect(varnm, info.var)
        }
        # for-loop
        for (i in seq_along(varnm))
        {
            idx <- 0L
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "annotation/info", silent=TRUE)
                if (!is.null(n))
                {
                    if (varnm[i] %in% ls.gdsn(n))
                    {
                        idx <- j
                        break
                    }
                }
            }
            if (idx < 1L) stop("internal error, info field.")

            need <- FALSE
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "annotation/info", silent=TRUE)
                if (is.null(n))
                {
                    need <- TRUE
                    break
                }
            }

            n <- index.gdsn(flist[[idx]], paste0("annotation/info/", varnm[i]))
            n1 <- index.gdsn(flist[[idx]],
                paste0("annotation/info/@", varnm[i]), silent=TRUE)

            dp <- objdesp.gdsn(n)
            dp$dim[length(dp$dim)] <- 0L
            n2 <- .AddVar(storage.option, varInfo, varnm[i], storage=dp$storage,
                valdim=dp$dim)
            .MergeNodeAttr(n2, flist[idx], paste0("annotation/info/", varnm[i]))

            n3 <- n1
            if (!is.null(n1) | need)
            {
                n3 <- .AddVar(storage.option, varInfo,
                    paste("@", varnm[i], sep=""), storage="int32",
                    visible=FALSE)
            }

            for (j in seq_along(flist))
            {
                f <- flist[[j]]
                n4 <- index.gdsn(f, paste0("annotation/info/", varnm[i]),
                    silent=TRUE)
                n5 <- index.gdsn(f, paste0("annotation/info/@", varnm[i]),
                    silent=TRUE)

                if (!is.null(n4))
                {
                    append.gdsn(n2, n4)
                    if (!is.null(n5))
                    {
                        append.gdsn(n3, n5)
                    } else {
                        if (!is.null(n3))
                        {
                            cnt <- objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
                            .repeat_gds(n3, 1L, cnt)
                        }
                    }
                } else {
                    if (is.null(n3))
                        stop(paste0("annotation/info/@", varnm[i]), " error.")
                    cnt <- objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
                    .repeat_gds(n3, 0L, cnt)
                }
            }

            readmode.gdsn(n2)
            if (!is.null(n3)) readmode.gdsn(n3)
        }


        ## VCF FORMAT
        if (verbose) cat("annotation/format\n")
        varFormat <- addfolder.gdsn(varAnnot, "format")
        varnm <- NULL
        for (i in seq_along(flist))
        {
            n <- index.gdsn(flist[[i]], "annotation/format", silent=TRUE)
            if (!is.null(n))
                varnm <- unique(c(varnm, ls.gdsn(n)))
        }
        if (!is.null(fmt.var))
        {
            s <- setdiff(fmt.var, varnm)
            if (length(s) > 0L)
            {
                warning("No FORMAT variable(s): ", paste(s, collapse=", "),
                    immediate.=TRUE)
            }
            varnm <- intersect(varnm, fmt.var)
        }
        # for-loop
        for (i in seq_along(varnm))
        {
            if (verbose) cat("\t", varnm[i], " [", sep="")
            idx <- 0L
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "annotation/format", silent=TRUE)
                if (!is.null(n))
                {
                    if (varnm[i] %in% ls.gdsn(n))
                    {
                        idx <- j
                        break
                    }
                }
            }
            if (idx < 1L) stop("internal error, format field.")

            n <- index.gdsn(flist[[idx]],
                paste0("annotation/format/", varnm[i]))
            n1 <- index.gdsn(flist[[idx]],
                paste0("annotation/format/", varnm[i], "/data"))
            n2 <- index.gdsn(flist[[idx]],
                paste0("annotation/format/", varnm[i], "/@data"))

            n3 <- addfolder.gdsn(varFormat, varnm[i])
            .MergeNodeAttr(n3, flist[idx],
                paste0("annotation/format/", varnm[i]))
            dp <- objdesp.gdsn(n1)
            dp$dim[length(dp$dim)] <- 0L
            n4 <- .AddVar(storage.option, n3, "data", storage=dp$storage,
                valdim=dp$dim)
            n5 <- .AddVar(storage.option, n3, "@data", storage="int32",
                visible=FALSE)

            for (j in seq_along(flist))
            {
                if (verbose)
                    cat(ifelse(j > 1L, ", ", ""), j, sep="")
                f <- flist[[j]]
                n6 <- index.gdsn(f, paste0("annotation/format/", varnm[i]),
                    silent=TRUE)
                if (!is.null(n6))
                {
                    sid <- seqGetData(f, "sample.id")
                    n7 <- index.gdsn(n6, "data")
                    if (identical(sid, samp.id))
                    {
                        append.gdsn(n4, n7)
                    } else {
                        d <- objdesp.gdsn(n7)$dim
                        if (d[length(d)-1L] != length(sid))
                            stop("Invalid FORMAT variable.")
                        s <- vector("list", length(d))
                        s[length(s)-1L] <- match(samp.id, sid)
                        assign.gdsn(n4, n7, seldim=s, append=TRUE)
                    }
                    append.gdsn(n5, index.gdsn(n6, "@data"))
                } else {
                    cnt <- objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
                    .repeat_gds(n5, 0L, cnt)
                }
            }

            readmode.gdsn(n4)
            readmode.gdsn(n5)
            if (verbose) cat("]\n")
        }


        ## annotation folder
        if (verbose) cat("sample.annotation\n")
        varSamp <- addfolder.gdsn(gfile, "sample.annotation")
        varnm <- NULL
        for (i in seq_along(flist))
        {
            n <- index.gdsn(flist[[i]], "sample.annotation", silent=TRUE)
            if (!is.null(n))
                varnm <- unique(c(varnm, ls.gdsn(n)))
        }
        if (!is.null(samp.var))
        {
            s <- setdiff(samp.var, varnm)
            if (length(s) > 0L)
            {
                warning("No SAMPLE variable(s): ", paste(s, collapse=", "),
                    immediate.=TRUE)
            }
            varnm <- intersect(varnm, samp.var)
        }
        # for-loop
        for (i in seq_along(varnm))
        {
            idx <- 0L
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "annotation/format", silent=TRUE)
                if (!is.null(n))
                {
                    if (varnm[i] %in% ls.gdsn(n))
                    {
                        idx <- j
                        break
                    }
                }
            }
            if (idx < 1L) stop("internal error, format field.")

            copyto.gdsn(varSamp, index.gdsn(flist[[idx]],
                paste0("sample.annotation/", varnm[i])))
        }

    } else {
        stop("not implemented yet!")
    }

    for (f in flist) seqClose(f)
    flist <- list()
    closefn.gds(gfile)
    on.exit()

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }


    ##################################################
    # optimize access efficiency

    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.fn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # return
    invisible(normalizePath(out.fn))
}



#######################################################################
# Storage options for the SeqArray GDS file
#
seqStorage.Option <- function(compression=c("ZIP_RA", "ZIP_RA.max", "LZ4_RA",
    "LZ4_RA.max", "none"), float.mode="float32",
    geno.compress=NULL, info.compress=NULL, format.compress=NULL,
    index.compress=NULL, ...)
{
    # check
    compression <- match.arg(compression)
    if (compression == "none") compression <- ""

    stopifnot(is.character(float.mode), length(float.mode) > 0L)
    if (!is.null(geno.compress))
        stopifnot(is.character(geno.compress), length(geno.compress)==1L)
    if (!is.null(info.compress))
        stopifnot(is.character(info.compress), length(info.compress)==1L)
    if (!is.null(format.compress))
        stopifnot(is.character(format.compress), length(format.compress)==1L)
    if (!is.null(index.compress))
        stopifnot(is.character(index.compress), length(index.compress)==1L)

    rv <- list(compression = compression, float.mode = float.mode,
        geno.compress   = ifelse(is.null(geno.compress), compression,
            geno.compress),
        info.compress   = ifelse(is.null(info.compress), compression,
            info.compress),
        format.compress = ifelse(is.null(format.compress),
            ifelse(compression=="", "", paste(compression, ":8M", sep="")),
            format.compress),
        index.compress = ifelse(is.null(index.compress), compression,
            index.compress),
        ...)
    class(rv) <- "SeqGDSStorageClass"
    return(rv)
}



#######################################################################
# Apply functions in parallel
#
seqParallel <- function(cl=getOption("seqarray.parallel", FALSE),
    gdsfile, FUN, split=c("by.variant", "by.sample", "none"),
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

    if (verbose) cat("Delete INFO variable(s):")
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

    if (verbose) cat("Delete FORMAT variable(s):")
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
