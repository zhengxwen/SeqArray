#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Data Management of Large-scale Whole-Genome Sequence Variant Calls
#


#######################################################################
# Export to a GDS file
#
seqExport <- function(gdsfile, out.fn, info.var=NULL, fmt.var=NULL,
    samp.var=NULL, optimize=TRUE, digest=TRUE, verbose=TRUE, verbose.clean=NA)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.fn), length(out.fn)==1L)

    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))
    stopifnot(is.null(samp.var) | is.character(samp.var))

    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.clean), length(verbose.clean)==1L)
    if (is.na(verbose.clean)) verbose.clean <- verbose

    #######################################################################

    cp <- function(folder, sel, name, name2=NULL, show=verbose)
    {
        if (is.null(name2)) name2 <- name
        if (show)
        {
            if (name %in% c("sample.id", "variant.id"))
            {
                ss <- .seldim(gdsfile)
                n <- ifelse(name=="sample.id", ss[2L], ss[3L])
                cat("    ", name2, " (", .pretty(n), ")", sep="")
            } else
                cat("   ", name2)
        }
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
        readmode.gdsn(dst)
        .DigestCode(dst, digest, show)
    }

    cp2 <- function(folder, samp.sel, var.sel, name, show=verbose)
    {
        if (show) cat("   ", name)

        dat1 <- index.gdsn(gdsfile, paste(name, "data", sep="/"))
        dat2 <- index.gdsn(gdsfile, paste(name, "~data", sep="/"), silent=TRUE)

        idx <- index.gdsn(gdsfile, paste(name, "@data", sep="/"), silent=TRUE)
        if (is.null(idx))
            idx <- index.gdsn(gdsfile, "genotype/@data")
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
        flag <- rep(var.sel, num)

        sel[[length(dm)]] <- flag
        sel[[length(dm)-1L]] <- samp.sel
        assign.gdsn(dst1, dat1, append=FALSE, seldim=sel)
        readmode.gdsn(dst1)
        .DigestCode(dst1, digest, show)

        if (!is.null(dat2))
        {
            sel[[length(dm)]] <- samp.sel
            sel[[length(dm)-1L]] <- flag
            assign.gdsn(dst2, dat2, append=FALSE, seldim=sel)
            readmode.gdsn(dst2)
            .DigestCode(dst2, digest, show)
        }

        assign.gdsn(dstidx, idx, append=FALSE, seldim=var.sel)
        readmode.gdsn(dstidx)
        .DigestCode(dstidx, digest, FALSE)
    }

    cp.phase <- function(folder, samp.sel, var.sel, show=verbose)
    {
        if (show) cat("    phase")
        dat1 <- index.gdsn(gdsfile, "phase/data")
        dat2 <- index.gdsn(gdsfile, "phase/~data", silent=TRUE)

        dst1 <- add.gdsn(folder, "data", storage=dat1)
        put.attr.gdsn(dst1, val=dat1)
        if (!is.null(dat2))
        {
            dst2 <- add.gdsn(folder, "~data", storage=dat2)
            put.attr.gdsn(dst2, val=dat2)
        }

        dm <- objdesp.gdsn(dat1)$dim
        sel <- vector("list", length(dm))
        sel[[length(dm)]] <- var.sel
        sel[[length(dm)-1L]] <- samp.sel
        assign.gdsn(dst1, dat1, append=FALSE, seldim=sel)
        readmode.gdsn(dst1)
        .DigestCode(dst1, digest, show)

        if (!is.null(dat2))
        {
            sel[[length(dm)]] <- samp.sel
            sel[[length(dm)-1L]] <- var.sel
            assign.gdsn(dst2, dat2, append=FALSE, seldim=sel)
            readmode.gdsn(dst2)
            .DigestCode(dst2, digest, show)
        }
    }

    cp.info <- function(folder, sel, name, name2, show=verbose)
    {
        if (show) cat("   ", name2)
        src <- index.gdsn(gdsfile, name2)
        idx <- index.gdsn(gdsfile, .var_path(name2, "@"), silent=TRUE)

        while (grepl("/", name, fixed=TRUE))
        {
            ss <- unlist(strsplit(name, "/", fixed=TRUE))
            name <- paste(ss[-1L], collapse="/")
            if (is.null(index.gdsn(folder, ss[1L], silent=TRUE)))
                folder <- addfolder.gdsn(folder, ss[1L])
            else
                folder <- index.gdsn(folder, ss[1L])
        }

        dst <- add.gdsn(folder, name, storage=src)
        put.attr.gdsn(dst, val=src)

        if (is.null(idx))
        {
            dm <- objdesp.gdsn(src)$dim
            ss <- vector("list", length(dm))
            ss[[length(dm)]] <- sel
            assign.gdsn(dst, src, append=FALSE, seldim=ss)
            readmode.gdsn(dst)
            .DigestCode(dst, digest, show)
        } else {
            dstidx <- add.gdsn(folder, paste("@", name, sep=""), storage=idx)
            put.attr.gdsn(dstidx, val=idx)
            dm <- objdesp.gdsn(src)$dim
            ss <- vector("list", length(dm))
            ss[[length(dm)]] <- .Call(SEQ_SelectFlag, sel, read.gdsn(idx))
            assign.gdsn(dst, src, append=FALSE, seldim=ss)
            readmode.gdsn(dst)
            .DigestCode(dst, digest, show)
            assign.gdsn(dstidx, idx, append=FALSE, seldim=sel)
            readmode.gdsn(dstidx)
            .DigestCode(dstidx, digest, FALSE)
        }
    }

    #######################################################################

    # create the GDS file
    outfile <- createfn.gds(out.fn)
    on.exit({ closefn.gds(outfile) })

    if (verbose)
        cat("Export to ", sQuote(out.fn), ":\n", sep="")

    # copy folders and attributes
    put.attr.gdsn(outfile$root, val=gdsfile$root)
    copyto.gdsn(outfile, index.gdsn(gdsfile, "description"))

    # the selection
    S <- seqGetFilter(gdsfile)
    nsamp <- sum(S$sample.sel)
    nsnp <- sum(S$variant.sel)

    ## sample.id, etc
    cp(outfile, S$sample.sel, "sample.id")
    cp(outfile, S$variant.sel, "variant.id")
    cp(outfile, S$variant.sel, "position")
    cp(outfile, S$variant.sel, "chromosome")
    cp(outfile, S$variant.sel, "allele")
    sync.gds(outfile)

    ## genotype
    node <- addfolder.gdsn(outfile, "genotype")
    put.attr.gdsn(node, val=index.gdsn(gdsfile, "genotype"))
    if (exist.gdsn(gdsfile, "genotype/data") && nsamp && nsnp)
    {
        cp2(node, S$sample.sel, S$variant.sel, "genotype")
        if (prod(objdesp.gdsn(index.gdsn(gdsfile, "genotype/extra.index"))$dim) <= 0)
        {
            copyto.gdsn(node, index.gdsn(gdsfile, "genotype/extra.index"))
            copyto.gdsn(node, index.gdsn(gdsfile, "genotype/extra"))
        } else  # TODO
            stop("Not implemented in 'genotype/extra.index', please contact the author.")
        sync.gds(outfile)
    } else {
        if (verbose) cat("    genotype\n")
    }

    ## phase
    node <- addfolder.gdsn(outfile, "phase")
    put.attr.gdsn(node, val=index.gdsn(gdsfile, "phase"))
    if (exist.gdsn(gdsfile, "phase/data") && nsamp && nsnp)
    {
        cp.phase(node, S$sample.sel, S$variant.sel)
        if (prod(objdesp.gdsn(index.gdsn(gdsfile, "phase/extra.index"))$dim) <= 0)
        {
            copyto.gdsn(node, index.gdsn(gdsfile, "phase/extra.index"))
            copyto.gdsn(node, index.gdsn(gdsfile, "phase/extra"))
        } else  # TODO
            stop("Not implemented in 'phase/extra.index', please contact the author.")
        sync.gds(outfile)
    } else {
        if (verbose) cat("    phase\n")
    }

    ## annotation
    node <- addfolder.gdsn(outfile, "annotation")
    put.attr.gdsn(node, val=index.gdsn(gdsfile, "annotation"))
    node.info <- NULL
    node.fmt <- NULL
    lst <- ls.gdsn(index.gdsn(gdsfile, "annotation"), recursive=FALSE)
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
            lst.info <- ls.gdsn(index.gdsn(gdsfile, "annotation/info"),
                recursive=TRUE, include.dirs=FALSE)
            if (!is.null(info.var))
            {
                s <- setdiff(info.var, lst.info)
                if (length(s) > 0L)
                    warning("No INFO variable: ", paste(s, collapse=","))
                lst.info <- intersect(lst.info, info.var)
            }
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
            if (nsamp > 0L)
            {
                lst.fmt <- ls.gdsn(index.gdsn(gdsfile, "annotation/format"))
                if (!is.null(fmt.var))
                {
                    s <- setdiff(fmt.var, lst.fmt)
                    if (length(s) > 0L)
                        warning("No FORMAT variable: ", paste(s, collapse=","))
                    lst.fmt <- intersect(lst.fmt, fmt.var)
                }
            } else {
                lst.fmt <- character()
            }
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

    sync.gds(outfile)

    ## sample.annotation
    node <- addfolder.gdsn(outfile, "sample.annotation")
    n2 <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
    if (!is.null(n2) && nsamp>0L)
    {
        put.attr.gdsn(node, val=n2)
        lst <- ls.gdsn(n2)
        if (!is.null(samp.var))
        {
            s <- setdiff(samp.var, lst)
            if (length(s) > 0L)
                warning("No sample variable: ", paste(s, collapse=","))
            lst <- intersect(lst, samp.var)
        }
        for (nm in lst)
            cp(node, S$sample.sel, nm, paste("sample.annotation", nm, sep="/"))
    }

    # add RLE of chromosome
    .optim_chrom(outfile)

    # close the file
    on.exit()
    closefn.gds(outfile)
    if (verbose) cat("Done.\n")

    # optimize access efficiency
    if (optimize)
    {
        if (verbose.clean)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.fn, verbose=verbose.clean)
    }

    # output
    invisible(normalizePath(out.fn))
}


#######################################################################
# Recompress the GDS file
#
seqRecompress <- function(gds.fn, compress=c("ZIP", "LZ4", "LZMA", "Ultra",
    "UltraMax", "none"), exclude=character(), optimize=TRUE, verbose=TRUE)
{
    stopifnot(is.character(gds.fn), length(gds.fn)==1L)
    stopifnot(is.character(exclude))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    compress <- match.arg(compress)
    if (compress == "ZIP")
    {
        node_compress <- fmt_compress <- idx_compress <- "ZIP_RA"
    } else if (compress == "LZ4")
    {
        node_compress <- fmt_compress <- idx_compress <- "LZ4_RA"
    } else if (compress == "LZMA")
    {
        node_compress <- fmt_compress <- idx_compress <- "LZMA_RA"
    } else if (compress == "Ultra")
    {
        node_compress <- "LZMA_RA.ultra:4M"
        fmt_compress <- "LZMA_RA.ultra:8M"
        idx_compress <- "LZMA.max"
    } else if (compress == "UltraMax")
    {
        node_compress <- "LZMA_RA.ultra_max:8M"
        fmt_compress <- "LZMA_RA.ultra_max:8M"
        idx_compress <- "LZMA.max"
    } else if (compress == "none")
    {
        node_compress <- fmt_compress <- idx_compress <- ""
    } else {
        stop("Not implemented.")
    }

    # open the file
    f <- openfn.gds(gds.fn, readonly=FALSE)
    on.exit({ closefn.gds(f) })
    if (verbose)
        cat("Open", sQuote(gds.fn), "...\n")

    nm_lst <- ls.gdsn(f, include.hidden=TRUE, recursive=TRUE)
    nm_lst <- nm_lst[!grepl("^description", nm_lst)]
    nm_lst <- setdiff(nm_lst, c("@chrom_rle_val", "@chrom_rle_len"))
    nm_lst <- setdiff(nm_lst, exclude)
    for (nm in nm_lst)
    {
        n <- index.gdsn(f, nm)
        dp <- objdesp.gdsn(n)
        if (dp$is.array & prod(dp$dim)>0L)
        {
            if (verbose) cat("   ", nm)
            # compress
            sz <- dp$size
            if (grepl("@", basename(nm), fixed=TRUE))
            {
                compression.gdsn(n, idx_compress)
            } else if (grepl("^annotation/format/", nm))
            {
                compression.gdsn(n, fmt_compress)
            } else {
                compression.gdsn(n, node_compress)
            }
            sz2 <- objdesp.gdsn(n)$size
            if (verbose)
            {
                v <- (1 - sz2/sz) * 100
                if (v >= 0)
                    cat(sprintf("\t(deflated %.1f%%)", v))
                else
                    cat(sprintf("\t(inflated %.1f%%)", -v))
            }
            # digest
            if (!is.null(get.attr.gdsn(n)$md5))
            {
                digest.gdsn(n, algo="md5", action="add")
                if (verbose)
                    cat("  MD5:", get.attr.gdsn(n)$md5, "\n", sep="")
            } else if (verbose)
                cat("\n")
        }
    }

    # close the file
    on.exit()
    closefn.gds(f)
    if (verbose) cat("Done.\n")

    # optimize access efficiency
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(gds.fn, verbose=verbose)
    }

    invisible()
}
