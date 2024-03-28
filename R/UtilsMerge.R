#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Data Management of Large-scale Whole-Genome Sequence Variant Calls
#


#######################################################################
# Merge multiple GDS files
#

# check whether the chromosomes are overlapping or not (for faster checking)
.is_variant_overlap_bychr <- function(flist)
{
    lst <- lapply(flist, function(f) unique(seqGetData(f, "chromosome")))
    for (i in seq_len(length(lst))[-1L])
    {
        for (j in seq_len(i-1L))
        {
            s <- intersect(lst[[i]], lst[[j]])
            if (length(s)) return(TRUE)
        }
    }
    FALSE
}

# append values in the nodes
.append_node_variant <- function(varFolder, varnm, storage, storage.option,
    flist, nVariant, digest, verbose)
{
    if (verbose)
        cat("    ", varnm, " ", sep="")
    n <- .AddVar(storage.option, varFolder, basename(varnm), storage=storage)
    .append_gds(n, flist, varnm, verbose)
    .DigestCode(n, digest, verbose)
    if (nVariant != objdesp.gdsn(n)$dim)
        stop(sprintf("Invalid number of variants in '%s'.", varnm))
    n
}


# merge different variants for INFO
.append_info_variant <- function(flist, varnm, storage.option, varInfo,
    digest, verbose)
{
    h <- "."
    if (.crayon()) h <- crayon::blurred(h)

    for (i in seq_along(varnm))
    {
        if (verbose) cat("        ", varnm[i], " ", sep="")
        idx <- 0L
        for (j in seq_along(flist))
        {
            n <- index.gdsn(flist[[j]], "annotation/info", silent=TRUE)
            if (!is.null(n))
            {
                if (varnm[i] %in% ls.gdsn(n, recursive=TRUE, include.dirs=FALSE))
                {
                    idx <- j
                    break
                }
            }
        }
        if (is.na(idx) || (idx<1L)) stop("internal error, info field.")

        need <- FALSE
        for (j in seq_along(flist))
        {
            n <- index.gdsn(flist[[j]],
                paste("annotation/info/", varnm[i], sep=""), silent=TRUE)
            if (is.null(n))
            {
                need <- TRUE
                break
            }
        }

        s <- paste0("annotation/info/", varnm[i])
        n <- index.gdsn(flist[[idx]], s)
        n1 <- index.gdsn(flist[[idx]], .var_path(s, "@"), silent=TRUE)

        dp <- objdesp.gdsn(n)
        dp$dim[length(dp$dim)] <- 0L
        n2 <- .AddVar(storage.option, varInfo, varnm[i], storage=dp$storage,
            valdim=dp$dim)
        .MergeNodeAttr(n2, flist[idx], paste0("annotation/info/", varnm[i]))

        n3 <- n1
        if (!is.null(n1) | need)
        {
            n3 <- .AddVar(storage.option, varInfo, paste0("@", varnm[i]),
                storage="int32", visible=FALSE)
        }

        for (j in seq_along(flist))
        {
            f <- flist[[j]]
            s <- paste0("annotation/info/", varnm[i])
            n4 <- index.gdsn(f, s, silent=TRUE)
            n5 <- index.gdsn(f, .var_path(s, "@"), silent=TRUE)

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
                        .append_rep_gds(n3, as.raw(1L), cnt)
                    }
                }
            } else {
                if (is.null(n3))
                {
                    stop(.var_path(paste0("annotation/info/", varnm[i]), "@"),
                        " error.")
                }
                cnt <- objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
                .append_rep_gds(n3, as.raw(0L), cnt)
            }
            if (verbose)
            {
                cat(h)
                flush(stdout()); flush.console()
            }
        }

        readmode.gdsn(n2)
        .DigestCode(n2, digest, verbose)

        if (!is.null(n3))
        {
            readmode.gdsn(n3)
            .DigestCode(n3, digest, FALSE)
        }
    }
    invisible()
}


# merge for FORMAT
.append_format <- function(flist, varnm, samp.id, gfile, diffVar, nVariant,
    varidx, storage.option, digest, verbose)
{
    if (verbose)
        .cat("    annotation/format (", paste(varnm, collapse=","), ")")
    nSamp <- length(samp.id)
    varFormat <- index.gdsn(gfile, "annotation/format")

    for (i in seq_along(varnm))
    {
        if (verbose)
        {
            cat("        ", varnm[i], " ", sep="")
            cat(if (.crayon()) crayon::blurred("[") else "[")
        }
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
            paste("annotation/format/", varnm[i], sep=""))
        dp <- objdesp.gdsn(n1)
        n4 <- .AddVar(storage.option, n3, "data", storage=dp$storage,
            valdim=c(nSamp, 0L))
        n5 <- .AddVar(storage.option, n3, "@data", storage="int32",
            visible=FALSE)

        if (diffVar)
        {
            # merge different variants
            for (j in seq_along(flist))
            {
                if (verbose)
                {
                    h <- paste0(ifelse(j > 1L, ",", ""), j)
                    if (.crayon()) h <- crayon::blurred(h)
                    cat(h)
                    flush(stdout()); flush.console()
                }
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
                    .append_rep_gds(n5, as.raw(0L), cnt)
                }
            }

        } else {
            # merge different samples
            .Call(SEQ_MergeFormat, nVariant, varidx, flist,
                paste("annotation/format/", varnm[i], sep=""),
                gfile, list(verbose=verbose, na=rep(NA_integer_, nSamp)))
        }

        readmode.gdsn(n4)
        readmode.gdsn(n5)
        if (verbose)
            cat(if (.crayon()) crayon::blurred("]") else "]")
        .DigestCode(n4, digest, verbose)
        .DigestCode(n5, digest, FALSE)

        sync.gds(gfile)
    }
    invisible()
}


# Merge files
seqMerge <- function(gds.fn, out.fn, storage.option="LZMA_RA",
    info.var=NULL, fmt.var=NULL, samp.var=NULL, optimize=TRUE, digest=TRUE,
    geno.pad=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gds.fn))
    if (length(gds.fn) < 1L)
        stop("'gds.fn' should have at least one file.")
    stopifnot(is.character(out.fn), length(out.fn)==1L)

    if (is.character(storage.option))
        storage.option <- seqStorageOption(storage.option)
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))

    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))
    stopifnot(is.null(samp.var) | is.character(samp.var))

    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(geno.pad), length(geno.pad)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        .cat(sprintf("Preparing merging %d GDS files (%s):", length(gds.fn),
            date()))
    }

    # open all GDS files
    flist <- vector("list", length(gds.fn))
    on.exit({ for (f in flist) seqClose(f) })
    for (i in seq_along(gds.fn))
    {
        if (verbose)
            cat("    opening ", sQuote(gds.fn[i]), "\n", sep="")
        flist[[i]] <- seqOpen(gds.fn[i], allow.duplicate=TRUE)
    }
    if (verbose)
    {
        s <- sum(file.size(gds.fn), na.rm=TRUE)
        cat("   ", .pretty(s), "bytes in total\n")
    }

    # samples
    samp.id <- samp2.id <- seqGetData(flist[[1L]], "sample.id")
    for (f in flist[-1L])
    {
        s <- seqGetData(f, "sample.id")
        samp.id <- unique(c(samp.id, s))
        samp2.id <- intersect(samp2.id, s)
    }

    nSamp <- length(samp.id)
    if (verbose)
    {
        cat(sprintf("    %d sample%s in total (%d sample%s in common)\n",
            nSamp, .plural(nSamp), length(samp2.id), .plural(length(samp2.id))))
    }

    # variants
    if (!.is_variant_overlap_bychr(flist))
    {
        variant2.id <- NULL
        nVariant <- 0L
        for (f in flist)
        {
            nVariant <- nVariant + objdesp.gdsn(index.gdsn(f, "variant.id"))$dim
        }
    } else {
        variant.id <- variant2.id <- seqGetData(flist[[1L]], "$chrom_pos2")
        if (verbose)
        {
            cat(sprintf("    [%-2d] %s (%s variant%s)\n", 1L, basename(gds.fn[1L]),
                .pretty(length(variant.id)), .plural(length(variant.id))))
        }
        for (i in seq_along(flist)[-1L])
        {
            s <- seqGetData(flist[[i]], "$chrom_pos2")
            if (verbose)
            {
                cat(sprintf("    [%-2d] %s (%s variant%s)\n", i,
                    basename(gds.fn[i]), .pretty(length(s)), .plural(length(s))))
            }

            # variant id maybe not unique
            s1 <- intersect(variant.id, s)
            if (length(s1) <= 0L)
                variant.id <- c(variant.id, s)
            else
                variant.id <- c(variant.id, setdiff(s, s1))
            variant2.id <- intersect(variant2.id, s)
            remove(s, s1)
        }
        nVariant <- length(variant.id)
    }

    if (verbose)
    {
        cat(sprintf("    %s variant%s in total, %s variant%s in common\n",
            .pretty(nVariant), .plural(nVariant),
            .pretty(length(variant2.id)), .plural(length(variant2.id))))
    }

    # common samples
    is_merge_variant <- length(samp2.id)>0L || length(samp.id)==0L
    if (is_merge_variant)
    {
        if (length(variant2.id) > 0L)
        {
            stop("There are overlapping on both samples and variants, ",
                "please merge different samples and variants respectively.")
        }
    } else {
        varidx <- vector("list", length(flist))
        for (i in seq_along(flist))
        {
            varidx[[i]] <- match(seqGetData(flist[[i]], "$chrom_pos2"),
                variant.id)
            if (is.unsorted(varidx[[i]], strictly=TRUE))
                stop("File ", i, ": chromosomes and positions are unsorted.")
        }
    }


    ## create a GDS file
    gfile <- createfn.gds(out.fn)
    on.exit({ if (!is.null(gfile)) closefn.gds(gfile) }, add=TRUE)

    if (verbose)
        cat("Output:\n    ", normalizePath(out.fn), "\n", sep="")

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
        .AddVar(storage.option, n, "vcf.alt", unique(alt), closezip=TRUE,
            visible=FALSE)
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
        .AddVar(storage.option, n, "vcf.contig", unique(contig), closezip=TRUE,
            visible=FALSE)
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
        .AddVar(storage.option, n, "vcf.header", unique(header), closezip=TRUE,
            visible=FALSE)
    }
        

    ## add sample.id
    if (verbose) cat("Variables:\n    sample.id")
    n <- .AddVar(storage.option, gfile, "sample.id", samp.id, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    ## add variant.id
    if (verbose) cat("    variant.id")
    n <- .AddVar(storage.option, gfile, "variant.id", seq_len(nVariant),
        storage="int32", closezip=TRUE)
    .DigestCode(n, digest, verbose)

    if (is_merge_variant)
    {
        ## merge different variants

        ## add position, chromsome, allele
        # TODO: need to check whether position can be stored in 'int32'
        .append_node_variant(gfile, "position", "int32", storage.option,
            flist, nVariant, digest, verbose)
        .append_node_variant(gfile, "chromosome", "string", storage.option,
            flist, nVariant, digest, verbose)
        .append_node_variant(gfile, "allele", "string", storage.option,
            flist, nVariant, digest, verbose)
        sync.gds(gfile)

        ## add a folder for genotypes and phasing info
        if (verbose) cat("    genotype, phase [")
        varGeno <- addfolder.gdsn(gfile, "genotype")
        .MergeNodeAttr(varGeno, flist, "genotype")
        varPhase <- addfolder.gdsn(gfile, "phase")

        if (length(samp.id) > 0L)
        {
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
            if (ploidy > 2L)
            {
                n2 <- .AddVar(storage.option, varPhase, "data", storage="bit1",
                    valdim=c(ploidy-1L, nSamp, 0L))
            } else {
                n2 <- .AddVar(storage.option, varPhase, "data", storage="bit1",
                    valdim=c(nSamp, 0L))
            }

            geno_idx <- NULL
            for (i in seq_along(flist))
            {
                if (verbose)
                {
                    cat(ifelse(i > 1L, ",", ""), i, sep="")
                    flush.console()
                }
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
                geno_idx <- c(geno_idx, read.gdsn(index.gdsn(flist[[i]],
                    "genotype/@data")))
                if (geno.pad && i<length(flist) && ploidy==2L)
                {
                    if (prod(objdesp.gdsn(n1)$dim) %% 4L)
                    {
                        f <- flist[[i]]
                        seqSetFilter(f, variant.sel=.dim(f)[3L], verbose=FALSE)
                        v <- seqGetData(f, "genotype")
                        seqResetFilter(f, verbose=FALSE)
                        dim(v) <- NULL
                        vv <- rep(0L, ploidy*nSamp)
                        vv[is.na(v)] <- -1L
                        append.gdsn(n1, vv)
                        k <- length(geno_idx)
                        geno_idx[k] <- geno_idx[k] + 1L
                        if (verbose) cat(".")
                    }
                }
            }

            readmode.gdsn(n1)
            readmode.gdsn(n2)

            n <- .AddVar(storage.option, varGeno, "@data", geno_idx,
                storage="uint8", visible=FALSE)
            .DigestCode(n, digest, FALSE)

            # TODO
            n <- .AddVar(storage.option, varGeno, "extra.index",
                storage="int32", valdim=c(3L,0L))
            put.attr.gdsn(n, "R.colnames",
                c("sample.index", "variant.index", "length"))
            .AddVar(storage.option, varGeno, "extra", storage="int16",
                closezip=TRUE)

            n <- .AddVar(storage.option, varPhase, "extra.index",
                storage="int32", valdim=c(3L,0L))
            put.attr.gdsn(n, "R.colnames",
                c("sample.index", "variant.index", "length"))
            .AddVar(storage.option, varPhase, "extra", storage="bit1",
                closezip=TRUE)

            # sync file
            sync.gds(gfile)
            if (verbose) cat("]\n          ")
            .DigestCode(index.gdsn(gfile, "genotype/data"), digest, verbose)
            if (verbose) cat("          ")
            .DigestCode(index.gdsn(gfile, "phase/data"), digest, verbose)
        } else {
            if (verbose) cat("]\n")
        }

        sync.gds(gfile)

        ## add annotation folder
        varAnnot <- addfolder.gdsn(gfile, "annotation")

        # add id
        if (verbose) cat("    annotation/id ")
        n <- .AddVar(storage.option, varAnnot, "id", storage="string")
        .append_gds(n, flist, "annotation/id", verbose)
        .DigestCode(n, digest, verbose)

        # add qual
        if (verbose) cat("    annotation/qual ")
        n <- .AddVar(storage.option, varAnnot, "qual", storage="float")
        .append_gds(n, flist, "annotation/qual", verbose)
        .DigestCode(n, digest, verbose)

        # add filter
        nm <- "annotation/filter"
        if (verbose) cat("   ", nm)
        dp <- NULL; v <- NULL
        for (i in seq_along(flist))
        {
            dp <- rbind(dp, seqSummary(flist[[i]], "$filter", check="none",
                verbose=FALSE))
            a <- seqGetData(flist[[i]], nm)
            if (is.factor(a)) a <- as.character(a)
            v <- c(v, a)
        }
        dp <- unique(dp)
        if (is.factor(v) || is.character(v))
            v <- factor(v, dp$ID)
        n <- .AddVar(storage.option, varAnnot, "filter", v)
        put.attr.gdsn(n, "Description", dp$Description)
        .DigestCode(n, digest, verbose)

    } else {

        ## merge different samples

        v <- strsplit(variant.id, ":|_")

        ## add position, chromsome, allele
        # TODO: need to check whether position can be stored in 'int32'
        if (verbose) cat("    position")
        n <- .AddVar(storage.option, gfile, "position", sapply(v, `[`, i=2L),
            storage="int32")
        .DigestCode(n, digest, verbose)

        if (verbose) cat("    chromosome")
        n <- .AddVar(storage.option, gfile, "chromosome", sapply(v, `[`, i=1L),
            storage="string")
        .DigestCode(n, digest, verbose)

        if (verbose) cat("    allele")
        n <- .AddVar(storage.option, gfile, "allele", storage="string")
        .Call(SEQ_MergeAllele, nVariant, varidx, flist, n)
        readmode.gdsn(n)
        .DigestCode(n, digest, verbose)

        sync.gds(gfile)

        ## add a folder for genotypes
        if (verbose) cat("    genotype [")
        varGeno <- addfolder.gdsn(gfile, "genotype")
        .MergeNodeAttr(varGeno, flist, "genotype")

        ploidy <- seqSummary(flist[[1L]], check="none", verbose=FALSE)$ploidy
        for (i in seq_along(gds.fn))
        {
            if (!identical(ploidy, seqSummary(flist[[1L]], check="none",
                    verbose=FALSE)$ploidy))
                stop("ploidy should be the same!")
        }
        .AddVar(storage.option, varGeno, "data", storage="bit2",
            valdim=c(ploidy, nSamp, 0L))
        .AddVar(storage.option, varGeno, "@data", storage="uint8",
            visible=FALSE)

        # writing
        .Call(SEQ_MergeGeno, c(nVariant, nSamp, ploidy), varidx, flist,
            gfile, list(verbose=verbose))
        .DigestCode(readmode.gdsn(index.gdsn(varGeno, "data")), digest, verbose)
        .DigestCode(readmode.gdsn(index.gdsn(varGeno, "@data")), digest, FALSE)

        n <- .AddVar(storage.option, varGeno, "extra.index", storage="int32",
            valdim=c(3L,0L))
        put.attr.gdsn(n, "R.colnames",
            c("sample.index", "variant.index", "length"))
        readmode.gdsn(n)
        n <- .AddVar(storage.option, varGeno, "extra", storage="int16")
        readmode.gdsn(n)


        ## add phase folder
        if (verbose) cat("    phase [")
        varPhase <- addfolder.gdsn(gfile, "phase")
        if (ploidy > 2L)
        {
            .AddVar(storage.option, varPhase, "data", storage="bit1",
                valdim=c(ploidy-1L, nSamp, 0L))
        } else {
            .AddVar(storage.option, varPhase, "data", storage="bit1",
                valdim=c(nSamp, 0L))
        }

        # writing
        .Call(SEQ_MergePhase, c(nVariant, nSamp, ploidy), varidx, flist,
            gfile, list(verbose=verbose))
        .DigestCode(readmode.gdsn(index.gdsn(varPhase, "data")), digest, verbose)

        n <- .AddVar(storage.option, varPhase, "extra.index",
            storage="int32", valdim=c(3L,0L))
        put.attr.gdsn(n, "R.colnames",
            c("sample.index", "variant.index", "length"))
        readmode.gdsn(n)
        n <- .AddVar(storage.option, varPhase, "extra", storage="bit1")
        readmode.gdsn(n)

        sync.gds(gfile)

        ## add annotation folder
        varAnnot <- addfolder.gdsn(gfile, "annotation")

        # add id
        if (verbose) cat("    annotation/id")
        s <- seqGetData(flist[[1L]], "annotation/id")
        for (f in flist[-1L])
        {
            s1 <- seqGetData(f, "annotation/id")
            if (!identical(s, s1))
            {
                warning("'annotation/id' are not identical, ",
                    "and the first existing 'id' is taken.", immediate.=TRUE)
            }
        }
        n <- .AddVar(storage.option, varAnnot, "id", storage="string")
        .Call(SEQ_MergeInfo, nVariant, varidx, flist, "annotation/id",
            gfile, list(verbose=verbose))
        .DigestCode(n, digest, verbose)

        # add qual
        if (verbose) cat("    annotation/qual")
        s <- seqGetData(flist[[1L]], "annotation/qual")
        for (f in flist[-1L])
        {
            s1 <- seqGetData(f, "annotation/qual")
            if (!identical(s, s1))
            {
                warning("'annotation/qual' are not identical, ",
                    "and the first existing 'qual' is taken.", immediate.=TRUE)
            }
        }
        n <- .AddVar(storage.option, varAnnot, "qual", storage="float")
        .Call(SEQ_MergeInfo, nVariant, varidx, flist, "annotation/qual",
            gfile, list(verbose=verbose))
        .DigestCode(n, digest, verbose)

        # add filter
        if (verbose) cat("    annotation/filter")
        s <- seqGetData(flist[[1L]], "annotation/filter")
        for (f in flist[-1L])
        {
            s1 <- seqGetData(f, "annotation/filter")
            if (!identical(s, s1))
            {
                warning("'annotation/filter' are not identical, ",
                    "and the first existing 'filter' is taken.", immediate.=TRUE)
            }
        }
        n <- .AddVar(storage.option, varAnnot, "filter", storage="int32")
        .MergeNodeAttr(n, flist[1L], "annotation/filter")
        .Call(SEQ_MergeInfo, nVariant, varidx, flist, "annotation/filter",
            gfile, list(verbose=verbose))
        .DigestCode(n, digest, verbose)
    }

    sync.gds(gfile)

    ####  VCF INFO  ####

    varInfo <- addfolder.gdsn(varAnnot, "info")
    varnm <- NULL
    for (i in seq_along(flist))
    {
        n <- index.gdsn(flist[[i]], "annotation/info", silent=TRUE)
        if (!is.null(n))
        {
            varnm <- unique(c(varnm, ls.gdsn(n, recursive=TRUE,
                include.dirs=FALSE)))
        }
    }
    if (!is.null(info.var))
    {
        s <- setdiff(info.var, varnm)
        if (length(s) > 0L)
        {
            warning("No INFO variable(s): ", paste(s, collapse=", "),
                immediate.=TRUE)
        }
        varnm <- unique(intersect(info.var, varnm))
    }
    if (verbose)
        .cat("    annotation/info (", paste(varnm, collapse=","), ")")

    if (is_merge_variant)
    {
        # merge different variants
        .append_info_variant(flist, varnm, storage.option, varInfo, digest,
            verbose)
    } else {
        # merge different samples
        for (i in seq_along(varnm))
        {
            vnm <- paste0("annotation/info/", varnm[i])
            idx <- 0L
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "annotation/info", silent=TRUE)
                if (!is.null(n))
                {
                    if (varnm[i] %in% ls.gdsn(n, recursive=TRUE, include.dirs=FALSE))
                    {
                        idx <- j
                        break
                    }
                }
            }
            if (idx < 1L) stop("internal error, info field.")

            warnflag <- TRUE
            s <- seqGetData(flist[[idx]], vnm, .padNA=FALSE)
            for (f in flist[-idx])
            {
                n <- index.gdsn(f, vnm, silent=TRUE)
                if (!is.null(n))
                {
                    s1 <- seqGetData(f, vnm, .padNA=FALSE)
                    if (!identical(s, s1))
                    {
                        if (warnflag)
                        {
                            warning("'", vnm, "' are not identical, and ",
                                "the first existing '", varnm[i], "' is taken.",
                                call.=FALSE, immediate.=TRUE)
                            warnflag <- FALSE
                        }
                    }
                }
            }

            need <- FALSE
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]],
                    paste("annotation/info/", varnm[i], sep=""), silent=TRUE)
                if (is.null(n))
                {
                    need <- TRUE
                    break
                }
            }

            if (verbose) cat("        ", varnm[i], sep="")

            s <- paste0("annotation/info/", varnm[i])
            n <- index.gdsn(flist[[idx]], s)
            n1 <- index.gdsn(flist[[idx]], .var_path(s, "@"), silent=TRUE)

            dp <- objdesp.gdsn(n)
            dp$dim[length(dp$dim)] <- 0L
            n2 <- .AddVar(storage.option, varInfo, varnm[i], storage=dp$storage,
                valdim=dp$dim)
            .MergeNodeAttr(n2, flist[idx], paste0("annotation/info/", varnm[i]))

            n3 <- n1
            if (!is.null(n1) | need)
            {
                n3 <- .AddVar(storage.option, varInfo, paste0("@", varnm[i]),
                    storage="int32", visible=FALSE)
            }

            .Call(SEQ_MergeInfo, nVariant, varidx, flist,
                paste0("annotation/info/", varnm[i]),
                gfile, list(verbose=verbose))
            readmode.gdsn(n2)
            .DigestCode(n2, digest, verbose)

            if (!is.null(n3))
            {
                readmode.gdsn(n3)
                .DigestCode(n3, digest, FALSE)
            }
        }
    }

    sync.gds(gfile)

    ####  VCF FORMAT  ####

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
        varnm <- unique(intersect(fmt.var, varnm))
    }

    .append_format(flist, varnm, samp.id, gfile, is_merge_variant, nVariant,
        varidx, storage.option, digest, verbose)


    ####  sample annotation  ####

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
    if (verbose)
    {
        cat("    sample.annotation (",
            paste(varnm, collapse=","), ")\n", sep="")
    }
    if (is_merge_variant)
    {
        ## merge different variants
        for (i in seq_along(varnm))
        {
            idx <- 0L
            for (j in seq_along(flist))
            {
                n <- index.gdsn(flist[[j]], "sample.annotation", silent=TRUE)
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

            if (verbose) cat("       ", varnm[i])
            copyto.gdsn(varSamp, index.gdsn(flist[[idx]],
                paste0("sample.annotation/", varnm[i])))
            .DigestCode(index.gdsn(varSamp, varnm[i]), digest, verbose)
        }
    } else {
        ## merge different samples
        for (i in seq_along(varnm))
        {
            if (verbose) cat("       ", varnm[i])
            s <- NULL
            for (j in seq_along(flist))
            {
                f <- flist[[j]]
                vnm <- paste("sample.annotation/", varnm[i], sep="")
                n <- index.gdsn(f, vnm, silent=TRUE)
                if (!is.null(n))
                    s <- c(s, read.gdsn(n))
                else {
                    s <- c(s, rep(NA_integer_,
                        length(read.gdsn(index.gdsn(f, "sample.id")))))
                }
            }

            n <- .AddVar(storage.option, varSamp, varnm[i], s, closezip=TRUE)
            .DigestCode(n, digest, verbose)
        }
    }

    # add RLE of chromosome
    .optim_chrom(gfile)

    ####  close files  ####

    for (f in flist) seqClose(f)
    flist <- list()
    closefn.gds(gfile)
    on.exit()
    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }


    ####  optimize access efficiency  ####

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
# Reset the variant IDs in multiple GDS files
#
seqResetVariantID <- function(gds.fn, set=NULL, digest=TRUE, optimize=TRUE,
    verbose=TRUE)
{
    stopifnot(is.character(gds.fn), length(gds.fn)>0L)
    stopifnot(is.null(set) || is.logical(set))
    if (is.logical(set))
        stopifnot(length(gds.fn) == length(set))
    stopifnot(is.logical(digest), length(digest)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    n <- 0L
    i <- 1L
    for (k in seq_along(gds.fn))
    {
        fn <- gds.fn[k]
        if (verbose)
            cat(sprintf("[%2d] %s ...", i, fn))
        i <- i + 1L
        f <- seqOpen(fn, readonly=FALSE)
        dp <- objdesp.gdsn(index.gdsn(f, "variant.id"))
        len <- prod(dp$dim)
        v <- seq_len(len) + n
        if (is.null(set) || isTRUE(set[k]))
        {
            nd <- add.gdsn(f, "variant.id", v, storage="int", compress=dp$compress,
                replace=TRUE, closezip=TRUE)
            n <- n + len
            if (digest)
                .DigestCode(nd, TRUE, verbose)
            seqClose(f)
            if (optimize)
                cleanup.gds(fn, verbose=FALSE)
            if (verbose)
            {
                cat("    set new variant id: ", v[1L], " ... ", v[length(v)],
                    " [", length(v), "]\n", sep="")
            }
        } else {
            n <- n + len
            seqClose(f)
            if (verbose)
            {
                cat("\n    skip variant id: ", v[1L], " ... ", v[length(v)],
                    " [", length(v), "]\n", sep="")
            }
        }
    }

    invisible()
}
