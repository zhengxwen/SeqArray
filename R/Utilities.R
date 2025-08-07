#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Data Management of Large-scale Whole-Genome Sequence Variant Calls
#


#######################################################################

# R expressions called in C internally
.lang_eval <- list(
    new_rle = quote(new("Rle", values=values, lengths=lengths)),
    as_List = quote(as(object, "List"))
)

# Compressed Lists in IRanges
.List_IRanges_value <- list(
    logical   = LogicalList(),
    integer   = IntegerList(),
    numeric   = NumericList(),
    character = CharacterList(),
    raw       = RawList(),
    factor    = FactorList(),
    "unlistData", "partitioning", "end"
)

.onLoad <- function(lib, pkg)
{
    .Call(SEQ_Pkg_Init, .dim_name, process_count, process_index,
        .lang_eval, .List_IRanges_value)
    TRUE
}

.Last.lib <- function(libpath)
{
    # parallel object
    cl <- getOption("seqarray.parallel")
    if (inherits(cl, "cluster"))
        stopCluster(cl)
    options(seqarray.parallel=NULL)

    # multicore parallel object if available
    cl <- getOption("seqarray.multicore")
    if (inherits(cl, "cluster"))
        stopCluster(cl)
    options(seqarray.multicore=NULL)
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
            stopCluster(opt)
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
            cl <- makeCluster(num.cores)
            clusterCall(cl, function() {
                library(SeqArray, quietly=TRUE, verbose=FALSE)
                TRUE
            })
            cl
        }

        if (is.logical(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster)
            {
                cl <- detectCores() - 1L
                if (cl <= 1L) cl <- 2L
                cluster <- setup(cl)
                if (verbose)
                    cat("Enable the computing cluster with", cl, "R processes.\n")
            } else {
                if (verbose) cat("No computing cluster.\n")
            }
        } else if (is.numeric(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster > 1L)
            {
                cl <- cluster
                cluster <- setup(cluster)
                if (verbose)
                    cat("Enable the computing cluster with", cl, "R processes.\n")
            }
        }
    } else {
        # unix forking technique
        if (identical(cluster, TRUE))
        {
            n <- detectCores() - 1L
            if (n <= 1L) n <- 2L
            if (verbose)
                cat("Enable the computing cluster with", n, "forked R processes.\n")
        } else if (is.numeric(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (cluster > 1L)
            {
                if (verbose)
                    cat("Enable the computing cluster with", cluster, "forked R processes.\n")
            }
        }
    }

    options(seqarray.parallel=cluster)
    invisible()
}


seqMulticoreSetup <- function(num, type=c("psock", "fork"), verbose=TRUE)
{
    # check
    stopifnot(is.logical(num) || is.numeric(num), length(num)==1L)
    type <- match.arg(type)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # multicore parallel object if available
    cl <- getOption("seqarray.multicore")
    if (inherits(cl, "cluster"))
    {
        stopCluster(cl)
        options(seqarray.multicore=NULL)
    }

    # setup
    if (isTRUE(num))
        num <- detectCores()
    if (!is.na(num) && (num > 1L))
    {
        if ((type == "psock") || (.Platform$OS.type == "windows"))
        {
            cl <- makeCluster(num, outfile="")
            if (isTRUE(verbose))
                cat("Enable the multicore cluster: ")
        } else {
            cl <- makeForkCluster(num)
            if (isTRUE(verbose))
                cat("Enable the multicore cluster via forking: ")
        }
        if (isTRUE(verbose)) print(cl)
        options(seqarray.multicore=cl)
    } else {
        if (isTRUE(verbose))
            cat("Disable the pre-defined multicore cluster.\n")
    }

    # return nothing
    invisible()
}



#######################################################################
# Get the parallel parameters in SeqArray
#
seqGetParallel <- function()
{
    getOption("seqarray.parallel", FALSE)
}



#######################################################################
# Storage options for the SeqArray GDS file
#
seqStorageOption <- function(compression=c("ZIP_RA", "ZIP_RA.fast",
    "ZIP_RA.max", "LZ4_RA", "LZ4_RA.fast", "LZ4_RA.max",
    "LZMA_RA", "LZMA_RA.fast", "LZMA_RA.max", "Ultra", "UltraMax", "none"),
    mode=NULL, float.mode="float32", geno.compress=NULL, info.compress=NULL,
    format.compress=NULL, index.compress=NULL, ...)
{
    # check
    compression <- match.arg(compression)
    z <- switch(compression,
        none="", Ultra="LZMA_RA.ultra", UltraMax="LZMA_RA.ultra_max")
    if (!is.null(z)) compression <- z

    stopifnot(is.null(mode) | is.character(mode))
    stopifnot(is.character(float.mode), length(float.mode) > 0L)

    if (!is.null(geno.compress))
        stopifnot(is.character(geno.compress), length(geno.compress)==1L)
    if (!is.null(info.compress))
        stopifnot(is.character(info.compress), length(info.compress)==1L)
    if (!is.null(format.compress))
        stopifnot(is.character(format.compress), length(format.compress)==1L)
    if (!is.null(index.compress))
        stopifnot(is.character(index.compress), length(index.compress)==1L)

    if (compression %in% c("ZIP_RA.max", "LZ4_RA.max", "LZMA_RA.max"))
    {
        suf_b <- ":1M"; suf_i <- ":1M"; suf_f <- ":4M"
        if (is.null(index.compress)) index.compress <- compression
    } else if (compression=="LZMA_RA.ultra")
    {
        suf_b <- ":4M"; suf_i <- ":4M"; suf_f <- ":8M"
        if (is.null(index.compress)) index.compress <- "LZMA.max"
    } else if (compression=="LZMA_RA.ultra_max")
    {
        suf_b <- suf_i <- suf_f <- ":8M"
        if (is.null(index.compress)) index.compress <- "LZMA.max"
    } else {
        suf_b <- suf_i <- ""
        suf_f <- ":1M"
        if (is.null(index.compress)) index.compress <- compression
    }

    if (is.null(geno.compress))
    {
        if (compression == "LZMA_RA.ultra")
            geno.compress <- "LZMA_RA.ultra:1M"
        else if (compression == "LZMA_RA.ultra_max")
            geno.compress <- "LZMA_RA.ultra_max:8M"
    }

    rv <- list(compression = paste0(compression, suf_b),
        mode = mode, float.mode = float.mode,
        geno.compress = ifelse(is.null(geno.compress), compression,
            geno.compress),
        info.compress = ifelse(is.null(info.compress),
            ifelse(compression=="", "", paste0(compression, suf_i)),
            info.compress),
        format.compress = ifelse(is.null(format.compress),
            ifelse(compression=="", "", paste0(compression, suf_f)),
            format.compress),
        index.compress = index.compress,
        ...)
    class(rv) <- "SeqGDSStorageClass"
    return(rv)
}



#######################################################################
# Apply functions in parallel
#

# get the temporary file names for showing the status of child processes
.get_status_filename <- function(njobs)
{
    i <- as.integer(Sys.time())
    if (is.na(i)) i <- 0L
    i <- (i %% 90000L) + 10000L
    basename(tempfile(
        pattern=sprintf("_parallel_%d_n%d-%d_", i, njobs, seq_len(njobs)),
        tmpdir=".", fileext=".txt"))
}

# initialize the internal variables for child processes
.init_proc <- function(idx=1L, cnt=1L, fname=FALSE)
{
    .Call(SEQ_SetProcess, idx, cnt, fname)
}

# set the block current index and total count
.set_proc_block <- function(idx=0L, cnt=0L)
{
    .Call(SEQ_SetProcessBlock, idx, cnt)
}

# finalize the child process by resetting the status
.done_proc <- function()
{
    .init_proc()
    .set_proc_block()
    .packageEnv$process_status_fname <- NULL
}

# get the multiple for balancing jobs, the default is 4
.get_bl_multiple <- function()
{
    m <- getOption("seqarray.balancing_multiple", process_balancing_multiple)
    if (!is.numeric(m) || length(m)!=1L || !is.finite(m))
        stop("'seqarray.balancing_multiple' should be a positive number.")
    m
}

# auto set the block size
.auto_bl_size <- function(cnt, njobs)
{
    m <- .get_bl_multiple()
    n <- as.integer(ceiling(cnt / (njobs * m)))
    if (n <= 0L) n <- 1L
    n
}

# get the number of blocks according to the block size
.get_nblock_bl_size <- function(njobs, cnt, bl_size)
{
    if (is.na(bl_size) || (bl_size<=0L))
        bl_size <- .auto_bl_size(cnt, njobs)
    nblock <- cnt %/% bl_size
    if (cnt %% bl_size) nblock <- nblock + 1L
    if (nblock < njobs)
    {
        bl_size <- cnt %/% njobs
        if (cnt %% njobs) bl_size <- bl_size + 1L
        nblock <- cnt %/% bl_size
        if (cnt %% bl_size) nblock <- nblock + 1L
    }
    c(nblock, bl_size)  # return
}

# prepare data for load balancing with a GDS file
.prepare_balancing <- function(gdsfile, njobs, dm, split, .bl_size)
{
    if (split == "by.variant")
    {
        totcnt <- dm[3L]
        v <- .get_nblock_bl_size(njobs, totcnt, .bl_size)
        nblock <- v[1L]
        .bl_size <- v[2L]
        ii <- seqGetData(gdsfile, "$variant_index")
    } else {
        totcnt <- dm[2L]
        v <- .get_nblock_bl_size(njobs, totcnt, .bl_size)
        nblock <- v[1L]
        .bl_size <- v[2L]
        ii <- seqGetData(gdsfile, "$sample_index")
    }
    # return
    list(totcnt = totcnt, nblock = nblock, bl_size = .bl_size,
        sel_idx = ii[seq.int(1L, by=.bl_size, length.out=nblock)])
}


# run with a single core
.run_single_core <- function(gdsfile, FUN, split, .combine, .selection.flag,
    .initialize, .finalize, .initparam, ...)
{
    # initialize
    if (is.function(.initialize)) .initialize(1L, .initparam)
    on.exit({
        if (is.function(.finalize)) .finalize(1L, .initparam)        
    })
    # call the user-defined function
    if (!is.numeric(gdsfile))
    {
        if (.selection.flag)
        {
            # gdsfile should be GDS file
            stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
            dm <- .seldim(gdsfile)
            if (split == "by.variant")
                ans <- FUN(gdsfile, rep(TRUE, dm[3L]), ...)
            else if (split == "by.sample")
                ans <- FUN(gdsfile, rep(TRUE, dm[2L]), ...)
            else
                ans <- FUN(gdsfile, NULL, ...)
        } else {
            if (is.null(gdsfile))
                ans <- FUN(...)
            else
                ans <- FUN(gdsfile, ...)
        }
        # finalize
        if (is.function(.combine))
        {
            if (length(formals(args(.combine))) == 1L)
                .combine(ans)
        }
    } else {
        # initial values
        argone <- FALSE
        if (is.function(.combine))
            argone <- length(formals(args(.combine))) == 1L
        is_c_list <-
            identical(.combine, "unlist") || identical(.combine, "list")
        ans <- if (is_c_list) vector("list", gdsfile) else NULL
        # for-loop
        for (i in seq_len(gdsfile))
        {
            v <- FUN(i, ...)  # call
            if (is.function(.combine))
            {
                if (argone)
                {
                    .combine(v)
                } else {
                    ans <- if (is.null(ans)) v else .combine(ans, v)
                }
            } else if (is_c_list)
            {
                # assignment NULL would remove it from the list
                if (!is.null(v)) ans[[i]] <- v
            }
        }
    }
    # return
    ans
}


# run with a cluster object
.run_parallel_cluster <- function(cl, gdsfile, FUN, split, .combine,
    .selection.flag, .initialize, .finalize, .initparam,
    .balancing, .bl_size, .bl_progress, .status_file, ...)
{
    # number of nodes
    njobs <- length(cl)
    st_fname <- NULL
    if (isTRUE(.status_file))
    {
        st_fname <- .get_status_filename(njobs)
        file.create(st_fname, showWarnings=FALSE)
        on.exit(unlink(st_fname, force=TRUE))
    }

    if (split %in% c("by.variant", "by.sample"))
    {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
        dm <- .seldim(gdsfile)
        if (split == "by.variant")
        {
            if (njobs > dm[3L]) cl <- cl[1L:dm[3L]]
        } else {
            if (njobs > dm[2L]) cl <- cl[1L:dm[2L]]
        }
        sel <- seqGetFilter(gdsfile, .useraw=TRUE)
        sel$sample.sel <- memCompress(sel$sample.sel, type="gzip")
        sel$variant.sel <- memCompress(sel$variant.sel, type="gzip")
    } else {
        sel <- list(sample.sel=raw(), variant.sel=raw())
    }

    # call the user-defined function
    if (!isTRUE(.balancing) || split=="none")
    {
        totnum <- njobs
        if (is.numeric(gdsfile)) totnum <- as.integer(gdsfile)
        # initialize the cluster
        clusterApply(cl, seq_len(njobs),
            fun = function(i, njobs, st_fname, totnum, init, param)
        {
            # load the package
            library(SeqArray, quietly=TRUE, verbose=FALSE)
            # export to global variables
            .init_proc(i, njobs, st_fname)
            .set_proc_block(0L, totnum)
            .packageEnv$process_status_fname <- st_fname
            if (is.function(init)) init(i, param)
            # return none
            NULL
        },  njobs=njobs, st_fname=st_fname,
            totnum = if (is.numeric(gdsfile)) totnum else 0L,
            init=.initialize, param=.initparam)
        # exit call
        on.exit({
            clusterApply(cl, seq_len(njobs), fun = function(i, fc, param)
            {
                if (is.function(fc)) fc(i, param)
                .done_proc()
            }, fc=.finalize, param=.initparam)
        }, add=TRUE)
        # call in parallel
        if (is.null(gdsfile))
        {
            ans <- .DynamicClusterCall(cl, totnum, .fun = function(i, FUN, ...)
            {
                # export to global variables
                .set_proc_block(i, NULL)
                # call the user-defined function
                FUN(...)
            }, .combinefun=.combine, FUN = FUN, ...)
        } else if (inherits(gdsfile, "SeqVarGDSClass"))
        {
            ans <- .DynamicClusterCall(cl, totnum, .fun =
                function(i, .cnt, .gds.fn, .sel, FUN, .split, .selection.flag, ...)
            {
                # export to global variables
                .set_proc_block(i, NULL)
                # open the file
                f <- seqOpen(.gds.fn, allow.duplicate=TRUE)
                on.exit(seqClose(f))
                # set filter
                seqSetFilter(f,
                    sample.sel = memDecompress(.sel[[1L]], type="gzip"),
                    variant.sel = memDecompress(.sel[[2L]], type="gzip"),
                    verbose=FALSE)
                s <- .Call(SEQ_SplitSelection, f, .split, i, .cnt,
                    .selection.flag)
                # call the user-defined function
                if (.selection.flag) FUN(f, s, ...) else FUN(f, ...)
            },  .combinefun=.combine, .cnt=totnum, .gds.fn = gdsfile$filename,
                .sel = sel, FUN = FUN, .split = split,
                .selection.flag = .selection.flag, ...
            )
        } else {
            stopifnot(is.numeric(gdsfile))
            ans <- .DynamicClusterCall(cl, totnum, .fun = function(i, FUN, ...)
            {
                # export to global variables
                .set_proc_block(i, NULL)
                # call the user-defined function
                FUN(i, ...)
            },  .combinefun=.combine, FUN = FUN, ...)
        }

    } else {
        ## load balancing

        # selection indexing
        v <- .prepare_balancing(gdsfile, njobs, dm, split, .bl_size)
        .bl_size <- v$bl_size
        progress <- if (.bl_progress) .seqProgress(v$totcnt, njobs) else NULL
        updatefun <- function(i) .seqProgForward(progress, .bl_size)

        # initialize the cluster & distribute data
        clusterApply(cl, seq_len(njobs),
            fun = function(i, njobs, st_fname, nblock, gdsfn, sel, totcnt,
                init, param)
        {
            # load the package
            library(SeqArray, quietly=TRUE, verbose=FALSE)
            # export to global variables
            .init_proc(i, njobs, st_fname)
            .set_proc_block(1L, nblock)
            .packageEnv$process_status_fname <- st_fname
            # save internally
            .packageEnv$gfile <- seqOpen(gdsfn, allow.duplicate=TRUE)
            .packageEnv$sample.sel <- memDecompress(sel[[1L]], type="gzip")
            .packageEnv$variant.sel <- memDecompress(sel[[2L]], type="gzip")
            .packageEnv$totcnt <- totcnt
            # set filter
            seqSetFilter(.packageEnv$gfile,
                sample.sel = .packageEnv$sample.sel,
                variant.sel = .packageEnv$variant.sel,
                verbose=FALSE)
            # call
            if (is.function(init)) init(i, param)
            NULL
        },  njobs = njobs, st_fname = st_fname, nblock = v$nblock,
            gdsfn = gdsfile$filename, sel = sel, totcnt = v$totcnt,
            init = .initialize, param = .initparam
        )
        remove(sel)

        # finalize
        on.exit({
            clusterApply(cl, seq_len(njobs), fun = function(i, fc, param)
            {
                if (is.function(fc)) fc(i, param)
                if (inherits(.packageEnv$gfile, "SeqVarGDSClass"))
                    seqClose(.packageEnv$gfile)
                .packageEnv$gfile <- NULL
                .done_proc()
            }, fc=.finalize, param=.initparam)
        })

        # do parallel
        ans <- .DynamicClusterCall(cl, v$nblock, .fun =
            function(i, FUN, .split, .sel_idx, .bl_size, .selection.flag, ...)
        {
            # set the block index
            .set_proc_block(i, NULL)
            # set filter
            s <- .Call(SEQ_SplitSelectionX, .packageEnv$gfile, i, .split,
                .sel_idx, .packageEnv$variant.sel, .packageEnv$sample.sel,
                .bl_size, .selection.flag, .packageEnv$totcnt)
            # call the user-defined function
            if (.selection.flag)
                FUN(.packageEnv$gfile, s, ...)
            else
                FUN(.packageEnv$gfile, ...)
        }, .combinefun=.combine, .updatefun=updatefun, FUN = FUN,
            .split = (split=="by.variant"), .sel_idx = v$sel_idx,
            .bl_size = .bl_size, .selection.flag = .selection.flag, ...
        )
        # close progress object
        remove(progress)
    }

    if (is.list(ans) & identical(.combine, "unlist"))
        ans <- unlist(ans, recursive=FALSE)
    # return
    ans
}


# call the user-defined function in parallel
seqParallel <- function(cl=seqGetParallel(), gdsfile, FUN,
    split=c("by.variant", "by.sample", "none"), .combine="unlist",
    .selection.flag=FALSE, .initialize=NULL, .finalize=NULL, .initparam=NULL,
    .balancing=FALSE, .bl_size=NA_integer_, .bl_progress=FALSE,
    .status_file=FALSE, ...)
{
    # check
    stopifnot(is.null(cl) | is.logical(cl) | is.numeric(cl) |
        inherits(cl, "cluster") | inherits(cl, "BiocParallelParam"))
    stopifnot(is.null(gdsfile) | inherits(gdsfile, "SeqVarGDSClass") |
        is.numeric(gdsfile))
    if (is.numeric(gdsfile))
    {
        if (!is.finite(gdsfile) || gdsfile<1L)
        {
            stop("'gdsfile' should be NULL, a SeqVarGDSClass object ",
                "or a positive integer.")
        }
    }
    stopifnot(is.function(FUN))
    split <- match.arg(split)
    stopifnot(is.character(.combine) | is.function(.combine))
    stopifnot(is.null(.initialize) | is.function(.initialize))
    stopifnot(is.null(.finalize) | is.function(.finalize))
    stopifnot(is.logical(.selection.flag), length(.selection.flag)==1L)
    stopifnot(is.logical(.status_file), length(.status_file)==1L)
    stopifnot(is.logical(.balancing), length(.balancing)==1L)
    if (isTRUE(.balancing))
    {
        stopifnot(is.numeric(.bl_size), length(.bl_size)==1L)
        stopifnot(is.logical(.bl_progress), length(.bl_progress)==1L)
    }

    if (is.character(.combine))
    {
        stopifnot(length(.combine) == 1L)
        if (!(.combine %in% c("unlist", "list", "none")))
            .combine <- match.fun(.combine)
    }

    # check dimension
    if (is.null(gdsfile) || is.numeric(gdsfile))
    {
        if (split != "none")
            stop("'split' should be 'none' if 'gdsfile' is NULL or an integer.")
    } else {
        dm <- .seldim(gdsfile)
        if (split == "by.variant")
        {
            if (dm[3L] <= 0L) stop("No variants selected.")
        } else if (split == "by.sample") {
            if (dm[2L] <= 0L) stop("No samples selected.")
        }
    }

    # initialize internally
    .init_proc()

    # get the number of workers
    njobs <- .NumParallel(cl)
    cl <- .McoreParallel(cl)
    if (njobs <= 1L)
    {
        # run with a single core
        ans <- .run_single_core(gdsfile, FUN, split, .combine, .selection.flag,
            .initialize, .finalize, .initparam, ...)

    } else if (inherits(cl, "cluster"))
    {
        # multiple processes with a cluster object
        ans <- .run_parallel_cluster(cl, gdsfile, FUN, split, .combine,
            .selection.flag, .initialize, .finalize, .initparam,
            .balancing, .bl_size, .bl_progress, .status_file, ...)

    } else if (inherits(cl, "BiocParallelParam"))
    {
        ## multiple processes with a predefined cluster from BiocParallel

        if (is.function(.initialize))
        {
            BiocParallel::bplapply(seq_len(njobs), function(i) .initialize(i),
                BPPARAM=cl)
        }
        on.exit({
            if (is.function(.finalize))
            {
                BiocParallel::bplapply(seq_len(njobs), function(i) .finalize(i),
                    BPPARAM=cl)
            }
        })

        if (split %in% c("by.variant", "by.sample"))
        {
            sel <- seqGetFilter(gdsfile, .useraw=TRUE)
        } else {
            sel <- list(sample.sel=raw(), variant.sel=raw())
        }

        ans <- BiocParallel::bplapply(seq_len(njobs), FUN =
            function(.proc_idx, .proc_cnt, .gds.fn, .sel_sample, .sel_variant,
                .FUN, .split, .selection.flag, ...)
        {
            # load the package
            library(SeqArray, quietly=TRUE, verbose=FALSE)
            # export to global variables
            .init_proc(.proc_idx, .proc_cnt, NULL)

            if (is.null(.gds.fn))
            {
                # call the user-defined function
                .FUN(...)
            } else {
                # open the file
                .file <- seqOpen(.gds.fn, allow.duplicate=TRUE)
                on.exit({ seqClose(.file) })

                # set filter
                seqSetFilter(.file,
                    sample.sel = memDecompress(.sel_sample, type="gzip"),
                    variant.sel = memDecompress(.sel_variant, type="gzip"),
                    verbose=FALSE)
                .ss <- .Call(SEQ_SplitSelection, .file, .split, .proc_idx,
                    .proc_cnt, .selection.flag)

                # call the user-defined function
                if (.selection.flag) FUN(.file, .ss, ...) else FUN(.file, ...)
            }

        },  .proc_cnt = njobs, .gds.fn = gdsfile$filename,
            .sel_sample = memCompress(sel$sample.sel, type="gzip"),
            .sel_variant = memCompress(sel$variant.sel, type="gzip"),
            .FUN = FUN, .split = split, .selection.flag=.selection.flag,
            ..., BPPARAM=cl)

        if (is.list(ans))
        {
            if (identical(.combine, "unlist"))
            {
                ans <- unlist(ans, recursive=FALSE)
            } else if (is.function(.combine))
            {
                rv <- ans[[1L]]
                for (i in seq_len(length(ans)-1L) + 1L)
                    rv <- .combine(rv, ans[[i]])
                ans <- rv
            }
        }

    } else {
        ## forking processes

        if (getOption("seqarray.nofork", FALSE) || .Platform$OS.type=="windows")
        {
            # no forking on windows
            capture.output({ cl <- makeCluster(njobs, outfile="") })
            on.exit(stopCluster(cl))
            return(seqParallel(cl, gdsfile, FUN, split, .combine, .selection.flag,
                .initialize, .finalize, .initparam,
                .balancing, .bl_size, .bl_progress, ...))
        }

        st_fname <- NULL
        if (isTRUE(.status_file))
        {
            st_fname <- .get_status_filename(njobs)
            file.create(st_fname, showWarnings=FALSE)
            on.exit(unlink(st_fname, force=TRUE), add=TRUE)
        }

        if (is.function(.initialize))
            .initialize(NA_integer_, .initparam)
        if (is.function(.finalize))
            on.exit(.finalize(NA_integer_, .initparam), add=TRUE)

        if (split %in% c("by.variant", "by.sample"))
        {
            if (split == "by.variant")
            {
                if (njobs > dm[3L]) njobs <- dm[3L]
            } else {
                if (njobs > dm[2L]) njobs <- dm[2L]
            }
        }

        if (!isTRUE(.balancing) || split=="none")
        {
            ans <- .DynamicForkCall(njobs, njobs,
                .fun = function(.ji, .idx, .st_fname, ...)
                {
                    # export to global variables
                    .init_proc(.ji, njobs, .st_fname)
                    .set_proc_block()
                    .packageEnv$process_status_fname <- .st_fname
                    # set exit
                    on.exit({
                        .init_proc()
                        .packageEnv$process_status_fname <- NULL
                    })
                    if (!is.null(gdsfile))
                    {
                        # set filter
                        .ss <- .Call(SEQ_SplitSelection, gdsfile, split,
                            .ji, njobs, .selection.flag)
                        # call the user-defined function
                        if (.selection.flag)
                            FUN(gdsfile, .ss, ...)
                        else
                            FUN(gdsfile, ...)
                    } else {
                        FUN(...)
                    }
                },  .combinefun=.combine, .updatefun=NULL,
                    .st_fname=st_fname, ...)
        } else {
            ## load balancing

            # selection indexing
            .sel <- seqGetFilter(gdsfile)
            if (split == "by.variant")
            {
                if (is.na(.bl_size) || (.bl_size<=0L))
                    .bl_size <- .auto_bl_size(dm[3L], njobs)
                totnum <- dm[3L] %/% .bl_size
                if (dm[3L] %% .bl_size) totnum <- totnum + 1L
                if (totnum < njobs)
                {
                    .bl_size <- dm[3L] %/% njobs
                    if (dm[3L] %% njobs) .bl_size <- .bl_size + 1L
                    totnum <- dm[3L] %/% .bl_size
                    if (dm[3L] %% .bl_size) totnum <- totnum + 1L
                }
                .sel_idx <- which(.sel$variant.sel)
            } else {
                if (is.na(.bl_size) || (.bl_size<=0L))
                    .bl_size <- .auto_bl_size(dm[2L], njobs)
                totnum <- dm[2L] %/% .bl_size
                if (dm[2L] %% .bl_size) totnum <- totnum + 1L
                if (totnum < njobs)
                {
                    .bl_size <- dm[2L] %/% njobs
                    if (dm[2L] %% njobs) .bl_size <- .bl_size + 1L
                    totnum <- dm[2L] %/% .bl_size
                    if (dm[2L] %% .bl_size) totnum <- totnum + 1L
                }
                .sel_idx <- which(.sel$sample.sel)
            }
            .proglen <- length(.sel_idx)
            progress <- if (.bl_progress) .seqProgress(.proglen, njobs) else NULL

            .sel_idx <- .sel_idx[seq.int(1L, by=.bl_size, length.out=totnum)]
            .sel <- seqGetFilter(gdsfile, .useraw=TRUE)
            split <- split == "by.variant"

            # export to global variables
            .packageEnv$process_status_fname <- st_fname
            on.exit({ .packageEnv$process_status_fname <- NULL }, add=TRUE)

            # do parallel
            ans <- .DynamicForkCall(njobs, totnum, .fun = function(ji, i, ...)
            {
                # set the process & block index
                .init_proc(ji, njobs, NULL)
                .set_proc_block(i, totnum)
                # set filter
                .ss <- .Call(SEQ_SplitSelectionX, gdsfile, i, split,
                    .sel_idx, .sel$variant.sel, .sel$sample.sel,
                    .bl_size, .selection.flag, .proglen)
                # call the user-defined function
                if (.selection.flag)
                    FUN(gdsfile, .ss, ...)
                else
                    FUN(gdsfile, ...)
            },  .combinefun=.combine,
                .updatefun=function(i) .seqProgForward(progress, .bl_size), ...)

            remove(progress)
        }

        if (is.list(ans) & identical(.combine, "unlist"))
            ans <- unlist(ans, recursive=FALSE)
    }

    # output
    if (identical(.combine, "none"))
        invisible()
    else
        ans
}


seqParApply <- function(cl=seqGetParallel(), x, FUN, load.balancing=TRUE, ...)
{
    # check
    njobs <- .NumParallel(cl, "cl")
    parallel <- .McoreParallel(cl)
    stopifnot(is.logical(load.balancing), length(load.balancing)==1L)

    if (njobs <= 1L)
    {
        ## a single process
        ans <- lapply(x, FUN, ...)
    } else if (inherits(cl, "BiocParallelParam"))
    {
        ans <- BiocParallel::bplapply(x, FUN, ..., BPPARAM=cl)
    } else {
        if (!inherits(cl, "cluster"))
        {
            if (.Platform$OS.type == "windows")
                cl <- makeCluster(njobs)
            else
                cl <- makeForkCluster(njobs)
            on.exit({ stopCluster(cl) })
        }

        # a load balancing version
        if (isTRUE(load.balancing))
            ans <- clusterApplyLB(cl, x, FUN, ...)
        else
            ans <- clusterApply(cl, x, FUN, ...)
    }

    ans
}



#######################################################################
# Modify the SeqArray data structure
#######################################################################

#######################################################################
# Delete data variables
#

seqDelete <- function(gdsfile, info.var=character(), fmt.var=character(),
    samp.var=character(), verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    if (gdsfile$readonly)
        stop("The GDS file is read-only.")

    stopifnot(is.character(info.var))
    stopifnot(is.character(fmt.var))
    stopifnot(is.character(samp.var))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose) cat("Delete INFO variable(s):")
    if (length(info.var) == 0L) cat(" <None>")
    for (i in seq_along(info.var))
    {
        nm <- info.var[i]
        n <- index.gdsn(gdsfile, paste0("annotation/info/", nm))
        delete.gdsn(n, force=TRUE)
        n <- index.gdsn(gdsfile, paste0("annotation/info/@", nm), silent=TRUE)
        if (!is.null(n))
            delete.gdsn(n, force=TRUE)
        if (verbose) cat(ifelse(i>1L, ",", ""), nm)
    }
    if (verbose) cat("\n")

    if (verbose) cat("Delete FORMAT variable(s):")
    if (length(fmt.var) == 0L) cat(" <None>")
    for (i in seq_along(fmt.var))
    {
        nm <- fmt.var[i]
        n <- index.gdsn(gdsfile, paste0("annotation/format/", nm))
        delete.gdsn(n, force=TRUE)
        if (verbose) cat(ifelse(i>1L, ",", ""), nm)
    }
    if (verbose) cat("\n")

    if (verbose) cat("Delete Sample Annotation variable(s):")
    if (length(samp.var) == 0L) cat(" <None>")
    for (i in seq_along(samp.var))
    {
        nm <- samp.var[i]
        n <- index.gdsn(gdsfile, paste0("sample.annotation/", nm))
        delete.gdsn(n, force=TRUE)
        if (verbose) cat(ifelse(i>1L, ",", ""), nm)
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


seqTranspose <- function(gdsfile, var.name, compress=NULL, digest=TRUE, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & is.vector(var.name))
    stopifnot(length(var.name) == 1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

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
        moveto.gdsn(newnode, node, relpos="after")

        # write data
        apply.gdsn(node, margin=length(dm)-1L, as.is="none", FUN=function(g) {
            append.gdsn(newnode, g)
        }, .useraw=TRUE)

        readmode.gdsn(newnode)
        .DigestCode(newnode, digest, FALSE)
        if (verbose)
            print(newnode, attribute=TRUE, attribute.trim=TRUE)
    } else
        warning("It is a vector.")

    invisible()
}



#######################################################################
# Optimize data by transposing
#

.optim_chrom <- function(gdsfile)
{
    n <- index.gdsn(gdsfile, "chromosome")
    readmode.gdsn(n)
    chr <- read.gdsn(n)
    s <- rle(chr)
    n1 <- add.gdsn(gdsfile, "@chrom_rle_val", s$values, replace=TRUE,
        visible=FALSE)
    n2 <- add.gdsn(gdsfile, "@chrom_rle_len", s$lengths, replace=TRUE,
        visible=FALSE)
    moveto.gdsn(n2, n)
    moveto.gdsn(n1, n)
    invisible()
}

seqOptimize <- function(gdsfn, target=c("chromosome", "by.sample"),
    format.var=TRUE, cleanup=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfn), length(gdsfn)==1L)
    target <- match.arg(target)
    stopifnot(is.logical(format.var) || is.character(format.var))
    stopifnot(is.logical(cleanup), length(cleanup)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    gdsfile <- seqOpen(gdsfn, readonly=FALSE)
    on.exit({ seqClose(gdsfile) })

    if ("by.sample" %in% target)
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
    } else if ("chromosome" %in% target)
    {
        if (verbose)
            cat("Adding run-length encoding for chromosome coding ...")
        .optim_chrom(gdsfile)
        if (verbose)
            cat(" [Done]\n")
    }

    if (cleanup)
    {
        on.exit()
        seqClose(gdsfile)
        cleanup.gds(gdsfn, verbose=verbose)
    }

    invisible()
}
