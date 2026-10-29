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
        process_block_count, process_block_index,
        .lang_eval, .List_IRanges_value)
    TRUE
}



#######################################################################
# Get the file name of an example
#
seqExampleFileName <- function(type=c("gds", "vcf", "KG_Phase1", "dosage"))
{
    type <- match.arg(type)
    fn <- switch(type,
        gds = "CEU_Exon.gds",
        vcf = "CEU_Exon.vcf.gz",
        KG_Phase1 = "1KG_phase1_release_v3_chr22.gds",
        dosage = "Imputed_dosages.gds"
    )
    system.file("extdata", fn, package="SeqArray")
}



#######################################################################
# Setup the parallel parameters in SeqArray
#

.release_pkg_func <- function(env)
{
    # release seqarray.parallel
    cl <- .PkgEnv$seqarray.parallel
    if (inherits(cl, "cluster"))
    {
        .PkgEnv$seqarray.parallel <- NULL
        parallel::stopCluster(cl)
    }
    # release seqarray.multicore
    cl <- .PkgEnv$seqarray.multicore
    if (inherits(cl, "cluster"))
    {
        .PkgEnv$seqarray.multicore <- NULL
        parallel::stopCluster(cl)
    }
    # return none
    invisible()
}

.reg_finalizer <- function()
{
    if (is.null(.PkgEnv$.finalizer))
    {
        .PkgEnv$.finalizer <- new.env()
        reg.finalizer(.PkgEnv$.finalizer, .release_pkg_func, onexit=TRUE)
    }
    invisible()
}

seqParallelSetup <- function(cluster=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.null(cluster) || is.logical(cluster) ||
        is.numeric(cluster) || inherits(cluster, "cluster"))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (is.na(verbose)) verbose <- FALSE

    if (is.null(cluster) || isFALSE(cluster))
    {
        cl <- .PkgEnv$seqarray.parallel
        if (inherits(cl, "cluster"))
        {
            .PkgEnv$seqarray.parallel <- NULL
            stopCluster(cl)
            if (verbose)
                cat("Stop the computing cluster.\n")
        }
        return(invisible())
    }

    # Windows platform does not support forking, we have to setup a cluster
    if (.Platform$OS.type == "windows")
    {
        setup <- function(num.cores)
        {
            cl <- makeCluster(num.cores, outfile="")
            clusterCall(cl, function() {
                library(SeqArray, quietly=TRUE, verbose=FALSE)
                TRUE
            })
            cl
        }

        if (is.logical(cluster))
        {
            stopifnot(length(cluster) == 1L)
            if (isTRUE(cluster))
            {
                cl <- detectCores() - 1L
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
        if (isTRUE(cluster))
        {
            n <- detectCores() - 1L
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
                if (verbose)
                {
                    cat("Enable the computing cluster with", cluster,
                        "forked R processes.\n")
                }
            }
        }
    }

    .PkgEnv$seqarray.parallel <- cluster
    # need a finalizer when end the R session
    .reg_finalizer()
    # return
    invisible()
}


seqMulticoreSetup <- function(num=TRUE, type=c("psock", "fork"), verbose=TRUE)
{
    # check
    stopifnot(is.logical(num) || is.numeric(num), length(num)==1L)
    type <- match.arg(type)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # multicore parallel object if available
    cl <- .PkgEnv$seqarray.multicore
    if (inherits(cl, "cluster"))
    {
        stopCluster(cl)
        .PkgEnv$seqarray.multicore <- NULL
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
        .PkgEnv$seqarray.multicore <- cl
        # need a finalizer when end the R session
        .reg_finalizer()
        # return
        return(invisible(cl))

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
    ans <- .PkgEnv$seqarray.parallel
    if (is.null(ans)) ans <- FALSE
    ans
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
.init_proc <- function(idx=0L, cnt=0L, fname=FALSE)
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
    .PkgEnv$process_status_fname <- NULL
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


####  Run with a single core  ####

# run with a single core
.run_single_core <- function(gdsfile, FUN, split, .combine, .selection.flag,
    .initialize, .finalize, .initparam, ...)
{
    # initialize
    .init_proc(idx=1L, cnt=1L)
    if (is.function(.initialize)) .initialize(1L, .initparam)
    on.exit({
        .init_proc()
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


####  run with a cluster object  ####

# load the packages (e.g., '.package' in seqParallel())
.load_package <- function(package)
{
    for (pkg in package)
    {
        # 'quietly=TRUE' does not suppress the messages of masked objects
        suppressPackageStartupMessages(library(pkg, character.only=TRUE,
            quietly=TRUE, verbose=FALSE))
    }
    invisible()
}

.fc_parallel_initializer <- function(i, njobs, st_fname, gdsfn, totnum,
    init, param, package)
{
    # load the package
    .load_package(attr(gdsfn, "pkgname"))
    library(SeqArray, quietly=TRUE, verbose=FALSE)
    .load_package(package)
    # export to global variables
    .init_proc(i, njobs, st_fname)
    .set_proc_block(0L, totnum)
    .PkgEnv$process_status_fname <- st_fname
    if (is.function(init)) init(i, param)
    NULL
}

.fc_parallel_initializer_bl <- function(i, njobs, st_fname, nblock, gdsfn,
    sel, totcnt, init, param, package)
{
    # load the package
    .load_package(attr(gdsfn, "pkgname"))
    library(SeqArray, quietly=TRUE, verbose=FALSE)
    .load_package(package)
    # export to global variables
    .init_proc(i, njobs, st_fname)
    .set_proc_block(1L, nblock)
    .PkgEnv$process_status_fname <- st_fname
    # save internally
    .PkgEnv$gdsfile <- seqOpen(gdsfn, allow.duplicate=TRUE)
    .PkgEnv$sample.sel <- .decompress(sel[[1L]])
    .PkgEnv$variant.sel <- .decompress(sel[[2L]])
    .PkgEnv$totcnt <- totcnt
    # set filter
    seqSetFilter(.PkgEnv$gdsfile,
        sample.sel = .PkgEnv$sample.sel,
        variant.sel = .PkgEnv$variant.sel,
        verbose=FALSE)
    # call
    if (is.function(init)) init(i, param)
    NULL
}

.fc_parallel_finalizer <- function(i, fc, param)
{
    if (is.function(fc)) fc(i, param)
    if (inherits(.PkgEnv$gdsfile, "SeqVarGDSClass"))
        seqClose(.PkgEnv$gdsfile)
    .PkgEnv$gdsfile <- NULL
    remove("gdsfile", envir=.PkgEnv)
    .done_proc()
}

.fc_parallel_do <- function(i, .cnt, .gds.fn, .sel, FUN, .split,
    .selection.flag, ...)
{
    # export to global variables
    .set_proc_block(i, NULL)
    # open the file
    f <- seqOpen(.gds.fn, allow.duplicate=TRUE)
    on.exit(seqClose(f))
    # set filter
    seqSetFilter(f,
        sample.sel = .decompress(.sel[[1L]]),
        variant.sel = .decompress(.sel[[2L]]),
        verbose=FALSE)
    s <- .Call(SEQ_SplitSelection, f, .split, i, .cnt, .selection.flag)
    # call the user-defined function
    if (.selection.flag) FUN(f, s, ...) else FUN(f, ...)
}

.fc_parallel_do_bl <- function(i, FUN, .split, .sel_idx, .bl_size,
    .selection.flag, ...)
{
    # set the block index
    .set_proc_block(i, NULL)
    # set filter
    .s <- .Call(SEQ_SplitSelectionX, .PkgEnv$gdsfile, i, .split,
        .sel_idx, .PkgEnv$variant.sel, .PkgEnv$sample.sel,
        .bl_size, .selection.flag, .PkgEnv$totcnt)
    # call the user-defined function
    if (.selection.flag)
        FUN(.PkgEnv$gdsfile, .s, ...)
    else
        FUN(.PkgEnv$gdsfile, ...)
}

.fc_parallel_do_0 <- function(i, FUN, ...)
{
    # export to global variables
    .set_proc_block(i, NULL)
    # call the user-defined function
    FUN(...)
}

.fc_parallel_do_1 <- function(i, FUN, ...)
{
    # export to global variables
    .set_proc_block(i, NULL)
    # call the user-defined function
    FUN(i, ...)
}

# run with a cluster object
.run_parallel_cluster <- function(cl, gdsfile, FUN, split, .combine,
    .selection.flag, .initialize, .finalize, .initparam,
    .balancing, .bl_size, .bl_progress, .status_file, .package, ...)
{
    # number of nodes
    stopifnot(inherits(cl, "cluster"))
    njobs <- length(cl)
    # status files
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
        sel$sample.sel <- .compress(sel$sample.sel)
        sel$variant.sel <- .compress(sel$variant.sel)
    } else {
        sel <- list(sample.sel=raw(), variant.sel=raw())
    }

    # call the user-defined function
    if (!isTRUE(.balancing) || split=="none")
    {
        totnum <- njobs
        if (is.numeric(gdsfile)) totnum <- as.integer(gdsfile)

        # initialize the cluster
        clusterApply(cl, seq_len(njobs), fun = .fc_parallel_initializer,
            njobs=njobs, st_fname=st_fname,
            gdsfn = if (inherits(gdsfile, "SeqVarGDSClass"))
                gdsfile$filename else NULL,
            totnum = if (is.numeric(gdsfile)) totnum else 0L,
            init=.initialize, param=.initparam, package=.package)
        # exit call
        on.exit({
            clusterApply(cl, seq_len(njobs), fun=.fc_parallel_finalizer,
                fc=.finalize, param=.initparam)
        }, add=TRUE)

        # call in parallel
        if (is.null(gdsfile))
        {
            ans <- .DynamicClusterCall(cl, totnum, .fun = .fc_parallel_do_0,
                .combinefun=.combine, FUN = FUN, ...)
        } else if (inherits(gdsfile, "SeqVarGDSClass"))
        {
            ans <- .DynamicClusterCall(cl, totnum, .fun = .fc_parallel_do,
                .combinefun=.combine, .cnt=totnum, .gds.fn = gdsfile$filename,
                .sel = sel, FUN = FUN, .split = split,
                .selection.flag = .selection.flag, ...
            )
        } else {
            stopifnot(is.numeric(gdsfile))
            ans <- .DynamicClusterCall(cl, totnum, .fun = .fc_parallel_do_1,
                .combinefun=.combine, FUN = FUN, ...)
        }

    } else {
        ## load balancing

        # selection indexing
        v <- .prepare_balancing(gdsfile, njobs, dm, split, .bl_size)
        .bl_size <- v$bl_size
        progress <- if (.bl_progress) .seqProgress(v$totcnt) else NULL
        updatefun <- if (.bl_progress)
            function(i) .seqProgForward(progress, .bl_size) else NULL

        # initialize the cluster & distribute data
        clusterApply(cl, seq_len(njobs),
            fun = .fc_parallel_initializer_bl,
            njobs = njobs, st_fname = st_fname, nblock = v$nblock,
            gdsfn = gdsfile$filename, sel = sel, totcnt = v$totcnt,
            init = .initialize, param = .initparam, package = .package
        )
        remove(sel)
        # exit call
        on.exit({
            clusterApply(cl, seq_len(njobs), fun=.fc_parallel_finalizer,
                fc=.finalize, param=.initparam)
        })

        # do parallel
        ans <- .DynamicClusterCall(cl, v$nblock,
            .fun = .fc_parallel_do_bl,
            .combinefun=.combine, .updatefun=updatefun, FUN = FUN,
            .split = (split=="by.variant"), .sel_idx = v$sel_idx,
            .bl_size = .bl_size, .selection.flag = .selection.flag, ...
        )
        # close progress object
        remove(progress)
    }

    # return
    if (is.list(ans) & identical(.combine, "unlist"))
        ans <- unlist(ans, recursive=FALSE)
    ans
}


####  Run using forking in Unix/Linux    ####

# run using forking in Unix/Linux
.run_forking <- function(njobs, gdsfile, FUN, split, .combine, .selection.flag,
    .initialize, .finalize, .initparam, .balancing, .bl_size, .bl_progress,
    .status_file, ...)
{
    # check # of nodes
    stopifnot(is.numeric(njobs))
    # status files
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
            if (njobs > dm[3L]) njobs <- dm[3L]
        } else {
            if (njobs > dm[2L]) njobs <- dm[2L]
        }
    }

    # export to global variables
    .init_proc(0L, 0L, st_fname)
    .set_proc_block()
    .PkgEnv$process_status_fname <- st_fname
    on.exit({
        .done_proc()
        .PkgEnv$process_status_fname <- NULL
    }, add=TRUE)

    # the first argument is NA for forking
    if (is.function(.initialize))
        .initialize(NA_integer_, .initparam)
    if (is.function(.finalize))
        on.exit(.finalize(NA_integer_, .initparam), add=TRUE)

    # call the user-defined function
    if (!isTRUE(.balancing) || split=="none")
    {
        # run
        if (is.null(gdsfile))
        {
            ans <- .DynamicForkCall(njobs, njobs, .fun = function(.ji, .i, ...)
            {
                .init_proc(.ji, njobs, NULL)  # export to global variables
                FUN(...)
            }, .combinefun=.combine, .updatefun=NULL, ...)
        } else if (inherits(gdsfile, "SeqVarGDSClass"))
        {
            ans <- .DynamicForkCall(njobs, njobs, .fun = function(.ji, .i, ...)
            {
                # actually .ji = .i
                # export to global variables
                .init_proc(.ji, njobs, NULL)
                # set filter
                .s <- .Call(SEQ_SplitSelection, gdsfile, split, .i, njobs,
                    .selection.flag)
                # call the user-defined function
                if (.selection.flag)
                    FUN(gdsfile, .s, ...)
                else
                    FUN(gdsfile, ...)
            }, .combinefun=.combine, .updatefun=NULL, ...)
        } else {
            stopifnot(is.numeric(gdsfile))
            totnum <- as.integer(gdsfile)
            # call forking
            ans <- .DynamicForkCall(njobs, totnum, .fun = function(.ji, .i, ...)
            {
                # export to global variables
                .init_proc(.ji, njobs, NULL)
                .set_proc_block(.i, totnum)
                FUN(.i, ...)
            }, .combinefun=.combine, .updatefun=NULL, ...)
        }
    } else {
        ## load balancing

        # selection indexing
        v <- .prepare_balancing(gdsfile, njobs, dm, split, .bl_size)
        .bl_size <- v$bl_size
        .sel_idx <- v$sel_idx
        totcnt <- v$totcnt
        nblock <- v$nblock
        remove(v)
        progress <- if (.bl_progress) .seqProgress(totcnt) else NULL
        updatefun <- if (.bl_progress)
            function(i) .seqProgForward(progress, .bl_size) else NULL
        .sel <- seqGetFilter(gdsfile, .useraw=TRUE)
        split <- split == "by.variant"

        # do parallel
        ans <- .DynamicForkCall(njobs, nblock, .fun = function(.ji, .i, ...)
        {
            # set the process & block index
            .init_proc(.ji, njobs, NULL)
            .set_proc_block(.i, nblock)
            # set filter
            .s <- .Call(SEQ_SplitSelectionX, gdsfile, .i, split,
                .sel_idx, .sel$variant.sel, .sel$sample.sel, .bl_size,
                .selection.flag, totcnt)
            # call the user-defined function
            if (.selection.flag)
                FUN(gdsfile, .s, ...)
            else
                FUN(gdsfile, ...)
        }, .combinefun=.combine,
            .updatefun=function(i) .seqProgForward(progress, .bl_size), ...)

        remove(progress)
    }

    # return
    if (is.list(ans) & identical(.combine, "unlist"))
        ans <- unlist(ans, recursive=FALSE)
    ans
}


####  Run with BiocParallel  ####

# call .dofun(i, ...) for i in 1:.num with a BiocParallelParam object, like
#   .DynamicClusterCall(): the jobs are sent to the workers one by one, and
#   the results are combined by .combinefun once they are received
.bp_dynamic_call <- function(cl, .num, .dofun, .combinefun, .updatefun=NULL,
    ...)
{
    argone <- FALSE
    if (is.function(.combinefun))
        argone <- length(formals(args(.combinefun))) == 1L
    if (identical(.combinefun, "unlist") | identical(.combinefun, "list"))
        ans <- vector("list", .num)
    else
        ans <- NULL

    # the number of received results
    k <- 0L
    # the error in .combinefun or .updatefun, which is raised after all
    #   running jobs are done, since an error in bpiterate(REDUCE=) would leave
    #   their results to the next call on 'cl'
    err <- NULL
    # combine a result d = list(job index, value)
    reduce <- function(x, d)
    {
        if (!is.null(err)) return(NULL)
        err <<- tryCatch({
            dv <- d[[2L]]
            if (is.function(.combinefun))
            {
                if (argone)
                {
                    .combinefun(dv)
                } else {
                    if (is.null(ans))
                        ans <<- dv
                    else
                        ans <<- .combinefun(ans, dv)
                }
            } else if (.combinefun %in% c("unlist", "list"))
            {
                # assignment NULL would remove it from the list
                if (!is.null(dv)) ans[[d[[1L]]]] <<- dv
            }
            k <<- k + 1L
            if (is.function(.updatefun)) .updatefun(k)
            NULL
        }, error = function(e) e)
        NULL
    }

    if (inherits(cl, "DoparParam"))
    {
        # bpiterate() is not supported, so combine after all jobs are done
        v <- BiocParallel::bplapply(seq_len(.num), .bp_tag, .dofun=.dofun,
            ..., BPPARAM=cl)
        for (d in v) reduce(NULL, d)
    } else {
        # the job index, no new job after an error
        i <- 0L
        iter <- function() if (is.null(err) && i < .num) (i <<- i + 1L)
        BiocParallel::bpiterate(iter, .bp_tag, .dofun=.dofun, ...,
            REDUCE=reduce, init=NULL, BPPARAM=cl)
    }
    if (!is.null(err)) stop(err)

    # output
    ans
}

# call .dofun(i, ...) in a worker, and return list(i, value)
.bp_tag <- function(i, .dofun, ...) list(i, .dofun(i, ...))

# the functions below are called in the workers of BiocParallel; they are not
#   closures, so that the caller's environment is not sent to the workers;
#   the user-defined function and its arguments are passed in '.fun' and
#   '.args' (a list), so that the arguments of the user-defined function
#   (e.g., '.FUN' in seqApply()) do not match the arguments of bplapply() or
#   these functions

# initialize a worker once per call (identified by '.init$token'): load the
#   packages, open the GDS file for load balancing and call .initialize();
#   'i' is the first job run by the worker in this call, and it is used as the
#   process id, which is in 1..njobs for SnowParam and MulticoreParam since
#   the first njobs jobs are sent to different workers
.bp_init_worker <- function(i, .init)
{
    if (identical(.PkgEnv$bp_worker$token, .init$token)) return(invisible())
    # the previous call is not finalized in this worker
    .bp_clear_worker()
    # load the packages
    .load_package(.init$pkgname)
    library(SeqArray, quietly=TRUE, verbose=FALSE)
    .load_package(.init$package)
    .PkgEnv$bp_worker <- list(id=i)
    # export to global variables; no status file if 'i' is out of range, which
    #   is possible with the backends other than SnowParam and MulticoreParam
    st <- if (i <= .init$cnt) .init$st_fname
    .init_proc(i, .init$cnt, if (is.null(st)) FALSE else st)
    .set_proc_block(0L, .init$blcnt)
    .PkgEnv$process_status_fname <- st
    # open the file for load balancing
    if (!is.null(.init$gds.fn))
    {
        f <- seqOpen(.init$gds.fn, allow.duplicate=TRUE)
        .PkgEnv$bp_worker$gdsfile <- f
        ss <- .decompress(.init$sel[[1L]])
        vs <- .decompress(.init$sel[[2L]])
        seqSetFilter(f, sample.sel=ss, variant.sel=vs, verbose=FALSE)
        .PkgEnv$bp_worker$sample.sel <- ss
        .PkgEnv$bp_worker$variant.sel <- vs
    }
    # call the user-defined function
    if (is.function(.init$initialize)) .init$initialize(i, .init$param)
    .PkgEnv$bp_worker$token <- .init$token
    invisible()
}

# close the GDS file opened in the worker, and remove the state of the call
.bp_clear_worker <- function()
{
    f <- .PkgEnv$bp_worker$gdsfile
    .PkgEnv$bp_worker <- NULL
    if (inherits(f, "SeqVarGDSClass")) seqClose(f)
}

# call .finalize() if the worker is initialized in the call '.token'
.bp_final <- function(i, .token, .fun, .param)
{
    on.exit({ .bp_clear_worker(); .done_proc() })
    s <- .PkgEnv$bp_worker
    if (identical(s$token, .token) && is.function(.fun))
        .fun(s$id, .param)
    NULL
}

# call FUN(...) if 'gdsfile' is NULL
.bp_do_0 <- function(i, .cnt, .fun, .args, .init)
{
    .bp_init_worker(i, .init)
    .init_proc(i, .cnt, NULL)
    .set_proc_block(i, NULL)
    do.call(.fun, .args, quote=TRUE)
}

# call FUN(i, ...) if 'gdsfile' is numeric
.bp_do_1 <- function(i, .fun, .args, .init)
{
    .bp_init_worker(i, .init)
    .set_proc_block(i, NULL)
    do.call(.fun, c(list(i), .args), quote=TRUE)
}

# call FUN(gdsfile, ...) on the i-th part of the selection
.bp_do <- function(i, .cnt, .gds.fn, .sel, .fun, .args, .split,
    .selection.flag, .init)
{
    .bp_init_worker(i, .init)
    .init_proc(i, .cnt, NULL)
    .set_proc_block(i, NULL)
    # open the file
    f <- seqOpen(.gds.fn, allow.duplicate=TRUE)
    on.exit(seqClose(f))
    # set filter
    seqSetFilter(f,
        sample.sel = .decompress(.sel[[1L]]),
        variant.sel = .decompress(.sel[[2L]]),
        verbose=FALSE)
    s <- .Call(SEQ_SplitSelection, f, .split, i, .cnt, .selection.flag)
    # call the user-defined function
    do.call(.fun, c(if (.selection.flag) list(f, s) else list(f), .args),
        quote=TRUE)
}

# call FUN(gdsfile, ...) on the i-th block of the selection, using the file
#   opened in .bp_init_worker()
.bp_do_bl <- function(i, .fun, .args, .split, .sel_idx, .bl_size,
    .selection.flag, .totcnt, .init)
{
    .bp_init_worker(i, .init)
    .set_proc_block(i, NULL)
    w <- .PkgEnv$bp_worker
    f <- w$gdsfile
    # set filter
    s <- .Call(SEQ_SplitSelectionX, f, i, .split, .sel_idx, w$variant.sel,
        w$sample.sel, .bl_size, .selection.flag, .totcnt)
    # call the user-defined function
    do.call(.fun, c(if (.selection.flag) list(f, s) else list(f), .args),
        quote=TRUE)
}

# run with a BiocParallelParam object
.run_biocparallel <- function(cl, gdsfile, FUN, split, .combine,
    .selection.flag, .initialize, .finalize, .initparam,
    .balancing, .bl_size, .bl_progress, .status_file, .package, ...)
{
    # number of workers
    nworker <- BiocParallel::bpnworkers(cl)
    njobs <- nworker
    # the arguments passed to FUN
    args <- list(...)

    if (split %in% c("by.variant", "by.sample"))
    {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
        dm <- .seldim(gdsfile)
        # no more jobs than the selected variants or samples
        njobs <- min(njobs, if (split == "by.variant") dm[3L] else dm[2L])
        sel <- seqGetFilter(gdsfile, .useraw=TRUE)
        sel$sample.sel <- .compress(sel$sample.sel)
        sel$variant.sel <- .compress(sel$variant.sel)
    } else {
        dm <- NULL
        sel <- list(sample.sel=raw(), variant.sel=raw())
    }
    # status files
    st_fname <- NULL
    if (isTRUE(.status_file))
    {
        st_fname <- .get_status_filename(njobs)
        file.create(st_fname, showWarnings=FALSE)
        on.exit(unlink(st_fname, force=TRUE))
    }

    # start the workers if needed, so that they are kept during this call
    isup <- BiocParallel::bpisup(cl)
    if (!isup) BiocParallel::bpstart(cl)
    # the identifier of this call, so that each worker is initialized once
    .PkgEnv$bp_ncall <- if (is.null(.PkgEnv$bp_ncall)) 1L else
        .PkgEnv$bp_ncall + 1L
    token <- paste(Sys.getpid(), .PkgEnv$bp_ncall)
    init <- list(token=token,
        pkgname = if (inherits(gdsfile, "SeqVarGDSClass"))
            attr(gdsfile$filename, "pkgname"),
        package=.package, cnt=njobs, blcnt=0L, st_fname=st_fname,
        initialize=.initialize, param=.initparam)
    # finalize the workers on exit, and stop the workers if started here;
    #   each worker receives one of the jobs with SnowParam and
    #   MulticoreParam, but it is not guaranteed with the other backends
    on.exit({
        tryCatch({
            if (is.function(.finalize) || isup)
            {
                BiocParallel::bplapply(seq_len(nworker), .bp_final,
                    .token=token, .fun=.finalize, .param=.initparam,
                    BPPARAM=cl)
            }
        }, finally = if (!isup) BiocParallel::bpstop(cl))
    }, add=TRUE)

    if (!isTRUE(.balancing) || split == "none")
    {
        ####  non-balancing path  ####

        if (is.null(gdsfile))
        {
            # no GDS file: call FUN(...) on each worker
            ans <- .bp_dynamic_call(cl, njobs, .dofun = .bp_do_0,
                .combinefun=.combine, .cnt=njobs, .fun=FUN, .args=args,
                .init=init)

        } else if (inherits(gdsfile, "SeqVarGDSClass"))
        {
            # GDS file: open, set filter, split, call FUN
            ans <- .bp_dynamic_call(cl, njobs, .dofun = .bp_do,
                .combinefun=.combine, .cnt=njobs, .gds.fn=gdsfile$filename,
                .sel=sel, .fun=FUN, .args=args, .split=split,
                .selection.flag=.selection.flag, .init=init)

        } else {
            # numeric gdsfile: call FUN(i, ...) for i in 1:gdsfile
            stopifnot(is.numeric(gdsfile))
            init$blcnt <- as.integer(gdsfile)
            ans <- .bp_dynamic_call(cl, as.integer(gdsfile),
                .dofun = .bp_do_1, .combinefun=.combine, .fun=FUN,
                .args=args, .init=init)
        }

    } else {
        ####  load-balancing path  ####

        # selection indexing
        v <- .prepare_balancing(gdsfile, njobs, dm, split, .bl_size)
        progress <- if (.bl_progress) .seqProgress(v$totcnt) else NULL
        updatefun <- if (.bl_progress)
            function(i) .seqProgForward(progress, v$bl_size) else NULL

        # the file is opened once in each worker
        init$gds.fn <- gdsfile$filename
        init$sel <- sel
        init$blcnt <- v$nblock
        remove(sel)

        # distribute blocks across workers
        ans <- .bp_dynamic_call(cl, v$nblock, .dofun = .bp_do_bl,
            .combinefun=.combine, .updatefun=updatefun, .fun=FUN, .args=args,
            .split=(split=="by.variant"), .sel_idx=v$sel_idx,
            .bl_size=v$bl_size, .selection.flag=.selection.flag,
            .totcnt=v$totcnt, .init=init)

        # close progress object
        remove(progress)
    }

    # return
    if (is.list(ans) & identical(.combine, "unlist"))
        ans <- unlist(ans, recursive=FALSE)
    ans
}


# call the user-defined function in parallel
seqParallel <- function(cl=seqGetParallel(), gdsfile, FUN,
    split=c("by.variant", "by.sample", "none"), .combine="unlist",
    .selection.flag=FALSE, .initialize=NULL, .finalize=NULL, .initparam=NULL,
    .balancing=FALSE, .bl_size=NA_integer_, .bl_progress=FALSE,
    .status_file=FALSE, .proc_time=FALSE, .package=NULL, ...)
{
    # check
    stopifnot(is.null(cl) || is.logical(cl) || is.numeric(cl) ||
        inherits(cl, "cluster") || inherits(cl, "BiocParallelParam"))
    stopifnot(is.null(gdsfile) || inherits(gdsfile, "SeqVarGDSClass") ||
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
    stopifnot(is.logical(.proc_time), length(.proc_time)==1L)
    stopifnot(is.null(.package) || is.character(.package))
    stopifnot(is.logical(.balancing), length(.balancing)==1L)
    if (is.na(.balancing))
        .balancing <- isTRUE(getOption("seqarray.balancing", TRUE))
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
    tm <- as.double(Sys.time())
    .init_proc()

    # get the number of workers
    cl <- .McoreParallel(cl)
    njobs <- .NumParallel(cl)
    if (njobs <= 1L)
    {
        # run with a single core
        .load_package(.package)
        ans <- .run_single_core(gdsfile, FUN, split, .combine, .selection.flag,
            .initialize, .finalize, .initparam, ...)

    } else if (inherits(cl, "cluster"))
    {
        # multiple processes with a cluster object
        ans <- .run_parallel_cluster(cl, gdsfile, FUN, split, .combine,
            .selection.flag, .initialize, .finalize, .initparam,
            .balancing, .bl_size, .bl_progress, .status_file, .package, ...)

    } else if (inherits(cl, "BiocParallelParam"))
    {
        ## multiple processes with a predefined cluster from BiocParallel
        ans <- .run_biocparallel(cl, gdsfile, FUN, split, .combine,
            .selection.flag, .initialize, .finalize, .initparam,
            .balancing, .bl_size, .bl_progress, .status_file, .package, ...)

    } else {
        # forking processes, when cl is an integer
        if (getOption("seqarray.nofork", FALSE) || .Platform$OS.type=="windows")
        {
            # no forking on windows
            capture.output({ cl <- makeCluster(njobs, outfile="") })
            on.exit(stopCluster(cl))
            # call
            return(seqParallel(cl, gdsfile, FUN, split, .combine,
                .selection.flag, .initialize, .finalize, .initparam,
                .balancing, .bl_size, .bl_progress, .status_file,
                .package=.package, ...))
        }
        # load the packages before forking
        .load_package(.package)
        # call
        ans <- .run_forking(njobs, gdsfile, FUN, split, .combine,
            .selection.flag, .initialize, .finalize, .initparam,
            .balancing, .bl_size, .bl_progress, .status_file, ...)
    }

    # show the elapsed time
    if (isTRUE(.proc_time) && isFALSE(.balancing && .bl_progress))
    {
        tm <- as.double(Sys.time()) - tm
        .cat("\r[==================================================] 100%",
            ", complete, ", .Call(SEQ_SecToTime, tm), "    ")
    }

    # output
    if (identical(.combine, "none"))
        invisible()
    else
        ans
}


seqParApply <- function(cl=seqGetParallel(), x, FUN, .balancing=TRUE, ...)
{
    # check
    cl <- .McoreParallel(cl)
    njobs <- .NumParallel(cl, "cl")
    stopifnot(is.logical(.balancing), length(.balancing)==1L)
    if (is.na(.balancing))
        .balancing <- isTRUE(getOption("seqarray.balancing", TRUE))

    if (njobs <= 1L)
    {
        # a single process
        ans <- lapply(x, FUN, ...)
    } else if (inherits(cl, "BiocParallelParam"))
    {
        # multiple cores
        ans <- BiocParallel::bplapply(x, FUN, ..., BPPARAM=cl)
    } else {
        # multiple cores
        if (!inherits(cl, "cluster"))
        {
            if (.Platform$OS.type == "windows")
                cl <- makeCluster(njobs)
            else
                cl <- makeForkCluster(njobs)
            on.exit({ stopCluster(cl) })
        }

        # a load balancing version
        if (.balancing)
            ans <- clusterApplyLB(cl, x, FUN, ...)
        else
            ans <- clusterApply(cl, x, FUN, ...)
    }

    # output
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


seqTranspose <- function(gdsfile, var.name, compress=NULL, digest=TRUE,
    verbose=TRUE)
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

# Run-length encoding on chromosome codes
.optim_chrom <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "gds.class"))
    # read
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
    # return
    invisible()
}

seqOptimize <- function(gdsfile, target=c("chromosome", "by.sample"),
    format.var=TRUE, cleanup=TRUE, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "gds.class") || is.character(gdsfile))
    target <- match.arg(target)
    stopifnot(is.logical(format.var) || is.character(format.var))
    stopifnot(is.logical(cleanup), length(cleanup)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (is.na(verbose)) verbose <- FALSE

    # the open file
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1L)
        if (verbose)
            .cat("Open ", sQuote(basename(gdsfile)))
        gdsfn <- gdsfile
        gdsfile <- seqOpen(gdsfn, readonly=FALSE)
        on.exit(seqClose(gdsfile))
    } else {
        # do not call cleanup.gds
        cleanup <- FALSE
    }
    if (gdsfile$readonly)
        stop("'gdsfile' should not be read-only.")

    # process
    if ("by.sample" %in% target)
    {
        # genotype
        if (verbose) cat("Working on 'genotype' ...\n")
        .Transpose(gdsfile, "genotype/data", "~")

        # phase
        if (verbose) cat("Working on 'phase' ...\n")
        .Transpose(gdsfile, "phase/data", "~")

        # annotation - format
        if (isTRUE(format.var) || is.character(format.var))
        {
            n <- index.gdsn(gdsfile, "annotation/format", silent=TRUE)
            if (!is.null(n))
            {
                nm <- ls.gdsn(n)
                if (isTRUE(format.var)) format.var <- nm
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
            cat("Adding run-length encoding for chromosomes ...")
        # run-length encoding on chromosome codes
        .optim_chrom(gdsfile)
        if (verbose) cat("  [Done]\n")
    }

    if (isTRUE(cleanup))
    {
        on.exit()
        seqClose(gdsfile)
        cleanup.gds(gdsfn, verbose=verbose)
    }
    invisible()
}
