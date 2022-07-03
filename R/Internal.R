#######################################################################
#
# Package Name: SeqArray
#
# Description: Internal functions
#


#######################################################################
# Package-wide variables

.packageEnv <- new.env()

# the index of current children process
process_index <- 1L
# the number of children processes
process_count <- 1L

## R objects for class, dimnames ...
.dim_name <- list(
    geno_dim2 = list(allele=NULL, sample=NULL),
    geno_dim3 = list(allele=NULL, sample=NULL, variant=NULL),
    dosage_dim = list(sample=NULL, variant=NULL),
    data_dim  = c("length", "data"),
    data_dim2 = list(sample=NULL, variant=NULL),
    data_class = "SeqVarDataList"
)


#######################################################################
# Internal C functions
#
.cfunction0 <- function(name)
{
    fn <- function() NULL
    f <- quote(.Call(SEQ_ExternalName0))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction <- function(name)
{
    fn <- function(x) NULL
    f <- quote(.Call(SEQ_ExternalName1, x))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction2 <- function(name)
{
    fn <- function(x, y) NULL
    f <- quote(.Call(SEQ_ExternalName2, x, y))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction3 <- function(name)
{
    fn <- function(x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName3, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction4 <- function(name)
{
    fn <- function(w, x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName4, w, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction5 <- function(name)
{
    fn <- function(v, w, x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName5, v, w, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}


#######################################################################
# Get the numbers of selected samples and variants
#
.seldim <- function(gdsfile)
{
    # seldim[1L] -- ploidy
    # seldim[2L] -- # of selected samples
    # seldim[3L] -- # of selected variants
    .Call(SEQ_Summary, gdsfile, "genotype")$seldim
}

.dim <- function(gdsfile)
{
    # dim[1L] -- ploidy
    # dim[2L] -- # of total samples
    # dim[3L] -- # of total variants
    .Call(SEQ_Summary, gdsfile, "genotype")$dim
}


#######################################################################
# Detect whether there are integer genotypes or numeric dosages
#
.has_geno <- function(gdsfile)
{
    exist.gdsn(gdsfile, "genotype/data")
}

.has_dosage <- function(gdsfile)
{
    nm <- getOption("seqarray.node_ds", "annotation/format/DS")
    nm2 <- paste0(nm, "/data")
    if (!exist.gdsn(gdsfile, nm2))
        stop("No ", sQuote("genotype/data"), " or ", sQuote(nm2), ".")
    nm
}


#######################################################################
# Display functions
#
.cat <- function(...) cat(..., "\n", sep="")

.crayon <- function()
{
    flag <- getOption("gds.crayon", TRUE)
    isTRUE(flag) && requireNamespace("crayon", quietly=TRUE)
}



#######################################################################
# get pretty number with big mark
#
.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

.plural <- function(num) if (num > 1L) "s" else ""

.pretty_size <- function(s)
{
    size <- function(x)
    {
        if (x >= 1024^4)
            sprintf("%.1fT", x / 1024^4)
        else if (x >= 1024^3)
            sprintf("%.1fG", x / 1024^3)
        else if (x >= 1024^2)
            sprintf("%.1fM", x / 1024^2)
        else if (x >= 1024)
            sprintf("%.1fK", x / 1024)
        else
            sprintf("%gB", x)
    }
    sapply(s, size)
}



#######################################################################
# Variable path
#
.var_path <- function(var.name, prefix)
{
    nm <- unlist(strsplit(var.name, "/"))
    i <- length(nm)
    if (i > 0L)
    {
        nm[i] <- paste(prefix, nm[i], sep="")
        nm <- paste(nm, collapse="/")
    } else
        nm <- character(0)
    nm
}



#######################################################################
# append values
#
.append_gds <- function(target.node, gdslist, varname, verbose)
{
    s <- "."
    if (.crayon()) s <- crayon::blurred(s)
    .MergeNodeAttr(target.node, gdslist, varname)
    for (i in seq_along(gdslist))
    {
        append.gdsn(target.node, index.gdsn(gdslist[[i]], varname))
        if (verbose)
        {
            cat(s)
            flush(stdout()); flush.console()
        }
    }
    readmode.gdsn(target.node)
    invisible()
}

.append_rep_gds <- function(node, elm, count)
{
    .Call(SEQ_AppendFill, node, elm, count)
}



#######################################################################
# Data Type define VCF Format
#
.vcf_type <- function(node)
{
    a <- get.attr.gdsn(node)
    if (is.null(a$Type))
    {
        v <- switch(as.character(objdesp.gdsn(node)$type[1L]),
            Raw="Integer", Integer="Integer", Real="Float", Factor="String",
            String="String", Logical="Flag", ".")
    } else
        v <- as.character(a$Type[1L])
    v
}



#######################################################################
# Clear the internal variable map for better performance in seqGetData()
# if only a few variable accessed in the following steps (e.g.,
# seqBlockApply(), seqUnitApply())
#
.clear_varmap <- function(gdsfile)
{
    .Call(SEQ_ClearVarMap, gdsfile)
}



#######################################################################
# Open and close a connection,
# Please always call '.close_conn' after '.open_bin' and '.open_text'

.last_str <- function(s, len)
{
    substring(s, nchar(s)-len+1L, nchar(s))
}

.open_bin <- function(filename)
{
    stopifnot(is.character(filename))
    con2 <- NULL
    fmt <- ""

    if ((substr(filename, 1L, 6L) == "ftp://") |
        (substr(filename, 1L, 7L) == "http://"))
    {
        if (.last_str(filename, 3L) == ".gz")
        {
            con <- gzcon(url(filename, "rb"))
            fmt <- "gz"
        } else
            con <- url(filename, "rb")
    } else {
        if (.last_str(filename, 3L) == ".gz")
        {
            con <- gzfile(filename, "rb")
            fmt <- "gz"
        } else if (.last_str(filename, 3L) == ".xz")
        {
            con <- xzfile(filename, "rb")
            fmt <- "xz"
        } else 
            con <- file(filename, "rb")
    }

    # open(con)
    list(filename=filename, con=con, con2=con2, fmt=fmt)
}

.open_text <- function(filename, require.txtmode=FALSE)
{
    stopifnot(is.character(filename))
    con2 <- NULL

    if ((substr(filename, 1L, 6L) == "ftp://") |
        (substr(filename, 1L, 7L) == "http://"))
    {
        if (.last_str(filename, 3L) == ".gz")
        {
            con <- gzcon(url(filename, "rb"))
            if (require.txtmode)
            {
                fn <- tempfile(fileext=".tmpfile")
                write(readLines(con), file=fn, sep="\n")
                close(con)
                con <- fn
                con2 <- NULL
            }
        } else
            con <- url(filename, "rt")
    } else
        con <- file(filename, "rt")

    # open(con)
    list(filename=filename, con=con, con2=con2)
}

.close_conn <- function(conn)
{
    if (is.character(conn$con))
    {
        if (.last_str(conn$con, 8L) == ".tmpfile")
            unlink(conn$con, force=TRUE)
    } else if (inherits(conn$con, "connection"))
    {
        close(conn$con)
    }

    if (inherits(conn$con2, "connection"))
    {
        close(conn$con2)
    }
}



#######################################################################
# Parallel functions

# need parallel? how many? return 1 if no parallel
.NumParallel <- function(cl, nm="parallel")
{
    if (is.null(cl) | identical(cl, FALSE))
    {
        ans <- 1L
    } else if (is.numeric(cl))
    {
        if (length(cl) != 1L)
            stop("'parallel' should be length-one.")
        if (is.na(cl)) cl <- 1L
        if (cl < 1L) cl <- 1L
        mc <- getOption("seqarray.multicore")
        if (inherits(mc, "cluster"))
        {
            if (cl > length(mc)) cl <-  length(mc)
        }
        ans <- as.integer(cl)
    } else if (isTRUE(cl))
    {
        mc <- getOption("seqarray.multicore")
        if (inherits(mc, "cluster"))
        {
            ans <- length(mc)
        } else {
            ans <- detectCores() - 1L
            if (is.na(ans)) ans <- 2L
            if (ans <= 1L) ans <- 2L
        }
    } else if (inherits(cl, "cluster"))
    {
        ans <- length(cl)
        if (ans < 1L) ans <- 1L
    } else if (inherits(cl, "BiocParallelParam"))
    {
        ans <- BiocParallel::bpnworkers(cl)
        if (ans < 1L) ans <- 1L
    } else
        stop("Invalid '", nm, "'.")
    if (ans > 128L)  # limited by R itself
        stop("It is unable to allocate resources for more than 128 nodes.")
    ans
}

# check if the multicore cluster is specified
.McoreParallel <- function(parallel)
{
    if (is.numeric(parallel) || isTRUE(parallel))
    {
        mc <- getOption("seqarray.multicore")
        if (inherits(mc, "cluster"))
        {
            if (isTRUE(parallel) || (parallel >= length(mc)))
            {
                if (parallel > length(mc))
                {
                    warning("Using at most ", length(mc),
                        " cores in the user-defined multicore cluster",
                        call.=FALSE, immediate.=TRUE)
                }
                parallel <- mc
            } else {
                parallel <- mc[seq_len(parallel)]
            }
        }
    }
    parallel
}

# forking implements or not
.IsForking <- function(cl)
{
    if (.Platform$OS.type == "windows")
    {
        FALSE
    } else {
        isTRUE(cl) || is.numeric(cl)
    }
}



#######################################################################
# Parallel functions
#   cl -- a cluster object
#   .num -- the total number of segments
#   .fun -- a user-defined function
#   .combinefun -- a user-defined function for combining the returned values
#   .updatefun -- a user-defined function for updating progress (could be NULL)
#   ... -- other parameters passed to .fun
#

.parse_send_data <- quote(
    parallel:::sendData(con, list(type=type,data=value,tag=tag)))
.postNode <- function(con, type, value=NULL, tag=NULL) eval(.parse_send_data)

.sendCall <- function(con, fun, args, return=TRUE, tag=NULL)
{
    .postNode(con, "EXEC", list(fun=fun, args=args, return=return, tag=tag))
    NULL
}

.parse_recv_one_data <- quote(parallel:::recvOneData(cl))
.recvOneResult <- function(cl)
{
    v <- eval(.parse_recv_one_data)
    list(value=v$value$value, node=v$node, tag=v$value$tag)
}

.DynamicClusterCall <- function(cl, .num, .fun, .combinefun,
    .updatefun=NULL, ...)
{
    # in order to use the internal functions accessed by ':::'
    # the functions are all defined in 'parallel/R/snow.R'

    # check
    stopifnot(is.null(cl) | inherits(cl, "cluster"))
    stopifnot(is.numeric(.num))
    stopifnot(is.function(.fun))
    stopifnot(is.character(.combinefun) | is.function(.combinefun))
    stopifnot(is.null(.updatefun) | is.function(.updatefun))
    argone <- FALSE
    if (is.function(.combinefun))
        argone <- length(formals(args(.combinefun))) == 1L

    if (!is.null(cl))
    {
        ## parallel implementation

        if (identical(.combinefun, "unlist") | identical(.combinefun, "list"))
            ans <- vector("list", .num)
        else
            ans <- NULL

        p <- length(cl)
        if ((.num > 0L) && p)
        {
            ####  this closure is sending to all nodes

            argfun <- function(i) c(i, list(...))
            submit <- function(node, job)
                .sendCall(cl[[node]], .fun, argfun(job), tag = job)

            for (i in 1:min(.num, p)) submit(i, i)
            for (i in seq_len(.num))
            {
                d <- .recvOneResult(cl)
                j <- i + min(.num, p)
                if (j <= .num) submit(d$node, j)

                dv <- d$value
                if (inherits(dv, "try-error"))
                {
                    stop("One of the nodes produced an error: ",
                        as.character(dv))
                }

                if (is.function(.combinefun))
                {
                    if (argone)
                    {
                        .combinefun(dv)
                    } else {
                        if (is.null(ans))
                            ans <- dv
                        else
                            ans <- .combinefun(ans, dv)
                    }
                } else if (.combinefun %in% c("unlist", "list"))
                {
                    # assignment NULL would remove it from the list
                    if (!is.null(dv)) ans[[d$tag]] <- dv
                }

                if (!is.null(.updatefun)) .updatefun(i)
            }
        }
    } else {
        ####  serial implementation
        if (is.function(.combinefun))
        {
            ans <- NULL
            for (i in seq_len(.num))
            {
                dv <- .fun(i, ...)
                if (argone)
                {
                    .combinefun(dv)
                } else {
                    if (is.null(ans))
                        ans <- dv
                    else
                        ans <- .combinefun(ans, dv)
                }
            }
        } else if (identical(.combinefun, "none"))
        {
            for (i in seq_len(.num)) .fun(i, ...)
            ans <- NULL
        } else if (.combinefun %in% c("unlist", "list"))
        {
            ans <- vector("list", .num)
            for (i in seq_len(.num))
            {
                v <- .fun(i, ...)
                # assignment NULL would remove it from the list
                if (!is.null(v)) ans[[i]] <- v
            }
        }
    }

    # output
    ans
}



#######################################################################
# Parallel functions
#   ncore -- the number of cores
#   .num -- the total number of segments
#   .fun -- a user-defined function
#   .combinefun -- a user-defined function for combining the returned values
#   .updatefun -- a user-defined function for updating progress (could be NULL)
#   ... -- other parameters passed to .fun
#

.parse_prepare_cleanup <- quote(parallel:::prepareCleanup())
.mc_prepCleanup <- function() eval(.parse_prepare_cleanup)

.parse_cleanup <- quote(parallel:::cleanup(TRUE))
.mc_cleanup <- function() eval(.parse_cleanup)

.parse_process_id <- quote(parallel:::processID(jobs))
.mc_processID <- function(jobs) eval(.parse_process_id)

.parse_sel_child <- quote(parallel:::selectChildren(children, timeout))
.mc_selectChildren <- function(children, timeout) eval(.parse_sel_child)

.parse_read_child <- quote(parallel:::readChild(child))
.mc_readChild <- function(child) eval(.parse_read_child)

.DynamicForkCall <- function(ncore, .num, .fun, .combinefun, .updatefun, ...)
{
    # in order to use the internal functions accessed by ':::'
    # the functions are all defined in 'parallel/R/unix/mclapply.R'

    # check
    stopifnot(is.numeric(ncore), length(ncore)==1L)
    stopifnot(is.numeric(.num), length(.num)==1L)
    stopifnot(is.function(.fun))
    stopifnot(is.character(.combinefun) | is.function(.combinefun))
    stopifnot(is.null(.updatefun) | is.function(.updatefun))
    argone <- FALSE
    if (is.function(.combinefun))
        argone <- length(formals(args(.combinefun))) == 1L

    # all processes created from now on will be terminated by cleanup
    parallel::mc.reset.stream()
    .mc_prepCleanup()
    on.exit(.mc_cleanup())

    # initialize
    if (identical(.combinefun, "unlist") | identical(.combinefun, "list"))
        ans <- vector("list", .num)
    else
        ans <- NULL

    jobs <- lapply(seq_len(min(.num, ncore)), function(i)
        parallel::mcparallel(.fun(i, ...), name=NULL, mc.set.seed=TRUE, silent=FALSE))
    jobsp <- .mc_processID(jobs)
    jobid <- seq_along(jobsp)
    has.errors <- 0L

    finish <- rep(FALSE, .num)
    nexti <- length(jobid) + 1L
    while (!all(finish))
    {
        s <- .mc_selectChildren(jobs[!is.na(jobsp)], -1)
        if (is.null(s)) break  # no children, should not happen
        if (is.integer(s))
        {
            for (ch in s)
            {
                ji <- match(TRUE, jobsp==ch)
                ci <- jobid[ji]
                r <- .mc_readChild(ch)
                if (is.raw(r))
                {
                    child.res <- unserialize(r)
                    if (inherits(child.res, "try-error"))
                    {
                        has.errors <- has.errors + 1L
                        stop(child.res)
                    }
                    if (is.function(.combinefun))
                    {
                        if (inherits(child.res, "try-error"))
                            stop(child.res)
                        if (argone)
                        {
                            .combinefun(child.res)
                        } else {
                            if (is.null(ans))
                                ans <- child.res
                            else
                                ans <- .combinefun(ans, child.res)
                        }
                    } else if (.combinefun %in% c("unlist", "list"))
                    {
                        # assignment NULL would remove it from the list
                        if (!is.null(child.res)) ans[[ci]] <- child.res
                    }
                    if (!is.null(.updatefun)) .updatefun(ci)
                } else {
                    # the job has finished
                    finish[ci] <- TRUE
                    jobsp[ji] <- jobid[ji] <- NA_integer_
                    # still something to do
                    if (nexti <= .num)
                    {
                        jobid[ji] <- nexti
                        jobs[[ji]] <- parallel::mcparallel(.fun(nexti, ...),
                            name=NULL, mc.set.seed=TRUE, silent=FALSE)
                        jobsp[ji] <- .mc_processID(jobs[[ji]])
                        nexti <- nexti + 1L
                    }
                }
            }
        }
    }

    if (has.errors)
    {
        warning(sprintf("%d function calls resulted in an error", has.errors),
            immediate.=TRUE, domain=NA)
    }

    # output
    ans
}



#######################################################################
# Add a variable with corresponding compression mode
# See: seqStorageOption
#
.AddVar <- function(storage.option, node, varname, val=NULL,
    storage=storage.mode(val), valdim=NULL, closezip=FALSE, visible=TRUE)
{
    # check
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))
    if (inherits(node, "gds.class"))
        node <- node$root
    stopifnot(inherits(node, "gdsn.class"))

    # full variable name
    path <- name.gdsn(node, TRUE)
    idx_flag <- substr(varname, 1L, 1L) == "@"
    if (idx_flag && grepl("/", varname, fixed=TRUE))
        varname <- .var_path(substring(varname, 2L), "@")
    fullvarname <- paste(path, varname, sep="/")
    info_flag <- substr(fullvarname, 1L, 16L) == "annotation/info/"
    fmt_flag  <- substr(fullvarname, 1L, 18L) == "annotation/format/"
    varname2 <- fullvarname
    if (substr(fullvarname, 1L, 12L) == "description/")
    {
        varname2 <- "description"
    } else {
        if (path == "genotype")
            varname2 <- "genotype"
        else if (fmt_flag & !idx_flag)
            varname2 <- path
    }

    # compression flag
    compress.flag <- storage.option$compression
    if (is.null(compress.flag)) compress.flag <- ""
    nm <- names(storage.option)
    if (!is.null(nm))
    {
        i <- match(fullvarname, nm)
        j <- match(varname2, nm)
        if (is.na(i)) i <- j
        if (is.na(i))
        {
            if (path == "genotype")
                compress.flag <- storage.option$geno.compress
            else if (idx_flag)
                compress.flag <- storage.option$index.compress
            else if (info_flag)
                compress.flag <- storage.option$info.compress
            else if (fmt_flag)
                compress.flag <- storage.option$format.compress
        } else {
            compress.flag <- storage.option[[i]]
        }
    }

    # storage mode
    storage.param <- ""
    if (storage == "float")
    {
        stm <- "float32"
        nm <- names(storage.option$float.mode)
        if (is.null(nm))
        {
            if (length(storage.option$float.mode) > 0L)
                stm <- storage.option$float.mode[1L]
        } else {
            i <- match(fullvarname, nm)
            if (is.na(i))
            {
                i <- match("", nm)
                if (!is.na(i))
                    stm <- storage.option$float.mode[i]
            } else
                stm <- storage.option$float.mode[i]
        }
        s <- unlist(strsplit(stm, ":", fixed=TRUE))
        if (length(s) > 2L)
            stop("Invalid 'storage.option$float.mode'.")
        storage <- s[1L]
        if (length(s) == 2L)
            storage.param <- s[2L]
    }

    nm <- names(storage.option$mode)
    if (!is.null(nm))
    {
        i <- match(fullvarname, nm)
        j <- match(varname2, nm)
        if (is.na(i)) i <- j
        if (!is.na(i))
        {
            stm <- storage.option$mode[i]
            s <- unlist(strsplit(stm, ":", fixed=TRUE))
            if (length(s) > 2L)
                stop("Invalid 'storage.option$mode'.")
            storage <- s[1L]
            if (length(s) == 2L)
                storage.param <- s[2L]
            else
                storage.param <- ""
        }
    }

    # sub directory?
    if (grepl("/", varname, fixed=TRUE))
    {
        nm <- unlist(strsplit(varname, "/", fixed=TRUE))
        varname <- nm[length(nm)]
        for (m in nm[-length(nm)])
        {
            v <- index.gdsn(node, m, silent=TRUE)
            if (is.null(v))
                v <- addfolder.gdsn(node, m)
            node <- v
        }
    }

    # add.gdsn arguments
    args <- list(node=node, name=varname, val=val, storage=storage,
        valdim=valdim, compress=compress.flag, closezip=closezip,
        visible=visible)
    args <- c(args, eval(parse(text=paste("list(", storage.param, ")"))))
    do.call(add.gdsn, args)
}



#######################################################################
# Merge attributes
#
.MergeAttr <- function(val)
{
    stopifnot(is.list(val), length(val) > 0L)
    ans <- unlist(val, recursive=FALSE)
    if (!is.null(ans))
    {
        nm <- names(ans)
        if (is.null(nm)) stop("Invalid attributes.")
        if (any(nm == "")) stop("Invalid attributes.")
        while ((k = anyDuplicated(nm)) > 0L)
        {
            i <- which(nm[k] == nm)
            v <- unique(unlist(ans[i]))
            ans[[ i[1L] ]] <- v
            ans <- ans[ - i[-1L] ]
            nm <- nm[ - i[-1L] ]
        }
    }
    ans
}

.MergeNodeAttr <- function(target.node, srcfiles, varname)
{
    v <- vector("list", length(srcfiles))
    for (i in seq_along(srcfiles))
        v[[i]] <- get.attr.gdsn(index.gdsn(srcfiles[[i]], varname))
    v <- .MergeAttr(v)
    if (!is.null(v))
    {
        algo <- c("md5", "sha1", "sha256", "sha384", "sha512")
        algo <- c(algo, paste0(algo, "_r"))
        nm <- names(v)
        for (i in seq_along(v))
        {
            if (!(nm[i] %in% algo))
                put.attr.gdsn(target.node, nm[i], v[[i]])
        }
    }
    invisible()
}



#######################################################################
# Get the unique data frame
#
.DigestCode <- function(node, algo=TRUE, verbose=TRUE)
{
    if (!is.null(node))
    {
        stopifnot(inherits(node, "gdsn.class"))
        if (isTRUE(algo) | is.character(algo))
        {
            if (isTRUE(algo)) algo <- "md5"
            h <- digest.gdsn(node, algo=algo, action="add")
            if (verbose)
            {
                s <- paste0("  [", algo, ": ", h, "]")
                if (.crayon()) s <- crayon::blurred(s)
                cat(s, "\n", sep="")
            }
        } else if (verbose)
            cat("\n")
        flush.console()
    }
    invisible()
}


.DigestFile <- function(gfile, digest=TRUE, verbose=TRUE)
{
    ## digest
    flag <- verbose & (isTRUE(digest) | is.character(digest))
    if (flag) cat("Digests:\n")

    for (nm in c("sample.id", "variant.id", "position", "chromosome", "allele"))
    {
        if (flag) cat("   ", nm)
        .DigestCode(index.gdsn(gfile, nm), digest, verbose)
    }

    n <- index.gdsn(gfile, "genotype/data", silent=TRUE)
    if (!is.null(n))
    {
        if (flag) cat("    genotype")
        .DigestCode(n, digest, verbose)
        .DigestCode(index.gdsn(gfile, "genotype/@data", silent=TRUE), digest, FALSE)
    }

    n <- index.gdsn(gfile, "phase/data", silent=TRUE)
    if (!is.null(n))
    {
        if (flag) cat("    phase")
        .DigestCode(n, digest, verbose)
    }

    if (flag) cat("    annotation/id")
    .DigestCode(index.gdsn(gfile, "annotation/id"), digest, verbose)
    if (flag) cat("    annotation/qual")
    .DigestCode(index.gdsn(gfile, "annotation/qual"), digest, verbose)
    if (flag) cat("    annotation/filter")
    .DigestCode(index.gdsn(gfile, "annotation/filter"), digest, verbose)

    node <- index.gdsn(gfile, "annotation/info")
    for (n in ls.gdsn(node))
    {
        if (flag) cat("    annotation/info/", n, sep="")
        .DigestCode(index.gdsn(node, n), digest, verbose)
        n1 <- index.gdsn(node, paste0("@", n), silent=TRUE)
        if (!is.null(n1))
            .DigestCode(n1, digest, FALSE)
    }

    node <- index.gdsn(gfile, "annotation/format")
    for (n in ls.gdsn(node))
    {
        if (flag) cat("    annotation/format/", n, sep="")
        .DigestCode(index.gdsn(node, paste0(n, "/data"), silent=TRUE), digest, verbose)
        .DigestCode(index.gdsn(node, paste0(n, "/@data")), digest, FALSE)
    }

    node <- index.gdsn(gfile, "sample.annotation")
    for (n in ls.gdsn(node))
    {
        if (flag) cat("    sample.annotation/", n, sep="")
        .DigestCode(index.gdsn(node, n), digest, verbose)
    }

    invisible()
}



#######################################################################
# Convert to a VariantAnnotation object

.seqDebug <- function(gdsfile)
    invisible(.Call("SEQ_Debug", gdsfile))



#######################################################################
# Convert to a VariantAnnotation object

.seqProgress <- function(count, nproc)
    .Call(SEQ_Progress, count, nproc)

.seqProgForward <- function(progress, inc)
    invisible(.Call(SEQ_ProgressAdd, progress, inc))



#######################################################################
# Convert to a VariantAnnotation object
#

.variableLengthToList <- function(x)
{
    xl <- list()
    j <- 1L
    for (i in 1:length(x$length))
    {
        len <- x$length[i]
        if (len > 0L)
        {
            xl[[i]] <- x$data[j:(j+len-1)]
            j <- j+len
        } else {
            xl[[i]] <- NA
        }
    }
    xl
}

.toAtomicList <- function(x, type)
{
    switch(type,
        Character=CharacterList(x),
        String=CharacterList(x),
        Integer=IntegerList(x),
        Float=NumericList(x))
}

.variableLengthToMatrix <- function(x)
{
    xl <- list()
    i <- 1L
    for (j in seq_along(x))
    {
        for (k in seq_len(NROW(x[[j]])))
        {
            xl[[i]] <- x[[j]][k,]
            i <- i + 1L
        }
    }
    matrix(xl, nrow=NROW(x[[1L]]), ncol=length(x))
}
