#######################################################################
#
# Package Name: SeqArray
#
# Description: Internal functions
#



#######################################################################
# Get the numbers of selected samples and variants
#

.seldim <- function(gdsfile)
{
    z <- seqSummary(gdsfile, "genotype", check="none", verbose=FALSE)
    # z$seldim[1] -- # of selected samples
    # z$seldim[2] -- # of selected variants
    z$seldim
}



#######################################################################
# Internal C function
#

.cfunction0 <- function(name)
{
    fn <- function(x) NULL
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
# append the repeated values
#

.repeat_gds <- function(node, elm, count)
{
    val <- rep(elm, 65536L)
    while (count > 0L)
    {
        size <- ifelse(count <= 65536L, count, 65536L)
        append.gdsn(node, val[seq_len(size)], check=FALSE)
        count <- count - size
    }
    invisible()
}



#######################################################################
# Data Type define VCF Format
#

.vcf_type <- function(node)
{
    a <- get.attr.gdsn(node)
    if (is.null(a$Type))
    {
        i <- objdesp.gdsn(node)$type
        if (i %in% c("Raw", "Integer", "Logical"))
            v <- "Integer"
        else if (i %in% "Real")
            v <- "Float"
        else if (i %in% c("Factor", "String"))
            v <- "String"
        else
            v <- "."
    } else
        v <- a$Type
    v
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

    if ((substr(filename, 1L, 6L) == "ftp://") |
        (substr(filename, 1L, 7L) == "http://"))
    {
        if (.last_str(filename, 3L) == ".gz")
        {
            con <- gzcon(url(filename, "rb"))
        } else
            con <- url(filename, "rb")
    } else {
        if (.last_str(filename, 3L) == ".gz")
        {
            con <- gzfile(filename, "rb")
        } else if (.last_str(filename, 3L) == ".xz")
        {
            con <- xzfile(filename, "rb")
        } else 
            con <- file(filename, "rb")
    }

    # open(con)
    list(filename=filename, con=con, con2=con2)
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
# Load Parallel package
#

.loadparallel <- function()
{
    if (!requireNamespace("parallel"))
        stop("The 'parallel' package should be installed.")
    invisible()
}


#######################################################################
# Parallel functions
#

.DynamicClusterCall <- function(cl, .num, .fun, .combinefun,
    .stopcluster, ...)
{
    # in order to use the internal functions accessed by ':::'
    # the functions are all defined in 'parallel/R/snow.R'

    .SendData <- parse(text=
        "parallel:::sendData(con, list(type=type,data=value,tag=tag))")
    .RecvOneData <- parse(text="parallel:::recvOneData(cl)")

    postNode <- function(con, type, value = NULL, tag = NULL)
    {
        eval(.SendData)
    }
    sendCall <- function(con, fun, args, return = TRUE, tag = NULL)
    {
        postNode(con, "EXEC",
            list(fun = fun, args = args, return = return, tag = tag))
        NULL
    }
    recvOneResult <- function(cl)
    {
        v <- eval(.RecvOneData)
        list(value = v$value$value, node = v$node, tag = v$value$tag)
    }


    #################################################################
    # check
    stopifnot(is.null(cl) | inherits(cl, "cluster"))
    stopifnot(is.numeric(.num))
    stopifnot(is.function(.fun))
    stopifnot(is.character(.combinefun) | is.function(.combinefun))
    stopifnot(is.logical(.stopcluster))

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
                sendCall(cl[[node]], .fun, argfun(job), tag = job)

            for (i in 1:min(.num, p)) submit(i, i)
            for (i in seq_len(.num))
            {
                d <- recvOneResult(cl)
                j <- i + min(.num, p)

                stopflag <- FALSE
                if (j <= .num)
                {
                    submit(d$node, j)
                } else {
                    if (.stopcluster)
                    {
                        parallel::stopCluster(cl[d$node])
                        cl <- cl[-d$node]
                        stopflag <- TRUE
                    }
                }

                dv <- d$value
                if (inherits(dv, "try-error"))
                {
                    if (.stopcluster)
                        parallel::stopCluster(cl)
                    stop("One of the nodes produced an error: ", as.character(dv))
                }

                if (is.function(.combinefun))
                {
                    if (is.null(ans))
                        ans <- dv
                    else
                        ans <- .combinefun(ans, dv)
                } else if (.combinefun %in% c("unlist", "list"))
                {
                    ans[[d$tag]] <- dv
                }

                if (stopflag)
                    message(sprintf("Stop \"job %d\".", d$node))
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
                if (is.null(ans))
                    ans <- dv
                else
                    ans <- .combinefun(ans, dv)
            }
        } else if (identical(.combinefun, "none"))
        {
            for (i in seq_len(.num)) .fun(i, ...)
            ans <- NULL
        } else if (.combinefun %in% c("unlist", "list"))
        {
            ans <- vector("list", .num)
            for (i in seq_len(.num))
                ans[[i]] <- .fun(i, ...)
        }
    }

    ans
}
