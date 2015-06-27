#######################################################################
#
# Package Name: SeqArray
#
# Description: Internal functions
#



#######################################################################
# Internal C function
#

.cfunction <- function(name)
{
    fn <- function(arg) { NULL }
    f <- quote(.Call(EXTERNALNAME, arg))
    f[[2L]] <- getNativeSymbolInfo(name, "SeqArray")$address
    body(fn) <- f
    fn
}

.cfunction2 <- function(name)
{
    fn <- function(arg1, arg2) { NULL }
    f <- quote(.Call(EXTERNALNAME, arg1, arg2))
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
    stopifnot(is.null(.combinefun) | is.character(.combinefun) |
        is.function(.combinefun))
    stopifnot(is.logical(.stopcluster))

    if (!is.null(cl))
    {
        ############################################################
        # parallel implementation

        if (is.null(.combinefun))
            ans <- vector("list", .num)
        else
            ans <- NULL

        p <- length(cl)
        if ((.num > 0L) && p)
        {
            ## **** this closure is sending to all nodes
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
                    stop("One node produced an error: ", as.character(dv))
                }

                if (!is.character(.combinefun))
                {
                    if (is.null(.combinefun))
                    {
                        ans[[d$tag]] <- dv
                    } else {
                        if (is.null(ans))
                            ans <- dv
                        else
                            ans <- .combinefun(ans, dv)
                    }
                }

                if (stopflag)
                    message(sprintf("Stop \"job %d\".", d$node))
            }
        }
    } else {

        ############################################################
        # serial implementation

        if (is.null(.combinefun))
        {
            ans <- vector("list", .num)
            for (i in seq_len(.num))
                ans[[i]] <- .fun(i, ...)
        } else if (is.character(.combinefun))
        {
            for (i in seq_len(.num)) .fun(i, ...)
            ans <- NULL
        } else {
            ans <- NULL
            for (i in seq_len(.num))
            {
                dv <- .fun(i, ...)
                if (is.null(ans))
                    ans <- dv
                else
                    ans <- .combinefun(ans, dv)
            }
        }
    }

    ans
}
