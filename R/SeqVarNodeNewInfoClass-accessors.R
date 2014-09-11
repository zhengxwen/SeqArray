#######################################################################
#######################################################################
#
#
#

#######################################################################
# Add a new variable to the INFO field
#
seqInfoNewVar <- function(gdsfile, var.name, variant.id, val,
    description="", compress=c("ZIP.MAX", ""), no.data.index=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    stopifnot(is.character(var.name) & is.vector(var.name))
    stopifnot(length(var.name) == 1)

    stopifnot(is.vector(variant.id))
    stopifnot(length(variant.id) > 0)

    if (is.vector(val) | is.factor(val))
    {
        if (length(variant.id) != length(val))
            stop("`val' should have the same length as `variant.id'.")
    } else if (is.matrix(val))
    {
        if (length(variant.id) != ncol(val))
            stop("`val' should have the same number of columns as the length of `variant.id'.")
        stopifnot(nrow(val) > 0)
    } else
        stop("`val' should be a vector or matrix.")

    stopifnot(is.character(description) & is.vector(description))
    stopifnot(length(description) == 1)

    compress <- match.arg(compress)

	stopifnot(is.logical(no.data.index) & is.vector(no.data.index))
    stopifnot(length(no.data.index) == 1)

    # determine the storage mode -- type
    if (is.factor(val))
    {
        type <- "String"
        stm <- "int32"
    } else {
        s <- storage.mode(val)
        if (s == "integer")
        {
            type <- "Integer"
            stm <- "int32"
        } else if (s == "logical")
        {
            type <- "Flag"
            if (!is.vector(val))
                stop("`val' should be a vector if it is `logical'.")
            stm <- "bit1"
        } else if (s %in% c("double", "numeric"))
        {
            type <- "Float"
            stm <- "float64"
        } else if (s == "character")
        {
            type <- "String"
            stm <- "string"
        } else
            stop("`val' should be numeric or character-type.")
    }

    # determine number
    if (is.vector(val) | is.factor(val))
        number <- if (type != "Flag") 1 else 0
    else
        number <- nrow(val)

    ####  align variant.id  ####

    seqSetFilter(gdsfile, action="push", verbose=FALSE)
    on.exit({ seqSetFilter(gdsfile, action="pop", verbose=FALSE) })
    vid <- seqGetData(gdsfile, "variant.id")

    map <- match(variant.id, vid)
    if (any(is.na(map)))
        stop("`variant.id' should be the variant IDs in the specified GDS file.")
    if (is.unsorted(map, strictly=TRUE))
        stop("`variant.id' should have the same order as the IDs in the specified GDS file.")

    if (no.data.index & (length(variant.id) < length(vid)) & !is.matrix(val))
    {
        if (!(is.logical(val) & any(is.na(val))) & !is.character(val))
        {
            val <- val[match(vid, variant.id)]
            if (is.double(val))
                val[is.na(val)] <- NaN
            variant.id <- vid
        }
    }

    ####  add variable(s)  ####

    node <- add.gdsn(index.gdsn(gdsfile, "annotation/info"), var.name,
        val=val, storage=stm, compress=compress, closezip=TRUE)
    put.attr.gdsn(node, "Number", number)
    put.attr.gdsn(node, "Type", type)
    put.attr.gdsn(node, "Description", description)

    if (length(variant.id) < length(vid))
    {
        # need an index variable
        idx <- integer(length(vid))
        idx[map] <- 1L
        node <- add.gdsn(index.gdsn(gdsfile, "annotation/info"),
            paste("@", var.name, sep=""),
            val=idx, compress=compress, closezip=TRUE)
        put.attr.gdsn(node, "R.invisible")
    }

    sync.gds(gdsfile)

    invisible()
}



#######################################################################
# Add a new variable to the INFO field
#
.seqInfoNewVarEx <- function(gdsfile, var.name,
    number=c(".", "A", "G"),
    type=c("Integer", "Float", "Flag", "Character", "String"),
    description="", compress=c("ZIP.MAX", ""))
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name)==1))
    stopifnot(is.character(description))
    type <- match.arg(type)
    compress <- match.arg(compress)

    stopifnot(is.vector(number))
    stopifnot(length(number) == 1)
    if (is.character(number))
    {
    	number <- match.arg(number)
        vdim <- c(0)
        num <- 0
    } else if (is.numeric(number))
    {
        number <- as.integer(number)
        stopifnot(is.finite(number))
        if (number < 1)
            stopifnot("`number' should be greater than zero.")
        if (type == "Flag")
            number <- 0
        num <- number

        if (number > 1)
            vdim <- c(number, 0)
        else
            vdim <- c(0)
    }

    # storage mode
    stm <- c("int32", "float64", "bit1", "string", "string")
    names(stm) <- c("Integer", "Float", "Flag", "Character", "String")

    # add position
    node <- add.gdsn(index.gdsn(gdsfile, "annotation/info"), var.name,
        storage=stm[type], valdim=vdim, compress=compress)
    put.attr.gdsn(node, "Number", number)
    put.attr.gdsn(node, "Type", type)
    put.attr.gdsn(node, "Description", description)


    # output
    new("SeqVarNodeNewInfoClass", gdsn=node, index=NULL,
    	number=num, type=type, ext=NULL)
}



#######################################################################
# 
#
# setMethod("seqAppend", "SeqVarNodeNewInfoClass",
#     function(object, variant.id=NULL, val)
#     {
#     }
# )



#######################################################################
# 
#
# setMethod("seqClose", "SeqVarNodeNewInfoClass",
#     function(object)
#     {
#         print(object)
#     }
# )
