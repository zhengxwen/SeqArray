#######################################################################
#
# Package Name: SeqArray
#
# The interface to the class "SeqVarGDSClass"
#


#######################################################################
# Open a sequencing-variant GDS file
#
seqOpen <- function(gds.fn, readonly=TRUE)
{
    # check
    stopifnot(is.character(gds.fn) & (length(gds.fn)==1))
    stopifnot(is.logical(readonly))

    # open the file
    ans <- openfn.gds(gds.fn, readonly=readonly)

    # check header
    n <- index.gdsn(ans, "description", silent=TRUE)
    if (is.null(n))
    {
        closefn.gds(ans)
        stop(sprintf("'%s' is not a sequencing-variant GDS file.", gds.fn))
    }

    tmp <- get.attr.gdsn(n)
    if (!("sequence.variant.format" %in% names(tmp)))
    {
        closefn.gds(ans)
        stop(sprintf("'%s' is not a sequencing-variant GDS file.", gds.fn))
    }

    .Call("seq_Open_Init", ans, PACKAGE="SeqArray")

    new("SeqVarGDSClass", ans)
}



#######################################################################
# Close a sequencing-variant GDS file
#
setMethod("seqClose", "SeqVarGDSClass",
    function(object)
    {
        .Call("seq_File_Done", object, PACKAGE="SeqArray")
        closefn.gds(object)
        invisible()
    }
)



#######################################################################
# To set a working space with selected samples and variants
#
seqSetFilter <- function(gdsfile, sample.id=NULL, variant.id=NULL,
    samp.sel=NULL, variant.sel=NULL, action=c("set", "push", "pop"),
    verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(verbose))
    action <- match.arg(action)

    if (action == "push")
    {
        .Call("seq_FilterPush", gdsfile, PACKAGE="SeqArray")
    } else if (action == "pop")
    {
        .Call("seq_FilterPop", gdsfile, PACKAGE="SeqArray")
        return(invisible())
    }

    if (!is.null(sample.id))
    {
        stopifnot(is.vector(sample.id))
        stopifnot(is.numeric(sample.id) | is.character(sample.id))
        .Call("seq_SetSpaceSample", gdsfile, sample.id, verbose,
            PACKAGE="SeqArray")
    } else if (!is.null(samp.sel))
    {
        stopifnot(is.vector(samp.sel) & is.logical(samp.sel))
        .Call("seq_SetSpaceSample", gdsfile, samp.sel, verbose,
            PACKAGE="SeqArray")
    }

    if (!is.null(variant.id))
    {
        stopifnot(is.vector(variant.id))
        stopifnot(is.numeric(variant.id) | is.character(variant.id))
        .Call("seq_SetSpaceVariant", gdsfile, variant.id, verbose,
            PACKAGE="SeqArray")
    } else if (!is.null(variant.sel))
    {
        stopifnot(is.vector(variant.sel) & is.logical(variant.sel))
        .Call("seq_SetSpaceVariant", gdsfile, variant.sel, verbose,
            PACKAGE="SeqArray")
    } else {
        if (is.null(sample.id) & is.null(samp.sel))
        {
            .Call("seq_SetSpaceSample", gdsfile, NULL, verbose,
                PACKAGE="SeqArray")
            .Call("seq_SetSpaceVariant", gdsfile, NULL, verbose,
                PACKAGE="SeqArray")
        }
    }

    invisible()
}



#######################################################################
# To get a working space
#
seqGetFilter <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))

    .Call("seq_GetSpace", gdsfile, PACKAGE="SeqArray")
}



#######################################################################
# Get data from a working space with selected samples and variants
#
seqGetData <- function(gdsfile, var.name)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name)==1))

    .Call("seq_GetData", gdsfile, var.name, PACKAGE="SeqArray")
}



#######################################################################
# Apply functions over margins on a working space with selected samples and variants
#
seqApply <- function(gdsfile, var.name, FUN,
    margin = c("by.variant"),
    as.is = c("list", "integer", "double", "character", "none"),
    var.index = c("none", "relative", "absolute"), ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name) > 0))
    var.index <- match.arg(var.index)
    var.index <- match(var.index, c("none", "relative", "absolute"))

    FUN <- match.fun(FUN)
    margin <- match.arg(margin)
    as.is <- match.arg(as.is)

    if (margin == "by.variant")
    {
        # C call
        rv <- .Call("seq_Apply_Variant", gdsfile, var.name, FUN, as.is,
            var.index, new.env(), PACKAGE="SeqArray")
        if (as.is == "none") return(invisible())
    }
    rv
}



#######################################################################
# Apply functions via a sliding window over variants
#
seqSlidingWindow <- function(gdsfile, var.name, win.size, shift=1, FUN,
    as.is = c("list", "integer", "double", "character", "none"),
    var.index = c("none", "relative", "absolute"), ...)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & (length(var.name) > 0))

    stopifnot(is.numeric(win.size) & (length(win.size)==1))
    win.size <- as.integer(win.size)
    stopifnot(is.finite(win.size))
    if (win.size < 1)
        stop("`win.size' should be greater than 0.")

    stopifnot(is.numeric(shift) & (length(shift)==1))
    shift <- as.integer(shift)
    stopifnot(is.finite(shift))
    if (shift < 1)
        stop("`shift' should be greater than 0.")

    as.is <- match.arg(as.is)
    var.index <- match.arg(var.index)
    var.index <- match(var.index, c("none", "relative", "absolute"))

    FUN <- match.fun(FUN)

    # C call
    rv <- .Call("seq_SlidingWindow", gdsfile, var.name, win.size, shift,
    	FUN, as.is, var.index, new.env(), PACKAGE="SeqArray")

    if (as.is == "none") return(invisible())
    rv
}



#######################################################################
# Summarize the GDS file
#
seqSummary <- function(gdsfile, varname=NULL,
    check=c("check", "full.check", "none"), verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(varname) | is.character(varname))
    if (is.character(varname))
        stopifnot(length(varname) == 1)
    check <- match.arg(check)

    # initialize a GDS object
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1)
        gds <- seqOpen(gdsfile)
        on.exit(seqClose(gds))   
    } else {
        gds <- gdsfile
    }


    # whether specify the variable name or not
    if (is.null(varname))
    {
        ########################################################
        check.dim <- function(node, dim.len)
        {
            dm <- objdesp.gdsn(node)$dim
            if (check != "")
            {
                if (length(dm) != dim.len)
                {
                    print(node)
                    stop("Error in dimension!")
                }
            }
            dm
        }


        ########################################################
        # initialize
        ans <- list()
        ans$filename <- gds$filename
        if (verbose)
            cat("File: ", ans$filename, "\n", sep="")


        ########################################################
        # check header
        n <- index.gdsn(gds, "description")
        tmp <- get.attr.gdsn(n)
        if ("sequence.variant.format" %in% names(tmp))
        {
            ans$sequence.variant.format <- tmp$sequence.variant.format
            if (verbose)
            {
                cat("Sequence Variant Format: ", ans$sequence.variant.format,
                    "\n", sep="")
            }
        }

        ########################################################
        # the number of samples
        n <- index.gdsn(gds, "sample.id")
        dm <- check.dim(n, 1)
        n.samp <- dm[1]
        ans$num.of.sample <- n.samp

        ########################################################
        # the number of variants
        n <- index.gdsn(gds, "variant.id")
        dm <- check.dim(n, 1)
        n.variant <- dm[1]
        ans$num.of.variant <- n.variant

        if (verbose)
        {
            cat("The number of samples: ", n.samp, "\n", sep="")
            cat("The number of variants: ", n.variant, "\n", sep="")
        }


        ########################################################
        # position
        n <- index.gdsn(gds, "position")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'position' variable.")

        ########################################################
        # chromosome
        n <- index.gdsn(gds, "chromosome")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'chromosome' variable.")
        if (verbose)
        {
            chr <- seqGetData(gds, "chromosome")
            tab <- table(chr, exclude=NULL)
            names(dimnames(tab)) <- "The chromosomes:"
            print(tab)
            rm(chr)
        }

        ########################################################
        # allele
        n <- index.gdsn(gds, "allele")
        dm <- check.dim(n, 1)
        if (dm != n.variant)
            stop("Error: the length of the 'allele' variable.")
        if ((check != "") | verbose)
            nallele <- .Call("seq_NumOfAllele", n, PACKAGE="SeqArray")
        if (verbose)
        {
            tab <- table(nallele)
            names(dimnames(tab)) <- "The number of alleles per site:"
            print(tab)
        }


        ########################################################
        # genotype
        n <- index.gdsn(gds, "genotype/data")
        dm <- check.dim(n, 3)
        if (dm[2] != n.samp)
        {
            print(n)
            stop("Invalid number of samples.")
        }


        ########################################################
        # phase
        n <- index.gdsn(gds, "phase/data")
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) > 2)
            dm <- dm[-c(length(dm)-1, length(dm))]
        if (dm[1] != n.samp)
        {
            print(n)
            stop("Invalid number of samples.")
        }
        if (dm[2] != n.variant)
        {
            print(n)
            stop("Invalid number of variants.")
        }


        ########################################################
        # annotation/id
        n <- index.gdsn(gds, "annotation/id", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }

        ########################################################
        # annotation/qual
        n <- index.gdsn(gds, "annotation/qual", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }

        ########################################################
        # annotation/filter
        n <- index.gdsn(gds, "annotation/filter", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check.dim(n, 1)
            if (dm != n.variant)
            {
                print(n)
                stop("Invalid number of variants.")
            }
        }


        ########################################################
        # annotation/info
        n <- index.gdsn(gds, "annotation/info", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            vars <- unlist(strsplit(vars, "@", fixed=TRUE))
            vars <- unique(vars[vars != ""])
            dat <- NULL
            if (verbose)
                cat("Annotation, information variables:\n")
            for (nm in vars)
            {
                n <- index.gdsn(gds, paste("annotation/info/", nm, sep=""))
                a <- get.attr.gdsn(n)
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ", a$Description,
                        "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name = nm,
                    number = a$Number, type = a$Type, description = a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$info <- dat
        }


        ########################################################
        # annotation/format
        n <- index.gdsn(gds, "annotation/format", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- NULL
            if (verbose)
                cat("Annotation, format variables:\n")
            for (nm in vars)
            {
                n <- index.gdsn(gds, paste("annotation/format/", nm, sep=""))
                a <- get.attr.gdsn(n)
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ", a$Description,
                        "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name = nm,
                    number = a$Number, type = a$Type, description = a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$format <- dat
        }

        # output
        return(invisible(ans))

    } else {

        # get a description of variable
        ans <- .Call("seq_VarSummary", gds, varname, PACKAGE="SeqArray")
        ans
    }
}



#######################################################################
# Delete data variables
#

seqDelete <- function(gdsfile, info.varname=NULL, format.varname=NULL)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    if (gdsfile$readonly)
        stop("The GDS file is read-only.")
    stopifnot(is.character(info.varname) | is.null(info.varname))
    stopifnot(is.character(format.varname) | is.null(format.varname))
    if (is.null(info.varname) & is.null(format.varname))
        stop("There is no variable.")

    # call C function
    .Call("seq_Delete", gdsfile, info.varname, format.varname,
    	PACKAGE="SeqArray")

    # return
    invisible()
}




#######################################################################
# Merge multiple GDS files
#

# seqMissing <- function(gds.obj, by.samp=TRUE, by.snp=TRUE)
# {
#   # check
#   stopifnot(inherits(gds.obj, "SeqVarGDSClass"))
#   stopifnot(is.logical(by.samp) & (length(by.samp)==1))
#   stopifnot(is.logical(by.snp) & (length(by.snp)==1))
# 
#   rv <- list()
#   if (by.snp)
#   {
#       rv[["by.snp"]] <- apply.gdsn(index.gdsn(gds.obj, "genotype/data"),
#           margin=3, as.is = "double", FUN = function(g)
#               { .Call("seq_missing_snp", g, PACKAGE="SeqArray") }
#       )
#   }
#   if (by.samp)
#   {
#       node <- index.gdsn(gds.obj, "genotype/data")
#       m <- objdesp.gdsn(node)$dim[2]
#       n <- integer(m)
#       apply.gdsn(node, margin=3, as.is = "none", FUN = function(g)
#           { .Call("seq_missing_samp", g, n, PACKAGE="SeqArray") })
#       rv[["by.samp"]] <- n / m
#   }
# 
#   rv
# }

# seqAlleleFreq <- function(gds.obj)
# {
#   # check
#   stopifnot(inherits(gds.obj, "SeqVarGDSClass"))
# 
#   apply.gdsn(index.gdsn(gds.obj, "genotype/data"),
#       margin=3, as.is = "double", FUN = function(g)
#           { .Call("seq_allele_freq", g, PACKAGE="SeqArray") })
# }




#######################################################################
# Modify the SeqArray object
#######################################################################

#######################################################################
# Transpose data variable(s)
#

seqTranspose <- function(gdsfile, var.name, compress=NULL, verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(var.name) & is.vector(var.name))
    stopifnot(length(var.name) == 1)

    node <- index.gdsn(gdsfile, var.name)
    desp <- objdesp.gdsn(node)
    dm <- desp$dim
    if (length(dm) > 1)
    {
        # dimension
        dm <- c(dm[-(length(dm)-1)], 0)
        # folder
        index <- unlist(strsplit(var.name, "/"))
        if (length(index) <= 1)
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
        apply.gdsn(node, margin=length(dm)-1, as.is="none", FUN=function(g) {
            append.gdsn(newnode, g)
        })

        readmode.gdsn(newnode)
    } else
        warning("It is a vector.")

    invisible()
}
