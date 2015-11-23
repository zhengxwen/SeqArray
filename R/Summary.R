#######################################################################
#
# Package Name: SeqArray
#


#######################################################################
# Summarize the GDS file
#

.summary_filter <- function(gdsfile, action, verbose)
{
    at <- get.attr.gdsn(index.gdsn(gdsfile, "annotation/filter"))
    id <- at$R.levels
    dp <- at$Description
    if (is.null(dp))
    {
        # old version
        n <- index.gdsn(gdsfile, "description/vcf.filter", silent=TRUE)
        if (!is.null(n))
        {
            at <- read.gdsn(n)
            dp <- at$Description[match(id, at$ID)]
        } else
            dp <- rep(NA_character_, length(id))
    }
    ans <- list(id=id, description=dp)
    if (action %in% c("check", "full.check"))
    {
        ans$tab <- table(seqGetData(gdsfile, "annotation/filter"))
    }
    ans
}

.summary_info <- function(gdsfile, action, verbose)
{
    info <- index.gdsn(gdsfile, "annotation/info")
    ans <- NULL
    for (nm in ls.gdsn(info))
    {
        a <- get.attr.gdsn(index.gdsn(info, nm))
        if (is.null(a$Number)) a$Number <- "."
        if (is.null(a$Type)) a$Type <- .vcf_type(index.gdsn(info, nm))
        if (is.null(a$Description)) a$Description <- ""
        if (is.null(a$Source)) a$Source <- NA_character_
        if (is.null(a$Version)) a$Version <- NA_character_

        d <- data.frame(ID=nm, Number=a$Number, Type=a$Type,
            Description=a$Description, Source=a$Source, Version=a$Version,
            stringsAsFactors=FALSE)
        ans <- rbind(ans, d)
    }
    if (is.null(ans))
    {
        s <- character()
        ans <- data.frame(ID=s, Number=s, Type=s, Description=s,
            Source=s, Version=s, stringsAsFactors=FALSE)
    }
    if (verbose) print(ans)
    invisible(ans)
}

.summary_format <- function(gdsfile, action, verbose)
{
    # genotype
    a <- get.attr.gdsn(index.gdsn(gdsfile, "genotype"))
    if (is.null(a$VariableName)) a$VariableName <- "GT"
    if (is.null(a$Description)) a$Description <- ""
    ans <- data.frame(ID=a$VariableName, Number=1, Type="String",
        Description=a$Description, stringsAsFactors=FALSE)

    # the FORMAT field
    fmt <- index.gdsn(gdsfile, "annotation/format")
    for (nm in ls.gdsn(fmt))
    {
        a <- get.attr.gdsn(index.gdsn(fmt, nm))
        if (is.null(a$Number)) a$Number <- "."
        if (is.null(a$Type)) a$Type <- .vcf_type(index.gdsn(info, nm))
        if (is.null(a$Description)) a$Description <- ""

        d <- data.frame(ID=nm, Number=a$Number, Type=a$Type,
            Description=a$Description, stringsAsFactors=FALSE)
        ans <- rbind(ans, d)
    }

    if (verbose) print(ans)
    invisible(ans)
}


seqSummary <- function(gdsfile, varname=NULL,
    check=c("check", "full.check", "none"), verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.null(varname) | is.character(varname))
    if (is.character(varname))
        stopifnot(length(varname) == 1L)
    check <- match.arg(check)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # initialize a GDS object
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile) == 1L)
        gds <- seqOpen(gdsfile)
        on.exit(seqClose(gds))   
    } else {
        gds <- gdsfile
    }


    # whether specify the variable name or not
    if (is.null(varname))
    {
        ########
        check_dim <- function(node, dimlen)
        {
            dm <- objdesp.gdsn(node)$dim
            if (check != "none")
            {
                if (length(dm) != dimlen)
                {
                    print(node)
                    stop("Error in dimension!")
                }
            }
            dm
        }

        STOP <- function(...)
        {
            if (check != "none") stop(...)
        }
        STOP2 <- function(node, ...)
        {
            if (check != "none")
                { print(node); stop(...) }
        }


        ########
        # initialize
        ans <- list()
        ans$filename <- gds$filename
        if (verbose)
            cat("File: ", ans$filename, "\n", sep="")


        ########
        # check header
        tmp <- get.attr.gdsn(gds$root)
        if ("FileVersion" %in% names(tmp))
        {
            ans$format.version <- tmp$FileVersion
            if (verbose)
                cat("Format Version: ", ans$format.version, "\n", sep="")
        }

        ########
        n <- index.gdsn(gds, "description/reference", silent=TRUE)
        if (is.null(n))
            ans$reference <- character()
        else
            ans$reference <- as.character(read.gdsn(n))
        n <- index.gdsn(gds, "description/vcf.header", silent=TRUE)
        if (!is.null(n))
        {
            tmp <- read.gdsn(n)
            if (is.data.frame(tmp))
            {
                ans$reference <- c(ans$reference,
                    tmp$value[tmp$id == "reference"])
            } else
                STOP("Invalid 'description/vcf.header'.")
        }
        ans$reference <- unique(ans$reference)
        if (verbose)
        {
            if (length(ans$reference) > 0)
            {
                cat("Reference: ", paste(ans$reference, collapse=", "),
                    "\n", sep="")
            } else
                cat("Reference: unknown\n")
        }

        ########
        n <- index.gdsn(gds, "genotype/data")
        ans$ploidy <- objdesp.gdsn(n)$dim[1L]

        ########
        # the number of samples
        n <- index.gdsn(gds, "sample.id")
        dm <- check_dim(n, 1L)
        n.samp <- dm[1L]
        ans$num.of.sample <- n.samp

        ########
        # the number of variants
        n <- index.gdsn(gds, "variant.id")
        dm <- check_dim(n, 1L)
        n.variant <- dm[1L]
        ans$num.of.variant <- n.variant

        if (verbose)
        {
            cat("Ploidy: ", ans$ploidy, "\n", sep="")
            cat("Number of samples: ", n.samp, "\n", sep="")
            cat("Number of variants: ", n.variant, "\n", sep="")
        }


        ########
        # position
        n <- index.gdsn(gds, "position")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'position' variable.")

        ########
        # chromosome
        n <- index.gdsn(gds, "chromosome")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'chromosome' variable.")
        if (verbose & (check != "none"))
        {
            chr <- seqGetData(gds, "chromosome")
            tab <- table(chr, exclude=NULL)
            names(dimnames(tab)) <- "Chromosomes:"
            print(tab)
            rm(chr)
        }

        ########
        # allele
        n <- index.gdsn(gds, "allele")
        dm <- check_dim(n, 1L)
        if (dm != n.variant)
            STOP("Error: the length of the 'allele' variable.")
        if (check != "none")
        {
            nallele <- .Call(SEQ_NumOfAllele, n)
            if (verbose)
            {
                tab <- table(nallele)
                names(dimnames(tab)) <- "Number of alleles per site:"
                print(tab)
            }
        }


        ########
        # genotype
        n <- index.gdsn(gds, "genotype/data")
        dm <- check_dim(n, 3L)
        if (dm[2L] != n.samp)
            STOP2(n, "Invalid dimension of sample in 'genotype/data'.")
        if (dm[3L] < n.variant)
            STOP2(n, "Invalid dimension of variant in 'genotype/data'.")

        n <- index.gdsn(gds, "genotype/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 3L)
            if (dm[2L] < n.variant)
                STOP(n, "Invalid dimension of variant in 'genotype/~data'.")
            if (dm[3L] != n.samp)
                STOP(n, "Invalid dimension of sample in 'genotype/~data'.")
        }


        ########
        # phase
        n <- index.gdsn(gds, "phase/data")
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) > 2L)
            dm <- dm[-c(length(dm)-1L, length(dm))]
        if (dm[1L] != n.samp)
            STOP2(n, "Invalid dimension of sample in 'phase/data'.")
        if (dm[2L] != n.variant)
            STOP2(n, "Invalid dimension of variant in 'phase/data'.")

        n <- index.gdsn(gds, "phase/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- objdesp.gdsn(n)$dim
            if (length(dm) > 2L)
                dm <- dm[-c(length(dm)-1L, length(dm))]
            if (dm[1L] != n.variant)
                STOP2(n, "Invalid dimension of variant in 'phase/~data'.")
            if (dm[2L] != n.samp)
                STOP2(n, "Invalid dimension of sample in 'phase/~data'.")
        }


        ########
        # annotation/id
        n <- index.gdsn(gds, "annotation/id", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }

        ########
        # annotation/qual
        n <- index.gdsn(gds, "annotation/qual", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }

        ########
        # annotation/filter
        n <- index.gdsn(gds, "annotation/filter", silent=TRUE)
        if (!is.null(n))
        {
            dm <- check_dim(n, 1L)
            if (dm != n.variant)
                STOP2(n, "Invalid number of variants.")
        }


        ########
        # annotation/info
        n <- index.gdsn(gds, "annotation/info", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            vars <- unlist(strsplit(vars, "@", fixed=TRUE))
            vars <- unique(vars[vars != ""])
            dat <- NULL
            if (verbose)
                cat("Annotation, INFO variable(s):\n")
            for (nm in vars)
            {
                n1 <- index.gdsn(n, nm)
                d <- objdesp.gdsn(n1)
                a <- get.attr.gdsn(n1)
                n2 <- index.gdsn(n, paste("@", nm, sep=""), silent=TRUE)
                if (is.null(n2))
                {
                    if (d$dim[length(d$dim)] != n.variant)
                        STOP2(n1, "Invalid dimension of variant.")
                } else {
                    dm <- check_dim(n2, 1L)
                    if (dm != n.variant)
                        STOP2(n1, "Invalid dimension of variant.")
                }
                if (is.null(a$Number))
                {
                    if (is.null(n2))
                        a$Number <- prod(d$dim[-length(d$dim)])
                    else
                        a$Number <- "A"
                }
                if (is.null(a$Type))
                {
                    a$Type <- .vcf_type(n1)
                }
                if (is.null(a$Description))
                {
                    a$Description <- "."
                }
                if (verbose)
                {
                    cat("\t", nm, ", ", a$Number, ", ", a$Type, ", ",
                        a$Description, "\n", sep="")
                }
                dat <- rbind(dat, data.frame(var.name=nm,
                    number=a$Number, type=a$Type, description=a$Description,
                    stringsAsFactors=FALSE))
            }
            ans$info <- dat
        }


        ########
        # annotation/format
        n <- index.gdsn(gds, "annotation/format", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- NULL
            if (verbose)
                cat("Annotation, FORMAT variable(s):\n")
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


        ########
        # sample.annotation
        n <- index.gdsn(gds, "sample.annotation", silent=TRUE)
        if (!is.null(n))
        {
            vars <- ls.gdsn(n)
            dat <- data.frame(var.name=vars, type=rep(".", length(vars)),
                stringsAsFactors=FALSE)
            for (i in seq_len(length(vars)))
                dat$type[i] <- .vcf_type(index.gdsn(n, vars[i]))
            if (verbose)
            {
                cat("Annotation, sample variable(s):\n")
                cat(paste("\t", vars, sep=""), sep="\n")
            }
            ans$sample.annot <- dat
        }

        # output
        invisible(ans)

    } else {
        if (varname == "annotation/filter")
        {
            .summary_filter(gds, check, verbose)
        } else if (varname == "annotation/info")
        {
            .summary_info(gds, check, verbose)
        } else if (varname == "annotation/format")
        {
            .summary_format(gds, check, verbose)
        } else {
            # get a description of variable
            .Call(SEQ_Summary, gds, varname)
        }
    }
}
