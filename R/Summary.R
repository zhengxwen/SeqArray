#######################################################################
#
# Package Name: SeqArray
#
# Description: Get the summary information
#


#######################################################################
# Summarize the GDS file
#

.str <- function(s)
{
    mode(s) <- "character"
    s[is.na(s)] <- "<NA>"
    s
}

.check_dim <- function(node, dimlen)
{
    dm <- objdesp.gdsn(node)$dim
    if (length(dm) != dimlen)
        stop(name.gdsn(node, fullname=TRUE), ": dimension error.")
    dm
}

.check_digest <- function(gdsfile, varname, verbose)
{
    ans <- NA
    n <- index.gdsn(gdsfile, varname, silent=TRUE)
    if (!is.null(n))
    {
        if (verbose)
        {
            cat("    ", name.gdsn(n, fullname=TRUE), ":", sep="")
            flush.console()
        }
        s <- ""
        v <- digest.gdsn(n, action="return")
        s1 <- paste(names(v)[which(v==TRUE)], collapse="/")
        s2 <- paste(names(v)[which(v==FALSE)], collapse="/")
        if (s1=="" & s2=="")
        {
            s <- paste0(s, "\tno digest")
        } else {
            if (s1 != "")
            {
                s <- paste0(s, "\t", sQuote(s1), " [OK]")
                ans <- TRUE
            }
            if (s2 != "")
            {
                s <- paste0(s, ifelse(s1!="", ", ", ""), sQuote(s2), " fails")
                ans <- FALSE
            }
        }
        s <- paste0(s, "\n")
        if (verbose)
        {
            if (.crayon())
            {
                if (identical(ans, NA))
                    s <- crayon::blurred(s)
                else if (identical(ans, FALSE))
                    s <- crayon::bold(s)
            }
            cat(s)
        }
    }
    ans
}

.summary_geno <- function(gdsfile, check, verbose, showsel=TRUE)
{
    ans <- .Call(SEQ_Summary, gdsfile, "genotype")
    if (check != "none")
    {
        nsamp <- ans$dim[2L]
        nvar  <- ans$dim[3L]

        n <- index.gdsn(gdsfile, "genotype/data")
        dm <- .check_dim(n, 3L)
        if (dm[2L] != nsamp)
            message("Invalid sample dimension in 'genotype/data'.")
        if (dm[3L] < nvar)
            message("Invalid variant dimension in 'genotype/data'.")

        n <- index.gdsn(gdsfile, "genotype/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- .check_dim(n, 3L)
            if (dm[2L] < nvar)
                message("Invalid variant dimension in 'genotype/~data'.")
            if (dm[3L] != nsamp)
                message("Invalid sample dimension in 'genotype/~data'.")
        }
    }
    if (verbose)
    {
        cat(paste(c("Ploidy:", "Number of samples:", "Number of variants:"),
            .pretty(ans$dim)), sep="\n")
        if (showsel)
        {
            cat(paste(c("Number of selected samples:",
                "Number of selected variants:"), .pretty(ans$seldim[-1L])),
                sep="\n")
        }
    }
    if (check == "full")
    {
        if (verbose) cat("Genotypes:\n")
        ans$digest <- c(
            "data"  = .check_digest(gdsfile, "genotype/data", verbose),
            "~data" = .check_digest(gdsfile, "genotype/~data", verbose),
            "@data" = .check_digest(gdsfile, "genotype/@data", verbose))
    }
    invisible(ans)
}

.summary_phase <- function(gdsfile, check, verbose)
{
    ans <- .Call(SEQ_Summary, gdsfile, "phase")
    if (verbose)
    {
        cat(paste(c("Ploidy:", "Number of samples:", "Number of variants:",
            "Number of selected samples:", "Number of selected variants:"),
            c(ans$dim, ans$seldim)), sep="\n")
    }

    if (check != "none")
    {
        nsamp <- ans$dim[2L]
        nvar  <- ans$dim[3L]

        n <- index.gdsn(gdsfile, "phase/data")
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) > 2L)
            dm <- dm[-c(length(dm)-1L, length(dm))]
        if (dm[1L] != nsamp)
            message("Invalid sample dimension in 'phase/data'.")
        if (dm[2L] != nvar)
            message("Invalid variant dimension in 'phase/data'.")

        n <- index.gdsn(gdsfile, "phase/~data", silent=TRUE)
        if (!is.null(n))
        {
            dm <- objdesp.gdsn(n)$dim
            if (length(dm) > 2L)
                dm <- dm[-c(length(dm)-1L, length(dm))]
            if (dm[1L] != nvar)
                message("Invalid variant dimension in 'phase/~data'.")
            if (dm[2L] != nsamp)
                message("Invalid sample dimension in 'phase/~data'.")
        }

        if (check == "full")
        {
            if (verbose) cat("Phase:\n")
            ans$digest <- c(
                "data"  = .check_digest(gdsfile, "phase/data", verbose),
                "~data" = .check_digest(gdsfile, "phase/~data", verbose))
        }
    }

    invisible(ans)
}

.summary_sample_id <- function(gdsfile, check, verbose)
{
    ans <- .Call(SEQ_Summary, gdsfile, "sample.id")
    if (check == "full")
    {
        ss <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
        if (is.vector(ss))
        {
            i <- anyDuplicated(ss)
            if (i > 0L)
                message("Duplicated sample id, e.g., ", ss[i])
            else
                cat("Sample ID: no duplicate.\n")
        } else
            message("Invalid dimension of 'sample.id'.")
        if (check == "full")
            attr(ans, "digest") <- .check_digest(gdsfile, "sample.id", verbose)
    }
    ans
}

.summary_variant_id <- function(gdsfile, check, verbose)
{
    ans <- .Call(SEQ_Summary, gdsfile, "variant.id")
    if (check == "full")
    {
        ss <- read.gdsn(index.gdsn(gdsfile, "variant.id"))
        if (is.vector(ss))
        {
            i <- anyDuplicated(ss)
            if (i > 0L)
                message("Duplicated variant id, e.g., ", ss[i])
            else
                cat("Variant ID: no duplicate.\n")
        } else
            message("Invalid dimension of 'variant.id'.")
        if (check == "full")
            attr(ans, "digest") <- .check_digest(gdsfile, "variant.id", verbose)
    }
    ans
}

.summary_position <- function(gdsfile, check, verbose)
{
    ans <- .Call(SEQ_Summary, gdsfile, "position")
    if (check != "none")
    {
        n <- index.gdsn(gdsfile, "position")
        dm <- .check_dim(n, 1L)
        if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
            message("Invalid dimension of 'position'.")
        if (check == "full")
            attr(ans, "digest") <- .check_digest(gdsfile, "position", verbose)
    }
    ans
}

.summary_chrom <- function(gdsfile, check, verbose)
{
    if (check != "none")
    {
        n <- index.gdsn(gdsfile, "chromosome")
        dm <- .check_dim(n, 1L)
        if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
            message("Invalid dimension of 'chromosome'.")
        if (verbose)
            cat("Chromosomes:\n")

        ans <- table(seqGetData(gdsfile, "chromosome"))
        class(ans) <- NULL
        i <- suppressWarnings(as.integer(names(ans)))
        ans <- ans[order(i)]

        if (verbose)
        {
            s <- format(paste("Chr", format(names(ans)), ": ", ans, sep=""))
            i <- 1L
            while (i <= length(s))
            {
                k <- i + 6L
                if (k > length(s)) k <- length(s) + 1L
                cat("    ", toString(s[seq.int(i, k-1L)]), "\n", sep="")
                i <- k
            }
        }
        if (check == "full")
            attr(ans, "digest") <- .check_digest(gdsfile, "chromosome", verbose)
        invisible(ans)
    } else
        .Call(SEQ_Summary, gdsfile, "chromosome")
}

.summary_allele <- function(gdsfile, check, verbose)
{
    if (verbose) cat("Alleles:\n")

    dm <- .check_dim(index.gdsn(gdsfile, "allele"), 1L)
    if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
        message("Invalid dimension of 'allele'.")
    n <- index.gdsn(gdsfile, "description/vcf.alt", silent=TRUE)
    if (!is.null(n))
        ans <- read.gdsn(n)
    else
        ans <- data.frame(ID=character(), Description=character())

    if (verbose)
    {
        if (nrow(ans) > 0L)
            cat(paste0("    ", ans$ID, ", ", ans$Description, "\n"), sep="")
        else
            cat("    ALT: <None>\n")
    }
    if (check != "none")
    {
        ans <- list(value=ans)
        ans$table <- table(seqNumAllele(gdsfile))
        if (verbose)
        {
            s <- sprintf("%s, %d(%.1f%%)", names(ans$table), ans$table,
                prop.table(ans$table)*100.0)
            cat("    tabulation: ", paste(s, collapse="; "), "\n", sep="")
        }
        if (check == "full")
            ans$digest <- .check_digest(gdsfile, "chromosome", verbose)
    }

    invisible(ans)
}

.summary_annot_id <- function(gdsfile, check, verbose)
{
    ans <- NULL
    vn <- "annotation/id"
    n <- index.gdsn(gdsfile, vn, silent=TRUE)
    if (!is.null(n))
    {
        ans <- .Call(SEQ_Summary, gdsfile, vn)
        if (check != "none")
        {
            dm <- .check_dim(n, 1L)
            if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
                message("Invalid dimension of 'annotation/id'.")
            if (check == "full")
            {
                if (verbose) cat("Annotation, ID:\n")
                ans <- .check_digest(gdsfile, vn, verbose)
                names(ans) <- "digest"
            }
        }
    }
    invisible(ans)
}

.summary_annot_qual <- function(gdsfile, check, verbose)
{
    ans <- NULL
    vn <- "annotation/qual"
    n <- index.gdsn(gdsfile, vn, silent=TRUE)
    if (!is.null(n))
    {
        ans <- .Call(SEQ_Summary, gdsfile, vn)
        if (check != "none")
        {
            dm <- .check_dim(n, 1L)
            if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
                message("Invalid dimension of 'annotation/qual'.")
            if (verbose)
                cat("Annotation, Quality:\n")
            ans <- summary(seqGetData(gdsfile, vn))
            class(ans) <- NULL
            if (verbose)
            {
                s <- names(ans)
                s[1L] <- "Min"; s[2L] <- "1st Qu"
                s[5L] <- "3rd Qu"; s[6L] <- "Max"
                cat("    ", paste(paste(s, ans, sep=": "), collapse=", "),
                    "\n", sep="")
            }
            if (check == "full")
                attr(ans, "digest") <- .check_digest(gdsfile, vn, verbose)
        }
    }
    invisible(ans)
}

.summary_filter <- function(gdsfile, check, verbose)
{
    if (verbose) cat("Annotation, FILTER:\n")

    n <- index.gdsn(gdsfile, "annotation/filter")
    if (check != "none")
    {
        dm <- .check_dim(n, 1L)
        if (dm != .Call(SEQ_Summary, gdsfile, "variant.id"))
            message("Invalid dimension of 'annotation/filter'.")
    }

    at <- get.attr.gdsn(n)
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
    ans <- data.frame(ID=id, Description=dp, stringsAsFactors=FALSE)

    if (check != "none")
    {
        ans <- list(value=ans)
        v <- seqGetData(gdsfile, "annotation/filter")
        a <- table(v)
        n <- length(v); remove(v)
        a <- unname(a[match(ans$value$ID, names(a))])
        a[is.na(a)] <- 0L
        ans$table <- a
        if (verbose)
        {
            if (nrow(ans$value) > 0L)
            {
                cat(paste0("    ", ans$value$ID, ", ",
                    .str(ans$value$Description), ", ", ans$table,
                    sprintf("(%0.1f%%)", a/n*100), "\n"),
                    sep="")
            } else
                cat("    <None>\n")
        }
    } else {
        if (verbose)
        {
            if (nrow(ans) > 0L)
            {
                cat(paste0("    ", ans$ID, ", ", .str(ans$Description), "\n"),
                    sep="")
            } else
                cat("    <None>\n")
        }
    }
    invisible(ans)
}

.summary_info <- function(gdsfile, check, verbose)
{
    if (verbose)
        cat("Annotation, INFO variable(s):\n")

    info <- index.gdsn(gdsfile, "annotation/info")
    ans <- NULL
    for (nm in ls.gdsn(info))
    {
        a <- get.attr.gdsn(index.gdsn(info, nm))
        if (is.null(a$Number)) a$Number <- "A"
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

    if (verbose)
    {
        if (nrow(ans) > 0L)
        {
            ss <- paste0("    ", ans$ID, ", ", ans$Number, ", ", ans$Type,
                ", ", ans$Description)
            for (i in seq_along(ss))
            {
                cat(ss[i])
                s <- ans$Source[i]; v <- ans$Version[i]
                if (!is.na(s) | !is.na(v))
                    cat(", ", .str(s), ", ", .str(v), sep="")
                cat("\n")
            }
        } else
            cat("    <None>\n")
    }
    invisible(ans)
}

.summary_format <- function(gdsfile, check, verbose, GT=TRUE)
{
    if (verbose)
        cat("Annotation, FORMAT variable(s):\n")

    # genotype
    if (GT)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, "genotype"))
        if (is.null(a$VariableName)) a$VariableName <- "GT"
        if (is.null(a$Description)) a$Description <- NA_character_
        ans <- data.frame(ID=a$VariableName, Number=1, Type="String",
            Description=a$Description[1L], stringsAsFactors=FALSE)
    } else {
        s <- character()
        ans <- data.frame(ID=s, Number=integer(), Type=s, Description=s,
            stringsAsFactors=FALSE)
    }

    # the FORMAT field
    fmt <- index.gdsn(gdsfile, "annotation/format")
    for (nm in ls.gdsn(fmt))
    {
        a <- get.attr.gdsn(index.gdsn(fmt, nm))
        if (is.null(a$Number)) a$Number <- "."
        if (is.null(a$Type)) a$Type <- .vcf_type(index.gdsn(fmt, nm))
        if (is.null(a$Description)) a$Description <- NA_character_

        d <- data.frame(ID=nm, Number=a$Number[1L], Type=a$Type[1L],
            Description=a$Description[1L], stringsAsFactors=FALSE)
        ans <- rbind(ans, d)
    }

    if (verbose)
    {
        if (nrow(ans) > 0L)
        {
            cat(paste0("    ", ans$ID, ", ", ans$Number, ", ", ans$Type,
                ", ", .str(ans$Description), "\n"), sep="")
        } else
            cat("    <None>\n")
    }
    invisible(ans)
}

.summary_samp_annot <- function(gdsfile, check, verbose)
{
    folder <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
    if (!is.null(folder))
    {
        if (verbose)
            cat("Annotation, sample variable(s):\n")
        ID <- ls.gdsn(folder)
        Types <- character(length(ID))
        Desp <- rep(NA_character_, length(ID))
        for (i in seq_len(length(ID)))
        {
            n <- index.gdsn(folder, ID[i])
            Types[i] <- .vcf_type(n)
            at <- get.attr.gdsn(n)
            if (!is.null(at$Description))
                Desp[i] <- as.character(at$Description)[1L]
        }
        if (verbose)
        {
            if (length(ID) > 0L)
            {
                cat(paste0("    ", ID, ", ", Types, ", ", .str(Desp), "\n"),
                    sep="")
            } else
                cat("    <None>\n")
        }
        ans <- data.frame(ID=ID, Type=Types, Description=Desp,
            stringsAsFactors=FALSE)
    } else {
        s <- character()
        ans <- data.frame(ID=s, Type=s, Description=s, stringsAsFactors=FALSE)
    }

    invisible(ans)
}

.summary_reference <- function(gdsfile, check, verbose)
{
    n <- index.gdsn(gdsfile, "description/reference", silent=TRUE)
    if (is.null(n))
        ans <- character()
    else
        ans <- as.character(read.gdsn(n))
    n <- index.gdsn(gdsfile, "description/vcf.header", silent=TRUE)
    if (!is.null(n))
    {
        tmp <- read.gdsn(n)
        if (is.data.frame(tmp))
            ans <- c(ans, tmp$value[tmp$id == "reference"])
        else
            message("Invalid 'description/vcf.header'.")
    }
    ans <- unique(ans)
    if (verbose)
    {
        if (length(ans) > 0)
            cat("Reference: ", paste(ans, collapse=", "), "\n", sep="")
        else
            cat("Reference: unknown\n")
    }
    invisible(ans)
}

.summary_contig <- function(gdsfile, check, verbose)
{
    n <- index.gdsn(gdsfile, "description/vcf.contig", silent=TRUE)
    if (!is.null(n))
    {
        ans <- read.gdsn(n)
        if (check != "none" & !is.data.frame(ans))
            message("'description/vcf.contig' should be data.frame.")
        if (verbose)
        {
            cat("Contigs:\n")
            if (is.data.frame(ans))
            {
                if (nrow(ans) > 0L)
                {
                    s <- paste0("    ", .str(ans[, 1L]))
                    for (i in seq_len(ncol(ans))[-1L])
                        s <- paste0(s, ", ", .str(ans[, i]))
                    s <- paste0(s, "\n")
                    cat(s, sep="")
                } else
                    cat("    <None>\n")
            }
        }
        invisible(ans)
    } else
        NULL
}

.summary_digest <- function(gdsfile)
{
    enum <- function(node)
    {
        ans <- NULL
        dp <- objdesp.gdsn(node)
        if (dp$type %in% c("Folder", "VFolder"))
        {
            for (nm in ls.gdsn(node, include.hidden=TRUE))
            {
                ans <- rbind(ans, enum(index.gdsn(node, nm)))
            }
        } else if (dp$is.array)
        {
            v <- digest.gdsn(node, action="return")
            if (any(!is.na(v)))
            {
                at <- get.attr.gdsn(node)
                at <- format(at[names(v)[!is.na(v)]])
                err <- any(v==FALSE, na.rm=TRUE)
                d <- data.frame(varname = dp$fullname,
                    digest = toString(paste(names(at), at, sep=": ")),
                    validation = !err, stringsAsFactors=FALSE)
                ans <- rbind(ans, d)
            }
        }
        ans
    }

    enum(gdsfile$root)
}


#######################################################################
# summarize
#
seqSummary <- function(gdsfile, varname=NULL,
    check=c("default", "none", "full"), verbose=TRUE)
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
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))   
    }

    if (is.null(varname))
    {
        ## initialize
        ans <- list()
        ans$filename <- gdsfile$filename
        if (verbose)
            cat("File: ", ans$filename, "\n", sep="")

        ## version
        at <- get.attr.gdsn(gdsfile$root)
        if ("FileVersion" %in% names(at))
        {
            ans$version <- at$FileVersion
            if (verbose)
                cat("Format Version: ", ans$version, "\n", sep="")
        }

        ## reference
        ans$reference <- .summary_reference(gdsfile, check, verbose)

        ## sample.id
        .summary_sample_id(gdsfile, check, verbose)

        ## variant.id
        .summary_variant_id(gdsfile, check, verbose)

        ## genotype
        at <- .summary_geno(gdsfile, check, verbose, FALSE)
        ans$ploidy <- at$dim[1L]
        ans$num.sample <- at$dim[2L]
        ans$num.variant <- at$dim[3L]

        ## phase
        .summary_phase(gdsfile, check, FALSE)

        ## position
        .summary_position(gdsfile, check, verbose)

        ## chromosome
        .summary_chrom(gdsfile, check, verbose)

        ## contigs
        .summary_contig(gdsfile, check, verbose)

        ## allele
        ans$allele <- .summary_allele(gdsfile, check, verbose)

        ## annotation/id
        .summary_annot_id(gdsfile, check, verbose)

        ## annotation/qual
        ans$annot_qual <- .summary_annot_qual(gdsfile, check, verbose)

        ## annotation/filter
        n <- index.gdsn(gdsfile, "annotation/filter", silent=TRUE)
        if (!is.null(n))
            ans$filter <- .summary_filter(gdsfile, check, verbose)

        ## annotation/info
        n <- index.gdsn(gdsfile, "annotation/info", silent=TRUE)
        if (!is.null(n))
            ans$info <- .summary_info(gdsfile, check, verbose)

        ## annotation/format
        n <- index.gdsn(gdsfile, "annotation/format", silent=TRUE)
        if (!is.null(n))
            ans$format <- .summary_format(gdsfile, check, verbose)

        ## sample.annotation
        n <- index.gdsn(gdsfile, "sample.annotation", silent=TRUE)
        if (!is.null(n))
            ans$sample.annot <- .summary_samp_annot(gdsfile, check, verbose)

        # output
        invisible(ans)

    } else {
        switch(varname,
            `genotype` = .summary_geno(gdsfile, check, verbose),
            `phase` = .summary_phase(gdsfile, check, verbose),
            `sample.id` = .summary_sample_id(gdsfile, check, verbose),
            `variant.id` = .summary_variant_id(gdsfile, check, verbose),
            `position` = .summary_position(gdsfile, check, verbose),
            `chromosome` = .summary_chrom(gdsfile, check, verbose),
            `allele` = .summary_allele(gdsfile, check, verbose),
            `annotation/id` = .summary_annot_id(gdsfile, check, verbose),
            `annotation/qual` = .summary_annot_qual(gdsfile, check, verbose),
            `annotation/filter` = .summary_filter(gdsfile, check, verbose),
            `annotation/info` = .summary_info(gdsfile, check, verbose),
            `annotation/format` = .summary_format(gdsfile, check, verbose, FALSE),
            `sample.annotation` = .summary_samp_annot(gdsfile, check, verbose),
            `$reference` = .summary_reference(gdsfile, check, verbose),
            `$contig` = .summary_contig(gdsfile, check, verbose),
            `$filter` = .summary_filter(gdsfile, check, verbose),
            `$alt` = .summary_allele(gdsfile, "none", verbose),
            `$info` = .summary_info(gdsfile, check, verbose),
            `$format` = .summary_format(gdsfile, check, verbose),
            `$digest` = .summary_digest(gdsfile),

            .Call(SEQ_Summary, gdsfile, varname)  # others
        )
    }
}


#######################################################################
# summarize
#
seqDigest <- function(gdsfile, varname, algo=c("md5"), verbose=FALSE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(varname), length(varname)==1L)
    algo <- match.arg(algo)

    if (requireNamespace("digest", quietly=TRUE))
    {
        .cfunction("FC_DigestInit")(algo)
        seqApply(gdsfile, varname, FUN=.cfunction("FC_DigestScan"),
            margin="by.variant", as.is="none", .useraw=NA, .progress=verbose)
        .cfunction("FC_DigestDone")(algo)
    } else
        stop("The digest package is not installed.")
}


#######################################################################
# summarize
#
seqSystem <- function()
{
    rv <- .Call(SEQ_System)
    rv$options <- list(
        seqarray.parallel = seqGetParallel()
    )
    rv
}


#######################################################################
# Perform the checking for the GDS file
#

.check_var_dim <- function(f, nm, type, cnt, verbose)
{
    if (verbose)
    {
        cat("    ", nm, sep="")
        flush.console()
    }
    dp <- objdesp.gdsn(index.gdsn(f, nm))
    if (isTRUE(dp$type != type))
        s <- "Invalid data type."
    else if (dp$is.array && (length(dp$dim)!=1L || dp$dim[1L]!=cnt))
        s <- "Invalid dimension."
    else
        s <- ""
    rv <- data.frame(ok=(s==""), info=s, stringsAsFactors=FALSE)
    rownames(rv) <- nm
    rv
}

seqCheck <- function(gdsfile, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.character(gdsfile) && length(gdsfile)==1L)
    {
        if (verbose)
            cat("Open '", gdsfile, "'\n", sep="")
        gdsfile <- seqOpen(gdsfile)
        on.exit({ seqClose(gdsfile) })
    }

    # the output value
    ans <- list()

    # hash check
    if (verbose)
        cat("Hash check:\n")
    nm_lst <- ls.gdsn(gdsfile, include.hidden=TRUE, recursive=TRUE)
    s <- sapply(nm_lst, function(nm) {
            ns <- names(get.attr.gdsn(index.gdsn(gdsfile, nm)))
            any(ns %in% c("md5", "sha1", "sha256", "sha384", "sha512"))
        })
    nm_lst <- nm_lst[s]
    ans$hash <- data.frame(
        algo = sapply(nm_lst, function(nm) {
            paste(intersect(c("md5", "sha1", "sha256", "sha384", "sha512"),
                names(get.attr.gdsn(index.gdsn(gdsfile, nm)))), collapse=",")
        }),
        ok = sapply(nm_lst, function(nm)
            all(.check_digest(gdsfile, nm, verbose))),
        stringsAsFactors=FALSE)

    # dimension check
    # gds_dm[1L] -- ploidy
    # gds_dm[2L] -- # of total samples
    # gds_dm[3L] -- # of total variants
    gds_dm <- .dim(gdsfile)
    if (verbose)
        cat("Dimension check:\n")

    # add to the result
    addtab <- function(dat) {
        tab <<- rbind(tab, dat)
        if (verbose) cat(ifelse(dat$ok, "\t[OK]\n", "\t[Error]\n"))
        invisible()
    }

    tab <- NULL
    addtab(
        .check_var_dim(gdsfile, "variant.id", NA, gds_dm[3L], verbose))
    addtab(
        .check_var_dim(gdsfile, "position", "Integer", gds_dm[3L], verbose))
    addtab(
        .check_var_dim(gdsfile, "chromosome", "String", gds_dm[3L], verbose))
    addtab(
        .check_var_dim(gdsfile, "allele", "String", gds_dm[3L], verbose))
    ans$dimension <- tab

    invisible(ans)
}
