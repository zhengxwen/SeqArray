#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Big Data Management of Genome-Wide Sequence Variants
#



#######################################################################
# Convert a VCF (sequence) file to a GDS file
#######################################################################

#######################################################################
# Parse the header of a VCF file
# http://www.1000genomes.org/wiki/analysis/variant-call-format
#

seqVCF.Header <- function(vcf.fn, getnum=FALSE)
{
    # check
    stopifnot(is.character(vcf.fn))
    stopifnot(is.logical(getnum), length(getnum)==1L)

    #########################################################
    # open the vcf file

    # parse and determine how many copies of genomes: haploid, diploid or others
    geno.text <- NULL
    nSample <- 0L
    nVariant <- 0L

    n <- 0L
    for (i in seq_along(vcf.fn))
    {
        opfile <- file(vcf.fn[i], open="rt")
        on.exit(close(opfile))

        # read header
        ans <- NULL
        while (length(s <- readLines(opfile, n=1L)) > 0L)
        {
            n <- n + 1L
            if (substr(s, 1L, 6L) != "#CHROM")
            {
                s <- substring(s, 3L)
                ans <- c(ans, s)
            } else {
                if (getnum)
                {
                    s <- scan(text=s, what=character(0), sep="\t",
                        quiet=TRUE)[-seq_len(9)]
                    if (length(s) > nSample) nSample <- length(s)
                }
                s <- readLines(opfile, n=1L)
                if (length(s) > 0L)
                {
                    ss <- scan(text=s, what=character(), sep="\t", quiet=TRUE)
                    geno.text <- c(geno.text, ss[-seq_len(9L)])
                }
                if (getnum)
                {
                    nVariant <- nVariant + length(s) +
                        length(count.fields(opfile))
                }
                break
            }
            if ((n %% 10000L) == 0L)
            {
                warning(sprintf(
                "There are too many lines in the header (>= %d). ", n),
                "In order not to slow down the conversion, please consider ",
                "deleting unnecessary annotations (like contig).",
                immediate.=TRUE)
            }
        }

        close(opfile)
        on.exit()
    }

    ans <- unique(ans)


    #########################################################
    # parse texts

    #########################################################
    # get names and values

    ValString <- function(txt)
    {
        # split by "="
        x <- as.integer(regexpr("=", txt, fixed=TRUE))
        rv <- matrix("", nrow=length(txt), ncol=2L)
        for (i in seq_along(txt))
        {
            if (x[i] > 0L)
            {
                rv[i,1L] <- substring(txt[i], 1L, x[i]-1L)
                rv[i,2L] <- substring(txt[i], x[i]+1L)
            }
        }
        rv
    }


    #########################################################
    # get names and values, length(txt) == 1L

    AsDataFrame <- function(txt, check.name)
    {
        if (substr(txt, 1L, 1L) == "<")
            txt <- substring(txt, 2L)
        if (substr(txt, nchar(txt), nchar(txt)) == ">")
            txt <- substr(txt, 1L, nchar(txt)-1L)

        v <- ValString(scan(text=txt, what=character(),
            sep=",", quote="\"", quiet=TRUE))

        # check
        ans <- list()
        if (!is.null(check.name))
        {
            if (is.null(names(check.name)))
                flag <- rep(TRUE, length(check.name))
            else
                flag <- (names(check.name) == "T")

            for (i in seq_along(check.name))
            {
                j <- match(check.name[i], v[,1L])
                if (is.na(j))
                {
                    if (flag[i])
                        stop("No '", check.name[i], "'.")
                    a <- NA_character_
                } else
                    a <- v[i,2L]

                ans[[ check.name[i] ]] <- a
            }
        } else {
            for (i in seq_len(nrow(v)))
                ans[[ v[i,1L] ]] <- v[i,2L]
        }
        as.data.frame(ans, stringsAsFactors=FALSE)
    }

    #########################################################
    # check number

    CheckNum <- function(number)
    {
        if (!(number %in% c("A", "G", ".")))
        {
            N <- suppressWarnings(as.integer(number))
            if (is.finite(N))
                N >= 0L
            else
                FALSE
        } else
            TRUE
    }


    #########################################################
    # start

    #########################################################
    # ploidy

    if (!is.null(geno.text))
    {
        txt <- unlist(sapply(geno.text, function(s) {
            scan(text=s, what=character(), sep=":", quiet=TRUE, nmax=1) },
            simplify=TRUE, USE.NAMES=FALSE))
        num <- sapply(strsplit(txt, "[|/]"), function(x) length(x) )
        tab <- table(num)
        ploidy <- as.integer(names(which.max(tab)))
    } else
        ploidy <- as.integer(NA)

    if (is.null(ans))
    {
        rv <- list(fileformat="unknown", info=NULL, filter=NULL, format=NULL,
            alt=NULL, contig=NULL, assembly=NULL, header=NULL,
            ploidy = ploidy)
        class(rv) <- "SeqVCFHeaderClass"
        return(rv)
    }

    ans <- ValString(ans)
    ans <- data.frame(id=ans[,1L], value=ans[,2L], stringsAsFactors=FALSE)


    #########################################################
    # fileformat=VCFv4.*

    s <- ans$value[ans$id == "fileformat"]
    if (length(s) < 1L)
        stop("no defined fileformat!")
    else if (length(s) > 1L)
        stop("Multiple fileformat!")
    else {
        if (substr(s, 1L, 4L) == "VCFv")
        {
            v <- suppressWarnings(as.double(substring(s, 5)))
            if (!is.finite(v)) v <- 0
            if (v < 4.0)
                stop("'fileformat' should be >= v4.0.")
        } else
            stop("Invalid 'fileformat'.")
    }
    fileformat <- s[1L]
    ans <- ans[ans$id != "fileformat", ]

    #########################################################
    # assembly=url

    s <- ans$value[ans$id == "assembly"]
    assembly <- if (length(s) > 0L) s else NULL
    ans <- ans[ans$id != "assembly", ]


    #########################################################
    # INFO=<ID=ID,Number=number,Type=type,Description="description">
    # INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">

    INFO <- NULL
    s <- ans$value[ans$id == "INFO"]
    for (i in seq_along(s))
    {
        v <- AsDataFrame(s[i], c(T="ID", T="Number", T="Type", T="Description",
            F="Source", F="Version"))
        if (!is.null(v))
        {
            if (!is.element(tolower(v$Type),
                    c("integer", "float", "flag", "character", "string")))
                stop("INFO=", s[i])
            if (!CheckNum(v$Number))
                stop("INFO=", s[i])
            INFO <- rbind(INFO, v)
        } else
            stop("INFO=", s[i])
    }
    ans <- ans[ans$id != "INFO", ]


    #########################################################
    # FILTER=<ID=ID,Description="description">

    FILTER <- NULL
    s <- ans$value[ans$id == "FILTER"]
    for (i in seq_along(s))
    {
        v <- AsDataFrame(s[i], c("ID", "Description"))
        if (!is.null(v))
            FILTER <- rbind(FILTER, v)
        else
            stop("FILTER=", s[i])
    }
    ans <- ans[ans$id != "FILTER", ]


    #########################################################
    # FORMAT=<ID=ID,Number=number,Type=type,Description="description">

    FORMAT <- NULL
    s <- ans$value[ans$id == "FORMAT"]
    for (i in seq_along(s))
    {
        v <- AsDataFrame(s[i], c("ID", "Number", "Type", "Description"))
        if (!is.null(v))
        {
            if (!is.element(tolower(v$Type),
                    c("integer", "float", "character", "string")))
                stop("FORMAT=", s[i])
            if (!CheckNum(v$Number))
                stop("FORMAT=", s[i])
            FORMAT <- rbind(FORMAT, v)
        } else
            stop("FORMAT=", s[i])
    }
    ans <- ans[ans$id != "FORMAT", ]


    #########################################################
    # ALT=<ID=type,Description=description>

    ALT <- NULL
    s <- ans$value[ans$id == "ALT"]
    for (i in seq_along(s))
    {
        v <- AsDataFrame(s[i], c("ID", "Description"))
        if (!is.null(v))
        {
            ALT <- rbind(ALT, v)
        } else
            stop("ALT=", s[i])
    }
    ans <- ans[ans$id != "ALT", ]


    #########################################################
    # contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

    contig <- NULL
    s <- ans$value[ans$id == "contig"]
    for (i in seq_along(s))
    {
        v <- AsDataFrame(s[i], NULL)
        if (!is.null(v))
        {
            contig <- rbind(contig, v)
        } else
            stop("contig=", s[i])
    }
    ans <- ans[ans$id != "contig", ]


    #########################################################
    # reference=""
    reference <- NULL
    s <- ans$value[ans$id == "reference"]
    if (length(s) > 0L) reference <- s
    ans <- ans[ans$id != "reference", ]


    #########################################################
    # output

    rv <- list(fileformat=fileformat, info=INFO, filter=FILTER, format=FORMAT,
        alt=ALT, contig=contig, assembly=assembly, reference=reference,
        header=ans, ploidy=ploidy)
    if (getnum)
    {
        rv$num.sample <- nSample
        rv$num.variant <- nVariant
    }
    class(rv) <- "SeqVCFHeaderClass"
    rv
}



#######################################################################
# get sample id from a VCF file
#

seqVCF.SampID <- function(vcf.fn)
{
    # check
    stopifnot(is.character(vcf.fn))
    stopifnot(length(vcf.fn) == 1L)

    # open the vcf file
    opfile <- file(vcf.fn[1L], open="rt")
    on.exit(close(opfile))

    # read header
    samp.id <- NULL
    while (length(s <- readLines(opfile, n=1L)) > 0L)
    {
        if (substr(s, 1L, 6L) == "#CHROM")
        {
            samp.id <- scan(text=s, what=character(0), sep="\t",
                quiet=TRUE)[-seq_len(9)]
            break
        }
    }
    if (is.null(samp.id))
        stop("Error VCF format: invalid sample id!")

    return(samp.id)
}



#######################################################################
# Convert a VCF file to a GDS file
#

seqVCF2GDS <- function(vcf.fn, out.fn, header=NULL,
    genotype.var.name="GT", genotype.storage=c("bit2", "bit4", "bit8"),
    storage.option=seqStorage.Option("ZIP_RA.default"),
    info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr",
    reference=NULL, start=1L, count=-1L, optimize=TRUE, raise.error=TRUE,
    digest=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf.fn), length(vcf.fn)>0L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)

    stopifnot(is.null(header) | inherits(header, "SeqVCFHeaderClass"))
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))

    stopifnot(is.character(genotype.var.name), length(genotype.var.name)==1L)
    stopifnot(!is.na(genotype.var.name))
    genotype.storage <- match.arg(genotype.storage)

    stopifnot(is.null(info.import) | is.character(info.import))
    stopifnot(is.null(fmt.import) | is.character(fmt.import))
    stopifnot(is.character(ignore.chr.prefix), length(ignore.chr.prefix)>0L)
    stopifnot(is.null(reference) | is.character(reference))

    stopifnot(is.numeric(start), length(start)==1L)
    stopifnot(is.numeric(count), length(count)==1L)

    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(raise.error), length(raise.error)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose) cat(date(), "\n", sep="")

    # check sample id
    samp.id <- NULL
    for (i in seq_along(vcf.fn))
    {
        if (is.null(samp.id))
        {
            samp.id <- seqVCF.SampID(vcf.fn[i])
            if (length(samp.id) <= 0L)
                stop("There is no sample in the VCF file.")
        } else {
            tmp <- seqVCF.SampID(vcf.fn[i])
            if (length(samp.id) != length(tmp))
            {
                stop(sprintf("The file '%s' has different sample id.",
                    vcf.fn[i]))
            }
            if (!all(samp.id == tmp))
            {
                stop(sprintf("The file '%s' has different sample id.",
                    vcf.fn[i]))
            }
        }
    }


    ##################################################
    # parse the header of VCF

    if (is.null(header))
        header <- seqVCF.Header(vcf.fn)
    if (!inherits(header, "SeqVCFHeaderClass"))
        stop("'header' should be NULL or returned from 'seqVCF.Header()'.")
    # 'seqVCF.Header'
    # returns list(fileformat=fileformat, info=INFO, filter=FILTER,
    #              format=FORMAT, alt=ALT, contig=contig, assembly=assembly,
    #              reference, header, ploidy)

    if (verbose)
    {
        cat("The Variant Call Format (VCF) header:\n")
        cat("    file format: ", header$fileformat, "\n", sep="")
        cat("    the number of sets of chromosomes (ploidy): ",
            header$ploidy, "\n", sep="")
        cat("    the number of samples: ", length(samp.id), "\n", sep="")
        cat("    GDS genotype storage: ", genotype.storage, "\n", sep="")
    }

    # check header
    tmp <- FALSE
    if (!is.null(header$format))
    {
        flag <- duplicated(header$format$ID)
        if (any(flag))
        {
            if (verbose)
            {
                cat(sprintf("Duplicated Format ID (%s) are removed.\n",
                    paste(header$format$ID[flag], collapse=",")))
            }
            header$format <- header$format[!flag, ]
        }
        if (nrow(header$format) > 0L)
            tmp <- TRUE
    }
    if (tmp)
    {
        geno_idx <- which(header$format$ID == genotype.var.name)
        if (length(geno_idx) <= 0L)
        {
            stop(sprintf(
            "There is no variable in the FORMAT field according to '%s'.",
                genotype.var.name))
        } else if (length(geno_idx) > 1L)
        {
            stop(sprintf(
            "There are more than one variable in the FORMAT field according to '%s'.",
                genotype.var.name))
        } else {
            if (tolower(header$format$Type[geno_idx]) != "string")
            {
                stop(sprintf(
                "'%s' in the FORMAT field should be string-type according to genotypes.",
                    genotype.var.name))
            }
            if (header$format$Number[geno_idx] != "1")
            {
                stop(sprintf(
                "'%s' in the FORMAT field should be one-length string-type.",
                    genotype.var.name))
            }
        }
        geno_format <- header$format[geno_idx, ]
        header$format <- header$format[-geno_idx, ]
    } else {
        message("\t",
            "variable id in the FORMAT field should be defined ahead, ",
            "and the undefined id is/are ignored during the conversion.")
        geno_format <- list(Description="Genotype")
    }



    #######################################################################
    # create a GDS file

    gfile <- createfn.gds(out.fn)
    on.exit(closefn.gds(gfile))

    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(gfile, "description")
    put.attr.gdsn(n, "vcf.fileformat", header$fileformat)
    if (!is.null(header$assembly))
        put.attr.gdsn(n, "vcf.assembly", header$assembly)
    reference <- as.character(unique(c(reference, header$reference)))
    .AddVar(storage.option, n, "reference", reference, closezip=TRUE,
        visible=FALSE)
    if (!is.null(header$alt))
    {
        if (nrow(header$alt) > 0L)
        {
            .AddVar(storage.option, n, "vcf.alt", header$alt, closezip=TRUE,
                visible=FALSE)
        }
    }
    if (!is.null(header$contig))
    {
        if (nrow(header$contig) > 0L)
        {
            .AddVar(storage.option, n, "vcf.contig", header$contig,
                closezip=TRUE, visible=FALSE)
        }
    }
    if (!is.null(header$header))
    {
        if (nrow(header$header) > 0L)
        {
            .AddVar(storage.option, n, "vcf.header", header$header,
                closezip=TRUE, visible=FALSE)
        }
    }


    ##################################################
    # add sample.id

    nSamp <- length(samp.id)
    .AddVar(storage.option, gfile, "sample.id", samp.id, closezip=TRUE)


    ##################################################
    # add basic site information

    # add variant.id
    .AddVar(storage.option, gfile, "variant.id", storage="int32")

    # add position
    # TODO: need to check whether position can be stored in 'int32'
    .AddVar(storage.option, gfile, "position", storage="int32")

    # add chromosome
    .AddVar(storage.option, gfile, "chromosome", storage="string")

    # add allele
    .AddVar(storage.option, gfile, "allele", storage="string")

    # add a folder for genotypes
    varGeno <- addfolder.gdsn(gfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", genotype.var.name[1L])
    put.attr.gdsn(varGeno, "Description", geno_format$Description[1L])

    # add data to the folder of genotype
    if (header$ploidy > 1L)
    {
        geno.node <- .AddVar(storage.option, varGeno, "data",
            storage=genotype.storage, valdim=c(header$ploidy, nSamp, 0L))
    } else {
        geno.node <- .AddVar(storage.option, varGeno, "data",
            storage=genotype.storage, valdim=c(1L, nSamp, 0L))
    }
    node <- .AddVar(storage.option, varGeno, "@data", storage="uint8",
        visible=FALSE)

    node <- .AddVar(storage.option, varGeno, "extra.index", storage="int32",
        valdim=c(3L,0L))
    put.attr.gdsn(node, "R.colnames",
        c("sample.index", "variant.index", "length"))
    .AddVar(storage.option, varGeno, "extra", storage="int16")


    # add phase folder
    varPhase <- addfolder.gdsn(gfile, "phase")
    if (header$ploidy > 1L)
    {
        # add data
        if (header$ploidy > 2L)
        {
            .AddVar(storage.option, varPhase, "data", storage="bit1",
                valdim=c(header$ploidy-1L, nSamp, 0L))
        } else {
            .AddVar(storage.option, varPhase, "data", storage="bit1",
                valdim=c(nSamp, 0L))
        }

        node <- .AddVar(storage.option, varPhase, "extra.index",
            storage="int32", valdim=c(3L,0L))
        put.attr.gdsn(node, "R.colnames",
            c("sample.index", "variant.index", "length"))
        .AddVar(storage.option, varPhase, "extra", storage="bit1")
    }


    ##################################################

    # add annotation folder
    varAnnot <- addfolder.gdsn(gfile, "annotation")

    # add id
    .AddVar(storage.option, varAnnot, "id", storage="string")
    # add qual
    .AddVar(storage.option, varAnnot, "qual", storage="float")
    # add filter
    varFilter <- .AddVar(storage.option, varAnnot, "filter", storage="int32")


    ##################################################
    # VCF INFO

    varInfo <- addfolder.gdsn(varAnnot, "info")

    if (!is.null(header$info))
    {
        flag <- duplicated(header$info$ID)
        if (any(flag))
        {
            if (verbose)
            {
                cat(sprintf("Duplicated Info ID (%s) are removed.\n",
                    paste(header$info$ID[flag], collapse=",")))
            }
            header$info <- header$info[!flag, ]
        }

        if (nrow(header$info) > 0L)
        {
            int_type <- integer(nrow(header$info))
            int_num  <- suppressWarnings(as.integer(header$info$Number))
            if (is.null(info.import))
                import.flag <- rep(TRUE, length(int_num))
            else
                import.flag <- header$info$ID %in% info.import

            for (i in 1:nrow(header$info))
            {
                # INFO Type
                switch(tolower(header$info$Type[i]),
                    integer = { mode <- "int32"; int_type[i] <- 1L },
                    float = { mode <- "float"; int_type[i] <- 2L },
                    flag = { mode <- "bit1"; int_type[i] <- 3L },
                    character = { mode <- "string"; int_type[i] <- 4L },
                    string = { mode <- "string"; int_type[i] <- 4L },
                    stop(sprintf("Unknown INFO Type: %s", header$info$Type[i]))
                )

                # INFO Number
                s <- header$info$Number[i]
                if (grepl("^[[:digit:]]", s))
                {
                    initdim <- as.integer(s)
                    if (mode == "bit1")
                    {
                        if (initdim != 0L)
                        {
                            print(header$info[i, ])
                            stop("The length of 'Flag' type should be ZERO!")
                        }
                    } else {
                        if (initdim <= 0L)
                        {
                            print(header$info[i, ])
                            stop("The length should be >0.")
                        } else if (initdim > 1L)
                            initdim <- c(initdim, 0L)
                        else
                            initdim <- 0L
                    }
                } else {
                    initdim <- 0L
                    if (s == ".")
                        int_num[i] <- -1L
                    else if (s == "A")
                        int_num[i] <- -2L
                    else if (s == "G")
                        int_num[i] <- -3L
                    else if (s == "R")
                        int_num[i] <- -4L
                    else {
                        stop(sprintf("Unknown INFO (%s) Number: %s",
                            header$info$ID[i], s))
                    }
                }

                # add
                if (import.flag[i])
                {
                    node <- .AddVar(storage.option, varInfo, header$info$ID[i],
                        storage=mode, valdim=initdim)
                    put.attr.gdsn(node, "Number", header$info$Number[i])
                    put.attr.gdsn(node, "Type", header$info$Type[i])
                    put.attr.gdsn(node, "Description", header$info$Description[i])
                    if (!is.na(header$info$Source[i]))
                        put.attr.gdsn(node, "Source", header$info$Source[i])
                    if (!is.na(header$info$Version[i]))
                        put.attr.gdsn(node, "Version", header$info$Version[i])

                    if (s %in% c(".", "A", "G"))
                    {
                        node <- .AddVar(storage.option, varInfo,
                            paste("@", header$info$ID[i], sep=""),
                            storage="int32", visible=FALSE)
                    }
                }
            }

            header$info$int_type <- as.integer(int_type)
            header$info$int_num  <- as.integer(int_num)
            header$info$import.flag <- import.flag
        }
    }


    ##################################################
    # VCF Format

    # add the FORMAT field
    varFormat <- addfolder.gdsn(varAnnot, "format")

    if (!is.null(header$format))
    {
        int_type <- integer(nrow(header$format))
        int_num  <- suppressWarnings(as.integer(header$format$Number))
        if (is.null(fmt.import))
            import.flag <- rep(TRUE, length(int_num))
        else
            import.flag <- header$format$ID %in% fmt.import
    } else {
        int_type <- integer()
        int_num <- integer()
        import.flag <- logical()
    }

    for (i in seq_along(int_type))
    {
        # FORMAT Type
        switch(tolower(header$format$Type[i]),
            integer = { mode <- "int32"; int_type[i] <- 1L },
            float = { mode <- "float"; int_type[i] <- 2L },
            character = { mode <- "string"; int_type[i] <- 4L },
            string = { mode <- "string"; int_type[i] <- 4L },
            stop(sprintf("Unknown FORMAT Type: %s", header$format$Type[i]))
        )

        # FORMAT Number
        s <- header$format$Number[i]
        if (grepl("^[[:digit:]]", s))
        {
            initdim <- as.integer(s)
            if (initdim <= 0)
            {
                print(header[, i])
                stop("The length should be >0.")
            } else if (initdim > 1)
                initdim <- c(initdim, nSamp, 0L)
            else
                initdim <- c(nSamp, 0L)
        } else {
            initdim <- c(nSamp, 0L)
            if (s == ".")
                int_num[i] <- -1L
            else if (s == "A")
                int_num[i] <- -2L
            else if (s == "G")
                int_num[i] <- -3L
            else if (s == "R")
                int_num[i] <- -4L
            else {
                stop(sprintf("Unknown FORMAT (%s) Number: %s",
                    header$format$ID[i], s))
            }
        }

        # add
        if (import.flag[i])
        {
            node <- addfolder.gdsn(varFormat, header$format$ID[i])
            put.attr.gdsn(node, "Number", header$format$Number[i])
            put.attr.gdsn(node, "Type", header$format$Type[i])
            put.attr.gdsn(node, "Description", header$format$Description[i])

            .AddVar(storage.option, node, "data", storage=mode, valdim=initdim)
            .AddVar(storage.option, node, "@data", storage="int32", visible=FALSE)
        }
    }

    if (!is.null(header$format))
    {
        header$format$int_type <- as.integer(int_type)
        header$format$int_num  <- as.integer(int_num)
        header$format$import.flag <- import.flag
    }


    ##################################################
    # add annotation folder

    addfolder.gdsn(gfile, "sample.annotation")


    ##################################################
    # sync file
    sync.gds(gfile)


    ##################################################
    # for-loop each file

    filterlevels <- header$filter$ID

    for (i in seq_along(vcf.fn))
    {
        opfile <- file(vcf.fn[i], open="rt")
        on.exit({ closefn.gds(gfile); close(opfile) })

        if (verbose)
        {
            cat(sprintf("Parsing '%s' (%s bytes)\n", vcf.fn[i],
                file.size(vcf.fn[i])))
            flush.console()
        }

        # call C function
        v <- .Call(SEQ_Parse_VCF4, vcf.fn[i], header, gfile$root,
            list(sample.num = as.integer(length(samp.id)),
                genotype.var.name = genotype.var.name,
                raise.error = raise.error, filter.levels = filterlevels,
                start = start, count = count,
                verbose = verbose),
            readLines, opfile, 512L,  # readLines(opfile, 512L)
            ignore.chr.prefix, new.env())

        filterlevels <- unique(c(filterlevels, v))
        if (verbose)
            print(geno.node)

        on.exit({ closefn.gds(gfile) })
        close(opfile)
    }

    if (length(filterlevels) > 0L)
    {
        put.attr.gdsn(varFilter, "R.class", "factor")
        put.attr.gdsn(varFilter, "R.levels", filterlevels)

        if (!is.null(header$filter))
        {
            if (nrow(header$filter) > 0L)
            {
                dp <- header$filter$Description[match(filterlevels,
                    header$filter$ID)]
                dp[is.na(dp)] <- ""
            }
        } else
            dp <- rep("", length(filterlevels))

        put.attr.gdsn(varFilter, "Description", dp)
    }

    .DigestFile(gfile, digest, verbose)

    on.exit()
    closefn.gds(gfile)


    ##################################################
    # optimize access efficiency

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.fn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a SeqArray GDS file to a VCF file
#

seqGDS2VCF <- function(gdsfile, vcf.fn, info.var=NULL, fmt.var=NULL,
    verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)
    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))

    # get a summary
    z <- seqSummary(gdsfile, check="none", verbose=FALSE)

    # the INFO field
    if (!is.null(info.var))
    {
        s <- z$info$ID
        if (is.null(s)) s <- character()
        if (length(setdiff(info.var, s)) > 0L)
            stop(paste("Not exist:", paste(setdiff(info.var, s), collapse=",")))
        if (length(info.var) > 0L)
            z$info <- z$info[match(info.var, z$info$ID), ]
        else
            z$info <- list()
    }

    # the FORMAT field
    if (!is.null(fmt.var))
    {
        s <- z$format$ID[-1L]
        if (is.null(s)) s <- character()
        if (length(setdiff(fmt.var, s)) > 0L)
            stop(paste("Not exist:", paste(setdiff(fmt.var, s), collapse=",")))
        if (length(fmt.var) > 0L)
            z$format <- z$format[match(fmt.var, z$format$ID), ]
        else
            z$format <- list()
    } else {
        z$format <- z$format[-1L, ]
    }


    ## double quote text if needed
    dq <- function(s, text=FALSE)
    {
        .Call(SEQ_Quote, s, text)
    }


    ######################################################
    # create an output text file

    vcf.fn <- vcf.fn[1L]
    ext <- substring(vcf.fn, nchar(vcf.fn)-2L)
    if (ext == ".gz")
    {
        ofile <- gzfile(vcf.fn, "wt")
    } else if (ext == ".bz")
    {
        ofile <- bzfile(vcf.fn, "wt")
    } else if (ext == ".xz")
    {
        ofile <- xzfile(vcf.fn, "wt")
    } else {
        ofile <- file(vcf.fn, open="wt")
    }
    op <- options("useFancyQuotes")
    options(useFancyQuotes = FALSE)
    on.exit({ options(op); close(ofile) })

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("Output: ", vcf.fn, "\n", sep="")
        cat("The INFO field: ", paste(z$info$ID, collapse=", "),
            "\n", sep="")
        cat("The FORMAT field: ", paste(z$format$ID, collapse=", "),
            "\n", sep="")
    }

    ######################################################
    # write the header
    a <- get.attr.gdsn(index.gdsn(gdsfile, "description"))

    # fileformat
    if (is.null(a$vcf.fileformat))
        a$vcf.fileformat <- "VCFv4.2"
    cat("##fileformat=", a$vcf.fileformat, "\n", sep="", file=ofile)

    # fileDate
    cat("##fileDate=", format(Sys.time(), "%Y%m%d"), "\n", sep="", file=ofile)

    # program, source
    aa <- get.attr.gdsn(gdsfile$root)
    if (is.null(aa$FileVersion))
        aa$FileVersion <- "v1.0"
    cat("##source=SeqArray_Format_", aa$FileVersion, "\n", sep="", file=ofile)

    # reference
    if (length(z$reference) > 0L)
        cat("##reference=", z$reference[1L], "\n", sep="", file=ofile)

    # assembly
    if (!is.null(a$vcf.assembly))
        cat("##assembly=", dq(a$vcf.assembly), "\n", sep="", file=ofile)

    # ALT=<ID=type,Description=description>
    for (i in seq_len(nrow(z$allele)))
    {
        s <- sprintf("##ALT=<ID=%s,Description=%s>",
            as.character(z$allele[i,1L]), dq(z$allele[i,2L]))
        writeLines(s, con=ofile)
    }

    # contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
    n <- index.gdsn(gdsfile, "description/vcf.contig", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        nm <- names(dat)
        for (i in seq_len(nrow(dat)))
        {
            s <- NULL
            for (j in seq_len(ncol(dat)))
                s[j] <- paste(nm[j], "=", dq(dat[i,j]), sep="")
            s <- paste(s, collapse=",")
            cat("##contig=<", s, ">\n", sep="", file=ofile)
        }
    }

    # the INFO field
    for (nm in z$info$ID)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile,
            paste("annotation/info/", nm, sep="")))
        s <- sprintf("ID=%s,Number=%s,Type=%s,Description=%s",
            nm, dq(a$Number), dq(a$Type), dq(a$Description, TRUE))
        if (!is.null(a$Source))
            s <- paste(s, ",Source=", dq(a$Source, TRUE), sep="")
        if (!is.null(a$Version))
            s <- paste(s, ",Version=", dq(a$Version, TRUE), sep="")
        cat("##INFO=<", s, ">\n", file=ofile, sep="")
    }

    # the FILTER field
    n <- index.gdsn(gdsfile, "annotation/filter", silent=TRUE)
    if (!is.null(n))
    {
        at <- get.attr.gdsn(n)
        id <- at$R.levels; dp <- at$Description
        for (i in seq_along(id))
        {
            cat(sprintf("##FILTER=<ID=%s,Description=%s>\n",
                dq(id[i]), dq(dp[i], TRUE)), file=ofile)
        }
    }

    # the FORMAT field
    a <- get.attr.gdsn(index.gdsn(gdsfile, "genotype"))
    cat(sprintf("##FORMAT=<ID=%s,Number=1,Type=String,Description=%s>\n",
        a$VariableName, dq(a$Description, TRUE)), file=ofile)
    for (nm in z$format$ID)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile,
            paste("annotation/format/", nm, sep="")))
        cat(sprintf("##FORMAT=<ID=%s,Number=%s,Type=%s,Description=%s>\n",
            nm, dq(a$Number), dq(a$Type), dq(a$Description, TRUE)), file=ofile)
    }

    # others ...
    n <- index.gdsn(gdsfile, "description/vcf.header", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        for (i in seq_len(nrow(dat)))
        {
            s <- dat[i,1L]
            if (!(s %in% c("fileDate", "source")))
                cat("##", s, "=", dq(dat[i,2L]), "\n", sep="", file=ofile)
        }
    }


    ######################################################
    # write the header -- samples

    cat(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
        seqGetData(gdsfile, "sample.id")), sep="\t", file=ofile)
    cat("\n", file=ofile)


    ######################################################
    # write the contents

    # the INFO field
    if (!is.null(z$info$ID))
        nm.info <- paste("annotation/info/", z$info$ID, sep="")
    else
        nm.info <- c()
    # the FORMAT field
    if (!is.null(z$format$ID))
        nm.format <- paste("annotation/format/", z$format$ID, sep="")
    else
        nm.format <- c()

    # initialize the variable length of INFO
    len.info <- NULL
    for (n in nm.info)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, n))
        a$Number <- if (is.null(a$Number)) "." else a$Number[1L]
        len.info <- c(len.info, a$Number)
    }
    len.info <- suppressWarnings(as.integer(len.info))

    # initialize the variable length of FORMAT
    len.fmt <- NULL
    for (n in nm.format)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, n))
        a$Number <- if (is.null(a$Number)) "." else a$Number[1L]
        len.fmt <- c(len.fmt, a$Number)
    }
    len.fmt <- suppressWarnings(as.integer(len.fmt))

    # call C function
    .Call(SEQ_InitOutVCF4, len.info, len.fmt)


    # variable names
    nm <- c("chromosome", "position", "annotation/id", "allele",
        "annotation/qual", "annotation/filter", "genotype", "phase")
    if (length(nm.info) > 0L) nm <- c(nm, nm.info)
    if (length(nm.format) > 0L) nm <- c(nm, nm.format)

    s <- c("chr", "pos", "id", "allele", "qual", "filter", "geno", "phase")
    # the INFO field
    if (length(nm.info) > 0L)
        s <- c(s, paste("info.", z$info$ID, sep=""))
    # the FORMAT field
    if (length(nm.format) > 0L)
        s <- c(s, paste("fmt.", z$format$ID, sep=""))
    names(nm) <- s

    # output lines variant by variant
    seqApply(gdsfile, nm, margin="by.variant", as.is="none",
        FUN = function(x)
        {
            cat(.Call(SEQ_OutVCF4, x), file=ofile)
        })

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(vcf.fn))
}



#######################################################################
# Convert a SeqArray GDS file to a SNP GDS file
#

seqGDS2SNP <- function(gdsfile, out.gdsfn,
    compress.geno="ZIP_RA", compress.annotation="ZIP_RA",
    optimize=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("SeqArray GDS to SNP GDS Format:\n")
    }

    # if it is a file name
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile)
        on.exit({ seqClose(gdsfile) })
    }

    # create GDS file
    gfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({ closefn.gds(gfile) }, add=TRUE)

    # add a flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")

    # add "sample.id"
    sampid <- seqGetData(gdsfile, "sample.id")
    add.gdsn(gfile, "sample.id", sampid,
        compress=compress.annotation, closezip=TRUE)

    # add "snp.id"
    add.gdsn(gfile, "snp.id", seqGetData(gdsfile, "variant.id"),
        compress=compress.annotation, closezip=TRUE)

    # add "snp.rs.id"
    if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
    {
        add.gdsn(gfile, "snp.rs.id", seqGetData(gdsfile, "annotation/id"),
            compress=compress.annotation, closezip=TRUE)
    }

    # add "snp.position"
    add.gdsn(gfile, "snp.position", seqGetData(gdsfile, "position"),
        compress=compress.annotation, closezip=TRUE)

    # add "snp.chromosome"
    add.gdsn(gfile, "snp.chromosome", seqGetData(gdsfile, "chromosome"),
        compress=compress.annotation, closezip=TRUE)

    # add "snp.allele"
    add.gdsn(gfile, "snp.allele",
        .cfunction("FC_AlleleStr")(seqGetData(gdsfile, "allele")),
        compress=compress.annotation, closezip=TRUE)

    # add "genotype"
    gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
        valdim=c(length(sampid), 0L), compress=compress.geno)
    put.attr.gdsn(gGeno, "sample.order")

    seqApply(gdsfile, "genotype", margin="by.variant", as.is="none",
        FUN = function(x) {
            g <- colSums(x==0L)
            g[is.na(g)] <- 3L
            append.gdsn(gGeno, g)
    })
    readmode.gdsn(gGeno)

    on.exit()
    seqClose(gdsfile)
    closefn.gds(gfile)

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.gdsfn))
}



#######################################################################
# Convert a SNP GDS file to a SeqArray GDS file
#

seqSNP2GDS <- function(gds.fn, out.gdsfn, compress.geno="ZIP_RA",
    compress.annotation="ZIP_RA", optimize=TRUE, digest=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gds.fn), length(gds.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("SNP GDS to SeqArray GDS Format:\n")
    }

    # open the existing SNP GDS
    srcfile <- openfn.gds(gds.fn)
    on.exit({ closefn.gds(srcfile) })

    nSamp <- prod(objdesp.gdsn(index.gdsn(srcfile, "sample.id"))$dim)
    nSNP <- prod(objdesp.gdsn(index.gdsn(srcfile, "snp.id"))$dim)

    n <- index.gdsn(srcfile, "genotype")
    dm <- objdesp.gdsn(n)$dim
    if (length(dm) != 2)
        stop("'genotype' of SNP GDS should be a matrix.")
    snpfirstdim <- TRUE
    rd <- names(get.attr.gdsn(n))
    if ("snp.order" %in% rd) snpfirstdim <- TRUE
    if ("sample.order" %in% rd) snpfirstdim <- FALSE
    if (snpfirstdim)
    {
        if ((dm[1]!=nSNP) || (dm[2]!=nSamp))
            stop("Invalid dimension of 'genotype'.")
    } else {
        if ((dm[1]!=nSamp) || (dm[2]!=nSNP))
            stop("Invalid dimension of 'genotype'.")
    }


    # create GDS file
    dstfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({ closefn.gds(dstfile) }, add=TRUE)

    put.attr.gdsn(dstfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(dstfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(dstfile, "description")
    put.attr.gdsn(n, "source.format", "SNPRelate GDS Format")

    # add sample.id
    if (verbose) cat("    sample.id")
    copyto.gdsn(dstfile, index.gdsn(srcfile, "sample.id"))
    n <- index.gdsn(dstfile, "sample.id")
    compression.gdsn(n, compress.annotation)
    .DigestCode(n, digest, verbose)

    # add variant.id
    if (verbose) cat("    variant.id")
    copyto.gdsn(dstfile, index.gdsn(srcfile, "snp.id"), "variant.id")
    n <- index.gdsn(dstfile, "variant.id")
    compression.gdsn(n, compress.annotation)
    .DigestCode(n, digest, verbose)

    # add position
    if (verbose) cat("    position")
    copyto.gdsn(dstfile, index.gdsn(srcfile, "snp.position"), "position")
    n <- index.gdsn(dstfile, "position")
    compression.gdsn(n, compress.annotation)
    .DigestCode(n, digest, verbose)

    # add chromosome
    if (verbose) cat("    chromosome")
    n <- add.gdsn(dstfile, "chromosome", storage="string",
        compress=compress.annotation)
    append.gdsn(n, index.gdsn(srcfile, "snp.chromosome"))
    readmode.gdsn(n)
    .DigestCode(n, digest, verbose)

    # add allele
    if (verbose) cat("    allele")
    n <- add.gdsn(dstfile, "allele", val = .cfunction("FC_AlleleStr2")(
        read.gdsn(index.gdsn(srcfile, "snp.allele"))),
        compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add a folder for genotypes
    if (verbose) cat("    genotype")
    n <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(n, "VariableName", "GT")
    put.attr.gdsn(n, "Description", "Genotype")


    # add genotypes to genotype/data
    n1 <- add.gdsn(n, "data", storage="bit2", valdim=c(2L, nSamp, 0L),
        compress=compress.geno)
    apply.gdsn(index.gdsn(srcfile, "genotype"), ifelse(snpfirstdim, 1L, 2L),
        FUN = .cfunction("FC_SNP2GDS"), as.is="gdsnode", target.node=n1)
    readmode.gdsn(n1)
    .DigestCode(n1, digest, verbose)

    n1 <- add.gdsn(n, "@data", storage="uint8", compress=compress.annotation,
        visible=FALSE)
    .repeat_gds(n1, 1L, nSNP)
    readmode.gdsn(n1)
    .DigestCode(n1, digest, FALSE)

    n1 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n1, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(n, "extra", storage="int16", compress=compress.geno, closezip=TRUE)


    # add a folder for phase information
    if (verbose) cat("    phase")
    n <- addfolder.gdsn(dstfile, "phase")

    n1 <- add.gdsn(n, "data", storage="bit1", valdim=c(nSamp, 0L),
        compress=compress.annotation)
    .repeat_gds(n1, 0L, as.double(nSNP)*nSamp)
    readmode.gdsn(n1)
    .DigestCode(n1, digest, verbose)

    n1 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.annotation, closezip=TRUE)
    put.attr.gdsn(n1, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(n, "extra", storage="bit1", compress=compress.annotation,
        closezip=TRUE)


    # add annotation folder
    n <- addfolder.gdsn(dstfile, "annotation")

    # add annotation/id
    if (verbose) cat("    annotation/id")
    n1 <- add.gdsn(n, "id", storage="string", compress=compress.annotation)
    if (is.null(index.gdsn(srcfile, "snp.rs.id", silent=TRUE)))
        assign.gdsn(n1, index.gdsn(srcfile, "snp.id"))
    else
        assign.gdsn(n1, index.gdsn(srcfile, "snp.rs.id"))
    readmode.gdsn(n1)
    .DigestCode(n1, digest, verbose)

    # add annotation/qual
    n1 <- add.gdsn(n, "qual", storage="float", compress=compress.annotation)
    .repeat_gds(n1, 100.0, nSNP)
    readmode.gdsn(n1)
    .DigestCode(n1, digest, FALSE)

    # add filter
    n1 <- add.gdsn(n, "filter", storage="int32", compress=compress.annotation)
    .repeat_gds(n1, 0L, nSNP)
    readmode.gdsn(n1)
    put.attr.gdsn(n1, "R.class", "factor")
    put.attr.gdsn(n1, "R.levels", c("PASS"))
    put.attr.gdsn(n1, "Description", c("All filters passed"))
    .DigestCode(n1, digest, FALSE)

    # add the INFO field
    n1 <- addfolder.gdsn(n, "info")
    n2 <- index.gdsn(srcfile, "snp.annot", silent=TRUE)
    if (!is.null(n2))
    {
        for (i in ls.gdsn(n2))
        {
            if (verbose) cat("    annotation/info/", i, sep="")
            copyto.gdsn(n1, index.gdsn(n2, i))
            n3 <- index.gdsn(n1, i)
            compression.gdsn(n3, compress.annotation)
           .DigestCode(n3, digest, verbose)
        }
    }

    # add the FORMAT field
    addfolder.gdsn(n, "format")

    # add sample annotation
    if (verbose) cat("    sample.annotation\n")
    n <- addfolder.gdsn(dstfile, "sample.annotation")
    n1 <- index.gdsn(srcfile, "sample.annot", silent=TRUE)
    if (!is.null(n1))
    {
        for (i in ls.gdsn(n1))
        {
            if (verbose) cat("    sample.annotation/", i, sep="")
            copyto.gdsn(n, index.gdsn(n1, i))
            n3 <- index.gdsn(n, i)
            compression.gdsn(n3, compress.annotation)
           .DigestCode(n3, digest, verbose)
        }
    }

    on.exit()
    closefn.gds(srcfile)
    closefn.gds(dstfile)

    ##################################################
    # optimize access efficiency

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.gdsfn))
}



#######################################################################
# Convert a PLINK BED file to a SeqArray GDS file
#

seqBED2GDS <- function(bed.fn, fam.fn, bim.fn, out.gdsfn,
    compress.geno="ZIP_RA", compress.annotation="ZIP_RA",
    optimize=TRUE, digest=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(bed.fn), length(bed.fn)==1L)
    stopifnot(is.character(fam.fn), length(fam.fn)==1L)
    stopifnot(is.character(bim.fn), length(bim.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("PLINK BED to SeqArray GDS Format:\n")
    }

    ##  open and detect bed.fn  ##

    bedfile <- .open_bin(bed.fn)
    on.exit({ .close_conn(bedfile) })
    bed_flag <- .Call(SEQ_ConvBEDFlag, bedfile$con, readBin, new.env())
    if (verbose)
    {
        cat("\tBED file: \"", bed.fn, "\"", sep="")
        if (bed_flag == 0)
            cat(" in the individual-major mode (SNP X Sample)\n")
        else
            cat(" in the SNP-major mode (Sample X SNP)\n")
    }

    ##  read fam.fn  ##

    f <- .open_text(fam.fn, TRUE)
    famD <- read.table(f$con, header=FALSE, stringsAsFactors=FALSE)
    .close_conn(f)

    names(famD) <- c("FamilyID", "InvID", "PatID", "MatID", "Sex", "Pheno")
    if (anyDuplicated(famD$InvID) == 0L)
    {
        sample.id <- famD$InvID
    } else {
        sample.id <- paste(famD$FamilyID, famD$InvID, sep="-")
        if (length(unique(sample.id)) != dim(famD)[1])
            stop("IDs in PLINK BED are not unique!")
    }
    if (verbose)
        cat("\tFAM file: \"", fam.fn, "\" (", nrow(famD), " samples)\n", sep="")

    ##  read bim.fn  ##

    f <- .open_text(bim.fn, TRUE)
    bimD <- read.table(f$con, header=FALSE, stringsAsFactors=FALSE)
    .close_conn(f)
    names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")
    if (verbose)
        cat("\tBIM file: \"", bim.fn, "\" (", nrow(bimD), " variants)\n", sep="")


    ##  create GDS file  ##

    dstfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({ closefn.gds(dstfile) }, add=TRUE)

    put.attr.gdsn(dstfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(dstfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(dstfile, "description")
    put.attr.gdsn(n, "source.format", "PLINK BED Format")

    # add sample.id
    if (verbose) cat("    sample.id")
    n <- add.gdsn(dstfile, "sample.id", sample.id, compress=compress.annotation,
        closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add variant.id
    if (verbose) cat("    variant.id")
    n <- add.gdsn(dstfile, "variant.id", seq_len(nrow(bimD)),
        compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add position
    if (verbose) cat("    position")
    n <- add.gdsn(dstfile, "position", bimD$pos, storage="int32",
        compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add chromosome
    if (verbose) cat("    chromosome")
    n <- add.gdsn(dstfile, "chromosome", bimD$chr, storage="string",
        compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add allele
    if (verbose) cat("    allele")
    n <- add.gdsn(dstfile, "allele", paste(bimD$allele1, bimD$allele2, sep=","),
        storage="string", compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add a folder for genotypes
    if (verbose) cat("    genotype")
    n <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(n, "VariableName", "GT")
    put.attr.gdsn(n, "Description", "Genotype")

    # add genotypes to genotype/data
    vg <- add.gdsn(n, "data", storage="bit2",
        valdim=c(2L, ifelse(bed_flag==0L, nrow(bimD), nrow(famD)), 0L),
        compress=compress.geno)
    # convert
    .Call(SEQ_ConvBED2GDS, vg, ifelse(bed_flag==0L, nrow(famD), nrow(bimD)),
        bedfile$con, readBin, new.env())
    readmode.gdsn(vg)

    n1 <- add.gdsn(n, "@data", storage="uint8",
        compress=ifelse(bed_flag==0L, "", compress.annotation),
        visible=FALSE)
    .repeat_gds(n1, 1L, nrow(bimD))
    readmode.gdsn(n1)
    .DigestCode(n1, digest, FALSE)

    n1 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n1, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(n, "extra", storage="int16", compress=compress.geno, closezip=TRUE)

    # sync file
    sync.gds(dstfile)

    # close the BED file
    on.exit({ closefn.gds(dstfile) })
    .close_conn(bedfile)

    if (bed_flag == 0L)
    {
        cat(" (transposed)")
        permdim.gdsn(vg, c(2L,1L))
    }
    .DigestCode(vg, digest, verbose)

    # add a folder for phase information
    if (verbose) cat("    phase")
    n <- addfolder.gdsn(dstfile, "phase")

    n1 <- add.gdsn(n, "data", storage="bit1", valdim=c(nrow(famD), 0L),
        compress=compress.annotation)
    .repeat_gds(n1, 0L, as.double(nrow(bimD))*nrow(famD))
    readmode.gdsn(n1)
    .DigestCode(n1, digest, TRUE)

    n1 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.annotation, closezip=TRUE)
    put.attr.gdsn(n1, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(n, "extra", storage="bit1", compress=compress.annotation,
        closezip=TRUE)


    # add annotation folder
    n <- addfolder.gdsn(dstfile, "annotation")

    # add annotation/id
    n1 <- add.gdsn(n, "id", bimD$snp.id, storage="string",
        compress=compress.annotation, closezip=TRUE)
    if (verbose) cat("    annotation/id")
    .DigestCode(n1, digest, verbose)

    # add annotation/qual
    n1 <- add.gdsn(n, "qual", storage="float", compress=compress.annotation)
    .repeat_gds(n1, 100.0, nrow(bimD))
    readmode.gdsn(n1)
    if (verbose) cat("    annotation/qual")
    .DigestCode(n1, digest, verbose)

    # add filter
    n1 <- add.gdsn(n, "filter", storage="int32", compress=compress.annotation)
    .repeat_gds(n1, 0L, nrow(bimD))
    readmode.gdsn(n1)
    put.attr.gdsn(n1, "R.class", "factor")
    put.attr.gdsn(n1, "R.levels", c("PASS"))
    put.attr.gdsn(n1, "Description", c("All filters passed"))
    if (verbose) cat("    annotation/filter")
    .DigestCode(n1, digest, verbose)

    # add the INFO field
    addfolder.gdsn(n, "info")
    # add the FORMAT field
    addfolder.gdsn(n, "format")

    # add sample annotation
    if (verbose) cat("    sample.annotation\n")
    n <- addfolder.gdsn(dstfile, "sample.annotation")

    n1 <- add.gdsn(n, "family", famD$FamilyID, compress=compress.annotation,
        closezip=TRUE)
    .DigestCode(n1, digest, FALSE)

    n1 <- add.gdsn(n, "father", famD$PatID, compress=compress.annotation,
        closezip=TRUE)
    .DigestCode(n1, digest, FALSE)

    n1 <- add.gdsn(n, "mother", famD$MatID, compress=compress.annotation,
        closezip=TRUE)
    .DigestCode(n1, digest, FALSE)

    sex <- rep("", length(sample.id))
    sex[famD$Sex==1L] <- "M"; sex[famD$Sex==2L] <- "F"
    n1 <- add.gdsn(n, "sex", sex, compress=compress.annotation, closezip=TRUE)
    .DigestCode(n1, digest, FALSE)

    n1 <- add.gdsn(n, "phenotype", famD$Pheno, compress=compress.annotation,
        closezip=TRUE)
    .DigestCode(n1, digest, FALSE)

    # sync file
    sync.gds(dstfile)


    ##################################################
    # optimize access efficiency

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    on.exit()
    closefn.gds(dstfile)
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.gdsfn))
}
