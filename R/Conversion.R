#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Big Data Management of Sequencing-based Genetic Variants
#



#######################################################################
# Convert a VCF (sequence) file to a GDS file
#######################################################################

#######################################################################
# Parse the header of a VCF file
# http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
#

seqVCF.Header <- function(vcf.fn)
{
    # check
    stopifnot(is.character(vcf.fn))

    #########################################################
    # open the vcf file

    # parse and determine how many copies of genomes: haploid, diploid or others
    geno.text <- NULL

    n <- 0L
    for (i in 1:length(vcf.fn))
    {
        opfile <- file(vcf.fn[i], open="rt")
        on.exit(close(opfile))

        # read header
        ans <- NULL
        while (length(s <- readLines(opfile, n=1)) > 0)
        {
            n <- n + 1L
            if (substr(s, 1, 6) != "#CHROM")
            {
                s <- substring(s, 3)
                ans <- c(ans, s)
            } else {
                s <- readLines(opfile, n=1)
                if (length(s) > 0)
                {
                    ss <- scan(text=s, what=character(), sep="\t", quiet=TRUE)
                    geno.text <- c(geno.text, ss[-(1:9)])
                }
                break
            }
            if ((n %% 10000L) == 0)
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

    value.string <- function(txt)
    {
        # split by "="
        x <- as.integer(regexpr("=", txt, fixed=TRUE))
        rv <- matrix("", nrow=length(txt), ncol=2)
        for (i in 1:length(txt))
        {
            if (x[i] > 0)
            {
                rv[i,1] <- substring(txt[i], 1, x[i]-1)
                rv[i,2] <- substring(txt[i], x[i]+1)
            }
        }
        rv
    }


    #########################################################
    # get names and values, length(txt) == 1

    val.data.frame <- function(txt, check.name)
    {
        if (substr(txt, 1, 1) == "<")
            txt <- substring(txt, 2)
        if (substr(txt, nchar(txt), nchar(txt)) == ">")
            txt <- substr(txt, 1, nchar(txt) - 1)

        rv <- value.string(scan(text=txt, what = character(),
            sep = ",", quote="\"", quiet=TRUE))

        # check
        ans <- list()
        if (!is.null(check.name))
        {
            if (nrow(rv) == length(check.name))
            {
                for (i in 1:nrow(rv))
                {
                    if (tolower(rv[i,1]) == tolower(check.name[i]))
                    {
                        ans[[ check.name[i] ]] <- rv[i,2]
                    } else
                        return(NULL)
                }
            }
        } else {
            for (i in 1:nrow(rv))
                ans[[ rv[i,1] ]] <- rv[i,2]
        }
        as.data.frame(ans, stringsAsFactors=FALSE)
    }

    #########################################################
    # check number

    check.num <- function(number)
    {
        if (!(number %in% c("A", "G", ".")))
        {
            N <- suppressWarnings(as.integer(number))
            if (is.finite(N))
                N >= 0
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
        num.ploidy <- as.integer(names(which.max(tab)))
    } else
        num.ploidy <- as.integer(NA)

    if (is.null(ans))
    {
        rv <- list(fileformat="unknown", info=NULL, filter=NULL, format=NULL,
            alt=NULL, contig=NULL, assembly=NULL, header=NULL,
            num.ploidy = num.ploidy)
        class(rv) <- "SeqVCFHeaderClass"
        return(rv)
    }

    ans <- value.string(ans)
    ans <- data.frame(id=ans[,1], value=ans[,2], stringsAsFactors=FALSE)

    #########################################################
    # fileformat=VCFv4.1

    s <- ans$value[ans$id == "fileformat"]
    if (length(s) < 1)
        stop("no defined fileformat!")
    else if (length(s) > 1)
        stop("Multiple fileformat!")
    fileformat <- s[1]
    ans <- ans[ans$id != "fileformat", ]

    #########################################################
    # assembly=url

    s <- ans$value[ans$id == "assembly"]
    assembly <- if (length(s) > 0) s else NULL
    ans <- ans[ans$id != "assembly", ]


    #########################################################
    # INFO=<ID=ID,Number=number,Type=type,Description="description">

    INFO <- NULL
    s <- ans$value[ans$id == "INFO"]
    for (i in seq_len(length(s)))
    {
        v <- val.data.frame(s[i], c("ID", "Number", "Type", "Description"))
        if (!is.null(v))
        {
            if (!is.element(tolower(v$Type), c("integer", "float", "flag", "character", "string")))
                stop(sprintf("INFO=%s", s[i]))
            if (!check.num(v$Number))
                stop(sprintf("INFO=%s", s[i]))
            INFO <- rbind(INFO, v)
        } else
            stop(sprintf("INFO=%s", s[i]))
    }
    ans <- ans[ans$id != "INFO", ]


    #########################################################
    # FILTER=<ID=ID,Description="description">

    FILTER <- NULL
    s <- ans$value[ans$id == "FILTER"]
    for (i in seq_len(length(s)))
    {
        v <- val.data.frame(s[i], c("ID", "Description"))
        if (!is.null(v))
            FILTER <- rbind(FILTER, v)
        else
            stop(sprintf("FILTER=%s", s[i]))
    }
    ans <- ans[ans$id != "FILTER", ]


    #########################################################
    # FORMAT=<ID=ID,Number=number,Type=type,Description="description">

    FORMAT <- NULL
    s <- ans$value[ans$id == "FORMAT"]
    for (i in seq_len(length(s)))
    {
        v <- val.data.frame(s[i], c("ID", "Number", "Type", "Description"))
        if (!is.null(v))
        {
            if (!is.element(tolower(v$Type), c("integer", "float", "character", "string")))
                stop(sprintf("FORMAT=%s", s[i]))
            if (!check.num(v$Number))
                stop(sprintf("FORMAT=%s", s[i]))
            FORMAT <- rbind(FORMAT, v)
        } else
            stop(sprintf("FORMAT=%s", s[i]))
    }
    ans <- ans[ans$id != "FORMAT", ]


    #########################################################
    # ALT=<ID=type,Description=description>

    ALT <- NULL
    s <- ans$value[ans$id == "ALT"]
    for (i in seq_len(length(s)))
    {
        v <- val.data.frame(s[i], c("ID", "Description"))
        if (!is.null(v))
        {
            ALT <- rbind(ALT, v)
        } else
            stop(sprintf("ALT=%s", s[i]))
    }
    ans <- ans[ans$id != "ALT", ]


    #########################################################
    # contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

    contig <- NULL
    s <- ans$value[ans$id == "contig"]
    for (i in seq_len(length(s)))
    {
        v <- val.data.frame(s[i], NULL)
        if (!is.null(v))
        {
            contig <- rbind(contig, v)
        } else
            stop(sprintf("contig=%s", s[i]))
    }
    ans <- ans[ans$id != "contig", ]


    #########################################################
    # output

    rv <- list(fileformat=fileformat, info=INFO, filter=FILTER, format=FORMAT,
        alt=ALT, contig=contig, assembly=assembly, header=ans,
        num.ploidy = num.ploidy)
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
    stopifnot(length(vcf.fn) == 1)

    # open the vcf file
    opfile <- file(vcf.fn[1], open="rt")
    on.exit(close(opfile))

    # read header
    samp.id <- NULL
    while (length(s <- readLines(opfile, n=1)) > 0)
    {
        if (substr(s, 1, 6) == "#CHROM")
        {
            samp.id <- scan(text=s, what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
            break
        }
    }
    if (is.null(samp.id))
        stop("Error VCF format: invalid sample id!")

    return(samp.id)
}



#######################################################################
# Convert a VCF (sequencing) file to a GDS file
#

seqVCF2GDS <- function(vcf.fn, out.fn, header=NULL,
    genotype.var.name="GT", compress.option=seqCompress.Option(),
    info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr",
    raise.error=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf.fn) & is.vector(vcf.fn))
    stopifnot(length(vcf.fn) > 0)

    stopifnot(is.character(out.fn) & is.vector(out.fn))
    stopifnot(length(out.fn) == 1)

    stopifnot(is.null(header) | inherits(header, "SeqVCFHeaderClass"))
    stopifnot(inherits(compress.option, "SeqGDSCompressFlagClass"))

    stopifnot(is.character(genotype.var.name) & is.vector(genotype.var.name))
    stopifnot(length(genotype.var.name) == 1)
    stopifnot(!is.na(genotype.var.name))

    stopifnot(is.null(info.import) |
        (is.character(info.import) & is.vector(info.import)))
    stopifnot(is.null(fmt.import) |
        (is.character(fmt.import) & is.vector(fmt.import)))

    stopifnot(is.character(ignore.chr.prefix))

    stopifnot(is.logical(raise.error) & is.vector(raise.error))
    stopifnot(length(raise.error) == 1)

    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose) == 1)


    # check sample id
    samp.id <- NULL
    for (i in 1:length(vcf.fn))
    {
        if (is.null(samp.id))
        {
            samp.id <- seqVCF.SampID(vcf.fn[i])
            if (length(samp.id) <= 0)
                stop("There is no sample in the VCF file.")
        } else {
            tmp <- seqVCF.SampID(vcf.fn[i])
            if (length(samp.id) != length(tmp))
                stop(sprintf("The file '%s' has different sample id.", vcf.fn[i]))
            if (!all(samp.id == tmp))
                stop(sprintf("The file '%s' has different sample id.", vcf.fn[i]))
        }
    }


    # define functions
    compress <- function(var.name)
    {
        if (var.name %in% names(compress.option))
            compress.option[[var.name]]
        else
            ""
    }


    ##################################################
    # parse the header of VCF

    if (is.null(header))
        header <- seqVCF.Header(vcf.fn)
    # 'seqVCF.Header'
    # returns list(fileformat=fileformat, info=INFO, filter=FILTER,
    #              format=FORMAT, alt=ALT, contig=contig, assembly=assembly,
    #              header=ans)

    if (verbose)
    {
        cat("The Variant Call Format (VCF) header:\n")
        cat("\tfile format: ", header$fileformat, "\n", sep="")
        cat("\tthe number of sets of chromosomes (ploidy): ",
            header$num.ploidy, "\n", sep="")
        cat("\tthe number of samples: ", length(samp.id), "\n", sep="")
    }

    # check header
    tmp <- FALSE
    if (!is.null(header$format))
    {
        if (nrow(header$format) > 0)
            tmp <- TRUE
    }
    if (tmp)
    {
        geno_idx <- which(header$format$ID == genotype.var.name)
        if (length(geno_idx) <= 0)
        {
            stop(sprintf(
            "There is no variable in the FORMAT field according to '%s'.",
                genotype.var.name))
        } else if (length(geno_idx) > 1)
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
    n <- addfolder.gdsn(gfile, "description")
    put.attr.gdsn(n, "sequence.variant.format", "v1.0")

    put.attr.gdsn(n, "vcf.fileformat", header$fileformat)
    if (!is.null(header$assembly))
        put.attr.gdsn(n, "vcf.assembly", header$assembly)
    if (!is.null(header$filter))
    {
        if (nrow(header$filter) > 0)
            put.attr.gdsn(add.gdsn(n, "vcf.filter", header$filter), "R.invisible")
    }
    if (!is.null(header$alt))
    {
        if (nrow(header$alt) > 0)
            put.attr.gdsn(add.gdsn(n, "vcf.alt", header$alt), "R.invisible")
    }
    if (!is.null(header$contig))
    {
        if (nrow(header$contig) > 0)
            put.attr.gdsn(add.gdsn(n, "vcf.contig", header$contig), "R.invisible")
    }
    if (!is.null(header$header))
    {
        if (nrow(header$header) > 0)
            put.attr.gdsn(add.gdsn(n, "vcf.header", header$header), "R.invisible")
    }


    ##################################################
    # add sample id

    nSamp <- length(samp.id)
    add.gdsn(gfile, "sample.id", samp.id, compress=compress("sample.id"),
        closezip=TRUE)


    ##################################################
    # add basic site information

    # add variant.id
    add.gdsn(gfile, "variant.id", storage="int32", valdim=c(0),
        compress=compress("variant.id"))

    # add position
    # TODO: need to check whether position can be stored in 'int32'
    add.gdsn(gfile, "position", storage="int32", valdim=c(0),
        compress=compress("position"))

    # add chromosome
    add.gdsn(gfile, "chromosome", storage="string", valdim=c(0),
        compress=compress("chromosome"))

    # add allele
    add.gdsn(gfile, "allele", storage="string", valdim=c(0),
        compress=compress("allele"))

    # add a folder for genotypes
    varGeno <- addfolder.gdsn(gfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", genotype.var.name[1])
    put.attr.gdsn(varGeno, "Description", geno_format$Description[1])

    # add data to the folder of genotype
    if (header$num.ploidy > 1)
    {
        geno.node <- add.gdsn(varGeno, "data", storage="bit2",
            valdim=c(header$num.ploidy, nSamp, 0),
            compress=compress("genotype"))
    } else {
        geno.node <- add.gdsn(varGeno, "data", storage="bit2",
            valdim=c(nSamp, 0), compress=compress("genotype"))
    }
    node <- add.gdsn(varGeno, "@data", storage="uint8", valdim=c(0),
        compress=compress("genotype"))
    put.attr.gdsn(node, "R.invisible")

    node <- add.gdsn(varGeno, "extra.index", storage="int32", valdim=c(3,0),
        compress=compress("genotype.extra"))
    put.attr.gdsn(node, "R.colnames", c("sample.index", "variant.index", "length"))
    add.gdsn(varGeno, "extra", storage="int16", valdim=c(0),
        compress=compress("genotype.extra"))


    # add phase folder
    if (header$num.ploidy > 1)
    {
        varPhase <- addfolder.gdsn(gfile, "phase")
        # add data
        if (header$num.ploidy > 2)
        {
            add.gdsn(varPhase, "data", storage="bit1",
                valdim=c(header$num.ploidy-1, nSamp, 0),
                compress=compress("phase"))
        } else {
            add.gdsn(varPhase, "data", storage="bit1", valdim=c(nSamp, 0),
                compress=compress("phase"))
        }

        node <- add.gdsn(varPhase, "extra.index", storage="int32",
            valdim=c(3,0), compress=compress("phase.extra"))
        put.attr.gdsn(node, "R.colnames",
            c("sample.index", "variant.index", "length"))
        add.gdsn(varPhase, "extra", storage="bit1", valdim=c(0),
            compress=compress("phase.extra"))
    }


    ##################################################

    # add annotation folder
    varAnnot <- addfolder.gdsn(gfile, "annotation")

    # add id
    add.gdsn(varAnnot, "id", storage="string", valdim=c(0),
        compress=compress("id"))

    # add qual
    add.gdsn(varAnnot, "qual", storage="float", valdim=c(0),
        compress=compress("qual"))

    # add filter
    varFilter <- add.gdsn(varAnnot, "filter", storage="int32", valdim=c(0),
        compress=compress("filter"))


    ##################################################
    # VCF INFO

    varInfo <- addfolder.gdsn(varAnnot, "info")

    if (!is.null(header$info))
    {
        if (nrow(header$info) > 0)
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
                    integer = { mode <- "int32"; int_type[i] <- 1 },
                    float = { mode <- "float32"; int_type[i] <- 2 },
                    flag = { mode <- "bit1"; int_type[i] <- 3 },
                    character = { mode <- "string"; int_type[i] <- 4 },
                    string = { mode <- "string"; int_type[i] <- 4 },
                    stop(sprintf("Unknown INFO Type: %s", header$info$Type[i]))
                )

                # INFO Number
                s <- header$info$Number[i]
                if (!(s %in% c(".", "A", "G")))
                {
                    initdim <- as.integer(s)
                    if (mode == "bit1")
                    {
                        if (initdim != 0)
                        {
                            print(header$info[i, ])
                            stop("The length of 'Flag' type should be ZERO!")
                        }
                    } else {
                        if (initdim <= 0)
                        {
                            print(header$info[i, ])
                            stop("The length should be >0.")
                        } else if (initdim > 1)
                            initdim <- c(initdim, 0)
                        else
                            initdim <- c(0)
                    }
                } else {
                    initdim <- c(0)
                    if (s == ".")
                        int_num[i] <- -1
                    else if (s == "A")
                        int_num[i] <- -2
                    else if (s == "G")
                        int_num[i] <- -3
                    else
                        stop("seqVCF2GDS: internal error!")
                }

                # add
                if (import.flag[i])
                {
                    node <- add.gdsn(varInfo, header$info$ID[i], storage=mode,
                        valdim=initdim, compress=compress("annotation"))
                    put.attr.gdsn(node, "Number", header$info$Number[i])
                    put.attr.gdsn(node, "Type", header$info$Type[i])
                    put.attr.gdsn(node, "Description", header$info$Description[i])

                    if (s %in% c(".", "A", "G"))
                    {
                        node <- add.gdsn(varInfo,
                            paste("@", header$info$ID[i], sep=""),
                            storage="int32", valdim=c(0),
                            compress=compress("info"))
                        put.attr.gdsn(node, "R.invisible")
                    }
                }
            }

            header$info$int_type <- as.integer(int_type)
            header$info$int_num  <- as.integer(int_num)
            header$info$import.flag <- import.flag
        }
    }


    ##################################################
    # VCF Variants

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
        int_type <- integer(0)
        int_num <- integer(0)
        import.flag <- logical(0)
    }

    for (i in seq_len(length(int_type)))
    {
        # FORMAT Type
        switch(tolower(header$format$Type[i]),
            integer = { mode <- "int32"; int_type[i] <- 1 },
            float = { mode <- "float32"; int_type[i] <- 2 },
            character = { mode <- "string"; int_type[i] <- 4 },
            string = { mode <- "string"; int_type[i] <- 4 },
            stop(sprintf("Unknown FORMAT Type: %s", header$format$Type[i]))
        )

        # FORMAT Number
        s <- header$format$Number[i]
        if (!(s %in% c(".", "A", "G")))
        {
            initdim <- as.integer(s)
            if (initdim <= 0)
            {
                print(header[, i])
                stop("The length should be >0.")
            } else if (initdim > 1)
                initdim <- c(initdim, nSamp, 0)
            else
                initdim <- c(nSamp, 0)
        } else {
            initdim <- c(nSamp, 0)
            if (s == ".")
                int_num[i] <- -1
            else if (s == "A")
                int_num[i] <- -2
            else if (s == "G")
                int_num[i] <- -3
            else
                stop("seqVCF2GDS: internal error!")
        }

        # add
        if (import.flag[i])
        {
            node <- addfolder.gdsn(varFormat, header$format$ID[i])
            put.attr.gdsn(node, "Number", header$format$Number[i])
            put.attr.gdsn(node, "Type", header$format$Type[i])
            put.attr.gdsn(node, "Description", header$format$Description[i])

            add.gdsn(node, "data", storage=mode, valdim=initdim,
                compress=compress("format"))
            tmp <- add.gdsn(node, "@data", storage="int32", valdim=c(0),
                compress=compress("format"))
            put.attr.gdsn(tmp, "R.invisible")
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
    for (i in 1:length(vcf.fn))
    {
        opfile <- file(vcf.fn[i], open="rt")
        on.exit({ closefn.gds(gfile); close(opfile) })

        if (verbose)
            cat("Parsing \"", vcf.fn[i], "\" ...\n", sep="")

        ##################################################
        # define wrapping R function for the C code 'SEQ_Parse_VCF4'
        #   because 'getConnection' in Rconnections.h is hidden
        # see https://stat.ethz.ch/pipermail/r-devel/2007-April/045264.html

        # call C function
        v <- .Call(SEQ_Parse_VCF4, vcf.fn[i], header, gfile$root,
            list(sample.num = as.integer(length(samp.id)),
                genotype.var.name = genotype.var.name,
                raise.error = raise.error, verbose = verbose),
            readLines, opfile, 512L,  # 'readLines(opfile, 512L)'
            ignore.chr.prefix, new.env())
        if (length(v) > 0)
        {
            put.attr.gdsn(varFilter, "R.class", "factor")
            put.attr.gdsn(varFilter, "R.levels", v)
        }

        if (verbose)
            print(geno.node)

        on.exit()
        close(opfile)
    }

    closefn.gds(gfile)


    ##################################################
    # optimize access efficiency

    if (verbose)
    {
        cat("Done.\n")
        cat("Optimize the access efficiency ...\n")
    }
    cleanup.gds(out.fn, verbose=verbose)

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a GDS file to a VCF (sequencing) file
#

seqGDS2VCF <- function(gdsfile, vcf.fn, info.var=NULL, fmt.var=NULL,
    verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(vcf.fn) & (length(vcf.fn)==1))
    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))

    # get a summary
    z <- seqSummary(gdsfile, check="none", verbose=FALSE)

    # the INFO field
    if (!is.null(info.var))
    {
        s <- z$info$var.name
        if (is.null(s)) s <- character()
        if (length(setdiff(info.var, s)) > 0)
            stop(paste("Not exist:", paste(setdiff(info.var, s), collapse=",")))
        if (length(info.var) > 0)
            z$info <- z$info[match(info.var, z$info$var.name), ]
        else
            z$info <- list()
    }

    # the FORMAT field
    if (!is.null(fmt.var))
    {
        s <- z$format$var.name
        if (is.null(s)) s <- character()
        if (length(setdiff(fmt.var, s)) > 0)
            stop(paste("Not exist:", paste(setdiff(fmt.var, s), collapse=",")))
        if (length(fmt.var) > 0)
            z$format <- z$format[match(fmt.var, z$format$var.name), ]
        else
            z$format <- list()
    }


    ## double quote text if needed
    dq <- function(s, text=FALSE)
    {
        .Call(SEQ_Quote, s, text)
    }


    ######################################################
    # create an output text file

    vcf.fn <- vcf.fn[1]
    ext <- substring(vcf.fn, nchar(vcf.fn)-2)
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
        cat("Output: ", vcf.fn, "\n", sep="")
        cat("The INFO field: ", paste(z$info$var.name, collapse=", "),
            "\n", sep="")
        cat("The FORMAT field: ", paste(z$format$var.name, collapse=", "),
            "\n", sep="")
    }

    ######################################################
    # write the header
    a <- get.attr.gdsn(index.gdsn(gdsfile, "description"))

    # fileformat
    if (is.null(a$vcf.fileformat))
        a$vcf.fileformat <- "VCFv4.1"
    cat("##fileformat=", a$vcf.fileformat, "\n", sep="", file=ofile)

    # fileDate
    cat("##fileDate=", format(Sys.time(), "%Y%m%d"), "\n", sep="", file=ofile)

    # program, source
    if (is.null(a$sequence.variant.format))
        a$sequence.variant.format <- "v1.0"
    cat("##source=SeqArray_RPackage_", a$sequence.variant.format, "\n", sep="", file=ofile)

    # assembly
    if (!is.null(a$vcf.assembly))
        cat("##assembly=", dq(a$vcf.assembly), "\n", sep="", file=ofile)

    # ALT=<ID=type,Description=description>
    n <- index.gdsn(gdsfile, "description/vcf.alt", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        s <- sprintf("##ALT=<ID=%s,Description=%s>",
            as.character(dat$ID), dq(dat$Description))
        writeLines(s, con=ofile)
    }

    # contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
    n <- index.gdsn(gdsfile, "description/vcf.contig", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        nm <- names(dat)
        for (i in 1:nrow(dat))
        {
            s <- NULL
            for (j in 1:ncol(dat))
                s[j] <- paste(nm[j], "=", dq(dat[i,j]), sep="")
            s <- paste(s, collapse=",")
            cat("##contig=<", s, ">\n", sep="", file=ofile)
        }
    }

    # the INFO field
    for (nm in z$info$var.name)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, paste("annotation/info/", nm, sep="")))
        cat(sprintf("##INFO=<ID=%s,Number=%s,Type=%s,Description=%s>\n",
            nm, dq(a$Number), dq(a$Type), dq(a$Description, TRUE)), file=ofile)
    }

    # the FILTER field
    n <- index.gdsn(gdsfile, "description/vcf.filter", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        for (i in 1:nrow(dat))
        {
            cat(sprintf("##FILTER=<ID=%s,Description=%s>\n",
                dq(dat$ID[i]), dq(dat$Description[i], TRUE)), file=ofile)
        }
    }

    # the FORMAT field
    a <- get.attr.gdsn(index.gdsn(gdsfile, "genotype"))
    cat(sprintf("##FORMAT=<ID=%s,Number=1,Type=String,Description=%s>\n",
        a$VariableName, dq(a$Description, TRUE)), file=ofile)
    for (nm in z$format$var.name)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, paste("annotation/format/", nm, sep="")))
        cat(sprintf("##FORMAT=<ID=%s,Number=%s,Type=%s,Description=%s>\n",
            nm, dq(a$Number), dq(a$Type), dq(a$Description, TRUE)), file=ofile)
    }

    # others ...
    n <- index.gdsn(gdsfile, "description/vcf.header", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        for (i in 1:nrow(dat))
        {
            s <- dat[i,1]
            if (!(s %in% c("fileDate", "source")))
                cat("##", s, "=", dq(dat[i,2]), "\n", sep="", file=ofile)
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
    if (!is.null(z$info$var.name))
        nm.info <- paste("annotation/info/", z$info$var.name, sep="")
    else
        nm.info <- c()
    # the FORMAT field
    if (!is.null(z$format$var.name))
        nm.format <- paste("annotation/format/", z$format$var.name, sep="")
    else
        nm.format <- c()

    # initialize the variable length of INFO
    len.info <- NULL
    for (n in nm.info)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, n))
        a$Number <- if (is.null(a$Number)) "." else a$Number[1]
        len.info <- c(len.info, a$Number)
    }
    len.info <- suppressWarnings(as.integer(len.info))

    # initialize the variable length of FORMAT
    len.fmt <- NULL
    for (n in nm.format)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile, n))
        a$Number <- if (is.null(a$Number)) "." else a$Number[1]
        len.fmt <- c(len.fmt, a$Number)
    }
    len.fmt <- suppressWarnings(as.integer(len.fmt))

    # call C function
    .Call(SEQ_InitOutVCF4, len.info, len.fmt)


    # variable names
    nm <- c("chromosome", "position", "annotation/id", "allele",
        "annotation/qual", "annotation/filter", "genotype", "phase")
    if (length(nm.info) > 0) nm <- c(nm, nm.info)
    if (length(nm.format) > 0) nm <- c(nm, nm.format)

    s <- c("chr", "pos", "id", "allele", "qual", "filter", "geno", "phase")
    # the INFO field
    if (length(nm.info) > 0)
        s <- c(s, paste("info.", z$info$var.name, sep=""))
    # the FORMAT field
    if (length(nm.format) > 0)
        s <- c(s, paste("fmt.", z$format$var.name, sep=""))
    names(nm) <- s

    # output lines variant by variant
    seqApply(gdsfile, nm, margin="by.variant", as.is="none",
        FUN = function(x) {
            s <- .Call(SEQ_OutVCF4, x)
            cat(s, file=ofile)
    })

    if (verbose)
        cat("Done.\n")

    # output
    invisible(normalizePath(vcf.fn))
}



#######################################################################
# Convert a Sequencing GDS file to a SNP GDS file
#

seqGDS2SNP <- function(gdsfile, out.gdsfn,
    compress.annotation="ZIP_RA.max", compress.geno="ZIP_RA.max",
    verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.gdsfn))
    stopifnot(is.vector(out.gdsfn) & (length(out.gdsfn)==1L))
    stopifnot(is.logical(verbose))

    if (verbose)
        cat("Sequencing GDS to SNP GDS Format:\n")

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
        valdim=c(length(sampid), 0), compress=compress.geno)
    put.attr.gdsn(gGeno, "sample.order")

    seqApply(gdsfile, "genotype", margin="by.variant", as.is="none",
        FUN = function(x) {
            g <- colSums(x==0)
            g[is.na(g)] <- 3L
            append.gdsn(gGeno, g)
    })

    readmode.gdsn(gGeno)

    if (verbose)
        cat("Done.\n")

    # output
    invisible(normalizePath(out.gdsfn))
}
