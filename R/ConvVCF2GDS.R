#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Data Management of Large-scale Whole-Genome Sequence Variant Calls
#



#######################################################################
# Format Conversion: VCF -> GDS
#######################################################################

# get split from the total count
.file_split <- function(count, pnum, start=1L, multiple=8L)
{
    .Call(SEQ_VCF_Split, start, count, pnum, multiple)
}

# need unique temporary file names
.get_temp_fn <- function(pnum, fn, tmpdir)
{
    ptmpfn <- character()
    while (length(ptmpfn) < pnum)
    {
        s <- tempfile(pattern=sprintf("%s_tmp%02d_", fn, length(ptmpfn)+1L),
            tmpdir=tmpdir)
        file.create(s)
        if (!(s %in% ptmpfn)) ptmpfn <- c(ptmpfn, s)
    }
    ptmpfn
}

# the interval (in the number of variants) of the file offsets recorded in
# .vcf_count_offset(), see the internal 'seek' attribute of 'start' in
# seqVCF2GDS()
.vcf_offset_step <- 1024L

# whether a VCF file is a plain text file (i.e., not compressed)
.is_text_file <- function(fn)
{
    if (grepl("^(ftp|http|https)://", fn, ignore.case=TRUE)) return(FALSE)
    f <- file(fn, "rb")
    on.exit(close(f))
    identical(readBin(f, "raw", 2L), charToRaw("##"))
}

# whether the file offsets of a VCF file can be used for seeking: TRUE if it
# is a plain text file (i.e., not compressed), or a BGZF file (i.e., bgzip)
.vcf_seekable <- function(fn)
{
    if (grepl("^(ftp|http|https)://", fn, ignore.case=TRUE)) return(FALSE)
    .is_text_file(fn) || .bgzf_is(fn)
}

# count the number of variants, and record the file offset of every
# '.vcf_offset_step' variants; a BGZF file is counted in parallel, since it
# can be split at the block boundaries without decompressing anything;
# return list(num, offset, uoffset, index, line)
.vcf_count_offset <- function(fn, parallel=FALSE, verbose=FALSE)
{
    if (!.bgzf_is(fn))
    {
        # not a BGZF file, no way to seek, so it has to be read sequentially
        infile <- file(fn, "rt")
        on.exit(close(infile))
        return(.Call(SEQ_VCF_NumLines, infile, TRUE, .vcf_offset_step, NULL,
            verbose))
    }

    # the BGZF file is read in the C code, to record the virtual offsets
    pnum <- .NumParallel(parallel)
    rg <- if (pnum > 1L) .Call(SEQ_bgzip_split, fn, pnum) else NULL
    if (NROW(rg) <= 1L)
    {
        return(.Call(SEQ_VCF_NumLines, fn, TRUE, .vcf_offset_step, NULL,
            verbose))
    }

    # count each part of the file in parallel
    lst <- seqParApply(parallel, seq_len(NROW(rg)),
        FUN = function(i, fn, rg, step)
        {
            # only the first part has the VCF header
            .Call(SEQ_VCF_NumLines, fn, i==1L, step, rg[i, ], FALSE)
        }, fn=fn, rg=rg, step=.vcf_offset_step)

    # combine: the variant indices of each part are shifted by the number of
    #     the variants in the previous parts
    num <- vapply(lst, function(z) z$num, 0)
    pre <- c(0, cumsum(num)[-length(num)])
    list(num = sum(num),
        offset  = unlist(lapply(lst, function(z) z$offset), use.names=FALSE),
        uoffset = unlist(lapply(lst, function(z) z$uoffset), use.names=FALSE),
        index   = unlist(lapply(seq_along(lst), function(i)
            lst[[i]]$index + pre[i]), use.names=FALSE),
        line = lst[[1L]]$line)
}

# for each starting variant index in 'start', return c(file index, file
# offset, within-block offset, the variant index at that offset, the VCF line
# number at that offset), or NULL if the file offset is not available;
# 'offset' is a list of the values returned by .vcf_count_offset() for each
# of the VCF files
.vcf_seek_list <- function(start, variant_count, offset)
{
    cum <- cumsum(variant_count)
    lapply(start, function(st)
    {
        i <- which(st <= cum)
        if (!length(i)) return(NULL)
        i <- i[1L]
        z <- offset[[i]]
        if (is.null(z)) return(NULL)
        # the number of variants in the previous files
        pre <- if (i > 1L) cum[i-1L] else 0
        # z$offset[k] is the file offset of the z$index[k]-th variant
        k <- findInterval(st - pre, z$index)
        if (k < 1L) return(NULL)
        m <- z$index[k]
        c(i, z$offset[k], z$uoffset[k], pre + m, z$line + m - 1)
    })
}


#######################################################################
# Parse the header of a VCF file
# http://www.1000genomes.org/wiki/analysis/variant-call-format
#

seqVCF_Header <- function(vcf.fn, getnum=FALSE, use_Rsamtools=NA,
    parallel=FALSE, verbose=TRUE)
{
    # note: 'use_Rsamtools' is deprecated and ignored, since counting the
    #   variants no longer needs the Rsamtools package
    # check
    if (!inherits(vcf.fn, "connection"))
    {
        stopifnot(is.character(vcf.fn))
        ilist <- seq_along(vcf.fn)
    } else {
        ilist <- 1L
    }
    stopifnot(is.logical(getnum), length(getnum)==1L)
    stopifnot(is.logical(use_Rsamtools), length(use_Rsamtools)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    #########################################################
    # open the vcf file

    # parse and determine how many copies of genomes: haploid, diploid or others
    geno.text <- NULL
    nSample <- nVariant <- 0L
    samp.id <- NULL

    # add a message if fails
    error_fn <- FALSE
    on.exit({
        if (!isFALSE(error_fn))
        {
            if (is.character(error_fn))
                message("Parsing the VCF header fails: ", error_fn)
            else
                message("Parsing the VCF header from a connection object fails.")
            error_fn <- FALSE
        }
    })

    n <- 0L
    for (i in ilist)
    {
        is_vcf_fn <- FALSE
        if (!inherits(vcf.fn, "connection"))
        {
            infile <- file(vcf.fn[i], open="rt")
            error_fn <- sQuote(vcf.fn[i])
            if (grepl("\\.bcf$", vcf.fn[i], ignore.case=TRUE))
            {
                s <- readChar(infile, 9L)
                if (substr(s, 1L, 4L) != "BCF\002")
                    stop(vcf.fn[i], " should be BCF2 format.")
            } else
                is_vcf_fn <- TRUE
        } else {
            infile <- vcf.fn
            error_fn <- TRUE
        }

        # read header
        ans <- NULL
        while (length(s <- readLines(infile, n=1L)) > 0L)
        {
            n <- n + 1L
            if (substr(s, 1L, 6L) != "#CHROM")
            {
                s <- substring(s, 3L)
                ans <- c(ans, s)
            } else {
                samp.id <- scan(text=s, what=character(), sep="\t",
                    quiet=TRUE)[-seq_len(9L)]
                nSample <- length(samp.id)
                if (is_vcf_fn)
                {
                    s <- readLines(infile, n=1L)
                    if (length(s) > 0L)
                    {
                        ss <- scan(text=s, what=character(), sep="\t",
                            quiet=TRUE)
                        geno.text <- c(geno.text, ss[-seq_len(9L)])
                    }
                    if (isTRUE(getnum))
                    {
                        nVariant <- nVariant + length(s) +
                            .Call(SEQ_VCF_NumLines, infile, FALSE, 0, NULL,
                                verbose)$num
                    }
                }
                break
            }
            if ((n %% 10000L) == 0L)
            {
                warning(sprintf(
                    "There are too many lines in the header (>= %d). ", n),
                    "In order not to slow down the conversion, please ",
                    "consider deleting unnecessary annotations (like contig).",
                    immediate.=TRUE)
            }
        }

        if (!inherits(vcf.fn, "connection")) close(infile)
    }

    # add a message if fails
    if (!inherits(vcf.fn, "connection"))
    {
        error_fn <- paste(sQuote(vcf.fn), collapse=", ")
    } else {
        error_fn <- TRUE
    }

    ans <- unique(ans)


    #########################################################
    # parse texts

    #########################################################
    # get names and values

    ValString <- function(txt)
    {
        # string splited by "="
        x <- as.integer(regexpr("=", txt, fixed=TRUE))
        rv <- matrix("", nrow=length(txt), ncol=2L)
        for (i in seq_along(txt))
        {
            if (x[i] > 0L)
            {
                rv[i,1L] <- trimws(substring(txt[i], 1L, x[i]-1L))
                rv[i,2L] <- trimws(substring(txt[i], x[i]+1L))
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
                    a <- v[j,2L]

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
        if (!(number %in% c("A", "G", ".", "R")))
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

    if (length(geno.text))
    {
        txt <- unlist(vapply(geno.text, function(s)
        {
            scan(text=s, what="", sep=":", quiet=TRUE, nmax=1L)
        }, "", USE.NAMES=FALSE))
        if (any(grepl(",", txt, fixed=TRUE)))
        {
            ploidy <- NA_integer_
        } else {
            num <- lengths(strsplit(txt, "[|/]"))
            ploidy <- max(num)
            if (ploidy > 2L)
            {
                if (sum(num<=2L) >= sum(num>2L))
                    ploidy <- 2L
                else
                    ploidy <- ((ploidy+1L) %/% 2L) * 2L
            }
        }
    } else
        ploidy <- NA_integer_

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
    # INFO=<ID=ID,Number=number,Type=type,Description="description",
    #       Source="source",Version="version">

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
            {
                stop("Invalid data type (", v$Type, ")\nINFO=", s[i])
            }
            if (!CheckNum(v$Number))
            {
                stop("Invalid number (", v$Number, ")\nINFO=", s[i])
            }
            INFO <- rbind(INFO, v)
        } else
            stop("INFO=", s[i])
    }
    INFO <- unique(INFO)
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
    FILTER <- unique(FILTER)
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
            {
                stop("Invalid data type (", v$Type, ")\nFORMAT=", s[i])
            }
            if (!CheckNum(v$Number))
            {
                stop("Invalid number (", v$Number, ")\nFORMAT=", s[i])
            }
            FORMAT <- rbind(FORMAT, v)
        } else
            stop("FORMAT=", s[i])
    }
    FORMAT <- unique(FORMAT)
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
    ALT <- unique(ALT)
    ans <- ans[ans$id != "ALT", ]


    #########################################################
    # contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

    contig <- NULL
    s <- ans$value[ans$id == "contig"]
    lst <- lapply(s, function(x) {
        v <- AsDataFrame(x, NULL)
        if (is.null(v)) stop("contig=", x)
        v
    })
    nmlst <- unique(unlist(lapply(lst, function(x) colnames(x))))
    xx <- lapply(lst, function(x) {
        for (nm in setdiff(nmlst, colnames(x)))
            x[[nm]] <- NA_character_
        x[, match(nmlst, colnames(x)), drop=FALSE]
    })
    contig <- unique(Reduce(rbind, xx))
    for (i in seq_len(NCOL(contig)))
    {
        s <- type.convert(contig[, i], as.is=TRUE, numerals="no.loss")
        if (is.logical(s) && all(is.na(s))) s <- as.integer(s)
        contig[, i] <- s
    }
    ans <- ans[ans$id != "contig", ]


    #########################################################
    # reference=""
    reference <- NULL
    s <- ans$value[ans$id == "reference"]
    if (length(s) > 0L) reference <- s
    reference <- unique(reference)
    ans <- ans[ans$id != "reference", ]


    #########################################################
    # output

    error_fn <- FALSE
    rv <- list(fileformat=fileformat, info=INFO, filter=FILTER, format=FORMAT,
        alt=ALT, contig=contig, assembly=assembly, reference=reference,
        header=ans, ploidy=ploidy)
    rv$num.sample <- nSample
    if (getnum)
        rv$num.variant <- nVariant
    rv$sample.id <- samp.id
    class(rv) <- "SeqVCFHeaderClass"
    rv
}

print.SeqVCFHeaderClass <- function(x, ...) str(x)



#######################################################################
# get sample id from a VCF file
#

seqVCF_SampID <- function(vcf.fn)
{
    if (!inherits(vcf.fn, "connection"))
    {
        # check
        stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)

        # open the vcf file
        infile <- file(vcf.fn[1L], open="rt")
        on.exit(close(infile))
    } else {
        infile <- vcf.fn
    }

    # read header
    samp.id <- NULL
    while (length(s <- readLines(infile, n=1L)) > 0L)
    {
        if (substr(s, 1L, 6L) == "#CHROM")
        {
            samp.id <- scan(text=s, what=character(0), sep="\t",
                quiet=TRUE)[-seq_len(9L)]
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
    storage.option="LZMA_RA", info.import=NULL, fmt.import=NULL,
    genotype.var.name="GT", ploidy=NA_integer_, ignore.chr.prefix="chr",
    scenario=c("general", "imputation"), reference=NULL, start=1L, count=-1L,
    variant_count=NA_integer_, split=c("block", "variant"), optimize=TRUE,
    raise.error=TRUE, digest=TRUE, use_Rsamtools=NA, parallel=FALSE,
    verbose=TRUE)
{
    # check
    if (!inherits(vcf.fn, "connection"))
        stopifnot(is.character(vcf.fn), length(vcf.fn)>0L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    stopifnot(is.null(header) | inherits(header, "SeqVCFHeaderClass"))

    scenario <- match.arg(scenario)
    split_given <- !missing(split)   # whether split="block" is asked for
    split <- match.arg(split)
    storage.tmp <- storage.option
    if (is.character(storage.option))
    {
        storage.option <- seqStorageOption(storage.option)
        if (scenario == "imputation")
        {
            storage.option$mode <- c(
                `annotation/format/DS`="packedreal16:offset=0,scale=0.0001",
                `annotation/format/GP`="packedreal16:offset=0,scale=0.0001"
            )
        }
    } else {
        scenario <- ""
    }

    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))
    stopifnot(is.character(genotype.var.name), length(genotype.var.name)==1L)
    if (is.na(genotype.var.name)) genotype.var.name <- "GT"

    stopifnot(is.null(info.import) | is.character(info.import))
    stopifnot(is.null(fmt.import) | is.character(fmt.import))
    stopifnot(is.numeric(ploidy), length(ploidy)==1L)
    if (!is.na(ploidy)) stopifnot(ploidy >= 1L)
    stopifnot(is.character(ignore.chr.prefix), length(ignore.chr.prefix)>0L)

    stopifnot(is.null(reference) | is.character(reference))
    stopifnot(is.numeric(start), length(start)==1L)
    stopifnot(is.numeric(count), length(count)==1L)
    # the internal 'seek' attribute set by the parallel jobs in seqVCF2GDS(),
    # c(file index, file offset in bytes, the variant index at that offset),
    # see .vcf_seek_list()
    seek_info <- attr(start, "seek")
    # the internal 'range' attribute set by the parallel jobs when
    # split="block", c(file index, open_addr, first_addr, end_addr), see
    # BGZF_Split() in the C code
    range_info <- attr(start, "range")
    attributes(start) <- NULL

    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(raise.error), length(raise.error)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(use_Rsamtools), length(use_Rsamtools)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    show_timeheader <- !isTRUE(attr(verbose, "header_no_time"))

    parallel <- .McoreParallel(parallel)
    pnum <- .NumParallel(parallel)

    # whether to split the input at the BGZF block boundaries, so that
    # counting the variants beforehand is not needed at all
    use_block <- FALSE
    if ((split == "block") && (pnum > 1L))
    {
        msg <- NULL
        if (inherits(vcf.fn, "connection"))
            msg <- "the input is a connection object"
        else if (!all(vapply(vcf.fn, .bgzf_is, TRUE)))
            msg <- "the input is not in the BGZF format (i.e., bgzip)"
        else if ((start != 1L) || (count >= 0L))
            msg <- "'start' or 'count' is used"
        else if (!identical(variant_count, NA_integer_))
            msg <- "'variant_count' is used"
        if (is.null(msg))
            use_block <- TRUE
        else if (verbose && split_given)
            .cat("    split=\"block\" is not used, since ", msg, ".")
    }

    if (inherits(vcf.fn, "connection"))
    {
        if (pnum > 1L)
            stop("No parallel support when the input is a connection object.")
        if (length(variant_count)!=1L || !is.na(variant_count))
        {
            warning("'variant_count' is not used in seqVCF2GDS() ",
                "when 'vcf.fn' is a connection object.")
        }
    } else if (!identical(variant_count, NA_integer_))
    {
        if (!is.numeric(variant_count))
            stop("'variant_count' should be a numeric vector.")
        if (length(variant_count) != length(vcf.fn))
            stop("'variant_count' and 'vcf.fn' should have the same length.")
    }
    if (verbose && show_timeheader) .cat("##< ", .tm())

    genotype.storage <- "bit2"

    # check sample id
    if (!inherits(vcf.fn, "connection"))
    {
        vcf.fn <- normalizePath(vcf.fn, mustWork=FALSE)
        samp.id <- NULL
        for (i in seq_along(vcf.fn))
        {
            if (is.null(samp.id))
            {
                samp.id <- seqVCF_SampID(vcf.fn[i])
                # if (length(samp.id) <= 0L)
                #     message("No sample in '", vcf.fn[i], "'")
            } else {
                tmp <- seqVCF_SampID(vcf.fn[i])
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
    }


    ##################################################
    # parse the header of VCF

    if (!inherits(vcf.fn, "connection"))
    {
        if (is.null(header))
            header <- seqVCF_Header(vcf.fn, verbose=FALSE)
    } else {
        if (is.null(header))
        {
            header <- seqVCF_Header(vcf.fn, verbose=FALSE)
            samp.id <- header$sample.id
        } else
            samp.id <- seqVCF_SampID(vcf.fn)
    }

    if (!inherits(header, "SeqVCFHeaderClass"))
        stop("'header' should be NULL or returned from 'seqVCF_Header()'.")
    # 'seqVCF_Header'
    # returns list(fileformat=fileformat, info=INFO, filter=FILTER,
    #              format=FORMAT, alt=ALT, contig=contig, assembly=assembly,
    #              reference, header, ploidy)

    # genome reference information
    reference <- as.character(unique(c(reference, header$reference)))
    # ploidy
    if (!is.na(ploidy)) header$ploidy <- ploidy

    if (verbose)
    {
        .cat("Variant Call Format (VCF) Import:")
        if (!inherits(vcf.fn, "connection"))
        {
            sz <- file.size(vcf.fn)
            if (length(vcf.fn) == 1L)
                .cat("    file:")
            else
                .cat("    files: [", .pretty_size(sum(sz)), " totally]")
            for (i in seq_along(vcf.fn))
            {
                .cat("        ", basename(vcf.fn[i]), " (",
                    .pretty_size(sz[i]), ")")
            }
        } else {
            .cat("    [connection object]")
        }
        .cat("    file format: ", header$fileformat)
        .cat("    genome reference: ", if (length(reference))
                paste(reference, collapse=", ") else "<unknown>")
        .cat("    # of sets of chromosomes (ploidy): ", header$ploidy)
        .cat("    # of samples: ", .pretty(length(samp.id)))
        .cat("    genotype field: ", genotype.var.name)
        .cat("    genotype storage: ", genotype.storage)

        if (!is.character(storage.tmp))
            storage.tmp <- "customized"
        .cat("    compression method: ", storage.tmp)
        if (identical(scenario, "imputation"))
        {
            cat("    scenario: imputation\n")
            if ("DS" %in% header$format$ID)
                cat("        annotation/format/DS: packedreal16\n")
            if ("GP" %in% header$format$ID)
                cat("        annotation/format/GP: packedreal16\n")
        }
        flush.console()
    }

    # check header
    oldheader <- header
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
            message(sprintf(
                "No variable '%s' in the FORMAT field.",
                genotype.var.name))
        } else if (length(geno_idx) > 1L)
        {
            stop(sprintf(
            "There are more than one variable in the FORMAT field according to '%s', please fix the header.",
                genotype.var.name))
        } else {
            if (tolower(header$format$Type[geno_idx]) != "string")
            {
                stop(sprintf(
                "'%s' in the FORMAT field should be string-type according to genotypes, please fix the header.",
                    genotype.var.name))
            }
            if (header$format$Number[geno_idx] != "1")
            {
                stop(sprintf(
                "'%s' in the FORMAT field should be one-length string-type, please fix the header.",
                    genotype.var.name))
            }
        }
        if (length(geno_idx) > 0L)
        {
            geno_format <- header$format[geno_idx, ]
            header$format <- header$format[-geno_idx, ]
        } else {
            geno_format <- list(Description="Genotype")
        }
    } else {
        if (length(samp.id) > 0L)
        {
            message("\t",
                "variable id in the FORMAT field should be defined ahead, ",
                "and the undefined id is/are ignored during the conversion.")
        }
        geno_format <- list(Description="Genotype")
    }


    #######################################################################
    # format conversion in parallel

    if (pnum > 1L)
    {
        if (verbose) .cat("    # of cores/jobs: ", pnum)

        if (use_block)
        {
            # =========================================================
            # split="block": split the BGZF file(s) at the block boundaries,
            #     so that the variants do not need to be counted beforehand;
            #     'variant.id' is not known in each job, and it is written
            #     after merging all of the temporary files

            # the number of the parts of each file, proportional to its size
            sz <- file.size(vcf.fn)
            k <- pmax(1L, as.integer(round(pnum * sz / sum(sz))))
            prg <- do.call(rbind, lapply(seq_along(vcf.fn), function(i)
                cbind(i, .Call(SEQ_bgzip_split, vcf.fn[i], k[i]))))
            nparts <- NROW(prg)

            if (nparts >= 2L)
            {
                # need unique temporary file names
                ptmpfn <- .get_temp_fn(nparts,
                    sub("^([^.]*).*", "\\1", basename(out.fn)), dirname(out.fn))
                if (verbose)
                {
                    .cat("    >>> writing to ", nparts, " files: <<<")
                    cat(sprintf("        %s\t[%s: %s .. %s]\n",
                        basename(ptmpfn), basename(vcf.fn)[prg[, 1L]],
                        .pretty_size(prg[, 3L]), .pretty_size(prg[, 4L])),
                        sep="")
                    flush.console()
                }

                # show information
                update_info <- function(i)
                {
                    .cat("        |> ", i, " [", .tm(), " done]")
                    flush.console()
                    NULL
                }

                # reset memory before calling parallel
                gc(FALSE, reset=TRUE, full=TRUE)

                # conversion in parallel
                seqParApply(parallel, seq_len(nparts), FUN = function(i,
                    vcf.fn, hdr, storage.option, info.import, fmt.import,
                    genotype.var.name, ignore.chr.prefix, scenario, optim,
                    raise.err, ptmpfn, prg, upd)
                {
                    # the range of the BGZF blocks assigned to this job
                    pstart <- 1L
                    attr(pstart, "range") <- prg[i, ]
                    tryCatch(
                    {
                        SeqArray::seqVCF2GDS(vcf.fn, ptmpfn[i], header=hdr,
                            storage.option=storage.option,
                            info.import=info.import, fmt.import=fmt.import,
                            genotype.var.name=genotype.var.name,
                            ignore.chr.prefix=ignore.chr.prefix,
                            start=pstart, count=-1L,
                            optimize=optim, scenario=scenario,
                            raise.error=raise.err,
                            digest=FALSE, parallel=FALSE, verbose=FALSE)
                        if (is.function(upd)) upd(i)
                        i
                    }, error = function(e) {
                        # capture full traceback
                        trace <- capture.output({
                            cat("Error: ", e$message, "\n", sep="")
                            traceback()
                        })
                        con <- file(paste0(ptmpfn[i], ".progress.txt"),
                            open="at")
                        writeLines(trace, con)
                        close(con)
                        stop(e$message)
                    })
                }, vcf.fn=vcf.fn, hdr=oldheader, storage.option=storage.option,
                    info.import=info.import, fmt.import=fmt.import,
                    genotype.var.name=genotype.var.name,
                    ignore.chr.prefix=ignore.chr.prefix, scenario=scenario,
                    optim=optimize, raise.err=raise.error, ptmpfn=ptmpfn,
                    prg=prg, upd=if (isTRUE(verbose)) update_info else NULL)

                if (verbose)
                    .cat("    >>> Done (", .tm(), ") <<<")

            } else {
                pnum <- 1L; use_block <- FALSE
                message("No use of parallel environment!")
            }

        } else {
            # get the number of variants in each VCF file
            # 'voffset[[i]]' stores the file offsets, if available, so that the
            #     parallel jobs can seek to the starting position directly
            voffset <- vector("list", length(vcf.fn))
            for (i in seq_along(vcf.fn))
            {
                v <- variant_count[i]
                if (is.na(v) || (v < 0L))
                {
                    fn <- vcf.fn[i]
                    if (.vcf_seekable(fn))
                    {
                        z <- .vcf_count_offset(fn, parallel)
                        variant_count[i] <- z$num
                        if (isTRUE(z$line >= 1)) voffset[[i]] <- z
                    } else {
                        variant_count[i] <- seqVCF_Header(fn, getnum=TRUE,
                            parallel=parallel, verbose=FALSE)$num.variant
                        if (verbose && !.is_text_file(fn))
                        {
                            cat("    Hint: '", basename(fn), "' is not in the ",
                                "BGZF format, so each of the parallel jobs has ",
                                "to decompress the file from the beginning; ",
                                "using bgzip instead of gzip is much faster.\n",
                                sep="")
                        }
                    }
                }
            }
            num_var <- sum(variant_count)
            if (anyNA(num_var)) stop("Getting invalid # of variants.")

            if (start < 1L)
                stop("'start' should be a positive integer if conversion in parallel.")
            else if (start > num_var)
                stop("'start' should not be greater than the total number of variants.")
            if (count < 0L)
                count <- num_var - start + 1L
            if (start+count > num_var+1L)
                stop("Invalid 'count'.")
            if (verbose)
                .cat("    # of variants: ", .pretty(count))

            if (count >= pnum)
            {
                # need unique temporary file names
                ptmpfn <- .get_temp_fn(pnum,
                    sub("^([^.]*).*", "\\1", basename(out.fn)), dirname(out.fn))
                psplit <- .file_split(count, pnum, start)
                if (verbose)
                {
                    .cat("    >>> writing to ", pnum, " files: <<<")
                    cat(sprintf("        %s\t[%s .. %s]\n", basename(ptmpfn),
                            .pretty(psplit[[1L]]),
                            .pretty(psplit[[1L]] + psplit[[2L]] - 1L)), sep="")
                    flush.console()
                }
                # the file offsets of the starting variants, if available
                pseek <- .vcf_seek_list(psplit[[1L]], variant_count, voffset)
                voffset <- NULL  # not needed by the parallel jobs

                # show information
                update_info <- function(i)
                {
                    .cat("        |> ", i, " [", .tm(), " done]")
                    flush.console()
                    NULL
                }
                if (!isTRUE(verbose)) update_info <- "none"

                # reset memory before calling parallel
                gc(FALSE, reset=TRUE, full=TRUE)

                # conversion in parallel
                seqParallel(parallel, NULL, FUN = function(
                    vcf.fn, header, storage.option, info.import, fmt.import,
                    genotype.var.name, ignore.chr.prefix, scenario, optim,
                    raise.err, ptmpfn, psplit, pseek, variant_count)
                {
                    i <- process_index  # the process id, starting from one
                    # the starting variant index, with the internal 'seek'
                    # attribute to avoid scanning the file from the beginning
                    pstart <- psplit[[1L]][i]
                    attr(pstart, "seek") <- pseek[[i]]
                    tryCatch(
                    {
                        SeqArray::seqVCF2GDS(vcf.fn, ptmpfn[i], header=oldheader,
                            storage.option=storage.option, info.import=info.import,
                            fmt.import=fmt.import,
                            genotype.var.name=genotype.var.name,
                            ignore.chr.prefix=ignore.chr.prefix,
                            start = pstart, count = psplit[[2L]][i],
                            variant_count=variant_count,
                            optimize=optim, scenario=scenario,
                            raise.error=raise.err,
                            digest=FALSE, parallel=FALSE, verbose=FALSE)
                        i  # return the process index
                    }, error = function(e) {
                        # capture full traceback
                        trace <- capture.output({
                            cat("Error: ", e$message, "\n", sep="")
                            traceback()
                        })
                        con <- file(paste0(ptmpfn[i], ".progress.txt"), open="at")
                        writeLines(trace, con)
                        close(con)
                        stop(e$message)
                    })
                }, split = "none", .combine = update_info,
                    vcf.fn=vcf.fn, header=header, storage.option=storage.option,
                    info.import=info.import, fmt.import=fmt.import,
                    genotype.var.name=genotype.var.name,
                    ignore.chr.prefix=ignore.chr.prefix, scenario=scenario,
                    optim=optimize, raise.err=raise.error,
                    ptmpfn=ptmpfn, psplit=psplit, pseek=pseek,
                    variant_count=variant_count)

                if (verbose)
                    .cat("    >>> Done (", .tm(), ") <<<")

            } else {
                pnum <- 1L
                message("No use of parallel environment!")
            }
        }
    }


    #######################################################################
    # create a new GDS file

    gfile <- createfn.gds(out.fn)
    on.exit({ if (!is.null(gfile)) closefn.gds(gfile) }, add=TRUE)

    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(gfile, "description")
    put.attr.gdsn(n, "vcf.fileformat", header$fileformat)
    if (!is.null(header$assembly))
        put.attr.gdsn(n, "vcf.assembly", header$assembly)
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

    # add chromosome
    .AddVar(storage.option, gfile, "chromosome", storage="string")

    # add position
    # TODO: need to check whether position can be stored in 'int32'
    .AddVar(storage.option, gfile, "position", storage="int32")

    # add allele
    .AddVar(storage.option, gfile, "allele", storage="string")

    # add a folder for genotypes
    varGeno <- addfolder.gdsn(gfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", genotype.var.name[1L])
    put.attr.gdsn(varGeno, "Description", geno_format$Description[1L])

    # add data to the folder of genotype
    if (is.na(header$ploidy)) header$ploidy <- 2L
    if (header$ploidy > 0L)
    {
        if (nSamp > 0L)
        {
            geno.node <- .AddVar(storage.option, varGeno, "data",
                storage=genotype.storage, valdim=c(header$ploidy, nSamp, 0L))
        } else
            geno.node <- NULL
    } else
        stop("Invalid ploidy.")
    node <- .AddVar(storage.option, varGeno, "@data", storage="uint8",
        visible=FALSE)

    node <- .AddVar(storage.option, varGeno, "extra.index", storage="int32",
        valdim=c(3L,0L))
    put.attr.gdsn(node, "R.colnames",
        c("sample.index", "variant.index", "length"))
    .AddVar(storage.option, varGeno, "extra", storage="int16")


    # add phase folder
    varPhase <- addfolder.gdsn(gfile, "phase")
    if (header$ploidy > 1L && nSamp > 0L)
    {
        # add data
        if (header$ploidy > 2L)
            dm <- c(header$ploidy-1L, nSamp, 0L)
        else if (header$ploidy > 1L)
            dm <- c(nSamp, 0L)
        else
            dm <- NULL

        if (!is.null(dm))
        {
            .AddVar(storage.option, varPhase, "data", storage="bit1", valdim=dm)
            node <- .AddVar(storage.option, varPhase, "extra.index",
                storage="int32", valdim=c(3L,0L))
            put.attr.gdsn(node, "R.colnames",
                c("sample.index", "variant.index", "length"))
            .AddVar(storage.option, varPhase, "extra", storage="bit1")
        }
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
    if (verbose) { prefix <- " "; cat("    INFO:") }

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

        if (NROW(header$info) > 0L)
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
                    if (verbose)
                    {
                        cat(prefix, header$info$ID[i], sep="")
                        prefix <- ","
                    }
                    node <- .AddVar(storage.option, varInfo, header$info$ID[i],
                        storage=mode, valdim=initdim)
                    put.attr.gdsn(node, "Number", header$info$Number[i])
                    put.attr.gdsn(node, "Type", header$info$Type[i])
                    put.attr.gdsn(node, "Description",
                        header$info$Description[i])
                    if (!is.na(header$info$Source[i]))
                        put.attr.gdsn(node, "Source", header$info$Source[i])
                    if (!is.na(header$info$Version[i]))
                        put.attr.gdsn(node, "Version", header$info$Version[i])

                    if (s %in% c(".", "A", "G", "R"))
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
    if (verbose) cat("\n")


    ##################################################
    # VCF Format

    # add the FORMAT field
    varFormat <- addfolder.gdsn(varAnnot, "format")
    if (verbose) { prefix <- " "; cat("    FORMAT:") }

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
            integer = { mode <- "vl_int"; int_type[i] <- 1L },
            float = { mode <- "float"; int_type[i] <- 2L },
            character = { mode <- "string"; int_type[i] <- 4L },
            string = { mode <- "string"; int_type[i] <- 4L },
            stop(sprintf("Unknown FORMAT Type: %s", header$format$Type[i]))
        )

        # FORMAT Number
        s <- header$format$Number[i]
        if (grepl("^[[:digit:]]", s))
        {
            if (as.integer(s) <= 0)
            {
                print(header[, i])
                stop("The length should be >0.")
            }
        } else {
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
            if (verbose)
            {
                cat(prefix, header$format$ID[i], sep="")
                prefix <- ","
            }
            node <- addfolder.gdsn(varFormat, header$format$ID[i])
            put.attr.gdsn(node, "Number", header$format$Number[i])
            put.attr.gdsn(node, "Type", header$format$Type[i])
            put.attr.gdsn(node, "Description", header$format$Description[i])

            if (nSamp > 0L)
            {
                .AddVar(storage.option, node, "data", storage=mode,
                    valdim=c(nSamp, 0L))
            }
            .AddVar(storage.option, node, "@data", storage="int32",
                visible=FALSE)
        }
    }

    if (!is.null(header$format))
    {
        header$format$int_type <- as.integer(int_type)
        header$format$int_num  <- as.integer(int_num)
        header$format$import.flag <- import.flag
    }
    if (verbose) cat("\n")


    ##################################################
    # add annotation folder

    addfolder.gdsn(gfile, "sample.annotation")


    ##################################################
    # sync file
    sync.gds(gfile)
    if (verbose)
        cat("Output:\n    ", out.fn, "\n", sep="")


    ##################################################
    # for-loop each file

    if (pnum <= 1L)
    {
        filterlevels <- header$filter$ID
        linecnt <- double(1L)

        # progress file
        prog_fn <- paste0(out.fn, ".progress.txt")
        progfile <- file(prog_fn, "wt")
        cat(">>> ", out.fn, " <<<\n", file=progfile, sep="")
        if (verbose)
            .cat("    [Progress Info: ", basename(prog_fn), "]")
        infile <- NULL
        on.exit({
            close(progfile)
            unlink(prog_fn, force=TRUE)
            if (inherits(infile, "connection")) close(infile)
        }, add=TRUE)

        if (!inherits(vcf.fn, "connection"))
        {

            for (i in seq_along(vcf.fn))
            {
                # split="block": only the assigned file is converted
                if (!is.null(range_info) && (range_info[1L] != i)) next
                if (!is.null(variant_count))
                {
                    cnt <- cumsum(variant_count)[i]
                    if (is.finite(cnt) & (start > cnt))
                    {
                        linecnt <- as.double(cnt)
                        next
                    }
                }

                if (verbose)
                {
                    cat(sprintf("Parsing '%s':\n", basename(vcf.fn[i])))
                    flush.console()
                }

                # seek to the starting position directly, if the file offset
                # is known, instead of scanning the file from the beginning
                skiphead <- TRUE; foffset <- brange <- NULL
                if (!is.null(range_info))
                {
                    # split="block": the range of the BGZF blocks, which is
                    # read in the C code; only the 1st part has the header
                    brange <- range_info[2:4]
                    skiphead <- (range_info[3L] <= 0)
                    infile <- vcf.fn[i]
                    linecnt <- 0
                } else {
                    if (!is.null(seek_info) && (seek_info[1L] == i))
                    {
                        linecnt <- as.double(seek_info[4L] - 1)
                        skiphead <- FALSE  # no VCF header at the file offset
                    }
                    if (skiphead || !.bgzf_is(vcf.fn[i]))
                    {
                        infile <- file(vcf.fn[i], open="rt")
                        if (!skiphead) seek(infile, where=seek_info[2L])
                    } else {
                        # let the C code read the BGZF file from the given
                        # virtual offset, since seek() does not work here
                        infile <- vcf.fn[i]
                        foffset <- seek_info[2:3]
                    }
                }

                # call C function
                v <- .Call(SEQ_VCF_Parse, vcf.fn[i], header, gfile$root,
                    list(sample.num = length(samp.id),
                        genotype.var.name = genotype.var.name,
                        infile = infile,
                        raise.error = raise.error, filter.levels = filterlevels,
                        start = start, count = count,
                        chr.prefix = ignore.chr.prefix,
                        progfile = progfile, use.file = skiphead,
                        line.base = if (skiphead || is.null(seek_info))
                            NULL else seek_info[5L],
                        file.offset = foffset, block.range = brange,
                        verbose = verbose),
                    linecnt, new.env())

                filterlevels <- unique(c(filterlevels, v))
                if (verbose && !is.null(geno.node))
                    print(geno.node)

                if (inherits(infile, "connection")) close(infile)
                infile <- NULL
            }

        } else {
            if (verbose)
            {
                cat("Parsing 'connection object':\n")
                flush.console()
            }

            # call C function
            v <- .Call(SEQ_VCF_Parse, "connection object", header, gfile$root,
                list(sample.num = length(samp.id),
                    genotype.var.name = genotype.var.name,
                    infile = vcf.fn,
                    raise.error = raise.error, filter.levels = filterlevels,
                    start = start, count = count,
                    chr.prefix = ignore.chr.prefix,
                    progfile = progfile, use.file = FALSE,
                    verbose = verbose),
                linecnt, new.env())

            filterlevels <- unique(c(filterlevels, v))
            if (verbose) print(geno.node)
        }

    } else {
        ## merge all temporary files

        # all GDS variables to be merged
        # with split="block", each job does not know the global variant
        #     index, so 'variant.id' is written after merging
        varnm <- c(if (!use_block) "variant.id", "position", "chromosome",
            "allele",
            "genotype/data", "genotype/@data",
            "genotype/extra", "genotype/extra.index",
            "phase/data", "phase/extra", "phase/extra.index",
            "annotation/id", "annotation/qual")

        if (is.null(index.gdsn(gfile, "phase/data", silent=TRUE)))
        {
            varnm <- setdiff(varnm,
                c("phase/data", "phase/extra", "phase/extra.index"))
        }

        s <- ls.gdsn(index.gdsn(gfile, "annotation/info"), include.hidden=TRUE)
        if (length(s) > 0L)
            varnm <- c(varnm, paste0("annotation/info/", s))

        s <- ls.gdsn(index.gdsn(gfile, "annotation/format"))
        if (length(s) > 0L)
        {
            varnm <- c(varnm, paste0("annotation/format/", rep(s, each=2L),
                    c("/data", "/@data")))
        }

        if (verbose) cat("Merging:\n")
        filtervar <- character()

        # open all temporary files
        for (fn in ptmpfn)
        {
            if (verbose)
                cat("    ", basename(fn), sep="")
            # open the gds file
            tmpgds <- seqOpen(fn, allow.duplicate=TRUE)
            # merge variables
            for (nm in varnm)
            {
                n <- index.gdsn(tmpgds, nm, silent=TRUE)
                if (!is.null(n))
                    append.gdsn(index.gdsn(gfile, nm), n)
            }
            # merge filter variable (a factor variable)
            filtervar <- c(filtervar, as.character(
                read.gdsn(index.gdsn(tmpgds, "annotation/filter"))))
            # close the file
            seqClose(tmpgds)
            if (verbose)
                .cat(" [done, ", .tm(), "]")
        }

        if (use_block)
        {
            # 'variant.id' is 1 .. the total number of variants
            n <- prod(objdesp.gdsn(index.gdsn(gfile, "position"))$dim)
            append.gdsn(index.gdsn(gfile, "variant.id"), seq_len(n))
        }

        filtervar <- as.factor(filtervar)
        append.gdsn(varFilter, filtervar)
        readmode.gdsn(varFilter)
        s <- levels(filtervar)
        filterlevels <- c(s, setdiff(header$filter$ID, s))

        # remove temporary files
        unlink(ptmpfn, force=TRUE)
    }

    if (length(filterlevels) > 0L)
    {
        put.attr.gdsn(varFilter, "R.class", "factor")
        put.attr.gdsn(varFilter, "R.levels", filterlevels)

        if (NROW(header$filter) > 0L)
        {
            dp <- header$filter$Description[
                match(filterlevels, header$filter$ID)]
            dp[is.na(dp)] <- ""
        } else
            dp <- rep("", length(filterlevels))

        put.attr.gdsn(varFilter, "Description", dp)
    }

    # RLE-coded chromosome
    .optim_chrom(gfile)

    # if there is no genotype
    n <- index.gdsn(gfile, "genotype/data", silent=TRUE)
    if (is.null(n) || prod(objdesp.gdsn(n)$dim) <= 0)
    {
        for (nm in c("genotype/data", "genotype/@data", "genotype/extra.index",
            "genotype/extra", "phase/data", "phase/extra.index", "phase/extra"))
        {
            n <- index.gdsn(gfile, nm, silent=TRUE)
            if (!is.null(n)) delete.gdsn(n)
        }
    }

    # create hash
    .DigestFile(gfile, digest, verbose)

    # close the GDS file
    closefn.gds(gfile)
    gfile <- NULL


    ##################################################
    # optimize access efficiency

    if (verbose)
        if (optimize) .cat("Done.  # ", .tm()) else cat("Done.\n")
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.fn, verbose=verbose)
    }
    if (verbose && show_timeheader) .cat("##> ", .tm())

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a BCF file to a GDS file
#

seqBCF2GDS <- function(bcf.fn, out.fn, header=NULL, storage.option="LZMA_RA",
    info.import=NULL, fmt.import=NULL, genotype.var.name="GT",
    ploidy=NA_integer_, ignore.chr.prefix="chr",
    scenario=c("general", "imputation"), reference=NULL, optimize=TRUE,
    raise.error=TRUE, digest=TRUE, bcftools="bcftools", verbose=TRUE)
{
    # check
    stopifnot(is.character(bcf.fn), length(bcf.fn)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    stopifnot(is.character(bcftools), length(bcftools)==1L)

    # command-line
    cmd <- paste(shQuote(bcftools), "view", shQuote(bcf.fn))
    if (verbose)
        cat("Running:\n    ", cmd, "\n", sep="")
    {
        f <- pipe(paste(shQuote(bcftools), "-v"))
        s <- readLines(f)
        close(f)
        if (length(s) <= 0L) stop("Please install bcftools.")
    }
    pipefile <- pipe(cmd, "rt")
    on.exit(close(pipefile))

    # run VCF to GDS
    seqVCF2GDS(pipefile, out.fn, header=header,
        storage.option=storage.option,
        info.import=info.import, fmt.import=fmt.import,
        genotype.var.name=genotype.var.name, ploidy=ploidy,
        ignore.chr.prefix=ignore.chr.prefix, scenario=scenario,
        reference=reference, optimize=optimize, raise.error=raise.error,
        digest=digest, verbose=verbose)

    # output
    invisible(normalizePath(out.fn))
}
