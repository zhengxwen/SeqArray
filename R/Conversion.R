#######################################################################
#
# Package Name: SeqArray
#
# Description:
#     Data Management of Large-scale Whole-Genome Sequence Variant Calls
#



#######################################################################
# Convert a SeqArray GDS file to a VCF file
#

seqGDS2VCF <- function(gdsfile, vcf.fn, info.var=NULL, fmt.var=NULL,
    use_Rsamtools=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    if (is.character(gdsfile))
        stopifnot(length(gdsfile)==1L)
    if (!inherits(vcf.fn, "connection"))
        stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)
    stopifnot(is.null(info.var) | is.character(info.var))
    stopifnot(is.null(fmt.var) | is.character(fmt.var))

    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))
    }

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

    bgzf <- FALSE
    if (!inherits(vcf.fn, "connection"))
    {
        ext <- substring(vcf.fn, nchar(vcf.fn)-2L)
        if (ext == ".gz")
        {
            if (.Platform$OS.type == "windows")
                use_Rsamtools <- FALSE
            if (isTRUE(use_Rsamtools) && requireNamespace("Rsamtools"))
            {
                ofile <- .Call(SEQ_bgzip_create, vcf.fn)
                bgzf <- TRUE
            } else {
                if (verbose)
                {
                    if (.Platform$OS.type != "windows")
                        message("Hint: install Rsamtools to enable the bgzf output.")
                }
                ofile <- gzfile(vcf.fn, "wb")
            }
        } else if (ext == ".bz")
        {
            ofile <- bzfile(vcf.fn, "wb")
        } else if (ext == ".xz")
        {
            ofile <- xzfile(vcf.fn, "wb")
        } else {
            ofile <- file(vcf.fn, open="wb")
        }
        on.exit(close(ofile), add=TRUE)
    } else {
        ofile <- vcf.fn
    }

    op <- options("useFancyQuotes")
    options(useFancyQuotes = FALSE)
    on.exit(options(op), add=TRUE)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("VCF Export: ", basename(vcf.fn), "\n", sep="")
        s <- .seldim(gdsfile)
        cat("    ", .pretty(s[2L]), " sample", .plural(s[2L]), ", ",
            .pretty(s[3L]), " variant", .plural(s[3L]), "\n", sep="")
        s <- paste(z$info$ID, collapse=", ")
        cat("    INFO Field: ", ifelse(s!="", s, "<none>"), "\n", sep="")
        s <- paste(z$format$ID, collapse=", ")
        cat("    FORMAT Field: ", ifelse(s!="", s, "<none>"), "\n", sep="")
        cat(ifelse(bgzf,
            "    output to BGZF format\n", "    output to a general gzip file\n"))
    }


    ######################################################
    # write the header

    txt <- character()
    a <- get.attr.gdsn(index.gdsn(gdsfile, "description"))

    # fileformat
    if (is.null(a$vcf.fileformat))
        a$vcf.fileformat <- "VCFv4.2"
    txt <- c(txt, paste0("##fileformat=", a$vcf.fileformat))

    # fileDate
    txt <- c(txt, paste0("##fileDate=", format(Sys.time(), "%Y%m%d")))

    # program, source
    aa <- get.attr.gdsn(gdsfile$root)
    if (is.null(aa$FileVersion))
        aa$FileVersion <- "v1.0"
    txt <- c(txt, paste0("##source=SeqArray_Format_", aa$FileVersion))

    # reference
    if (length(z$reference) > 0L)
        txt <- c(txt, paste0("##reference=", z$reference[1L]))

    # assembly
    if (!is.null(a$vcf.assembly))
        txt <- c(txt, paste0("##assembly=", dq(a$vcf.assembly)))

    # ALT=<ID=type,Description=description>
    for (i in seq_len(nrow(z$allele)))
    {
        s <- sprintf("##ALT=<ID=%s,Description=%s>",
            as.character(z$allele[i,1L]), dq(z$allele[i,2L]))
        txt <- c(txt, s)
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
            txt <- c(txt, paste0("##contig=<", s, ">"))
        }
    }

    # INFO field
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
        txt <- c(txt, paste0("##INFO=<", s, ">"))
    }

    # FILTER field
    n <- index.gdsn(gdsfile, "annotation/filter", silent=TRUE)
    if (!is.null(n))
    {
        at <- get.attr.gdsn(n)
        id <- at$R.levels; dp <- at$Description
        for (i in seq_along(id))
        {
            txt <- c(txt, sprintf("##FILTER=<ID=%s,Description=%s>",
                dq(id[i]), dq(dp[i], TRUE)))
        }
    }

    # FORMAT field
    a <- get.attr.gdsn(index.gdsn(gdsfile, "genotype"))
    txt <- c(txt, sprintf(
        "##FORMAT=<ID=%s,Number=1,Type=String,Description=%s>",
        a$VariableName, dq(a$Description, TRUE)))
    for (nm in z$format$ID)
    {
        a <- get.attr.gdsn(index.gdsn(gdsfile,
            paste("annotation/format/", nm, sep="")))
        txt <- c(txt, sprintf(
            "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=%s>",
            nm, dq(a$Number), dq(a$Type), dq(a$Description, TRUE)))
    }

    # others
    n <- index.gdsn(gdsfile, "description/vcf.header", silent=TRUE)
    if (!is.null(n))
    {
        dat <- read.gdsn(n)
        for (i in seq_len(nrow(dat)))
        {
            s <- dat[i,1L]
            if (!(s %in% c("fileDate", "source")))
                txt <- c(txt, paste0("##", s, "=", dq(dat[i,2L])))
        }
    }


    ######################################################
    # write the header -- samples

    txt <- c(txt, paste(
        c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", seqGetData(gdsfile, "sample.id")), collapse="\t"))
    writeLines(txt, ofile)


    ######################################################
    # write the contents

    # the INFO field
    nm.info <- character()
    if (!is.null(z$info$ID))
    {
        s <- z$info$ID
        if (length(s) > 0L)
            nm.info <- paste("annotation/info/", s, sep="")
    }

    # the FORMAT field
    nm.format <- character()
    if (!is.null(z$format$ID))
    {
        s <- z$format$ID
        if (length(s) > 0L)
            nm.format <- paste("annotation/format/", s, sep="")
    }

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

    # initialize
    dm <- .seldim(gdsfile)
    .Call(SEQ_ToVCF_Init, dm, len.info, len.fmt, ofile, verbose)
    on.exit({ .Call(SEQ_ToVCF_Done) }, add=TRUE)

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

    # C function name
    if (is.null(index.gdsn(gdsfile, "genotype/data", silent=TRUE)))
    {
        nm <- nm[!(nm %in% c("genotype", "phase"))]
        cfn <- "SEQ_ToVCF_NoGeno"
    } else if (is.null(index.gdsn(gdsfile, "phase/data", silent=TRUE)))
    {
        nm <- nm[nm != "phase"]
        cfn <- "SEQ_ToVCF_Haploid"
    } else {
        if (length(nm.format)>0L | dm[1L]!=2L)
            cfn <- "SEQ_ToVCF"
        else
            cfn <- "SEQ_ToVCF_Di_WrtFmt"
    }

    # output lines by variant
    seqApply(gdsfile, nm, margin="by.variant", as.is="none",
        FUN=.cfunction(cfn), .useraw=NA, .progress=verbose)

    # finalize
    .Call(SEQ_ToVCF_Done)
    on.exit({
        if (verbose)
            cat(date(), "    Done.\n", sep="")
    }, add=TRUE)

    # output
    if (!inherits(vcf.fn, "connection"))
        invisible(normalizePath(vcf.fn))
    else
        invisible()
}



#######################################################################
# Convert a SeqArray GDS file to a SNP GDS file
#

seqGDS2SNP <- function(gdsfile, out.gdsfn, dosage=FALSE,
    compress.geno="LZMA_RA", compress.annotation="LZMA_RA",
    ds.type=c("packedreal16", "float", "double"), optimize=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.logical(dosage) | is.character(dosage), length(dosage)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    ds.type <- match.arg(ds.type)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        if (isTRUE(dosage) | is.character(dosage))
            cat("SeqArray GDS to SNP GDS dosage format:\n")
        else
            cat("SeqArray GDS to SNP GDS format:\n")
    }

    # if it is a file name
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile)
        on.exit({ seqClose(gdsfile) })
    }

    # dosage variable name
    if (isTRUE(dosage) | is.character(dosage))
    {
        if (isTRUE(dosage))
            dosage <- "annotation/format/DS"
        else if (substr(dosage, 1L, 18L) != "annotation/format/")
            dosage <- paste0("annotation/format/", dosage)
    }

    if (verbose)
    {
        dm <- .seldim(gdsfile)
        cat("    # of samples: ", .pretty(dm[2L]), "\n", sep="")
        cat("    # of variants: ", .pretty(dm[3L]), "\n", sep="")
        if (is.character(dosage))
        {
            cat("    dosage compression: ", compress.geno, ", from [",
                dosage, "]\n", sep="")
        } else {
            cat("    genotype compression: ", compress.geno, "\n", sep="")
        }
        cat("    annotation compression: ", compress.annotation, "\n", sep="")
    }

    # create GDS file
    gfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({
        if (!is.null(gfile)) closefn.gds(gfile)
    }, add=TRUE)

    # add a flag
    if (!is.character(dosage))
    {
        put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    } else {
        put.attr.gdsn(gfile$root, "FileFormat", "IMPUTED_DOSAGE")
    }

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
    if (!is.character(dosage))
    {
        gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
            valdim=c(length(sampid), 0L), compress=compress.geno)
        put.attr.gdsn(gGeno, "sample.order")
        seqApply(gdsfile, "$dosage", as.is=gGeno, .useraw=TRUE, .progress=verbose,
            FUN = .cfunction("FC_GDS2SNP"))
    } else {
        .cfunction("FC_SetNumSamp")(length(sampid))
        gGeno <- add.gdsn(gfile, "genotype", storage=ds.type,
            valdim=c(length(sampid), 0L), compress=compress.geno)
        put.attr.gdsn(gGeno, "sample.order")
        seqApply(gdsfile, dosage, as.is=gGeno, .progress=verbose,
            FUN = .cfunction("FC_GDS2Dosage"))
    }

    readmode.gdsn(gGeno)
    closefn.gds(gfile)
    gfile <- NULL
    if (verbose)
        cat("Done.\n", date(), "\n", sep="")

    if (optimize)
    {
        if (verbose) cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.gdsfn))
}



#######################################################################
# Convert a SNP GDS file to a SeqArray GDS file
#

seqSNP2GDS <- function(gds.fn, out.fn, storage.option="LZMA_RA", major.ref=TRUE,
    ds.type=c("packedreal16", "float", "double"), optimize=TRUE, digest=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(gds.fn), length(gds.fn)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    ds.type <- match.arg(ds.type)
    stopifnot(is.logical(major.ref), length(major.ref)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.character(storage.option))
        storage.option <- seqStorageOption(storage.option)
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("SNP GDS to SeqArray GDS Format:\n")
    }

    # open the SNP GDS
    srcfile <- openfn.gds(gds.fn)
    on.exit({ closefn.gds(srcfile) })

    nSamp <- prod(objdesp.gdsn(index.gdsn(srcfile, "sample.id"))$dim)
    nSNP <- prod(objdesp.gdsn(index.gdsn(srcfile, "snp.id"))$dim)

    # check genotype type
    n <- index.gdsn(srcfile, "genotype")
    dm <- objdesp.gdsn(n)$dim
    if (length(dm) != 2L)
        stop("'genotype' of SNP GDS should be a matrix.")
    geno_type <- as.character(objdesp.gdsn(n)$type)
    if (!(geno_type %in% c("Integer", "Real")))
        stop("'genotype' should be integers or real numbers.")

    # determine how genotypes are stored
    snpfirstdim <- TRUE
    rd <- names(get.attr.gdsn(n))
    if ("snp.order" %in% rd) snpfirstdim <- TRUE
    if ("sample.order" %in% rd) snpfirstdim <- FALSE
    if (snpfirstdim)
    {
        if ((dm[1L]!=nSNP) || (dm[2L]!=nSamp))
            stop("Invalid dimension of 'genotype'.")
    } else {
        if ((dm[1L]!=nSamp) || (dm[2L]!=nSNP))
            stop("Invalid dimension of 'genotype'.")
    }


    # create GDS file
    dstfile <- createfn.gds(out.fn)
    # close the file at the end
    on.exit({ closefn.gds(dstfile) }, add=TRUE)

    put.attr.gdsn(dstfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(dstfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(dstfile, "description")
    put.attr.gdsn(n, "source.format", "SNPRelate GDS Format")

    # add sample.id
    if (verbose) cat("    sample.id")
    s <- read.gdsn(index.gdsn(srcfile, "sample.id"))
    n <- .AddVar(storage.option, dstfile, "sample.id", s, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add variant.id
    if (verbose) cat("    variant.id")
    s <- read.gdsn(index.gdsn(srcfile, "snp.id"))
    n <- .AddVar(storage.option, dstfile, "variant.id", s, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add position
    if (verbose) cat("    position")
    s <- read.gdsn(index.gdsn(srcfile, "snp.position"))
    n <- .AddVar(storage.option, dstfile, "position", s, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add chromosome
    if (verbose) cat("    chromosome")
    s <- read.gdsn(index.gdsn(srcfile, "snp.chromosome"))
    s <- as.character(s)
    n <- .AddVar(storage.option, dstfile, "chromosome", s, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add allele
    nd_allele <- .AddVar(storage.option, dstfile, "allele", storage="character")
    # add a folder for genotypes
    if (verbose) cat("    genotype")
    varGeno <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", "GT")
    put.attr.gdsn(varGeno, "Description", "Genotype")

    # major reference
    .cfunction("FC_SNP2GDS_Ref")(major.ref)

    nd_geno <- .AddVar(storage.option, varGeno, "data", storage="bit2",
        valdim=c(2L, nSamp, 0L))
    nd_geno_idx <- .AddVar(storage.option, varGeno, "@data",
        storage="uint8", visible=FALSE)
    n <- .AddVar(storage.option, varGeno, "extra.index", storage="int32",
        valdim=c(3L,0L), closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    .AddVar(storage.option, varGeno, "extra", storage="int16", closezip=TRUE)

    if (geno_type == "Integer")
    {
        # add genotypes to genotype/data
        apply.gdsn(list(index.gdsn(srcfile, "genotype"),
                index.gdsn(srcfile, "snp.allele")),
            c(ifelse(snpfirstdim, 1L, 2L), 1L),
            FUN=.cfunction("FC_SNP2GDS"), as.is="gdsnode",
            target.node=list(nd_geno, nd_allele))
        readmode.gdsn(nd_geno)
        .DigestCode(nd_geno, digest, verbose)

        if (verbose) cat("    allele")
        readmode.gdsn(nd_allele)
        .DigestCode(nd_allele, digest, verbose)

        .repeat_gds(nd_geno_idx, 1L, nSNP)
        readmode.gdsn(nd_geno_idx)
        .DigestCode(nd_geno_idx, digest, FALSE)

        sync.gds(dstfile)
    } else {
        # add dosages to annotation/format/DS
        readmode.gdsn(nd_geno)
        readmode.gdsn(nd_geno_idx)
    }

    # add a folder for phase information
    if (verbose) cat("    phase")
    varPhase <- addfolder.gdsn(dstfile, "phase")

    n <- .AddVar(storage.option, varPhase, "data", storage="bit1",
        valdim=c(nSamp, 0L))
    if (geno_type == "Integer")
        .repeat_gds(n, 0L, as.double(nSNP)*nSamp)
    readmode.gdsn(n)
    .DigestCode(n, digest, verbose)

    n <- .AddVar(storage.option, varPhase, "extra.index", storage="int32",
        valdim=c(3L,0L), closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    .AddVar(storage.option, varPhase, "extra", storage="bit1", closezip=TRUE)

    sync.gds(dstfile)


    # add annotation folder
    varAnnot <- addfolder.gdsn(dstfile, "annotation")

    # add annotation/id
    if (verbose) cat("    annotation/id")
    n <- .AddVar(storage.option, varAnnot, "id", storage="string")
    if (is.null(index.gdsn(srcfile, "snp.rs.id", silent=TRUE)))
        assign.gdsn(n, index.gdsn(srcfile, "snp.id"))
    else
        assign.gdsn(n, index.gdsn(srcfile, "snp.rs.id"))
    readmode.gdsn(n)
    .DigestCode(n, digest, verbose)

    # add annotation/qual
    n <- .AddVar(storage.option, varAnnot, "qual", storage="float")
    .repeat_gds(n, 100.0, nSNP)
    readmode.gdsn(n)
    .DigestCode(n, digest, FALSE)

    # add filter
    n <- .AddVar(storage.option, varAnnot, "filter", storage="int32")
    .repeat_gds(n, 1L, nSNP)
    readmode.gdsn(n)
    put.attr.gdsn(n, "R.class", "factor")
    put.attr.gdsn(n, "R.levels", c("PASS"))
    put.attr.gdsn(n, "Description", c("All filters passed"))
    .DigestCode(n, digest, FALSE)

    # add the INFO field
    n1 <- addfolder.gdsn(varAnnot, "info")
    n2 <- index.gdsn(srcfile, "snp.annot", silent=TRUE)
    if (!is.null(n2))
    {
        for (i in ls.gdsn(n2))
        {
            if (verbose) cat("    annotation/info/", i, sep="")
            copyto.gdsn(n1, index.gdsn(n2, i))
           .DigestCode(index.gdsn(n1, i), digest, verbose)
        }
    }

    # add the FORMAT field
    n <- addfolder.gdsn(varAnnot, "format")

    # add dosages to annotation/format/DS
    if (geno_type == "Real")
    {
        if (verbose) cat("    annotation/format/DS")
        varGeno <- addfolder.gdsn(n, "DS")
        put.attr.gdsn(varGeno, "Number", "1")
        put.attr.gdsn(varGeno, "Type", "Float")
        put.attr.gdsn(varGeno, "Description", "Estimated alternate allele dosage")

        nd_geno <- .AddVar(storage.option, varGeno, "data", storage=ds.type,
            valdim=c(nSamp, 0L))
        # add genotypes to genotype/data
        apply.gdsn(list(index.gdsn(srcfile, "genotype"),
                index.gdsn(srcfile, "snp.allele")),
            c(ifelse(snpfirstdim, 1L, 2L), 1L),
            FUN=.cfunction("FC_Dosage2GDS"), as.is="gdsnode",
            target.node=list(nd_geno, nd_allele))

        readmode.gdsn(nd_geno)
        .DigestCode(nd_geno, digest, verbose)

        if (verbose) cat("    allele")
        readmode.gdsn(nd_allele)
        .DigestCode(nd_allele, digest, verbose)

        nd_geno_idx <- .AddVar(storage.option, varGeno, "@data",
            storage="uint8", visible=FALSE)
        .repeat_gds(nd_geno_idx, 1L, nSNP)
        readmode.gdsn(nd_geno_idx)
        .DigestCode(nd_geno_idx, digest, FALSE)

        sync.gds(dstfile)
    }

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
           .DigestCode(index.gdsn(n, i), digest, verbose)
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
        cleanup.gds(out.fn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a PLINK BED file to a SeqArray GDS file
#

seqBED2GDS <- function(bed.fn, fam.fn, bim.fn, out.gdsfn, compress.geno="LZMA_RA",
    compress.annotation="LZMA_RA", optimize=TRUE, digest=TRUE, parallel=FALSE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(bed.fn), length(bed.fn)==1L)
    if (missing(fam.fn) && missing(bim.fn))
    {
        bed.fn <- gsub("\\.bed$", "", bed.fn, ignore.case=TRUE)
        fam.fn <- paste0(bed.fn, ".fam")
        bim.fn <- paste0(bed.fn, ".bim")
        bed.fn <- paste0(bed.fn, ".bed")
    }
    stopifnot(is.character(fam.fn), length(fam.fn)==1L)
    stopifnot(is.character(bim.fn), length(bim.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest) | is.character(digest), length(digest)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    pnum <- .NumParallel(parallel)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("PLINK BED to SeqArray GDS:\n")
    }

    ##  open and detect bed.fn  ##

    bedfile <- .open_bin(bed.fn)
    on.exit({ .close_conn(bedfile) })
    b <- as.integer(readBin(bedfile$con, raw(), 3L))
    if ((b[1L] != 0x6C) || (b[2L] != 0x1B))
        stop(sprintf("Invalid magic number (0x%02X,0x%02X).", b[1L], b[2L]))
    bed_flag <- b[3L] == 0L
    if (verbose)
    {
        cat("    bed file: '", bed.fn, "'", sep="")
        s <- ifelse(bed_flag, " (sample-major mode: [SNP, sample])",
            " (SNP-major mode: [sample, SNP])")
        if (.crayon()) s <- crayon::blurred(s)
        cat(s, "\n", sep="")
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
    {
        n <- nrow(famD)
        cat("    fam file: '", fam.fn, "' (", .pretty(n), " sample",
            .plural(n), ")\n", sep="")
    }

    ##  read bim.fn  ##

    f <- .open_text(bim.fn, TRUE)
    bimD <- read.table(f$con, header=FALSE, stringsAsFactors=FALSE)
    .close_conn(f)
    names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")
    if (verbose)
    {
        n <- nrow(bimD)
        cat("    bim file: '", bim.fn, "' (", .pretty(n), " variant",
            .plural(n), ")\n", sep="")
    }


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
    # RLE-coded chromosome
    .optim_chrom(dstfile)

    # add allele
    if (verbose) cat("    allele")
    n <- add.gdsn(dstfile, "allele", paste(bimD$allele2, bimD$allele1, sep=","),
        storage="string", compress=compress.annotation, closezip=TRUE)
    .DigestCode(n, digest, verbose)

    # add a folder for genotypes
    if (verbose)
        cat("    genotype", ifelse(verbose, ":\n", ""), sep="")
    n <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(n, "VariableName", "GT")
    put.attr.gdsn(n, "Description", "Genotype")

    # sync file
    sync.gds(dstfile)

    # add genotypes to genotype/data
    num4 <- ifelse(bed_flag, nrow(bimD), nrow(famD))
    cnt4 <- ifelse(bed_flag, nrow(famD), nrow(bimD))
    vg <- add.gdsn(n, "data", storage="bit2", valdim=c(2L, num4, 0L),
        compress=ifelse(bed_flag, "", compress.geno))

    if (pnum <= 1L)
    {
        # convert
        .Call(SEQ_ConvBED2GDS, vg, cnt4, bedfile$con, readBin, new.env(),
            if (verbose) stdout() else NULL)
        readmode.gdsn(vg)

    } else {
        # close bed file
        on.exit()
        .close_conn(bedfile)
        bedfile <- NULL
        # split to multiple files
        psplit <- .file_split(cnt4, pnum)
        # need unique temporary file names
        ptmpfn <- .get_temp_fn(pnum, basename(sub("^([^.]*).*", "\\1", out.gdsfn)),
            dirname(out.gdsfn))
        if (verbose)
        {
            cat(sprintf("    Writing to %d files:\n", pnum))
            cat(sprintf("        %s [%s..%s]\n", basename(ptmpfn),
                .pretty(psplit[[1L]]),
                .pretty(psplit[[1L]] + psplit[[2L]] - 1L)), sep="")
            flush.console()
        }

        # conversion in parallel
        seqParallel(parallel, NULL, FUN = function(bed.fn, tmp.fn, num4, psplit, cp)
        {
            library("SeqArray")
            # the process id, starting from one
            i <- process_index
            # open the bed file
            bedfile <- .open_bin(bed.fn)
            on.exit({ .close_conn(bedfile) })
            # create gds file
            f <- createfn.gds(tmp.fn[i])
            on.exit({ closefn.gds(f) }, add=TRUE)
            # progress file
            progfile <- file(paste0(tmp.fn[i], ".progress"), "wt")
            cat(tmp.fn[i], ":\n", file=progfile, sep="")
            on.exit({ close(progfile) }, add=TRUE)
            # new a gds node
            vg <- add.gdsn(f, "data", storage="bit2", valdim=c(2L, num4, 0L), compress=cp)
            # re-position the file
            cnt <- psplit[[2L]][i]
            if (cnt > 0L)
            {
                n4 <- (num4 %/% 4L) + (num4 %% 4L > 0L)
                n4 <- 3L + n4 * (psplit[[1L]][i] - 1L)
                if (bedfile$fmt == "")
                {
                    seek(bedfile$con, n4)
                } else {
                    # gz or xz
                    while (n4 > 0L)
                    {
                        m <- if (n4 <= 65536L) n4 else 65536L
                        readBin(bedfile$con, raw(), m)
                        n4 <- n4 - m
                    }
                }
                # convert
                .Call(SEQ_ConvBED2GDS, vg, cnt, bedfile$con, readBin, new.env(),
                    progfile)
            }
            readmode.gdsn(vg)
            invisible()
        }, split="none", bed.fn=bed.fn, tmp.fn=ptmpfn, num4=num4, psplit=psplit,
            cp=ifelse(bed_flag, "", compress.geno))

        if (verbose) cat("        merging files ...")
        for (i in seq_along(ptmpfn))
        {
            if (psplit[[2L]][i] > 0L)
            {
                f <- openfn.gds(ptmpfn[i])
                append.gdsn(vg, index.gdsn(f, "data"))
                closefn.gds(f)
            }
            unlink(ptmpfn[i], force=TRUE)
            unlink(paste0(ptmpfn[i], ".progress"), force=TRUE)
        }
        if (verbose) cat(" [Done]\n    ")
    }

    n1 <- add.gdsn(n, "@data", storage="uint8", compress=compress.annotation,
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
    if (!is.null(bedfile)) .close_conn(bedfile)

    if (bed_flag)
    {
        cat(" (transposed)")
        permdim.gdsn(vg, c(2L,1L))
    }
    if (verbose) cat("  ")
    .DigestCode(vg, digest, verbose)

    # add a folder for phase information
    if (verbose) cat("    phase")
    n <- addfolder.gdsn(dstfile, "phase")

    n1 <- add.gdsn(n, "data", storage="bit1", valdim=c(nrow(famD), 0L),
        compress=compress.annotation)
    .repeat_gds(n1, 0L, as.double(nrow(bimD))*nrow(famD))
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
    .repeat_gds(n1, 1L, nrow(bimD))
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

    ##################################################
    # optimize access efficiency

    if (verbose)
        cat("Done.\n", date(), "\n", sep="")
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



#######################################################################
# Convert a SeqArray GDS file to a PLINK BED file
#

seqGDS2BED <- function(gdsfile, out.fn, multi.row=FALSE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    if (is.character(gdsfile))
        stopifnot(length(gdsfile)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L, !is.na(out.fn))
    stopifnot(is.logical(multi.row), length(multi.row)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (multi.row) stop("'multi.row=TRUE' is not implemented yet.")

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("SeqArray GDS to PLINK BED:\n")
    }
    if (is.character(gdsfile))
    {
        if (verbose) cat("    open ", shQuote(gdsfile), "\n", sep="")
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))
    }
    if (verbose)
    {
        dm <- .seldim(gdsfile)
        cat("    # of samples: ", .pretty(dm[2L]), "\n", sep="")
        cat("    # of variants: ", .pretty(dm[3L]), "\n", sep="")
    }

    # fam file
    s <- seqGetData(gdsfile, "sample.id")
    n <- length(s)
    fam <- data.frame(FID=s, IID=s, FAT=rep(0L, n), MOT=rep(0L, n),
        sex=rep(0L, n), pheno=rep(-9L, n), stringsAsFactors=FALSE)
    famfn <- paste0(out.fn, ".fam")
    if (verbose) cat("    fam file: ", shQuote(famfn), "\n", sep="")
    write.table(fam, file=famfn, quote=FALSE, sep="\t", row.names=FALSE,
        col.names=FALSE)
    remove(fam)

    # bim file
    n <- .seldim(gdsfile)[3L]
    s <- seqGetData(gdsfile, "$alt")
    s[s == ""] <- "0"
    bim <- data.frame(
        chr = seqGetData(gdsfile, "chromosome"),
        var = seqGetData(gdsfile, "annotation/id"),
        mrg = rep(0L, n),
        pos = seqGetData(gdsfile, "position"),
        alt1 = gsub(",", "/", s, fixed=TRUE),
        alt2 = seqGetData(gdsfile, "$ref"),
        stringsAsFactors=FALSE)
    bimfn <- paste0(out.fn, ".bim")
    if (verbose) cat("    bim file: ", shQuote(bimfn), "\n", sep="")
    write.table(bim, file=bimfn, quote=FALSE, sep="\t", row.names=FALSE,
        col.names=FALSE)
    remove(bim)
    
    # bed file
    bedfn <- paste0(out.fn, ".bed")
    if (verbose) cat("    bed file: ", shQuote(bedfn), "\n", sep="")
    outf <- file(bedfn, "w+b")
    on.exit(close(outf), add=TRUE)
    writeBin(as.raw(c(0x6C, 0x1B, 0x01)), outf)
    seqApply(gdsfile, "$dosage", .cfunction("FC_GDS2BED"), as.is=outf,
        .useraw=TRUE, .progress=verbose)

    if (verbose)
        cat("Done.\n", date(), "\n", sep="")

    # output
    invisible(normalizePath(c(famfn, bimfn, bedfn)))
}
