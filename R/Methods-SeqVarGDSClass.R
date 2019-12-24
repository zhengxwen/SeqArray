# new a SeqVarGDSClass instance
SeqVarGDSClass <- function(gdsobj, ...)
{
    new("SeqVarGDSClass", gdsobj, ...)
}


# importFrom(GenomicRanges, granges)
setMethod("granges", "SeqVarGDSClass",
    function(x, ...)
    {
        vid <- seqGetData(x, "variant.id")
        chr <- seqGetData(x, "chromosome")
        pos <- seqGetData(x, "position")
        reflen <- nchar(seqGetData(x, "$ref"))
        reflen[reflen < 1L] <- 1L
        gr <- GRanges(seqnames=chr,
            ranges=IRanges(start=pos, end=pos+reflen-1L),
            ...)
        names(gr) <- vid
        gr
    })



####  from VariantAnnotation  ####

setMethod("ref", "SeqVarGDSClass", function(x)
{
    s <- seqGetData(x, "$ref")
    # remove invalid characters
    s <- gsub("[^ACGTMRWSYKVHDBNacgtmrwsykvhdbn\\-\\+\\.]", ".", s)
    DNAStringSet(s)
})


setMethod("alt", "SeqVarGDSClass", function(x)
{
    alt <- seqGetData(x, "$alt")
    s <- strsplit(alt, ",", fixed=TRUE)
    s[alt == ""] <- ""
    # remove invalid characters
    s <- sapply(s, function(x)
        gsub("[^ACGTMRWSYKVHDBNacgtmrwsykvhdbn\\-\\+\\.]", ".", x),
        simplify=FALSE)
    do.call(DNAStringSetList, s)
})


setMethod("qual", "SeqVarGDSClass", function(x)
{
    qual <- seqGetData(x, "annotation/qual")
    qual[is.na(qual)] <- NA	 # change NaN to NA
    qual
})


setMethod("filt", "SeqVarGDSClass", function(x)
{
    as.character(seqGetData(x, "annotation/filter"))
})


setMethod("fixed", "SeqVarGDSClass", function(x)
{
    DataFrame(REF=SeqArray::ref(x),
              ALT=SeqArray::alt(x),
              QUAL=SeqArray::qual(x),
              FILTER=SeqArray::filt(x))
})


setMethod("header", "SeqVarGDSClass", function(x)
{
    ## info
    seqsum <- seqSummary(x, check="none", verbose=FALSE)
    infoHd <- seqsum$info
    # names(infoHd)[2:4] <- c("Number", "Type", "Description")
    infoHd <- DataFrame(infoHd[,2:4], row.names=infoHd[,1L])

    ## geno
    genoHd <- seqsum$format
    # names(genoHd)[2:4] <- c("Number", "Type", "Description")
    genoHd <- DataFrame(genoHd[,2:4], row.names=genoHd[,1L])

    ## meta
    des <- get.attr.gdsn(index.gdsn(x, "description"))
    ff <- des$vcf.fileformat
    meta <- DataFrameList(fileformat=DataFrame(Value=ff, row.names="fileformat"))
    ref <- seqsum$reference
    if (length(ref) > 0L)
    {
        meta <- c(meta, DataFrameList(reference=DataFrame(Value=ref, row.names="reference")))
    }
    n <- index.gdsn(x, "description/vcf.header", silent=TRUE)
    if (!is.null(n))
    {
        des <- read.gdsn(n)
        ## ID=value header fields not parsed in GDS
        des <- des[!(des[,1L] %in% c("contig", "SAMPLE", "PEDIGREE")),]
        fields <- unique(des[,1L])
        for (f in fields) {
            des.f <- des[des[,1L] %in% f,,drop=FALSE]
            meta.des <- DataFrameList(DataFrame(Value=des.f[,2L], row.names=make.unique(des.f[,1L])))
            names(meta.des) <- f
            meta <- c(meta, meta.des)
        }
    }

    hdr <- c(meta, DataFrameList(INFO=infoHd, FORMAT=genoHd))

    ## fixed
    des <- seqsum$allele
    if (nrow(des) > 0L)
    {
        hdr[["ALT"]] <- DataFrame(Description=des[,2L], row.names=des[,1L])
    }
    des <- seqsum$filter
    des <- des[des$Description != "" & !is.na(des$Description),,drop=FALSE]
    if (nrow(des) > 0L)
    {
        hdr[["FILTER"]] <- DataFrame(Description=des[,2L], row.names=des[,1L])
    }

    hdr
})


setMethod("info", "SeqVarGDSClass", function(x, infovar=NULL)
{
    des <- seqSummary(x, "annotation/info", check="none", verbose=FALSE)
    if (!is.null(infovar))
        des <- des[des$ID %in% infovar, ]
    infoDf <- DataFrame(row.names=seqGetData(x, "variant.id"))
    if (nrow(des) > 0L)
    {
        for (i in seq_len(nrow(des)))
        {
            v <- seqGetData(x, paste0("annotation/info/", des$ID[i]), .padNA=FALSE)
            ## deal with variable length fields
            if (!is.null(names(v)))
            {
                vl <- .variableLengthToList(v)
                ## each element should have length number of alt alleles, even for NAs
                if (des$Number[i] == "A")
                {
                    nAlt <- lengths(SeqArray::alt(x))
                    addNA <- which(nAlt > 1L & is.na(vl))
                    for (ind in addNA)
                    {
                        vl[[ind]] <- rep(NA, nAlt[ind])
                    }
                }
                v <- .toAtomicList(vl, des$Type[i])
            } else if (!is.null(dim(v)))
            {
                ## v is a matrix with nrow="Number"
                vl <- list()
                for (j in 1:ncol(v))
                {
                    vl[[j]] <- v[,j]
                }
                v <- .toAtomicList(vl, des$Type[i])
            }
            if (is.atomic(v))
            {
                v[is.na(v)] <- NA  # change NaN to NA
                v[v %in% ""] <- NA
            }
            infoDf[[des$ID[i]]] <- v
        }
    }
    infoDf
})


setMethod("geno", "SeqVarGDSClass", function(x, genovar=NULL)
{
    ## genotype
    sample.id <- seqGetData(x, "sample.id")
    variant.id <- seqGetData(x, "variant.id")

    if (is.null(genovar) || "GT" %in% genovar)
    {
        gt <- seqApply(x, c(genovar="genotype", phase="phase"),
                       function(x) {
                           sep <- ifelse(x$phase, "|", "/")
                           paste(x$genovar[1L,], sep, x$genovar[2L,], sep="")
                       },
                       as.is="list", margin="by.variant")
        gt <- matrix(unlist(gt), ncol=length(gt),
                     dimnames=list(sample.id, variant.id))
        gt[gt == "NA/NA"] <- "."
        gt <- t(gt)

        genoSl <- SimpleList(GT=gt)
    } else {
        genoSl <- SimpleList()
    }

    ## all other fields
    des <- seqSummary(x, "annotation/format", check="none", verbose=FALSE)
    if (!is.null(genovar))
    {
        des <- des[des$ID %in% genovar,]
    }
    if (nrow(des) > 0L)
    {
        for (i in 1:nrow(des))
        {
            var.name <- paste("annotation/format/", des$ID[i], sep="")
            number <- suppressWarnings(as.integer(des$Number[i]))
            if (!is.na(number) && number > 1L)
            {
                v <- seqApply(x, var.name, function(v) v,
                              as.is="list", margin="by.variant")
                vm <- array(
                    dim=c(length(variant.id), length(sample.id), number),
                    dimnames=list(variant.id, sample.id, NULL))
                for (j in 1:length(v))
                {
                    if (is.null(v[[j]]))
                    {
                        vm[j,,] <- NA
                    } else {
                        vm[j,,] <- v[[j]]
                    }
                }
                v <- vm
            } else {
                v <- seqGetData(x, var.name, .padNA=FALSE)
                if (!is.null(names(v)))
                {
                    if (all(v$length == 1L) && !is.na(number) && number == 1L)
                    {
                        v <- v$data
                    } else {
                        v <- seqApply(x, var.name, function(v) v,
                                      as.is="list", margin="by.variant")
                        v <- .variableLengthToMatrix(v)
                    }
                }
                dimnames(v) <- list(sample.id, variant.id)
                v <- t(v)
            }
            genoSl[[des$ID[i]]] <- v
        }
    }
    if (is.null(names(genoSl))) names(genoSl) <- character()
    genoSl
})



####  from SummarizedExperiment	 ####

setMethod("rowRanges", "SeqVarGDSClass", function(x)
{
    granges(x,
        ID=seqGetData(x, "annotation/id"),
        REF=SeqArray::ref(x),
        ALT=SeqArray::alt(x),
        QUAL=SeqArray::qual(x),
        FILTER=SeqArray::filt(x))
})


setMethod("colData", "SeqVarGDSClass", function(x)
{
    sample.id <- seqGetData(x, "sample.id")
    node <- index.gdsn(x, "sample.annotation", silent=TRUE)
    if (!is.null(node))
        vars <- ls.gdsn(node)
    else
        vars <- character()
    if (length(vars) > 0)
    {
        annot <- lapply(vars, function(v) {
            seqGetData(x, paste0("sample.annotation/", v), .padNA=FALSE)
        })
        names(annot) <- vars
        DataFrame(Samples=seq_along(sample.id), annot, row.names=sample.id)
    } else {
        DataFrame(Samples=seq_along(sample.id), row.names=sample.id)
    }
})
