
#######################################################################
# Convert to a VariantAnnotation object
#

.variableLengthToList <- function(x)
{
    xl <- list()
    j <- 1L
    for (i in 1:length(x$length))
    {
        len <- x$length[i]
        if (len > 0L)
        {
            xl[[i]] <- x$data[j:(j+len-1)]
            j <- j+len
        } else {
            xl[[i]] <- NA
        }
    }
    xl
}

.toAtomicList <- function(x, type)
{
    switch(type,
        Character=CharacterList(x),
        String=CharacterList(x),
        Integer=IntegerList(x),
        Float=NumericList(x))
}

.variableLengthToMatrix <- function(x)
{
    xl <- list()
    i <- 1L
    for (j in seq_along(x))
    {
        for (k in seq_len(NROW(x[[j]])))
        {
            xl[[i]] <- x[[j]][k,]
            i <- i + 1L
        }
    }
    matrix(xl, nrow=NROW(x[[1L]]), ncol=length(x))
}


seqAsVCFInit <- function()
{
    # check flag
    if (.asVCF_Export()) return(invisible())

    # load packages to get generic methods
    requireNamespace("SummarizedExperiment", quietly=TRUE, verbose=FALSE)
    requireNamespace("VariantAnnotation", quietly=TRUE, verbose=FALSE)


    ####  from VariantAnnotation  ####

    setMethod("ref", "SeqVarGDSClass", function(x)
    {
	s <- seqGetData(x, "$ref")
	# remove invalid characters
	s <- gsub("[^ACGTMRWSYKVHDBNacgtmrwsykvhdbn\\-\\+\\.]", ".", s)
	Biostrings::DNAStringSet(s)
    }, where=globalenv())


    setMethod("alt", "SeqVarGDSClass", function(x)
    {
	alt <- seqGetData(x, "$alt")
	s <- strsplit(alt, ",", fixed=TRUE)
	s[alt == ""] <- ""
	# remove invalid characters
	s <- sapply(s, function(x)
            gsub("[^ACGTMRWSYKVHDBNacgtmrwsykvhdbn\\-\\+\\.]", ".", x),
            simplify=FALSE)
	do.call(Biostrings::DNAStringSetList, s)
    }, where=globalenv())


    setMethod("qual", "SeqVarGDSClass", function(x)
    {
	qual <- seqGetData(x, "annotation/qual")
	qual[is.na(qual)] <- NA	 # change NaN to NA
	qual
    }, where=globalenv())


    setMethod("filt", "SeqVarGDSClass", function(x)
    {
	as.character(seqGetData(x, "annotation/filter"))
    }, where=globalenv())


    setMethod("fixed", "SeqVarGDSClass", function(x)
    {
	S4Vectors::DataFrame(REF=ref(x),
                  ALT=alt(x),
                  QUAL=qual(x),
                  FILTER=filt(x))
    }, where=globalenv())

    
    setMethod("header", "SeqVarGDSClass", function(x)
    {
	sample.id <- seqGetData(x, "sample.id")

	## info
	seqsum <- seqSummary(x, check="none", verbose=FALSE)
	infoHd <- seqsum$info
	# names(infoHd)[2:4] <- c("Number", "Type", "Description")
	infoHd <- S4Vectors::DataFrame(infoHd[,2:4], row.names=infoHd[,1L])

	## geno
	genoHd <- seqsum$format
	# names(genoHd)[2:4] <- c("Number", "Type", "Description")
	genoHd <- S4Vectors::DataFrame(genoHd[,2:4], row.names=genoHd[,1L])

	## meta
	des <- get.attr.gdsn(index.gdsn(x, "description"))
	ff <- des$vcf.fileformat
	meta <- S4Vectors::DataFrame(Value=ff, row.names="fileformat")
	ref <- seqsum$reference
	if (length(ref) > 0L)
	{
            meta <- rbind(meta, S4Vectors::DataFrame(Value=ref, row.names="reference"))
	}
	n <- index.gdsn(x, "description/vcf.header", silent=TRUE)
	if (!is.null(n))
	{
            des <- read.gdsn(n)
            ## ID=value header fields not parsed in GDS
            des <- des[!(des[,1L] %in% c("contig", "SAMPLE", "PEDIGREE")),]
            ## deal with duplicate row names
            if (any(duplicated(des[,1L])))
            {
                des[,1] <- make.unique(des[,1L])
            }
            meta <- rbind(meta, S4Vectors::DataFrame(Value=des[,2L], row.names=des[,1L]))
	}

	hdr <- IRanges::DataFrameList(META=meta, INFO=infoHd, FORMAT=genoHd)

	## fixed
	des <- seqsum$allele
	if (nrow(des) > 0L)
	{
            hdr[["ALT"]] <- S4Vectors::DataFrame(Description=des[,2L], row.names=des[,1L])
	}
	des <- seqsum$filter
	des <- des[des$Description != "" & !is.na(des$Description),,drop=FALSE]
	if (nrow(des) > 0L)
	{
            hdr[["FILTER"]] <- S4Vectors::DataFrame(Description=des[,2L], row.names=des[,1L])
	}

	VariantAnnotation::VCFHeader(samples=sample.id, header=hdr)
    }, where=globalenv())


    setMethod("info", "SeqVarGDSClass", function(x, infovar=NULL)
    {
	des <- seqSummary(x, "annotation/info", check="none", verbose=FALSE)
	if (!is.null(infovar))
            des <- des[des$ID %in% infovar, ]
	infoDf <- S4Vectors::DataFrame(row.names=seqGetData(x, "variant.id"))
	if (nrow(des) > 0L)
	{
            for (i in seq_len(nrow(des)))
            {
                v <- seqGetData(x, paste0("annotation/info/", des$ID[i]))
                ## deal with variable length fields
                if (!is.null(names(v)))
                {
                    vl <- .variableLengthToList(v)
                    ## each element should have length number of alt alleles, even for NAs
                    if (des$Number[i] == "A")
                    {
                        nAlt <- lengths(alt(x))
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
    }, where=globalenv())


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

            genoSl <- S4Vectors::SimpleList(GT=gt)
	} else {
            genoSl <- S4Vectors::SimpleList()
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
                    v <- seqGetData(x, var.name)
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
    }, where=globalenv())


    setMethod("asVCF", "SeqVarGDSClass", function(x, chr.prefix="", info=NULL, geno=NULL)
    {
                  stopifnot(is.character(chr.prefix), length(chr.prefix)==1L)

                  seqsum <- seqSummary(x, check="none", verbose=FALSE)
                  if (!is.null(info))
                  {
                      validInfo <- seqsum$info$ID
                      infoDiff <- setdiff(info, c(validInfo, NA))
                      if (length(infoDiff) > 0)
                      {
                          warning(paste("info fields not present:", infoDiff))
                          info <- intersect(info, validInfo)
                      }
                  }
                  if (!is.null(geno))
                  {
                      validGeno <- seqsum$format$ID
                      genoDiff <- setdiff(geno, c(validGeno, NA))
                      if (length(genoDiff) > 0)
                      {
                          warning(paste("geno fields not present:", genoDiff))
                          geno <- intersect(geno, validGeno)
                      }
                  }
                  vcf <- VariantAnnotation::VCF(rowRanges=rowRanges(x),
                             colData=colData(x),
                             exptData=S4Vectors::SimpleList(header=header(x)),
                             fixed=fixed(x),
                             info=info(x, info=info),
                             geno=geno(x, geno=geno))
                  if (chr.prefix != "")
                      vcf <- GenomeInfoDb::renameSeqlevels(vcf, paste0(chr.prefix, GenomeInfoDb::seqlevels(vcf)))
                  vcf
    }, where=globalenv())



    ####  from SummarizedExperiment	 ####

    setMethod("rowRanges", "SeqVarGDSClass", function(x)
    {
	granges(x,
		ID=seqGetData(x, "annotation/id"),
		REF=ref(x),
		ALT=alt(x),
		QUAL=qual(x),
		FILTER=filt(x))
    }, where=globalenv())

    
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
                seqGetData(x, paste0("sample.annotation/", v))
            })
            names(annot) <- vars
            S4Vectors::DataFrame(Samples=seq_along(sample.id), annot, row.names=sample.id)
	} else {
            S4Vectors::DataFrame(Samples=seq_along(sample.id), row.names=sample.id)
	}
    }, where=globalenv())

    

    # set TRUE internally
    .asVCF_SetTRUE()
    invisible()
}


seqAsVCF <- function(gdsfile, ...)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    seqAsVCFInit()
    eval(parse(text="asVCF(gdsfile, ...)"))
}
