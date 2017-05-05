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
    if (!.asVCF_Export())
    {
        fn <- system.file("extdata", "ext_as_vcf.r", package="SeqArray")
        source(fn)
        # set TRUE internally
        .asVCF_SetTRUE()
    }
    invisible()
}


seqAsVCF <- function(gdsfile, ...)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    seqAsVCFInit()
    eval(parse(text="asVCF(gdsfile, ...)"))
}
