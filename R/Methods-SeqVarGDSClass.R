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


