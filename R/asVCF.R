#######################################################################
# Convert to a VariantAnnotation object
#

seqAsVCF <- function(x, chr.prefix="", info=NULL, geno=NULL)
{
    if (!requireNamespace("VariantAnnotation", quietly=TRUE, verbose=FALSE)) {
        stop("Please load VariantAnnotation to use this function")
    }
    
    stopifnot(is.character(chr.prefix), length(chr.prefix)==1L)

    sample.id <- seqGetData(x, "sample.id")
    hdr <- VariantAnnotation::VCFHeader(samples=sample.id, header=SeqArray::header(x))
    #hdr <- SeqArray::header(x)

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
    vcf <- VariantAnnotation::VCF(rowRanges=SeqArray::rowRanges(x),
                                  colData=SeqArray::colData(x),
                                  exptData=SimpleList(header=hdr),
                                  fixed=SeqArray::fixed(x),
                                  info=SeqArray::info(x, info=info),
                                  geno=SeqArray::geno(x, geno=geno))
    if (chr.prefix != "")
        vcf <- renameSeqlevels(vcf, paste0(chr.prefix, seqlevels(vcf)))
    vcf
}
