.test_fixed <- function(fv, fg) {
  checkIdentical(as.character(fv$REF), as.character(fg$REF))
  checkIdentical(fv$ALT, fg$ALT)
  checkIdentical(fv$QUAL, fg$QUAL)
  ## VCF sets filter==NA to "."
  fg$FILTER[is.na(fg$FILTER)] <- "."
  checkIdentical(fv$FILTER, fg$FILTER)
}

.test_rowData <- function(rdv, rdg) {
  checkIdentical(seqnames(rdv), seqnames(rdg))
  checkIdentical(start(rdv), start(rdg))
  checkIdentical(width(rdv), width(rdg))
  checkIdentical(strand(rdv), strand(rdg))
  fcols <- c("REF", "ALT", "QUAL", "FILTER")
  .test_fixed(mcols(rdv)[,fcols], mcols(rdg)[,fcols])
}

.test_colData <- function(cdv, cdg) {
  checkIdentical(cdv, cdg)
}

.test_header <- function(hdv, hdg) {
  checkIdentical(samples(hdv), samples(hdg))
  ## The meta(VCFHeader) getter now returns a DataFrameList, not DataFrame.
  ## The META DataFrame is now one element of the DataFrameList.
  ## tags with non-alphanumeric characters get ignored by scanBcfHeader
  meta.hdg <- meta(hdg)$META[grep("^[[:alnum:]]+$", row.names(meta(hdg)$META)),,drop=FALSE]
  checkIdentical(meta(hdv)$META, meta.hdg)
  checkIdentical(fixed(hdv), fixed(hdg))
  checkIdentical(info(hdv), info(hdg))
  checkIdentical(geno(hdv), geno(hdg))
}

.test_info <- function(iv, ig) {
  checkIdentical(sort(names(iv)), sort(names(ig)))
  for (n in names(iv)) {
    if (is(iv[[n]], "IntegerList")) {
      iv[[n]] <- IntegerList(lapply(iv[[n]], function(x) {if (length(x) > 0) x else NA}))
    } else if (is(iv[[n]], "CharacterList")) {
      iv[[n]] <- CharacterList(lapply(iv[[n]], function(x) {if (length(x) > 0) x else NA}))
    }
    checkEquals(iv[[n]], ig[[n]], paste(" ", n, "not identical"),
                tolerance=2*(.Machine$double.eps^0.5), checkNames=FALSE)
  }
}

.test_geno <- function(gv, gg) {
  for (n in names(gv)) {
    ## variant names are different
    checkEquals(nrow(gv[[n]]), nrow(gg[[n]]),
                paste(" ", n, "have different numbers of rows"))
    dimnames(gv[[n]])[1] <- dimnames(gg[[n]])[1]
    checkEquals(gv[[n]], gg[[n]], paste(" ", n, "not identical"),
                tolerance=10*(.Machine$double.eps^0.5), checkNames=FALSE)
  }
}

.test_asVCF <- function(vcffile, gdsfile) {
  vcf <- readVcf(vcffile, genome="hg19")
  gdsobj <- seqOpen(gdsfile)

##   .test_rowData(rowData(vcf), .rowData(gdsobj))
##   .test_colData(colData(vcf), .colData(gdsobj))
##   .test_header(header(vcf), .header(gdsobj))
##   .test_info(info(vcf), .info(gdsobj))
##   .test_geno(geno(vcf), .geno(gdsobj))

  vcfg <- asVCF(gdsobj)
  .test_rowData(rowData(vcf), rowData(vcfg))
  .test_colData(colData(vcf), colData(vcfg))
  .test_header(header(vcf), header(vcfg))
  .test_fixed(fixed(vcf), fixed(vcfg))
  .test_info(info(vcf), info(vcfg))
  .test_geno(geno(vcf), geno(vcfg))

  seqClose(gdsobj)
}

test_asVCF <- function() {
  vcffile <- seqExampleFileName("vcf")
  gdsfile <- seqExampleFileName("gds")
  .test_asVCF(vcffile, gdsfile)
}

test_asVCF_filterInHead <- function() {
  vcffile <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  gdsfile <- tempfile()
  seqVCF2GDS(vcffile, gdsfile)
  .test_asVCF(vcffile, gdsfile)
  unlink(gdsfile)
}

test_asVCF_altInHead <- function() {
  vcffile <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
  gdsfile <- tempfile()
  seqVCF2GDS(vcffile, gdsfile)
  .test_asVCF(vcffile, gdsfile)
  unlink(gdsfile)
}

## takes too long - use for development only
## test_asVCF_c22 <- function() {
##   vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
##   gdsfile <- tempfile()
##   seqVCF2GDS(vcffile, gdsfile)
##   .test_asVCF(vcffile, gdsfile)
##   unlink(gdsfile)
## }

test_info_geno <- function() {
  vcffile <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
  gdsfile <- tempfile()
  seqVCF2GDS(vcffile, gdsfile)

  info <- c("AN", "VT")
  geno <- "DS"

  vcf <- readVcf(vcffile, genome="hg19",
                 param=ScanVcfParam(info=info, geno=geno))
  gdsobj <- seqOpen(gdsfile)

  vcfg <- asVCF(gdsobj, info=info, geno=geno)
  .test_header(header(vcf), header(vcfg))
  .test_info(info(vcf), info(vcfg))
  .test_geno(geno(vcf), geno(vcfg))

  seqClose(gdsobj)
  unlink(gdsfile)
}

test_info_geno_na <- function() {
  vcffile <- seqExampleFileName("vcf")
  gdsfile <- seqExampleFileName("gds")
  info <- NA
  geno <- NA
  vcf <- readVcf(vcffile, genome="hg19",
                 param=ScanVcfParam(info=info, geno=geno))
  gdsobj <- seqOpen(gdsfile)

  vcfg <- asVCF(gdsobj, info=info, geno=geno)
  .test_header(header(vcf), header(vcfg))
  .test_info(info(vcf), info(vcfg))
  .test_geno(geno(vcf), geno(vcfg))

  seqClose(gdsobj)
}

test_filters <- function() {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  gdsfile <- tempfile()
  seqVCF2GDS(vcffile, gdsfile)

  gdsobj <- seqOpen(gdsfile)
  samples <- seqGetData(gdsobj, "sample.id")[1:5]
  variants <- seqGetData(gdsobj, "variant.id")[1:10]
  seqSetFilter(gdsobj, sample.id=samples, variant.id=variants)

  info <- c("AA", "AN")
  geno <- "GT"
  vcf <- readVcf(vcffile, genome="hg19",
                 param=ScanVcfParam(info=info, geno=geno,
                   samples=samples, which=granges(gdsobj)))
  vcfg <- asVCF(gdsobj, info=info, geno=geno)

  .test_rowData(rowData(vcf), rowData(vcfg))
  .test_colData(colData(vcf), colData(vcfg))
  .test_header(header(vcf), header(vcfg))
  .test_fixed(fixed(vcf), fixed(vcfg))
  .test_info(info(vcf), info(vcfg))
  .test_geno(geno(vcf), geno(vcfg))

  seqClose(gdsobj)
  unlink(gdsfile)
}

