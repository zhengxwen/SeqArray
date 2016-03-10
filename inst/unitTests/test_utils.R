test_refAllele <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T", "A", "AA")
  checkIdentical(c("A","AA","A","A", "A", "AA"), SeqArray:::.refAllele(x))
}

test_altAllele <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T", "A", "AA")
  checkIdentical(c("G","G","GG","G,T", "", ""), SeqArray:::.altAllele(x))
}

.test_variableLengthToList <- function() {
  x <- list(length=c(1,1,2,1,3,1),
            data=c(1,1,1:2,1,1:3,1))
  checkEquals(list(1,1,1:2,1,1:3,1), .variableLengthToList(x))
  x <- list(length=rep(1,5),
            data=1:5)
  checkEquals(list(1,2,3,4,5), .variableLengthToList(x))
}

.test_variableLengthToMatrix <- function() {
  x <- list(matrix(1:3, nrow=2, ncol=3, byrow=TRUE),
            matrix(4:6, nrow=2, ncol=3, byrow=TRUE),
            matrix(7:9, nrow=2, ncol=3, byrow=TRUE))
  xm <- matrix(list(1:3,1:3,4:6,4:6,7:9,7:9), nrow=2, ncol=3)
  checkIdentical(xm, .variableLengthToMatrix(x))
}

test_colData <- function() {
  ## test with no sample annotation  
  vcffile <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  gdsfile <- tempfile()
  seqVCF2GDS(vcffile, gdsfile, verbose=FALSE)
  gds <- seqOpen(gdsfile)
  annot <- colData(gds)
  checkTrue(ncol(annot) == 1)
  seqClose(gds)
  unlink(gdsfile)

  ## test with annotation
  gdsfile <- seqExampleFileName("gds")
  gds <- seqOpen(gdsfile)
  annot <- colData(gds)
  checkTrue(ncol(annot) > 1)
  seqClose(gds)
}
