# define new generic functions
setGeneric("seqClose", function(object, ...)
    standardGeneric("seqClose"))

setGeneric("seqSetFilter", function(object, variant.sel, ...)
    standardGeneric("seqSetFilter"))

# setGeneric("seqAppend", function(object, ...) standardGeneric("seqAppend"))

# redefine generics from VariantAnnotation and SummarizedExperiment to reduce package load time
setGeneric("ref", function(x, ...) standardGeneric("ref"))
setGeneric("alt", function(x, ...) standardGeneric("alt"))
setGeneric("qual", function(x, ...) standardGeneric("qual"))
setGeneric("filt", function(x, ...) standardGeneric("filt"))
setGeneric("fixed", function(x, ...) standardGeneric("fixed"))
setGeneric("header", function(x, ...) standardGeneric("header"))
setGeneric("info", function(x, ...) standardGeneric("info"))
setGeneric("geno", function(x, ...) standardGeneric("geno"))
setGeneric("rowRanges", function(x, ...) standardGeneric("rowRanges"))
setGeneric("colData", function(x, ...) standardGeneric("colData"))
