\name{seqAddValue}
\alias{seqAddValue}
\title{Add values to a GDS File}
\description{
    Add or modify the values in a GDS file with hash code
}
\usage{
seqAddValue(gdsfile, varnm, val, desp=character(), replace=FALSE, compress="LZMA_RA",
    packed=TRUE, packed.idx=TRUE, use_float32=TRUE, param=NULL,
    verbose=TRUE, verbose.attr=TRUE)
}
\arguments{
    \item{gdsfile}{character for file name, or a \code{\link{SeqVarGDSClass}}
        object}
    \item{varnm}{the variable name, e.g., "sample.id", "variant.id",
        "chromosome", "annotation/info/NEW_VARIABLE"}
    \item{val}{the R value can be integers, real numbers, characters,
        factor, logical, raw variable, \code{data.frame} or a list; a list of
        vectors is used for variable-length annotation data;
        \code{LogicalList}, \code{IntegerList}, \code{NumericList} and
        \code{CharacterList} from the \pkg{IRanges} package are also
        supported for variable-length data; or \code{NULL} for adding a
        new folder.
        For "annotation/info/VARIABLE", \code{val} can also be a function
        \code{val(k)} or \code{val(k, param)} returning the k-th block of
        values (any of the vector or list types above), or \code{NULL} if
        there is no more block. It is called with \code{k=1}, 2, 3, \ldots,
        i.e., the blocks are requested in order and each \code{k} is used
        exactly once, so \code{val} needs no internal state; it stops at the
        first \code{NULL}. The blocks are written to the GDS file one by one,
        so that the memory
        usage does not depend on the total number of variants. Note that
        \code{packed} and \code{packed.idx} are ignored in this case, since
        the data are unknown before writing}
    \item{desp}{variable description}
    \item{replace}{if \code{TRUE}, replace the existing variable silently if
        possible}
    \item{compress}{the compression method can be "" (no compression), see
        \code{\link[gdsfmt]{add.gdsn}}}
    \item{packed}{\code{TRUE}, pack data if there is any missing value}
    \item{packed.idx}{\code{TRUE}, store the index variable using integers with
        the fewest bits if possible}
    \item{use_float32}{if \code{TRUE}, use 32 bits instead of 64 bits to store
        floating numeric values (i.e., \code{double})}
    \item{param}{the second argument passed to \code{val(k, param)}, when
        \code{val} is a function taking two arguments; it is ignored
        otherwise}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.attr}{if \code{TRUE}, show attribute information in a GDS
        node}
}
\value{
    Return none.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{seqNewVarData}},
    \code{\link[gdsfmt]{add.gdsn}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
file.copy(gds.fn, "tmp.gds", overwrite=TRUE)

# display
(f <- seqOpen("tmp.gds", readonly=FALSE))

show(index.gdsn(f, "sample.id"))
seqAddValue(f, "sample.id", 1:90, replace=TRUE)
show(index.gdsn(f, "sample.id"))

show(index.gdsn(f, "chromosome"))
v <- seqGetData(f, "chromosome")
seqAddValue(f, "chromosome", paste0("chr", v), replace=TRUE)
show(index.gdsn(f, "chromosome"))
table(seqGetData(f, "chromosome"))

# annotation info
seqAddValue(f, "annotation/info/folder", NULL)  # add a new folder
seqAddValue(f, "annotation/info/folder/val", 1:1348, "random number")
seqAddValue(f, "annotation/info/folder/logical", rep(c(TRUE, FALSE), length.out=1348))
seqAddValue(f, "annotation/info/folder/packed", c(rep(2L, 1000), rep(NA, 348)))
seqAddValue(f, "annotation/info/newff",
    data.frame(x=1:1348, y=rep("s", 1348), stringsAsFactors=FALSE),
    desp=c("integer numbers", "character"))

# variable-length annotation info data
v <- lapply(1:1348, function(x) as.character(x))
v[[1]] <- 1:10
seqAddValue(f, "annotation/info/folder/val1", v)
head(seqGetData(f, "annotation/info/folder/val1", .tolist=TRUE))

# variable-length annotation info data using IntegerList
# v2 <- IRanges::IntegerList(lapply(1:1348, function(x) seq_len(x %% 5)))
# seqAddValue(f, "annotation/info/folder/val2", v2)
# head(seqGetData(f, "annotation/info/folder/val2", .tolist=TRUE))


# write the data block by block to reduce the memory usage, when 'val' is
#   a function returning the k-th block, or NULL if there is no more block
nvar <- 1348L; bsize <- 500L    # the number of variants in each block
nblock <- ceiling(nvar / bsize)
blk_st <- (seq_len(nblock) - 1L) * bsize + 1L
blk_cnt <- diff(c(blk_st, nvar + 1L))

# fixed-length data (one value per variant)
seqAddValue(f, "annotation/info/folder/val3", function(k)
    {
        if (k > nblock) return(NULL)    # no more block
        ii <- seq.int(blk_st[k], length.out=blk_cnt[k])
        ii * 1.5    # the values of the k-th block
    }, "a number per variant")
head(seqGetData(f, "annotation/info/folder/val3"))

# variable-length data, using the compact list(length, data) representation
#   instead of building a list of blk_cnt[k] elements
seqAddValue(f, "annotation/info/folder/val4", function(k)
    {
        if (k > nblock) return(NULL)
        ii <- seq.int(blk_st[k], length.out=blk_cnt[k])
        len <- ii %% 3L    # the number of values for each variant
        structure(list(length=len, data=rep(ii, len)), class="SeqVarDataList")
    }, "0..2 values per variant")
head(seqGetData(f, "annotation/info/folder/val4", .tolist=TRUE))

# 'param' is passed to val(k, param), if val() takes two arguments
seqAddValue(f, "annotation/info/folder/val5", function(k, param)
    {
        if (k > param$nblock) return(NULL)
        seq.int(param$st[k], length.out=param$cnt[k]) * 2L
    }, "a number per variant",
    param=list(nblock=nblock, st=blk_st, cnt=blk_cnt))
head(seqGetData(f, "annotation/info/folder/val5"))


# sample annotation
seqAddValue(f, "sample.annotation", data.frame(ii=1:90, y=rep("A", 90)),
    replace=TRUE)
seqAddValue(f, "sample.annotation/float", (1:90)/90)

# close the GDS file
seqClose(f)


# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
