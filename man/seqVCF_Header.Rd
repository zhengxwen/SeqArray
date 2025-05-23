\name{seqVCF_Header}
\alias{seqVCF_Header}
\title{Parse the Header of a VCF/BCF File}
\description{
    Parses the meta-information lines of a VCF or BCF file.
}
\usage{
seqVCF_Header(vcf.fn, getnum=FALSE, use_Rsamtools=NA, parallel=FALSE,
    verbose=TRUE)
}
\arguments{
    \item{vcf.fn}{the file name of VCF or BCF format;
        or a \code{\link{connection}} object for VCF format}
    \item{getnum}{if \code{TRUE}, return the total number of variants}
    \item{use_Rsamtools}{only applicable when \code{getnum=TRUE} and multiple
        cores are used;
        \code{NA}: using Rsamtools when it is installed; \code{FALSE}: not use
        the Rsamtools package; \code{TRUE}: to use Rsamtools, if it is not
        installed, the function fails}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{verbose}{when \code{getnum=TRUE} and \code{verbose=TRUE}, show
        the progress information for scanning the file}
}
\value{
    Return a list (with a class name "SeqVCFHeaderClass", S3 object):
    \item{fileformat}{the file format}
    \item{info}{the ID description in the INFO field}
    \item{filter}{the ID description in the FILTER field}
    \item{format}{the ID description in the FORMAT field}
    \item{alt}{the ID description in the ALT field}
    \item{contig}{the description in the contig field}
    \item{assembly}{the link of assembly}
    \item{reference}{genome reference, or \code{NULL} if unknown}
    \item{header}{the other header lines}
    \item{ploidy}{ploidy, two for humans}
    \item{num.sample}{the number of samples}
    \item{num.variant}{the number of variants, applicable only if
        \code{getnum=TRUE}}
    \item{sample.id}{a vector of sample IDs in the VCF/BCF file}
}
\details{
    The ID description contains four columns: \code{ID} -- variable name;
\code{Number} -- the number of elements, see the webpage of the 1000 Genomes
Project; \code{Type} -- data type; \code{Description} -- a variable
description.

    If multiple cores are used to get the total number of variants, the
Rsamtools package should be installed, and the indexing file (i.e., .csi or
.tbi) should be available along with \code{vcf.fn}.
}

\references{
    Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo,
M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., et al. (2011).
The variant call format and VCFtools. Bioinformatics 27, 2156-2158.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF_SampID}}, \code{\link{seqVCF2GDS}}
}

\examples{
# the VCF file
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# get sample id
seqVCF_Header(vcf.fn, getnum=TRUE)

# use a connection object
f <- file(vcf.fn, "r")
seqVCF_Header(f, getnum=TRUE)
close(f)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
