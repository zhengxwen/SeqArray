\name{seqVCF.Header}
\alias{seqVCF.Header}
\title{Parse the Header of a VCF File}
\description{
    Parses the header of a VCF file.
}
\usage{
seqVCF.Header(vcf.fn)
}
\arguments{
    \item{vcf.fn}{the file name of VCF}
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
    \item{header}{the other header lines}
    \item{num.ploidy}{the number of ploidy, two for humans}
}
\details{
    The ID description contains four columns:
    ID -- variable name;
    Number -- the number of elements, see the webpage of the 1000 Genomes Project;
    Type -- data type;
    Description -- a variable description.
}

\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
    Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF.SampID}}, \code{\link{seqVCF2GDS}}
}

\examples{
# the VCF file
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# get sample id
seqVCF.Header(vcf.fn)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
