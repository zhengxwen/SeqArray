\name{SeqVarGDSClass}
\docType{class}

\alias{SeqVarGDSClass-class}
\alias{SeqVarGDSClass}
\alias{granges,SeqVarGDSClass-method}
\alias{ref,SeqVarGDSClass-method}
\alias{alt,SeqVarGDSClass-method}
\alias{filt,SeqVarGDSClass-method}
\alias{qual,SeqVarGDSClass-method}
\alias{asVCF,SeqVarGDSClass-method}

\title{SeqVarGDSClass}

\description{
	A \code{SeqVarGDSClass} object provides access to a GDS file containing
Variant Call Format (VCF) data. It extends \code{\link[gdsfmt]{gds.class}}.
}

\details{
	A sequencing GDS file is created from a VCF file with
\code{\link{seqVCF2GDS}}. This file can be opened with \code{\link{seqOpen}}
to create a \code{SeqVarGDSClass} object.
}

\section{Accessors}{
	In the following code snippets \code{x} is a SeqVarGDSClass object.
	\describe{
		\item{}{
			\code{granges(x)}:
				Returns the chromosome and position of variants as a
				\code{GRanges} object.  Names correspond to the variant.id.
		}
		\item{}{
			\code{ref(x)}:
				Returns the reference alleles as a \code{\link{DNAStringSet}}.
		}
		\item{}{
			\code{alt(x)}:
			Returns the alternate alleles as a \code{\link{DNAStringSetList}}.
		}
		\item{}{
			\code{filt(x)}:
			Returns the filter data.
		}
		\item{}{
			\code{qual(x)}:
			Returns the quality scores.
		}
	}

	Other data can be accessed with \code{\link{seqGetData}}.
}

\section{Coercion methods}{
	In the following code snippets \code{x} is a SeqVarGDSClass object.
	\describe{
		\item{}{
			\code{asVCF(x, info=NULL, geno=NULL)}:
			Coerces a SeqVarGDSClass object to a
			\link{VCF-class} object.  Row names correspond
	to the variant.id.  \code{info}
			and \code{geno} specify the 'INFO' and
	'GENO' (FORMAT) fields to return, respectively.  If not
	specified, all fields are returned; if 'NA' no fields are
	returned.  Use \code{\link{seqSetFilter}} prior to calling
	\code{asVCF} to specify samples and variants to return.
		}
	}
}


\author{Xiuwen Zheng, Stephanie Gogarten}

\seealso{
  \code{\link[gdsfmt]{gds.class}},
  \code{\link{seqVCF2GDS}},
  \code{\link{seqOpen}},
  \code{\link{seqGetData}},
  \code{\link{seqSetFilter}},
  \code{\link{seqClose}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
gds

## sample ID
head(seqGetData(gds, "sample.id"))

## variants
granges(gds)

## alleles as comma-separated character strings
head(seqGetData(gds, "allele"))

## alleles as DNAStringSet or DNAStringSetList
ref(gds)
alt(gds)

## genotype
geno <- seqGetData(gds, "genotype")
dim(geno)
## dimensions are: allele, sample, variant
geno[1,1:10,1:5]

## rsID
head(seqGetData(gds, "annotation/id"))

## alternate allele count
head(seqGetData(gds, "annotation/info/AC"))

## individual read depth
depth <- seqGetData(gds, "annotation/format/DP")
names(depth)
## VCF header defined DP as variable-length data
table(depth$length)
## all length 1, so depth$data should be a sample by variant matrix
dim(depth$data)
depth$data[1:10,1:5]

seqClose(gds)
}
