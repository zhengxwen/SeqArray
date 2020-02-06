\name{SeqVarGDSClass}
\docType{class}

\alias{SeqVarGDSClass-class}
\alias{SeqVarGDSClass}
\alias{granges,SeqVarGDSClass-method}
\alias{ref,SeqVarGDSClass-method}
\alias{alt,SeqVarGDSClass-method}
\alias{filt,SeqVarGDSClass-method}
\alias{qual,SeqVarGDSClass-method}
\alias{fixed,SeqVarGDSClass-method}
\alias{header,SeqVarGDSClass-method}
\alias{info,SeqVarGDSClass-method}
\alias{geno,SeqVarGDSClass-method}
\alias{geno,SeqVarGDSClass,ANY-method}
\alias{rowRanges,SeqVarGDSClass-method}
\alias{colData,SeqVarGDSClass-method}
\alias{ref}
\alias{alt}
\alias{filt}
\alias{qual}
\alias{fixed}
\alias{header}
\alias{info}
\alias{geno}
\alias{rowRanges}
\alias{colData}

\title{SeqVarGDSClass}

\description{
	A \code{SeqVarGDSClass} object provides access to a GDS file containing
Variant Call Format (VCF) data. It extends \code{\link{gds.class}}.
}

\details{
	A SeqArray GDS file is created from a VCF file with
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
			\code{qual(x)}:
			Returns the quality scores.
		}
		\item{}{
			\code{filt(x)}:
			Returns the filter data.
		}
		\item{}{
			\code{fixed(x)}:
			Returns the fixed fields (ref, alt, qual, filt).
		}
		\item{}{
			\code{header(x)}:
			Returns the header as a \code{\link{DataFrameList}}.
		}
		\item{}{
			\code{rowRanges(x)}:
			Returns a \code{GRanges} object with metadata.
		}
		\item{}{
			\code{colData(x)}:
			Returns a \code{DataFrame} with sample identifiers and any
			information in the 'sample.annotation' node.
		}
		\item{}{
			\code{info(x, info=NULL)}:
			Returns the info fields as a \code{DataFrame}. \code{info} is a
			character vector with the names of fields to return (default is to
			return all).
		}
		\item{}{
			\code{geno(x, geno=NULL)}:
			Returns the geno (format) fields as a \code{SimpleList}. \code{geno}
			is a character vector with the names of fields to return (default
			is to return all).
		}
	}

	Other data can be accessed with \code{\link{seqGetData}}.
}

\section{Coercion methods}{
	In the following code snippets \code{x} is a SeqVarGDSClass object.
	\describe{
		\item{}{
			\code{seqAsVCF(x, chr.prefix="", info=NULL, geno=NULL)}:
			
		}
	}
}


\author{Stephanie Gogarten, Xiuwen Zheng}

\seealso{
    \code{\link{gds.class}}, \code{\link{seqOpen}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
gds

## sample ID
head(seqGetData(gds, "sample.id"))

## variants
granges(gds)

\dontrun{
## alleles as comma-separated character strings
head(seqGetData(gds, "allele"))

## alleles as DNAStringSet or DNAStringSetList
ref(gds)
v <- alt(gds)

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
}

seqClose(gds)
}
