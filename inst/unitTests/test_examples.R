#############################################################
#
# DESCRIPTION: run all examples in SeqArray
#

library(SeqArray)


#############################################################
#
# test functions
#

function.list <- c(
	"SeqArray-package",
	"SeqVarGDSClass-class",
	"seqAlleleFreq",
	"seqApply",
	"seqClose-methods",
	"seqCompress.Option",
	"seqDelete",
	"seqExampleFileName",
	"seqExport",
	"seqGDS2SNP",
	"seqGDS2VCF",
	"seqGetData",
	"seqGetFilter",
	"seqInfoNewVar",
	"seqMerge",
	"seqMissing",
	"seqNumAllele",
	"seqOpen",
	"seqOptimize",
	"seqParallel",
	"seqSetFilter",
	"seqSetFilterChrom",
	"seqSlidingWindow",
	"seqSummary",
	"seqTranspose",
	"seqVCF.Header",
	"seqVCF.SampID",
	"seqVCF2GDS"
)


test_examples <- function()
{
	sapply(function.list, FUN = function(func.name)
		{
			args <- list(
				topic=func.name,
				package="SeqArray",
				echo=FALSE, verbose=FALSE, ask=FALSE
			)
			suppressWarnings(do.call(example, args))
			NULL
		})
	invisible()
}
