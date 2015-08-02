#############################################################
#
# DESCRIPTION: run all examples in SeqArray
#

library(SeqArray)


#############################################################
#
# test functions
#

test_examples <- function()
{
	function.list <- readRDS(
		system.file("Meta", "Rd.rds", package="SeqArray"))$Name

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
