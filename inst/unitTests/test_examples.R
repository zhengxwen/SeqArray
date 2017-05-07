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
	if (Sys.info()[['sysname']] == "Windows")
		return(invisible())

	function.list <- readRDS(
		system.file("Meta", "Rd.rds", package="SeqArray"))$Name

	sapply(function.list, FUN = function(func.name)
		{
			args <- list(
				topic=func.name,
				package="SeqArray",
				echo=FALSE, verbose=FALSE, ask=FALSE
			)
			message("Running the examples in '", func.name, "()':")
			suppressWarnings(do.call(example, args))
			NULL
		})
	invisible()
}
