###########################################################################
#
# Test the parallel functions with a BiocParallelParam object
#

library(SeqArray)
library(RUnit)


# SnowParam started by the user, and MulticoreParam not started
.bp_list <- function()
{
	list(snow = BiocParallel::bpstart(BiocParallel::SnowParam(2)),
		mcore = BiocParallel::MulticoreParam(2))
}

.bp_stop <- function(lst)
{
	for (bp in lst)
		if (BiocParallel::bpisup(bp)) BiocParallel::bpstop(bp)
}


test.biocparallel_seqParallel <- function()
{
	if (!requireNamespace("BiocParallel", quietly=TRUE)) return(invisible())
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))
	lst <- .bp_list()
	on.exit(.bp_stop(lst), add=TRUE)

	vid <- seqGetData(f, "variant.id")
	sid <- seqGetData(f, "sample.id")
	for (nm in names(lst))
	{
		bp <- lst[[nm]]
		msg <- function(s) paste0(nm, ": ", s)
		checkIdentical(vid, seqParallel(bp, f, function(g)
			seqGetData(g, "variant.id")), msg("by.variant"))
		checkIdentical(sid, seqParallel(bp, f, function(g)
			seqGetData(g, "sample.id"), split="by.sample"), msg("by.sample"))
		checkIdentical(length(vid), seqParallel(bp, f, function(g, s) sum(s),
			.selection.flag=TRUE, .combine="+"), msg(".selection.flag"))
		checkIdentical(vid, seqParallel(bp, f, function(g)
			seqGetData(g, "variant.id"), .balancing=TRUE, .bl_size=100L),
			msg(".balancing=TRUE"))
		checkIdentical(length(vid), seqParallel(bp, f, function(g)
			length(seqGetData(g, "variant.id")), .combine="+",
			.balancing=TRUE, .bl_size=100L), msg(".balancing=TRUE, '+'"))
		checkIdentical(c(2L, 3L, 5L, 6L, 8L, 9L), seqParallel(bp, 9,
			function(i) if (i %% 3L) i + 1L, split="none"),
			msg("numeric gdsfile"))
		checkIdentical(list(1L, NULL, 3L), seqParallel(bp, 3,
			function(i) if (i != 2L) i, split="none", .combine="list"),
			msg("numeric gdsfile, 'list'"))
		checkIdentical(1:2, sort(seqParallel(bp, NULL,
			function() SeqArray:::process_index, split="none")),
			msg("NULL gdsfile"))
		# the arguments named like the internal ones
		checkIdentical(rep(list(1:5), 2L), seqParallel(bp, f,
			function(g, X, .FUN, .fun, .args, BPPARAM) c(X, .FUN, .fun, .args,
			BPPARAM), X=1L, .FUN=2L, .fun=3L, .args=4L, BPPARAM=5L,
			.combine="list"), msg("argument names"))
		# no more jobs than the selected variants
		seqSetFilter(f, variant.sel=3L, verbose=FALSE)
		checkIdentical(list(1L), seqParallel(bp, f, function(g)
			length(seqGetData(g, "variant.id")), .combine="list"),
			msg("1 variant selected"))
		seqResetFilter(f, verbose=FALSE)
	}

	invisible()
}


test.biocparallel_process <- function()
{
	if (!requireNamespace("BiocParallel", quietly=TRUE)) return(invisible())
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))
	lst <- .bp_list()
	on.exit(.bp_stop(lst), add=TRUE)

	pv <- function(...) c(SeqArray:::process_index, SeqArray:::process_count,
		SeqArray:::process_block_index, SeqArray:::process_block_count)
	for (nm in names(lst))
	{
		bp <- lst[[nm]]
		msg <- function(s) paste0(nm, ": ", s)
		# load balancing, 1348 variants in 27 blocks
		m <- do.call(rbind, seqParallel(bp, f, pv, .balancing=TRUE,
			.bl_size=50L, .combine="list"))
		checkTrue(all(m[,1L] %in% 1:2), msg("balancing, process_index"))
		checkTrue(all(m[,2L] == 2L), msg("balancing, process_count"))
		checkIdentical(1:27, m[,3L], msg("balancing, process_block_index"))
		checkTrue(all(m[,4L] == 27L), msg("balancing, process_block_count"))
		# numeric gdsfile
		m <- do.call(rbind, seqParallel(bp, 5, pv, split="none",
			.combine="list"))
		checkTrue(all(m[,1L] %in% 1:2), msg("numeric, process_index"))
		checkIdentical(1:5, m[,3L], msg("numeric, process_block_index"))
		checkTrue(all(m[,4L] == 5L), msg("numeric, process_block_count"))

		# .initialize and .finalize are called once in each worker
		d <- tempfile()
		dir.create(d)
		v <- seqParallel(bp, f, function(g) c(Sys.getpid(),
			get0("BP_INIT", envir=globalenv())),
			.initialize = function(i, p) assign("BP_INIT", i, envir=globalenv()),
			.finalize = function(i, p)
				file.create(file.path(p, paste(Sys.getpid(), i, sep="_"))),
			.initparam=d, .balancing=TRUE, .bl_size=50L, .combine="list")
		m <- unique(do.call(rbind, v))
		fn <- list.files(d)
		unlink(d, recursive=TRUE)
		checkTrue(all(m[,2L] %in% 1:2), msg(".initialize, process id"))
		checkIdentical(nrow(m), length(unique(m[,1L])),
			msg(".initialize, once in each worker"))
		checkIdentical(sort(paste(m[,1L], m[,2L], sep="_")), sort(fn),
			msg(".finalize, once in each initialized worker"))
	}

	# the workers not started are stopped after the call
	checkTrue(!BiocParallel::bpisup(lst$mcore), "bpisup()")
	checkTrue(BiocParallel::bpisup(lst$snow), "bpisup()")

	invisible()
}


test.biocparallel_functions <- function()
{
	if (!requireNamespace("BiocParallel", quietly=TRUE)) return(invisible())
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))
	lst <- .bp_list()
	on.exit(.bp_stop(lst), add=TRUE)

	gsum <- function(g) sum(g, na.rm=TRUE)
	# the sum per variant, independent of the blocks
	gsum2 <- function(g) colSums(g, dims=2L, na.rm=TRUE)
	units <- seqUnitSlidingWindows(f, win.size=50000L, win.shift=25000L,
		verbose=FALSE)
	v1 <- seqApply(f, "genotype", gsum, as.is="integer")
	v2 <- seqBlockApply(f, "genotype", gsum2, as.is="unlist", bsize=100L)
	v3 <- seqUnitApply(f, units, "genotype", gsum, as.is="unlist")
	v4 <- seqAlleleFreq(f, parallel=FALSE)
	v5 <- seqMissing(f, parallel=FALSE)
	v6 <- seqGet2bGeno(f, parallel=FALSE)
	for (nm in names(lst))
	{
		bp <- lst[[nm]]
		msg <- function(s) paste0(nm, ": ", s)
		checkIdentical(v1, seqApply(f, "genotype", gsum, as.is="integer",
			parallel=bp), msg("seqApply()"))
		checkIdentical(v2, seqBlockApply(f, "genotype", gsum2, as.is="unlist",
			bsize=100L, parallel=bp), msg("seqBlockApply()"))
		checkIdentical(v3, seqUnitApply(f, units, "genotype", gsum,
			as.is="unlist", parallel=bp), msg("seqUnitApply()"))
		checkIdentical(v4, seqAlleleFreq(f, parallel=bp), msg("seqAlleleFreq()"))
		checkIdentical(v5, seqMissing(f, parallel=bp), msg("seqMissing()"))
		checkIdentical(v6, seqGet2bGeno(f, parallel=bp), msg("seqGet2bGeno()"))
	}

	invisible()
}


test.biocparallel_error <- function()
{
	if (!requireNamespace("BiocParallel", quietly=TRUE)) return(invisible())
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))
	lst <- .bp_list()
	on.exit(.bp_stop(lst), add=TRUE)

	vid <- seqGetData(f, "variant.id")
	bp <- lst$snow
	# no result is left to the next call after an error
	checkException(seqParallel(bp, f, function(g) stop("FUN error"),
		.balancing=TRUE, .bl_size=100L), silent=TRUE)
	checkIdentical(vid, seqParallel(bp, f, function(g)
		seqGetData(g, "variant.id")), "after an error in FUN")
	checkException(seqParallel(bp, 6, function(i) { Sys.sleep(0.1*i); i },
		split="none", .combine=function(x, y) stop(".combine error")),
		silent=TRUE)
	checkIdentical(vid, seqParallel(bp, f, function(g)
		seqGetData(g, "variant.id")), "after an error in .combine")
	# no state is left in the workers
	v <- BiocParallel::bplapply(1:2, function(i)
		is.null(SeqArray:::.PkgEnv$bp_worker), BPPARAM=bp)
	checkTrue(all(unlist(v)), "no state left in the workers")

	invisible()
}


test.biocparallel_setup <- function()
{
	if (!requireNamespace("BiocParallel", quietly=TRUE)) return(invisible())
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit({ seqParallelSetup(FALSE, verbose=FALSE); seqClose(f) })

	v <- seqAlleleFreq(f, parallel=FALSE)
	bp <- BiocParallel::bpstart(BiocParallel::SnowParam(2))
	seqParallelSetup(bp, verbose=FALSE)
	checkIdentical(bp, seqGetParallel(), "seqGetParallel()")
	checkIdentical(v, seqAlleleFreq(f), "seqAlleleFreq() with seqGetParallel()")
	seqParallelSetup(FALSE, verbose=FALSE)
	checkIdentical(FALSE, seqGetParallel(), "seqParallelSetup(FALSE)")
	checkTrue(!BiocParallel::bpisup(bp), "seqParallelSetup(FALSE), bpisup()")

	# remove the setting of a number
	seqParallelSetup(2L, verbose=FALSE)
	seqParallelSetup(FALSE, verbose=FALSE)
	checkIdentical(FALSE, seqGetParallel(), "seqParallelSetup(2L), then FALSE")

	invisible()
}
