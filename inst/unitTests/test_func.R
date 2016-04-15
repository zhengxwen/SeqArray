###########################################################################
#
# Test SeqArray functions
#

library(SeqArray)
library(RUnit)


.create_valid_data <- function()
{
	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	set.seed(1000)
	samp.id <- seqGetData(f, "sample.id")
	variant.id <- seqGetData(f, "variant.id")

	# get results
	fcAlleleFreq <- list(
		d1=seqAlleleFreq(f, NULL),
		d2=seqAlleleFreq(f, 0L),
		d3=seqAlleleFreq(f, sample(c(0L,1L), length(variant.id), replace=TRUE)),
		d4=seqAlleleFreq(f, toupper(seqGetData(f, "annotation/info/AA")$data))
	)

	# validation data
	Valid <- list(
		fcAlleleFreq = fcAlleleFreq
	)
	save(Valid, file="Valid.RData", compress="xz")

	# close the GDS file
	seqClose(f)
}


# load the validation data
Valid <- get(load(system.file("unitTests", "data", "Valid.RData",
	package="SeqArray", mustWork=TRUE)))


test_allele_freq <- function()
{
	num.cores <- 2L

	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	samp.id <- seqGetData(f, "sample.id")
	variant.id <- seqGetData(f, "variant.id")

	# get results
	for (p in 1L:num.cores)
	{
		d <- seqAlleleFreq(f, NULL, parallel=p)
		checkEquals(Valid$fcAlleleFreq$d1, d, paste0("seqAlleleFreq 1:", p))
	}

	for (p in 1L:num.cores)
	{
		d <- seqAlleleFreq(f, 0L, parallel=p)
		checkEquals(Valid$fcAlleleFreq$d2, d, paste0("seqAlleleFreq 2:", p))
	}

	for (p in 1L:num.cores)
	{
		set.seed(1000)
		d <- seqAlleleFreq(f, sample(c(0L,1L), length(variant.id),
			replace=TRUE), parallel=p)
		checkEquals(Valid$fcAlleleFreq$d3, d, paste0("seqAlleleFreq 3:", p))
	}

	for (p in 1L:num.cores)
	{
		d <- seqAlleleFreq(f, toupper(seqGetData(f, "annotation/info/AA")$data),
			parallel=p)
		checkEquals(Valid$fcAlleleFreq$d4, d, paste0("seqAlleleFreq 4:", p))
	}

	# close the GDS file
	seqClose(f)
	invisible()
}


test_random_genotype <- function()
{
	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	seqResetFilter(f)
	gm <- seqGetData(f, "genotype")

	set.seed(100)
	for (i in 1:10)
	{
		x <- sample.int(dim(gm)[3L], round(dim(gm)[3L]/3))
		y <- sample.int(dim(gm)[2L], round(dim(gm)[2L]/3))
		seqSetFilter(f, variant.sel=x, sample.sel=y, verbose=FALSE)
		m1 <- seqGetData(f, "genotype")

		f1 <- rep(FALSE, dim(gm)[3L]); f1[x] <- TRUE
		f2 <- rep(FALSE, dim(gm)[2L]); f2[y] <- TRUE
		m2 <- gm[, f2, f1]
		checkEquals(m1, m2, "genotype: random access")
	}

	# close the GDS file
	seqClose(f)
	invisible()
}


test_dosage <- function()
{
	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	gm <- seqGetData(f, "genotype")
	dm <- seqGetData(f, "$dosage")
	dimnames(dm) <- NULL
	m <- (gm[1L,,]==0L) + (gm[2L,,]==0L)
	checkEquals(dm, m, "dosage, integer")

	gm <- seqGetData(f, "genotype", .useraw=TRUE)
	dm <- seqGetData(f, "$dosage", .useraw=TRUE)
	dimnames(dm) <- NULL
	m <- (gm[1L,,]==0L) + (gm[2L,,]==0L)
	m[gm[1L,,]==255L | gm[2L,,]==255L] <- 255L
	m <- as.raw(m)
	dim(m) <- dim(dm)
	checkEquals(dm, m, "dosage, raw")

	# close the GDS file
	seqClose(f)
	invisible()
}


test_random_dosage <- function()
{
	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	seqResetFilter(f)
	gm <- seqGetData(f, "$dosage")

	set.seed(200)
	for (i in 1:10)
	{
		x <- sample.int(dim(gm)[2L], round(dim(gm)[2L]/3))
		y <- sample.int(dim(gm)[1L], round(dim(gm)[1L]/3))
		seqSetFilter(f, variant.sel=x, sample.sel=y, verbose=FALSE)
		m1 <- seqGetData(f, "$dosage")

		f1 <- rep(FALSE, dim(gm)[2L]); f1[x] <- TRUE
		f2 <- rep(FALSE, dim(gm)[1L]); f2[y] <- TRUE
		m2 <- gm[f2, f1]
		checkEquals(m1, m2, "dosage: random access")
	}

	# close the GDS file
	seqClose(f)
	invisible()
}
