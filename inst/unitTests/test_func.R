
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
	# open the GDS file
	gds.fn <- seqExampleFileName("gds")
	f <- seqOpen(gds.fn)

	samp.id <- seqGetData(f, "sample.id")
	variant.id <- seqGetData(f, "variant.id")

	# get results
	for (p in 1:4)
	{
		d <- seqAlleleFreq(f, NULL, parallel=p)
		checkIdentical(Valid$fcAlleleFreq$d1, d, paste0("seqAlleleFreq 1:", i))
	}

	for (p in 1:4)
	{
		d <- seqAlleleFreq(f, 0L, parallel=p)
		checkIdentical(Valid$fcAlleleFreq$d2, d, paste0("seqAlleleFreq 2:", i))
	}

	for (p in 1:4)
	{
		set.seed(1000)
		d <- seqAlleleFreq(f, sample(c(0L,1L), length(variant.id),
			replace=TRUE), parallel=p)
		checkIdentical(Valid$fcAlleleFreq$d3, d, paste0("seqAlleleFreq 3:", i))
	}

	for (p in 1:4)
	{
		d <- seqAlleleFreq(f, toupper(seqGetData(f, "annotation/info/AA")$data),
			parallel=p)
		checkIdentical(Valid$fcAlleleFreq$d4, d, paste0("seqAlleleFreq 4:", i))
	}

	# close the GDS file
	seqClose(f)
	invisible()
}
