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
	on.exit(seqClose(f))

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

	invisible()
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
	on.exit(seqClose(f))

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
		# 'Rounding' was the default in versions prior to R_3.6.0
		# it is used for reproduction of the results created by R (v3.5.2)
		tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
			error=function(e) FALSE)
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

	invisible()
}


test_random_genotype <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

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

	invisible()
}


test_dosage <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

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

	invisible()
}


test_random_dosage <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

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

	invisible()
}


test_random_phase <- function()
{
	# open the GDS file
	file.copy(seqExampleFileName("gds"), "test.gds", overwrite=TRUE)

	f <- openfn.gds("test.gds", FALSE)
	m <- read.gdsn(index.gdsn(f, "phase/data"))
	s <- sample.int(2, length(m), TRUE) - 1L
	dim(s) <- dim(m)
	add.gdsn(index.gdsn(f, "phase"), "data", s, replace=TRUE)
	closefn.gds(f)

	f <- seqOpen("test.gds")
	seqResetFilter(f)
	gm <- seqGetData(f, "phase")

	set.seed(300)
	for (i in 1:10)
	{
		x <- sample.int(dim(gm)[2L], round(dim(gm)[2L]/3))
		y <- sample.int(dim(gm)[1L], round(dim(gm)[1L]/3))
		seqSetFilter(f, variant.sel=x, sample.sel=y, verbose=FALSE)
		m1 <- seqGetData(f, "phase")

		f1 <- rep(FALSE, dim(gm)[2L]); f1[x] <- TRUE
		f2 <- rep(FALSE, dim(gm)[1L]); f2[y] <- TRUE
		m2 <- gm[f2, f1]
		checkEquals(m1, m2, "phasing information: random access (1)")

		v <- seqApply(f, "phase", function(x) x, as.is="list")
		m3 <- matrix(unlist(v), nrow=length(v[[1L]]))
		checkEquals(m1, m3, "phasing information: random access (2)")
	}

	# close the GDS file
	seqClose(f)
	unlink("test.gds", force=TRUE)

	invisible()
}


test_random_info <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

	set.seed(200)
	num <- seqSummary(f, "genotype", check="none", verbose=FALSE)$dim[3L]

	for (nm in ls.gdsn(index.gdsn(f, "annotation/info")))
	{
		seqResetFilter(f, verbose=FALSE)
		dat <- seqGetData(f, paste0("annotation/info/", nm))
		if (is.list(dat)) dat <- dat$data
		dimnames(dat) <- NULL

		for (i in 1:5)
		{
			idx <- sample.int(num, num)
			rv <- vector("list", length(idx))
			for (k in idx)
			{
				seqSetFilter(f, variant.sel=k, verbose=FALSE)
				d <- unlist(seqApply(f, paste0("annotation/info/", nm),
					function(x) x, as.is="list"), recursive=FALSE)
				rv[[k]] <- d
			}

			m <- unlist(rv)
			checkEquals(dat, m, sprintf("INFO (%s): random access", nm))
		}
	}

	invisible()
}


test_random_format <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

	seqResetFilter(f)
	dat <- seqGetData(f, "annotation/format/DP")
	y <- dat$data
	dimnames(y) <- NULL

	set.seed(200)
	for (i in 1:5)
	{
		idx <- sample.int(length(dat$length), length(dat$length))
		rv <- vector("list", length(idx))
		for (k in idx)
		{
			seqSetFilter(f, variant.sel=k, verbose=FALSE)
			d <- unlist(seqApply(f, "annotation/format/DP", function(x) x,
				as.is="list"), recursive=FALSE)
			rv[[k]] <- d
		}

		m <- matrix(unlist(rv), nrow=length(d))
		checkEquals(y, m, "FORMAT: random access")
	}

	invisible()
}


test.apply_vs_blockapply <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

	n <- seqSummary(f, "genotype", verbose=FALSE)$dim[3L]

	# randoms variant set
	set.seed(1000)
	seqSetFilter(f, variant.sel=sample.int(n, 2/3*n))

	# genotype
	v1 <- seqApply(f, "genotype", function(x) mean(x, na.rm=TRUE),
		as.is="double")
	v2 <- seqBlockApply(f, "genotype", function(x)
		colMeans(x, na.rm=TRUE, dims=2L), as.is="unlist", bsize=128L)
	checkEquals(v1, v2, "Apply vs BlockApply: genotype")

	# phase
	v1 <- seqApply(f, "phase", function(x) mean(x, na.rm=TRUE),
		as.is="double")
	v2 <- seqBlockApply(f, "phase", function(x)
		colMeans(x, na.rm=TRUE), as.is="unlist", bsize=128L)
	checkEquals(v1, v2, "Apply vs BlockApply: phase")

	# annotation/info/AC
	v1 <- seqApply(f, "annotation/info/AC", function(x) x, as.is="double")
	v2 <- seqBlockApply(f, "annotation/info/AC", function(x) x,
		as.is="unlist", bsize=128L)
	checkEquals(v1, v2, "Apply vs BlockApply: AC")

	# annotation/info/BN
	v1 <- seqApply(f, "annotation/info/BN", function(x) x, as.is="double")
	v2 <- seqBlockApply(f, "annotation/info/BN", function(x) x,
		as.is="unlist", bsize=128L)
	checkEquals(v1, v2, "Apply vs BlockApply: BN")

	# annotation/format/DP
	v1 <- seqApply(f, "annotation/format/DP", function(x) mean(x, na.rm=TRUE),
		as.is="double")
	v2 <- seqBlockApply(f, "annotation/format/DP", function(x)
		colMeans(x$data, na.rm=TRUE), as.is="unlist", bsize=128L)
	checkEquals(v1, v2, "Apply vs BlockApply: DP")

	invisible()
}


test.dosage_alt <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

	# check integer
	g1 <- seqGetData(f, "$dosage")
	g2 <- seqGetData(f, "$dosage_alt")
	checkEquals(is.na(g1), is.na(g2), "GetData (int genotype): missing genotypes")
	g <- g1 + g2
	checkEquals(unique(c(g)), c(NA, 2L), "GetData (int genotype): sum of dosage and dosage_alt")

	# check RAW
	g1 <- seqGetData(f, "$dosage", .useraw=TRUE)
	g2 <- seqGetData(f, "$dosage_alt", .useraw=TRUE)
	checkEquals(g1==0xFF, g2==0xFF, "GetData (RAW genotype): missing genotypes")
	g <- as.integer(g1) + as.integer(g2)
	checkEquals(unique(g), c(510L, 2L), "GetData (RAW genotype): sum of dosage and dosage_alt")

	invisible()
}


test.parallel_balancing <- function()
{
	# open the GDS file
	f <- seqOpen(seqExampleFileName("gds"))
	on.exit(seqClose(f))

	# test 1
	p0 <- seqGetData(f, "position")

	p1 <- seqParallel(2, f, function(gds) seqGetData(gds, "position"))
	checkEquals(p0, p1, "Parallel load balancing: test 1 (1)")
	p2 <- seqParallel(2, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p2, "Parallel load balancing: test 1 (2)")
	p3 <- seqParallel(1, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p3, "Parallel load balancing: test 1 (3)")
	p4 <- seqParallel(2, f, function(gds, flag) p0[flag],
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE, .selection.flag=TRUE)
	checkEquals(p0, p4, "Parallel load balancing: test 1 (4)")

	cl <- makeCluster(2)
	p5 <- seqParallel(cl, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p5, "Parallel load balancing: test 1 (5)")
	p6 <- seqParallel(cl, f, function(gds, flag) p0[flag],
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE, .selection.flag=TRUE)
	checkEquals(p0, p6, "Parallel load balancing: test 1 (6)")
	stopCluster(cl)

	# test 2
	set.seed(1000)
	sel <- sample(c(FALSE, TRUE), length(p0), replace=TRUE)
	p0 <- p0[sel]
	seqSetFilter(f, variant.sel=sel)

	p1 <- seqParallel(2, f, function(gds) seqGetData(gds, "position"))
	checkEquals(p0, p1, "Parallel load balancing: test 2 (1)")
	p2 <- seqParallel(2, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p2, "Parallel load balancing: test 2 (2)")
	p3 <- seqParallel(1, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p3, "Parallel load balancing: test 2 (3)")
	p4 <- seqParallel(2, f, function(gds, flag) p0[flag],
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE, .selection.flag=TRUE)
	checkEquals(p0, p4, "Parallel load balancing: test 2 (4)")

	cl <- makeCluster(2)
	p5 <- seqParallel(cl, f, function(gds) seqGetData(gds, "position"),
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE)
	checkEquals(p0, p5, "Parallel load balancing: test 2 (5)")
	p6 <- seqParallel(cl, f, function(gds, flag) p0[flag],
		.balancing=TRUE, .bl_size=5, .bl_progress=TRUE, .selection.flag=TRUE)
	checkEquals(p0, p6, "Parallel load balancing: test 2 (6)")
	stopCluster(cl)

	invisible()
}
