library(SeqArray)
library(RUnit)


.popcnt <- function(x)
{
	ans <- 0L
	while(x > 0)
	{
		if (x %% 2) ans <- ans + 1L
		x <- x %/% 2
	}
	ans
}


test_popcnt <- function()
{
	set.seed(1000)
	v1 <- sample.int(2147483647L, 10000L, replace=TRUE)
	v2 <- sample.int(2147483647L, 10000L, replace=TRUE)

	c1 <- sapply(v1, .popcnt)
	c2 <- sapply(v2, .popcnt)

	t1 <- SeqArray:::.cfunction("test_array_popcnt32")(v1)
	t2 <- SeqArray:::.cfunction("test_array_popcnt32")(v2)
	t3 <- SeqArray:::.cfunction2("test_array_popcnt64")(v1, v2)

	checkEquals(c1, t1, "popcount u32")
	checkEquals(c2, t2, "popcount u32")
	checkEquals(c1+c2, t3, "popcount u64")

	invisible()
}


test_byte_count <- function()
{
	set.seed(1000)
	for (st in sample.int(1000L, 100L))
	{
		n <- 50000L + sample.int(64L, 1L) - 1L
		v <- sample.int(255L, n, replace=TRUE)
		v[sample.int(n, 25000L)] <- 0L
		v <- as.raw(v)

		n1 <- SeqArray:::.cfunction2("test_byte_count")(v, st)
		n2 <- sum(v[st:length(v)] != 0L)
		checkEquals(n1, n2, paste0("byte_count (start=", st, ")"))
	}

	invisible()
}


test_int_count <- function()
{
	set.seed(5000)
	for (st in sample.int(1000L, 100L))
	{
		n <- 50000L + sample.int(64L, 1L) - 1L
		v <- sample.int(256L, n, replace=TRUE)
		fd <- sample(v, 1L)

		n1 <- SeqArray:::.cfunction3("test_int32_count")(v, st, fd)
		n2 <- sum(v[st:length(v)] == fd)
		checkEquals(n1, n2, paste0("int_count (start=", st, ")"))
	}

	invisible()
}


test_int_count2 <- function()
{
	set.seed(5000)
	for (st in sample.int(1000L, 100L))
	{
		n <- 50000L + sample.int(64L, 1L) - 1L
		v <- sample.int(256L, n, replace=TRUE)
		fd1 <- sample(v, 1L)
		fd2 <- sample(v, 1L)

		n1 <- SeqArray:::.cfunction4("test_int32_count2")(v, st, fd1, fd2)
		n2 <- c(sum(v[st:length(v)] == fd1), sum(v[st:length(v)] == fd2))
		checkEquals(n1, n2, paste0("int_count2 (start=", st, ")"))
	}

	invisible()
}


test_int_replace <- function()
{
	set.seed(5000)
	for (st in sample.int(1000L, 100L))
	{
		n <- 50000L + sample.int(64L, 1L) - 1L
		v1 <- sample.int(256L, n, replace=TRUE)
		fd <- sample(v1, 1L)
		sub <- sample(v1, 1L)

		v2 <- SeqArray:::.cfunction4("test_int32_replace")(v1, st, fd, sub)
		a <- v1[st:length(v1)]
		a[a == fd] <- sub
		v1[st:length(v1)] <- a
		checkEquals(v1, v2, paste0("int_replace (start=", st, ")"))
	}

	invisible()
}
