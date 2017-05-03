.test_variableLengthToList <- function() {
  x <- list(length=c(1,1,2,1,3,1),
            data=c(1,1,1:2,1,1:3,1))
  checkEquals(list(1,1,1:2,1,1:3,1), .variableLengthToList(x))
  x <- list(length=rep(1,5),
            data=1:5)
  checkEquals(list(1,2,3,4,5), .variableLengthToList(x))
}

.test_variableLengthToMatrix <- function() {
  x <- list(matrix(1:3, nrow=2, ncol=3, byrow=TRUE),
            matrix(4:6, nrow=2, ncol=3, byrow=TRUE),
            matrix(7:9, nrow=2, ncol=3, byrow=TRUE))
  xm <- matrix(list(1:3,1:3,4:6,4:6,7:9,7:9), nrow=2, ncol=3)
  checkIdentical(xm, .variableLengthToMatrix(x))
}
