\name{seqOpen}
\alias{seqOpen}
\title{Open a Sequencing GDS File}
\description{
	Open a Sequencing GDS file
}
\usage{
seqOpen(gds.fn, readonly=TRUE)
}
\arguments{
	\item{gds.fn}{the file name}
	\item{readonly}{whether read-only or not}
}
\value{
	Return an object of class \code{\link[gdsfmt]{gds.class}}.
}
\details{
	It is strongly suggested to call \code{seqOpen} instead of
\code{\link[gdsfmt]{openfn.gds}}, since \code{seqOpen} will perform internal checking
for data integrality.
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# open the GDS file
gdsfile <- seqOpen(gds.fn)

# display the contents of the GDS file in a hierarchical structure
gdsfile

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
