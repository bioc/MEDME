\name{CGcount}
\alias{CGcount}
\title{ Determining the count of CpG dinucleotides for a set of genomic locations}
\description{
  The count of CpGs is determined in each window of size wsize, with or withouth weighting, for each probe according to its position, chromosome and genome realease
}
\usage{
CGcount(data, wsize = 1000, wFunction = "linear")
}
\arguments{
  \item{data}{ An object of class MEDMEset}
  \item{wsize}{ number; the size of the smoothing window, in bp }
  \item{wFunction}{string; the type of weighting function, to choose among linear, exp, log or none }
}
\details{Only human and mouse are currently supported. The respective genomic sequence metadata library needs to be downloaded from the Bioconductor website, installed and loaded (around 800Mb). Please
note that only the last genome release should be used. LiftOver UCSC tool could be used for batch conversion of old genomic position to the last genome release.}
\value{
An object of class MEDMEset is returned where the count of CpGs for each probe has been saved on the CGcount slot.
}
\seealso{
\code{\link{smooth}}
}
\examples{
data(testMEDMEset)
## just an example with the first 1000 probes
testMEDMEset = smooth(data = testMEDMEset[1:1000,])
library(BSgenome.Hsapiens.UCSC.hg18)
testMEDMEset = CGcount(data = testMEDMEset)
}
