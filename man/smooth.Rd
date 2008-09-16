\name{smooth}
\alias{smooth}
\title{ Determining weighted MeDIP data }
\description{
  MeDIP data from tiling arrays are smoothed by determining for each probe i the weighted average of the probes within a window of size wsize centered at i
}
\usage{
smooth(data, wsize=1000, wFunction='linear')
}
\arguments{
  \item{data}{An object of class MEDMEset}
  \item{wsize}{ number; the size of the smoothing window, in bp }
  \item{wFunction}{string; the type of weighting function, to choose among linear, exp, log or none }
}
\details{The un-smoothed data are read from the slot logR of the data MEDMEset and the resulting smoothed data are saved on the smoothed slot.}
\value{
An object of class MEDMEset. In particular, the smoothed data are saved on the smoothed slot.
}
\examples{
data(testMEDMEset)
# just an example with the first 1000 probes
testMEDMEset = smooth(data = testMEDMEset[1:1000,])
}
