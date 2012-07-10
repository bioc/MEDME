\name{MEDME.predict}
\alias{MEDME.predict}
\title{ Applying the logistic model on MeDIP enrichment data }
\description{
  This allows the probe-level determination of MeDIP smoothed data, as well as absolute and relative methylation levels (AMS and RMS respectively)
}
\usage{
MEDME.predict(data, MEDMEfit, MEDMEextremes = c(1,32), wsize = 1000, wFunction='linear')
}
\arguments{
  \item{data}{ An object of class MEDMEset}
  \item{MEDMEfit}{ the model obtained from the MEDME.model function }
  \item{MEDMEextremes}{ vector; the background and saturation values as determined by the fitting of the model on the calibration data}
  \item{wsize}{ number; the size of the smoothing window, in bp }
  \item{wFunction}{ string; the type of weighting function, to choose among linear, exp, log or none }
}
\value{
An object of class MEDMEset. The resulting smoothed data, the absolute  and relative methylation score (AMS and RMS) are saved in the smoothed, AMS and RMS slots, respectively.
}
\seealso{
\code{\link{smooth}}, \code{\link{CGcount}}, \code{\link{MEDME}}
}
\examples{
data(testMEDMEset)
## just an example with the first 1000 probes
testMEDMEset = smooth(data = testMEDMEset[1:1000, ])
library(BSgenome.Hsapiens.UCSC.hg18)
testMEDMEset = CGcount(data = testMEDMEset)
MEDMEmodel = MEDME(data = testMEDMEset, sample = 1, CGcountThr = 1, figName = NULL)
testMEDMEset = MEDME.predict(data = testMEDMEset, MEDMEfit = MEDMEmodel, MEDMEextremes = c(1,32), wsize = 1000, wFunction='linear')
}
\references{\url{http://genome.cshlp.org/cgi/content/abstract/gr.080721.108v1}}
