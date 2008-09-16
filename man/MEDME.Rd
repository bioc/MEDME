\name{MEDME}
\alias{MEDME}
\title{ Determining the logistic model of MeDIP enrichment in respect to the expected DNA methylation level }
\description{
  Probe-level MeDIP weighted enrichment is compared to the expected DNA methytlation level. The former is determined applying MeDIP protocol to a fully methylated DNA. The latter is determined as the count of CpGs for each probe. This is assumed to be the methylation level of each probe in a fully methylated sample.
}
\usage{
MEDME(data, sample, CGcountThr = 1, figName = NULL)
}
\arguments{
  \item{data}{ An object of class MEDMEset}
  \item{sample}{ Integer; the number of the sample to be used to fit the model, based on the order of samples in the smoothed slot}
  \item{CGcountThr}{ number; the threshold to avoid modelling probes with really low methylation level, i.e. CpG count }
  \item{figName}{ string; the name of the file reporting the model fitting }
}
\details{The model should be applied on calibration data containing MeDIP enrichment of fully methylated DNA, most likely artificially generated (see references). Nevertheless, in case chromosome or genome-wide human tiling arrays are used a regular sample could be used too. In fact, human genomic DNA is known to be hyper-methylated but in the promoter regions. Of course the performance of the method is expected to be somehow affected by this approximation.}
\value{
The logistic model as returned from the multdrc function from the drc R library
}
\seealso{
\code{\link{smooth}}, \code{\link{CGcount}}
}
\examples{
data(testMEDMEset)
## just an example with the first 1000 probes
testMEDMEset = smooth(data = testMEDMEset[1:1000, ])
testMEDMEset = CGcount(data = testMEDMEset)
MEDMEmodel = MEDME(data = testMEDMEset, sample = 1, CGcountThr = 1, figName = NULL)
}
\references{\url{http://genome.cshlp.org/cgi/content/abstract/gr.080721.108v1}}
