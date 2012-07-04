\name{MEDMEset-class}
\docType{class}
\alias{MEDMEset-class}
\alias{[,MEDMEset-method}
\alias{AMS,MEDMEset-method}
\alias{CG,MEDMEset-method}
\alias{chr,MEDMEset-method}
\alias{org,MEDMEset-method}
\alias{initialize,MEDMEset-method}
\alias{logR,MEDMEset-method}
\alias{pos,MEDMEset-method}
\alias{RMS,MEDMEset-method}
\alias{show,MEDMEset-method}
\alias{smoothed,MEDMEset-method}
\alias{AMS}
\alias{CG}
\alias{CGcounts}
\alias{chr}
\alias{org}
\alias{organism}
\alias{logR}
\alias{pos}
\alias{RMS}
\alias{smoothed}

\title{Class "MEDMEset"}
\description{ This class is used in MEDME library to store MeDIP derived DNA-methylation estimates and to save further elaboration of these, in association with chromosomal and positional probe information}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MEDMEset", ...)}.
This object could initially host the MeDIP normalized logRatio data, as returned by the \code{\link{MEDME.readFiles}} function. Afterwards,
the same obejct is returned by most of the MEDME library function. Each time, a new slot is filled with additional data, as smoothed
logR or Absolute/Relative Methylation Scores (AMS and RMS respectively). At the end of the analysis, usually after a call to the
\code{\link{MEDME.predict}} function, the \code{\link{MEDME.writeFiles}} function can be used to generate SGR or GFF files from this object.
}
\section{Slots}{
	 \describe{
    \item{\code{chr}:}{Object of class \code{"character"} : the probe-level chromosome asignments}
    \item{\code{pos}:}{Object of class \code{"numeric"} : the probe-level genomic position}
    \item{\code{logR}:}{Object of class \code{"matrix"} : the probe-level un-trasformed normalized MeDIP logRatios for each sample}
    \item{\code{smoothed}:}{Object of class \code{"matrix"} : the probe-level smoothed MeDIP logRatios for each sample}
    \item{\code{AMS}:}{Object of class \code{"matrix"}  : the probe-level Absolute Methylation Score for each sample}
    \item{\code{RMS}:}{Object of class \code{"matrix"}  : the probe-level Relative Methylation Score for each sample}
    \item{\code{CGcounts}:}{Object of class \code{"numeric"}  : the probe-level count of CpGs}
    \item{\code{organism}:}{Object of class \code{"character"} : the organism that the probe genomic positions are referring to, either hsa or mmu for homo sapiens or mus musculus respectively}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MEDMEset")}: subsets the object based on its probes and/or samples}
    \item{AMS}{\code{signature(object = "MEDMEset")}: extracts the Absolute Methylation Score from the AMS slot}
    \item{CG}{\code{signature(object = "MEDMEset")}: extracts the probe CpG count from the CGcounts slot}
    \item{chr}{\code{signature(object = "MEDMEset")}: extracts the probe chromosomal assignment}
    \item{org}{\code{signature(object = "MEDMEset")}: extracts the organism}
    \item{initialize}{\code{signature(.Object = "MEDMEset")}: automatically generates smoothed, AMS and RMS matrix when only the logR slot is filled}
    \item{logR}{\code{signature(object = "MEDMEset")}: extracts the matrix of MeDIP un-transformed logRatios }
    \item{pos}{\code{signature(object = "MEDMEset")}: extracts the probe genomic position}
    \item{RMS}{\code{signature(object = "MEDMEset")}: extracts the Relative Methylation Score from the RMS slot}
    \item{show}{\code{signature(object = "MEDMEset")}: prints a summary of the object content}
    \item{smoothed}{\code{signature(object = "MEDMEset")}: extracts the Absolute Methylation Score from the AMS slot}
	 }
}
\references{\url{http://genome.cshlp.org/cgi/content/abstract/gr.080721.108v1}}
\author{Mattia Pelizzola}
\seealso{
\code{\link{MEDME.readFiles}}, \code{\link{MEDME.writeFiles}}
}
\examples{
showClass("MEDMEset")
}
\keyword{classes}
