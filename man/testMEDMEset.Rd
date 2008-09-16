\name{testMEDMEset}
\alias{testMEDMEset}
\docType{data}
\title{Dataset of class MEDMEset for testing MEDME}
\description{
This dataset contains a subset of the data reported in references. It contains normalized un-smoothed probe-level MeDIP enrichment for 50000 probes. This is a random subset of a custom Nimblegen chromosome X tiling array. It is a two channels array with an resolution of 100bp and oligos of 60nt. The M value is reported only. The fullyMet column of the logR slot contains data from a calibration experiments where MeDIP has been applied to a fully methylated sample. The last two columns NBMEL and YUSAC2 contain DNA-methylation experimental data for two cell strains: NBMEL are newborn normal melanocytes cells and YUSAC2 a melanoma strain. Data was processed with within and between array normalization. The full dataset contains almost 380K probes. See references for details.
Chromosome, genomic position and logR of probes can be accessed with the methods chr, pos and logR respectively.
}
\usage{data(testMEDMEset)}
\format{MEDMEset}
\references{\url{http://genome.cshlp.org/cgi/content/abstract/gr.080721.108v1}}
