\name{MEDME.writeFiles}
\alias{MEDME.writeFiles}
\title{ writeFiles sgr or gff files from MEDME output }
\description{
  allows to write sgr or gff files after MEDME analysis
}
\usage{
MEDME.writeFiles(data, output, path = getwd(), format, featureLength  = NULL)
}
\arguments{
  \item{data}{ An object of class MEDMEset}
  \item{output}{ string; the name of the data slot to be written on the disk, either logR, smoothed, AMS or RMS}
  \item{path}{ string; the path where the files are stored; the current working directory is the default }
  \item{format}{ string; either sgr or gff to indicate the respective file formats }
  \item{featureLength}{ integer; in case of GFF file format the length of the features has to be provided to determine start and end positions }
}
\details{
One GFF or SGR file is provided for each sample of the data MEDMEset object.

In case of GFF files, tab-delimited files with header are provided with following fields for each probe:
chromosome, empty field, probe ids, start and stop chromosomal positions, and score and empty fields.

In case of sgr files, tab-delimited files with no header and chr, chr positions and score are provided.
}
