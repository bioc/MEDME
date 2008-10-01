\name{MEDME.readFiles}
\alias{MEDME.readFiles}
\title{ reading sgr or gff files for MEDME }
\description{
  allows to read sgr or gff files before submitting the data to MEDME analysis
}
\usage{
MEDME.readFiles(path = getwd(), files = NULL, format, organism)
}
\arguments{
  \item{path}{ string; the path where the files are stored; the current working directory is the default }
  \item{files}{ vector; optional vector of file names }
  \item{format}{ string; either sgr or gff to indicate the respective file formats }
  \item{organism}{ string; either hsa or mmu for homo sapiens and mus musculus respectively}
}
\details{
In case of GFF files (recommendend), tab-delimited files with header are expected with following fields:
chromosome, probe ids, start and stop chromosomal positions, and score are expected in columns 1, 3, 4, 5 and 6 repectively.
Multiple files are also expected to be in the same order of rows.

In case of sgr files (GFF is the preferred format), tab-delimited files with no header and chr, chr positions and score are expected in columns 1, 2 and 3 repectively.
Multiple files are also expected to be in the same order of rows.
}
\value{
An object of class MEDMEset. The column headers in the logR slot are determined from the file names. in case of SGR files the are not probe names and progressive numbers are used
in place of them. In case of GFF files the probe names are determined from the 3rd column.
}
