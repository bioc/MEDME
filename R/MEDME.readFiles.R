`MEDME.readFiles` <-
function(path = getwd(), files = NULL, format, hgRelease) {
    format = tolower(format)
    # checks
    if(format!='sgr' && format!='gff') {stop('format has to be either sgr or gff ..')}
    if(is.null(files)) {
        if(format=='sgr') {
            if(length(grep('sgr', dir(path)))==0) stop('there are not sgr files in the provided path ..')
        }
        if(format=='gff') {
            if(length(grep('gff', dir(path)))==0) stop('there are not gff files in the provided path ..')
        }
    }
    if (hgRelease != "hg17" && hgRelease != "hg18")
        stop("currently only hgReleases hg17 and hg18 are supported ..")


    chrs = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
    datain= NULL
    previousPos = NULL
    if(format=='sgr') {
        # Tab-delimited files with no header and chr, chr positions and score are expected in columns 1, 2 and 3 repectively.
        # They are also expected to be in the same order of rows.
        if(is.null(files)) {
            filenames = dir(path, full.names=TRUE)[grep('sgr', dir(path))]
            filenamesShort = dir(path, full.names=FALSE)[grep('sgr', dir(path))]
            samplenames = gsub('.sgr','', filenamesShort, fixed=TRUE)
        }
        else {
            filenames = files
            samplenames = files
        }
        for(filename in filenames) {
            filedata = read.table(filename, sep='\t', header=FALSE, row.names=NULL)
            datapos = filedata[,2]
            if(!is.null(previousPos) && !identical(datapos, previousPos)) stop('files are not concordant in order or assignement of probes positions ..')
            else previousPos = datapos
            datain = cbind(datain, filedata[,3])
        }
        datachr = filedata[,1]
        colnames(datain) = samplenames
        datain = data.frame(as.character(datachr), datapos, datain, check.names=FALSE, stringsAsFactors=FALSE)
        colnames(datain)[1:2] = c('chr','pos')
        warning('the resulting dataset lacks official probe-names as rownames ..')
    }
    if(format=='gff') {
        # Tab-delimited files with header and chr, start and stop chr positions, and score are expected in columns 1, 4, 5 and 6 repectively.
        # They are also expected to be in the same order of rows.
        if(is.null(files)) {
            filenames = dir(path, full.names=TRUE)[grep('gff', dir(path))]
            filenamesShort = dir(path, full.names=FALSE)[grep('gff', dir(path))]
            samplenames = gsub('.gff','', filenamesShort, fixed=TRUE)
        }
        else {
            filenames = files
            samplenames = files
        }
        for(filename in filenames) {
            filedata = read.table(filename, sep='\t', header=TRUE, row.names=NULL)
            datapos = round(rowMeans(cbind(as.numeric(filedata[,4]), as.numeric(filedata[,5]))))
            if(!is.null(previousPos) && !identical(datapos, previousPos)) stop('files are not concordant in order or assignement of probes positions ..')
            else previousPos = datapos
            datain = cbind(datain, filedata[,6])
        }
        datachr = filedata[,1]
        dataprobes = filedata[,3]
        colnames(datain) = samplenames
        datain = data.frame(as.character(datachr), datapos, datain, check.names=FALSE, stringsAsFactors=FALSE)
        colnames(datain)[1:2] = c('chr','pos')
        if(length(dataprobes)== length(unique(dataprobes))) rownames(datain) = dataprobes
        else {warning('no unique probe names provided on column 3; the resulting dataset lacks rownames ..')}
    }

    if(length(intersect(datachr, chrs))==0) warning('currently only human chromosome names formatted as chr1, chr2, .. , chrX, chrY are supported ..')
    
    # creating MEDMEset object
    MEDMEsetObj = new('MEDMEset', chr = datain$chr, pos = datain$pos, logR = as.matrix(datain[,3:ncol(datain)]), genomeRelease = hgRelease)
    return(MEDMEsetObj)
}
