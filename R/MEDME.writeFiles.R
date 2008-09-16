`MEDME.writeFiles` <-
function(data, output, path = getwd(), format, featureLength  = NULL) {
    format = tolower(format)
    # some checks
    if(class(data)!='MEDMEset') stop('data needs to be an object of class MEDMEset ..')
    if(output != 'logR' && output!= 'smoothed' && output!= 'AMS' && output!= 'RMS') stop('output must be one of logR, smoothed, AMS or RMS ..')
    if(format!='sgr' && format!='gff') {stop('format has to be either sgr or gff ..')}
    if(format == 'gff' && is.null(featureLength)) {stop('please provide feature Length ..')}

    if(output == 'logR') MEDMEout = logR(data)
    if(output == 'smoothed') MEDMEout = smoothed(data)
    if(output == 'AMS') MEDMEout = AMS(data)
    if(output == 'RMS') MEDMEout = RMS(data)
    
    if(format=='sgr') {
        for(i in 1:ncol(MEDMEout)) {
            filename = paste(path, '/', colnames(MEDMEout)[i], '.sgr', sep='')
            datamat = data.frame(chr(data), pos(data), MEDMEout[,i])
            datamat[,3] = signif(datamat[,3],4)
            write.table(datamat, file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
        }
    }
    if(format=='gff') {
        posLeft = round(pos(data) - featureLength/2)
        posRight = round(pos(data) + featureLength/2)
        colLabels = c('seqname','source','feature','start','end','score','strand','frame','group')
        for(i in 1:ncol(MEDMEout)) {
            filename = paste(path, '/', colnames(MEDMEout)[i], '.gff', sep='')
            datamat = data.frame(chr(data), '.', rownames(MEDMEout), posLeft, posRight, signif(MEDMEout[,i],4), '.', '.', '.')
            colnames(datamat) = colLabels
            write.table(datamat, file=filename, row.names = FALSE, col.names = TRUE, quote=FALSE)
        }
    }
}
