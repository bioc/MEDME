`smooth` <-
function(data, wsize=1000, wFunction='linear') {
    # checks
    if(class(data)!='MEDMEset') stop('data needs to be an object of class MEDMEset')
    if(class(wsize)!='numeric') stop('wsize needs to be a number ..')
    if(length(wsize)>1) stop('wsize needs to be 1 number ..')
    wsize = round(wsize)
    if(wFunction!='linear' && wFunction!='exp' && wFunction!='log' && wFunction!='none') stop('wFunction needs to be one of [linear, exp, log or none] ..')

    # weighting functions
    if(wFunction=='linear') wFun = 1
    if(wFunction=='exp') wFun = 2
    if(wFunction=='log') wFun = 3
    if(wFunction=='none') wFun = 0

    chrs = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
    probeChr = chr(data)
    probePos = pos(data)
    probeMeDIP = logR(data)
    rowN = rownames(probeMeDIP)
    probeMeDIP = as.data.frame(probeMeDIP)
    
    MeDIPw = NULL
    for(chr in chrs) {
        # extracting chr level data and sorting based on positions
        inds = which(probeChr == chr)
        if(length(inds) == 0) next
        cat(chr,' ')
        pos = probePos[inds]
        MeDIP = as.data.frame(probeMeDIP[inds,])
        sortedInds = sort(pos, index.return=TRUE)$ix
        pos = pos[sortedInds]
        MeDIP = as.data.frame(MeDIP[sortedInds,])
        posLength = length(pos)
        MeDIPwChr = MeDIP
        for(mcol in 1:ncol(MeDIP)) {
            res = .C("MEDMEweight", LENGTH=as.integer(posLength), POS=as.double(pos), MeDIP=as.double(MeDIP[,mcol]), WSIZE=as.double(wsize), WFUN=as.integer(wFun), BOOLEAN=as.integer(0), PACKAGE='MEDME')
            MeDIPwChr[,mcol] = res$MeDIP
        }
        rownames(MeDIPwChr) = rownames(MeDIP)
        MeDIPw = rbind(MeDIPw, MeDIPwChr)
    }
    MeDIPw = as.matrix(MeDIPw)[rownames(probeMeDIP),]
    MeDIPwset = new('MEDMEset', chr = probeChr, pos = probePos, logR = logR(data), smoothed = MeDIPw, AMS = AMS(data), RMS = RMS(data), CGcount = CG(data), organism=org(data))
    cat('\n')
    return(MeDIPwset)
}

