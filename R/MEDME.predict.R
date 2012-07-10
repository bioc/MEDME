`MEDME.predict` <-
function(data, MEDMEfit, MEDMEextremes = c(1,32), wsize = 1000, wFunction='linear') {
    # checks
    if(class(data)!='MEDMEset') stop('data needs to be an object of class MEDMEset ..')
    if(class(MEDMEfit)[1] != 'drc') stop('MEDMEfit needs to be an object of class drc ..')
    if(class(wsize)!='numeric') stop('wsize needs to be a number ..')
    if(length(wsize)>1) stop('wsize needs to be 1 number ..')
    wsize = round(wsize)
    if(wFunction!='linear' && wFunction!='exp' && wFunction!='log' && wFunction!='none') stop('wFunction needs to be one of [linear, exp, log or none] ..')

    probeChr = chr(data)
    probePos = pos(data)
    probeCGcounts = CG(data)
    probeMeDIP = logR(data)
    
    # weighting functions
    if(wFunction=='linear') wFun = 1
    if(wFunction=='exp') wFun = 2
    if(wFunction=='log') wFun = 3
    if(wFunction=='none') wFun = 0

    chrs = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
    probeMeDIP = as.data.frame(probeMeDIP)
    AMS = NULL
    RMS = NULL
    MeDIPw = NULL
    for(chr in chrs) {
        # extracting chr level data and sorting based on positions
        inds = which(probeChr == chr)
        if(length(inds) == 0) next
        print(chr)
        pos = probePos[inds]
        MeDIP = probeMeDIP[inds,]
        CGcount = probeCGcounts[inds]
        sortedInds = sort(pos, index.return=TRUE)$ix
        pos = pos[sortedInds]
        MeDIP = MeDIP[sortedInds,]
        CGcount = CGcount[sortedInds]
        MeDIPwChr = as.data.frame(matrix(, length(pos), ncol(MeDIP)), stringsAsFactors =FALSE, check.names=FALSE)
        AMSchr = as.data.frame(matrix(, length(pos), ncol(MeDIP)), stringsAsFactors =FALSE, check.names=FALSE)
        RMSchr = as.data.frame(matrix(, length(pos), ncol(MeDIP)), stringsAsFactors =FALSE, check.names=FALSE)
        MeDIPwChr[,1] = chr
        MeDIPwChr[,2] = pos
        AMSchr[,1] = chr
        AMSchr[,2] = pos
        RMSchr[,1] = chr
        RMSchr[,2] = pos
        for(mcol in 1:ncol(MeDIP)) {
            res = .C("MEDMEweight", LENGTH=as.integer(length(pos)), POS=as.double(pos), MeDIP=as.double(MeDIP[,mcol]), WSIZE=as.double(wsize), WFUN=as.integer(wFun), BOOLEAN=as.integer(1), CGcount=as.double(CGcount), AMS=as.double(MeDIP[,mcol]), RMS=as.double(MeDIP[,mcol]), as.double(c(coef(MEDMEfit),MEDMEextremes)), PACKAGE='MEDME')
            MeDIPwChr[,mcol] = res$MeDIP
            AMSchr[,mcol] = res$AMS
            RMStemp = res$RMS
            RMStemp[RMStemp>4]=4
            RMSchr[,mcol] = RMStemp
        }
        rownames(MeDIPwChr) = rownames(MeDIP)
        rownames(AMSchr) = rownames(MeDIP)
        rownames(RMSchr) = rownames(MeDIP)
        MeDIPw = rbind(MeDIPw, MeDIPwChr)
        AMS = rbind(AMS, AMSchr)
        RMS = rbind(RMS, RMSchr)
    }
    colnames(MeDIPw) = colnames(MeDIP)
    colnames(AMS) = colnames(MeDIP)
    colnames(RMS) = colnames(MeDIP)
    MeDIPw = as.matrix(MeDIPw)[rownames(probeMeDIP),]
    AMS = as.matrix(AMS)[rownames(probeMeDIP),]
    RMS = as.matrix(RMS)[rownames(probeMeDIP),]

    cat('done\n')

    MEDMEmsset = new('MEDMEset', chr = probeChr, pos = probePos, logR = logR(data), smoothed = as.matrix(MeDIPw), AMS = as.matrix(AMS), RMS= as.matrix(RMS), CGcounts = CG(data), organism=org(data))
    return(MEDMEmsset)
}

