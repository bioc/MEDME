CGcount<-
function (data, wsize = 1000, wFunction = "linear")
{
    if(class(data)!='MEDMEset') stop('data needs to be an object of class MEDMEset ..')
    if (class(wsize) != "numeric")
        stop("wsize needs to be a number ..")
    if (length(wsize) > 1) 
        stop("wsize needs to be 1 number ..")
    wsize = round(wsize)
    orgname = org(data)
    if (wFunction != "linear" && wFunction != "exp" && wFunction !=
        "log" && wFunction != "none")
        stop("wFunction needs to be one of [linear, exp, log or none] ..")
    if (wFunction == "linear")
        wFun = function(d, Wsize = wsize) {
            sum(1 - abs(d/(Wsize/2)))
        }
    if (wFunction == "exp") 
        wFun = function(d, Wsize = wsize) {
            sum(1 - d^2/(Wsize/2)^2)
        }
    if (wFunction == "log") 
        wFun = function(d, Wsize = wsize) {
            sum(1 - log(1 + abs(d)/(Wsize/18), 10))
        }
    if (wFunction == "none") 
        wFun = function(d, Wsize = wsize) {
            return(length(d))
        }

    probePos = pos(data)
    probeChr = chr(data)
    chrs = c(paste("chr", 1:22, sep = ""), "chrX", "chrY")
    CGwindow = array(NA, dim = length(probePos))
    pattern = DNAString("CG")
    halfw = round(wsize/2)
    for (chr in chrs) {
        chrinds = which(probeChr == chr)
        if (length(chrinds) == 0) 
            next
        cat(chr, " ")
        if(orgname == 'hsa') chrseq = Hsapiens[[chr]]
        else chrseq = Mmusculus[[chr]]
        allmatches<-start(matchPattern(pattern, chrseq))
        rm(chrseq)
        gc()
        pos = probePos[chrinds]
        mx<-max(c(pos,allmatches))+wsize
        hash<-logical(length=mx)
        hash[allmatches]<-TRUE  # hash the location of all CpGs
        for (i in 1:length(pos)) {
            ind<-(pos[i]-halfw+1):(pos[i]+halfw)
            ind<-ind[ind>0]
            reldist<- ind[hash[ind]] - pos[i]
            CGwindow[chrinds[i]] = wFun(reldist)
        }
        rm(allmatches,hash)
        gc()
    }
    cat("\n")
    MEDMEcgset = new('MEDMEset', chr = probeChr, pos = probePos, logR = logR(data), smoothed = smoothed(data), AMS = AMS(data), RMS = RMS(data), CGcount = as.numeric(CGwindow), organism=org(data))
    return(MEDMEcgset)
}

