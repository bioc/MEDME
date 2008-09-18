    #### class definition
    setClass(Class = 'MEDMEset', representation = representation(chr='character', pos='numeric', logR='matrix', smoothed='matrix', AMS='matrix', RMS='matrix', CGcount='numeric', genomeRelease='character'))
    # initialize method to generated empty smoothed, AMS and RMS matrix when an obect is created filling on ly the logR slot
    setMethod('initialize', 'MEDMEset', function(.Object, ...) {
        .Object <- callNextMethod()
        if(any(is.na(.Object@pos))) stop('NA are not allowed in pos slot ..')
        if(any(is.na(.Object@chr))) stop('NA are not allowed in chr slot..')
        if(length(.Object@pos)!=length(.Object@chr)) stop('chr and pos slots cannot contain data with different length ..')
        if(nrow(.Object@logR)!=length(.Object@pos)) stop('number of rows of logR is different from length of pos ..')
        if(is.null(rownames(.Object@logR))) stop('logR has to be a matrix with probeIds as rownames ..')

        Nrow = nrow(.Object@logR)
        Ncol = ncol(.Object@logR)
        NAmat = matrix(NA, Nrow, Ncol)
        rownames(NAmat) = rownames(.Object@logR)
        colnames(NAmat) = colnames(.Object@logR)
        if(nrow(.Object@smoothed) == 0) .Object@smoothed = NAmat
        if(nrow(.Object@AMS) == 0) .Object@AMS = NAmat
        if(nrow(.Object@RMS) == 0) .Object@RMS = NAmat
        
        # checking unexpected probe chromosomal assignments
        # eliminating probes with chromosomes not included in chrs
        chrs = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
        probeChr = .Object@chr
        extrachrInds = which(!(probeChr %in% chrs))
        if(length(extrachrInds)>0) {
            .Object@chr = .Object@chr[-extrachrInds]
            .Object@pos = .Object@pos[-extrachrInds]
            .Object@logR = .Object@logR[-extrachrInds, ]
            .Object@smoothed = .Object@smoothed[-extrachrInds, ]
            .Object@AMS = .Object@AMS[-extrachrInds, ]
            .Object@RMS = .Object@RMS[-extrachrInds, ]
            warning('probes assigned to chromosomes other than {chr1, .., chr22, chrX, chrY} have been excluded ..\n')
        }
        .Object
    })


    #### defining methods to extract chr, pos, logR, CGcount and genomeRelease (the latter with replacement too) as well as extension of [ and show methods
    # chr
    #chrfun = function(object) standardGeneric('chr')
    setGeneric('chr', function(object) standardGeneric('chr'))
    setMethod('chr','MEDMEset', function(object) object@chr)
    # pos
    #posfun = function(object) standardGeneric('pos')
    setGeneric('pos', function(object) standardGeneric('pos'))
    setMethod('pos','MEDMEset', function(object) object@pos)
    # logR
    #logRfun = function(object) standardGeneric('logR')
    setGeneric('logR', function(object) standardGeneric('logR'))
    setMethod('logR','MEDMEset', function(object) object@logR)
    # smoothed
    #smoothedfun = function(object) standardGeneric('smoothed')
    setGeneric('smoothed', function(object) standardGeneric('smoothed'))
    setMethod('smoothed','MEDMEset', function(object) object@smoothed)
    # AMS
    #AMSfun = function(object) standardGeneric('AMS')
    setGeneric('AMS', function(object) standardGeneric('AMS'))
    setMethod('AMS','MEDMEset', function(object) object@AMS)
    # RMS
    #RMSfun = function(object) standardGeneric('RMS')
    setGeneric('RMS', function(object) standardGeneric('RMS'))
    setMethod('RMS','MEDMEset', function(object) object@RMS)
    # CGcount
    #CGcountfun = function(object) standardGeneric('CG')
    setGeneric('CG', function(object) standardGeneric('CG'))
    setMethod('CG','MEDMEset', function(object) object@CGcount)
    # genomeRelease
    #gRfun = function(object) standardGeneric('gR')
    setGeneric('gR', function(object) standardGeneric('gR'))
    setMethod('gR','MEDMEset', function(object) object@genomeRelease)

    # subset
    setMethod('[','MEDMEset', function(x, i, j, drop) {
        if(missing(i)) i = 1:nrow(x@logR)
        if(missing(j)) j = 1:ncol(x@logR)
        x@chr = x@chr[i]
        x@pos = x@pos[i]
        x@logR = as.matrix(x@logR[i,j])
	x@smoothed = as.matrix(x@smoothed[i,j])
        x@AMS = as.matrix(x@AMS[i,j])
        x@RMS = as.matrix(x@RMS[i,j])
        x@CGcount = x@CGcount[i]
        x
    })

    # show
    setMethod('show','MEDMEset', function(object) {
            cat("S4 Object of class MEDMEset; ")
            cat(paste(nrow(logR(object)),'probes,', ncol(logR(object)),'samples\n'))
            maxP = min(3, length(object@chr))
            cat("\nchr : ")
            cat(object@chr[1:maxP])
            cat("\npos : ")
            cat(object@pos[1:maxP])
            cat("\nlogR :\n")
            show(object@logR[1:maxP,])
            cat("\nsmoothed :\n")
            show(object@smoothed[1:maxP,])
            cat("\nAMS :\n")
            show(object@AMS[1:maxP,])
            cat("\nRMS :\n")
            show(object@RMS[1:maxP,])
            cat("CGcount : ")
            cat(object@CGcount[1:maxP])
            cat("\ngenomeRelease : ")
            cat(object@genomeRelease)
            cat('\n')
        })

