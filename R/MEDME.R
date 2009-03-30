`MEDME` <-
function(data, sample, CGcountThr = 1, figName = NULL) {
    # checks
    if(class(data)!='MEDMEset') stop('data needs to be an object of class MEDMEset ..')
    if(sample < 0 || sample > ncol(logR(data))) stop('sample does not exist ..')
    if(CGcountThr <0 || class(CGcountThr)!='numeric') stop('CGcountThr has to be a positive number ..')
    probeCGcounts = CG(data)
    if(all(is.na(probeCGcounts))) stop('CGcount slot of data only contains NA ..')
    if(nrow(smoothed(data)) == 0) stop('please use the smooth function to determine smoothed data ..')

    probeMeDIPw = smoothed(data)[,sample]
    # filtering data
    x = probeCGcounts[probeCGcounts>CGcountThr]
    y = probeMeDIPw[probeCGcounts>CGcountThr]
    x = log2(x)

    # determining median for each CG bin
    dataMedian=NULL
    xRange = quantile(x, c(0.01, 0.99), na.rm=TRUE)
    bins= seq(xRange[1], xRange[2], 0.1)
    for(i in bins) {
        inds = which(x>=i & x <i+0.1)
        binM = median(y[inds], na.rm=TRUE)
        dataMedian = c(dataMedian, binM)
    }

    # plotting the MeDIPw vs probeCGcounts scatter plot and determining the logistic model
    # require(MASS)
    # require(drc)
    if(!is.null(figName)) png(paste(figName,'png',sep='.'), 600, 600)
        NAinds = which(is.na(x))
        if(length(NAinds)>0) meanSd.kde2d<-kde2d(x[-NAinds], y[-NAinds], n=25)
        else meanSd.kde2d<-kde2d(x, y, n=25)
    	contour(meanSd.kde2d,col=colors()[41:60],nlevels=20,cex=0.7, add=FALSE,ylim=c(-1.5,1.5), xlim=c(0,7), xlab='log2(mCG)', ylab='MeDIP log2R')
        points(bins, dataMedian, pch=19, cex=1.2, col='red')
	dmNAinds = which(is.na(dataMedian))
	if(length(dmNAinds)>0) {
		dataMedian = dataMedian[-dmNAinds]
		bins = bins[-dmNAinds]
	}
        MEDMEmod = drm(dataMedian~bins, data=data.frame(dataMedian=dataMedian, bins=bins), fct=LL.4())
        lines(bins, fitted(MEDMEmod), lwd=4, col='blue')
    if(!is.null(figName)) dev.off()
    return(MEDMEmod)
}

