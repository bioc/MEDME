// compile with:
// R CMD SHLIB MEDMEweight.c
// load in R with:
// dyn.load("MEDMEweight.so")
// to determine only MeDIPw run in R as:
// res = .C("MEDMEweight", LENGTH=as.integer(length), POS=as.double(pos), MeDIP=as.double(MeDIP), WSIZE=as.double(wsize), WFUN=as.integer(wfun), BOOLEAN=as.integer(0))
// to determine MeDIPw as well as AMS and RMS run in R as:
// res = .C("MEDMEweight", LENGTH=as.integer(length), POS=as.double(pos), MeDIP=as.double(MeDIP), WSIZE=as.double(wsize), WFUN=as.integer(wfun), BOOLEAN=as.integer(1), CGcount=as.double(CGcount), AMS=as.double(MeDIP), RMS=as.double(MeDIP), as.double(c(modelParameters, Min, Max)))

//length ids the length of pos and MeDIP arrays
//pos is the array of positions
//MeDIP is the array of rwa MeDIP logRatios
//wsize is the window size
//wfun is the weighting function (0->none, 1->linear, 2->exp, 3->log)
//MS is either 0 or 1 if AMS and RMS do not need /need to be determined
//CGcount is the array of CpGcounts
//AMS can be provided as MeDIP values, it will come back as AMS values
//RMS can be provided as MeDIP values, it will come back as RMS values
//modelParameters is an array with the B,C,D,E MEDME model parameters as well as the Min and Max X values of the model
//you get back MeDIPw in place of MeDIP as well as AMS and RMS values (if required)

# include <R.h>
# include <stdio.h>
# include <math.h>

void MEDMEweight(int *arraySize, double *pos, double *MeDIP, double *wsize, int *wfun, int *MS, double *CGcount, double *AMS, double *RMS, double *modelParameters) {
  double Wsize = floor(*wsize/2);
  double WsizeLog = *wsize/18;
  int i, j, leftPos, rightPos, leftInd, rightInd;
  double weight, relPos;
  double weightSum=0;
  double MeDIPsum=0;
  double MeDIPw[*arraySize];
  double AMSfinal[*arraySize];
  double RMSfinal[*arraySize];
  double B,C,D,E, Min, Max, Mlim;
  double AMSsum=0;
  double CGsum=0;


  if(*MS == 1) {
    // AMS and RMS are required, hence the model parameters are extracted
      B = modelParameters[0];
      C = modelParameters[1];
      D = modelParameters[2];
      E = modelParameters[3];
      Min = modelParameters[4];
      Max = modelParameters[5];
  }

  for(i=0; i<*arraySize; i++) {
    // for each probe i determining the range around it to be considered
    leftPos = pos[i] - Wsize;
    rightPos = pos[i] + Wsize;
    // the hundred up- and down-stream probes are considered to identify the probes in the (leftPos, rightPos) range
    leftInd = i - 100;
    if(leftInd < 0) {leftInd = 0;}
    rightInd = i + 100;
    if(rightInd > *arraySize -1) {rightInd = *arraySize -1;}

    for(j=leftInd; j<=rightInd; j++) {
      if(pos[j] > leftPos && pos[j] < rightPos) {
	// the j probe is in the range (leftPos, rightPos) 
	if(*wfun == 0) {weight = 1;} // no weighting is required
	if(*wfun == 1) { // linear weighting
	  relPos = fabs(pos[j] - pos[i]);
	  weight = 1 - relPos/Wsize;
	}
	if(*wfun == 2) { // exponential weighting
	  relPos = fabs(pos[j] - pos[i]);
	  weight = 1 - (relPos * relPos)/(Wsize * Wsize);
	}
	if(*wfun == 3) { // logarithmic weighting
	  relPos = fabs(pos[j] - pos[i]);
	  weight = 1 - log10(1 + (relPos/(WsizeLog)));
	}
	MeDIPsum = MeDIPsum + MeDIP[j] * weight;
	weightSum = weightSum + weight; 
      }
    }
    MeDIPw[i] = MeDIPsum / weightSum; // weighted MeDIP average
    MeDIPsum = 0;
    weightSum = 0;

    if(*MS == 1) {
      // determining AMS, saving it in AMSfinal
      Mlim = MeDIPw[i];
      if(Mlim > D) {Mlim = D;} // the model is only meaningful for MeDIPw in (C,D)
      if(Mlim < C) {Mlim = C;}
      AMSfinal[i] = pow(2, pow( pow(E,B) * (D - Mlim) / (Mlim - C), 1/B));
      if(AMSfinal[i] > Max) {AMSfinal[i] = Max;}
      if(AMSfinal[i] < Min) {AMSfinal[i] = Min;}
    }
  }

  if(*MS == 1) {
    // determining RMS, saving it in RMSfinal
    for(i=0; i<*arraySize; i++) {
      leftPos = pos[i] - Wsize;
      rightPos = pos[i] + Wsize;

      leftInd = i - 100;
      if(leftInd < 0) {leftInd = 0;}
      rightInd = i + 100;
      if(rightInd > *arraySize -1) {rightInd = *arraySize -1;}
      
      for(j=leftInd; j<=rightInd; j++) {
	if(pos[j] > leftPos && pos[j] < rightPos) {
	  AMSsum = AMSsum + AMSfinal[j];
	  CGsum = CGsum + CGcount[j];
	}
      }
      RMSfinal[i] = AMSsum / CGsum;
      AMSsum = 0;
      CGsum = 0;
    }
  }  


  //exporting results to input arguments
  for(i=0; i<*arraySize; i++) {
    MeDIP[i] = MeDIPw[i];
    if(*MS == 1) {
      AMS[i] = AMSfinal[i];
      RMS[i] = RMSfinal[i];
    }
  }
 
}
