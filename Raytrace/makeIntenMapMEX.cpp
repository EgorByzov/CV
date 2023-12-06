#include "mex.h"
#include <math.h>

void makeIntenMap(double* intenMap, double* rays, int numRays, double* pointsVector, 
                        int numPoints, double distance, bool isSmoothing){
    // unpack rays
    double* raysX = rays;
    double* raysY = rays + numRays;
    double* raysZ = rays + 2*numRays;
    double* raysEn = rays + 3*numRays;
    // resVector
    double* resX = pointsVector;
    double* resY = pointsVector + numPoints;
    double* resZ = pointsVector + 2*numPoints;
    // other definitions
    double tempDist;
    int i, j;
    
    // let's go!
    if (!isSmoothing){
        int numHits;
        int* hitInds = (int*)malloc(numPoints*sizeof(int));
                
        for (i = 0; i < numRays; i++) {
            
            numHits = 0;
            for (j = 0; j < numPoints; j++){
                tempDist = raysX[i] * resX[j] + raysY[i] * resY[j] + raysZ[i] * resZ[j];
                if (tempDist > distance){
                    hitInds[numHits] = j;
                    numHits ++;
                }
            }// for
            
            if (numHits > 0){
                for (j = 0; j < numHits; j++){
                    intenMap[hitInds[j]] += raysEn[i] / numHits;
                }
            }// if
            
        }// for
    }else{
        double sigDryan = sqrt(1 - 9*pow(distance,2.0));
        double multipl = 1.0 / 2.0 / 3.14159265358979323846 / pow(distance,2.0);
        
        for (i = 0; i < numRays; i++) {
            for (j = 0; j< numPoints; j++){
                tempDist = raysX[i] * resX[j] + raysY[i] * resY[j] + raysZ[i] * resZ[j];
                if (tempDist > sigDryan){
                    intenMap[j] += multipl * raysEn[i] * exp( -(1.0 - pow(tempDist, 2.0)) / 2.0 / pow(distance,2.0) );
                }
            }
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    ////////////////////////////
    // definitions
    double *rays, *intenMap, *pointsVector;
    double distance;
    int numRays, numPoints;
    bool* isSmoothing;
    
    ////////////////////////////
    // check of input arguments
    // number of arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:invalidNumInputs",
                "4 inputs required: rays, result vectors, distance and smoothing mode");
    } else if( nlhs > 1 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:maxlhs",
                "Too many output arguments");
    }
    // rays
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
            !(mxGetNumberOfDimensions(prhs[0]) == 2) ||
            !(mxGetN(prhs[0]) == 4) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongTraces",
                "rays must be double matrix Nx4");
    }
    // vectors
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            !(mxGetNumberOfDimensions(prhs[1]) == 2) ||
            !(mxGetN(prhs[1]) == 3) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongTraces",
                "Result vectors must be double matrix Nx3");
    }    
    // distance
    if ( (mxGetNumberOfDimensions(prhs[2]) != 2) ||
         (mxGetM(prhs[2]) != 1) || (mxGetN(prhs[2]) != 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "Distance must be scalars");
    }
    // isSmoothing
    if ( (mxGetNumberOfDimensions(prhs[3]) != 2) ||
         (mxGetM(prhs[3]) != 1) || (mxGetN(prhs[3]) != 1) ||
         (!mxIsLogical(prhs[3])) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "Smoothing mode must be logical");
    }
    
    //////////////////////////////////////////////
    // Prepare data for C call
    rays = mxGetPr(prhs[0]);    
    numRays = (int) mxGetM(prhs[0]);
    pointsVector = mxGetPr(prhs[1]);
    numPoints = (int) mxGetM(prhs[1]);
    distance = mxGetScalar(prhs[2]);
    isSmoothing = (bool*) mxGetLogicals(prhs[3]);
    plhs[0] = mxCreateDoubleMatrix(numPoints, 1, mxREAL);
    intenMap = mxGetPr(plhs[0]);    
    
    //////////////////////////////////////////////
    // Run computation
    makeIntenMap(intenMap, rays, numRays, pointsVector, numPoints, distance, isSmoothing[0]);
}
