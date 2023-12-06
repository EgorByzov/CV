#include "mex.h"
#include <math.h>

void calculatePointsInds(double* pointsInds, double* weight1, double* weight2,
                                                    double numPoints1, double numPoints2) {
    int i1 = 0, i2 = 0, ind = 0;
    
    while (i1 < numPoints1 - 1 || i2 < numPoints2 - 1){
        if (i1 < numPoints1 - 1 && i2 < numPoints2 - 1){
            if (weight1[i1 + 1] < weight2[i2 + 1]){
                // choose point on the 1st curve
                pointsInds[ind] = i1 + 1;
                pointsInds[ind + 1] = i2 + numPoints1 + 1;
                pointsInds[ind + 2] = i1 + 2;
                i1 = i1 + 1;
            }else {
                // choose point on the 2nd curve
                pointsInds[ind] = i1 + 1;
                pointsInds[ind + 1] = i2 + numPoints1 + 1;
                pointsInds[ind + 2] = i2 + numPoints1 + 2;
                i2 = i2 + 1;
            } //if
        }else if (i1 < numPoints1 - 1 && i2 == numPoints2 - 1) {
            // choose point on the 1st curve
            pointsInds[ind] = i1 + 1;
            pointsInds[ind + 1] = i2 + numPoints1 + 1;
            pointsInds[ind + 2] = i1 + 2;
            i1 = i1 + 1;                    
        }else if (i1 == numPoints1 - 1 && i2 < numPoints2 - 1){
            // choose point on the 2nd curve
            pointsInds[ind] = i1 + 1;
            pointsInds[ind + 1] = i2 + numPoints1 + 1;
            pointsInds[ind + 2] = i2 + numPoints1 + 2;
            i2 = i2 + 1;
        } // if
        ind = ind + 3;
    } // while
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    ////////////////////////////
    // definitions    
    double *pointsInds;
    double *weight1;
    double *weight2;
    double numPoints1, numPoints2;
    int numInds;
    
    ////////////////////////////
    // check of input arguments
    // number of arguments
    if (nrhs != 5) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:invalidNumInputs",
                "5 inputs required: pointsInds, weight1, weight2, numPoints1, numPoints2");
    } else if( nlhs > 1 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:maxlhs",
                "Too many output arguments");
    }
    // pointsInds
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
            !(mxGetNumberOfDimensions(prhs[0]) == 2) ||
            !(mxGetN(prhs[0]) == 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongTraces",
                "pointsInds must be double matrix Nx1");
    }
    // weight1, weight2
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            !(mxGetNumberOfDimensions(prhs[1]) == 2) ||
            !(mxGetN(prhs[1]) == 1) || 
         !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            !(mxGetNumberOfDimensions(prhs[2]) == 2) ||
            !(mxGetN(prhs[2]) == 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongTraces",
                "weight1, weight2 must be double matrices Nx1");
    }    
    // numPoints1, numPoints2
    if ( (mxGetNumberOfDimensions(prhs[3]) != 2) ||
            (mxGetNumberOfDimensions(prhs[4]) != 2) ||            
            (mxGetM(prhs[3]) != 1) || (mxGetN(prhs[3]) != 1) ||
            (mxGetM(prhs[4]) != 1) || (mxGetN(prhs[4]) != 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "numPoints1, numPoints2 must be scalars");
    }   
    
    //////////////////////////////////////////////
    // Prepare data for C call
    numPoints1 = mxGetScalar(prhs[3]);
    numPoints2 = mxGetScalar(prhs[4]);
    
//     pointsInds = mxGetPr(prhs[0]);    
//     numInds = (int) mxGetM(prhs[0]);    
//     if ( numInds != ( (numPoints1 + numPoints2 - 2) * 3 ) ) {
//         mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
//                 "Wrong length of pointsInds");
//     }
    
    weight1 = mxGetPr(prhs[1]);
    numInds = (int) mxGetM(prhs[1]);    
    if ( numInds != numPoints1 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "Wrong size of weight1");
    } 
    
    weight2 = mxGetPr(prhs[2]);
    numInds = (int) mxGetM(prhs[2]);    
    if ( numInds != numPoints2 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "Wrong size of weight1");
    }
    
    plhs[0] = mxCreateDoubleMatrix(numPoints1 + numPoints2 - 2, 3, mxREAL);   
    pointsInds = mxGetPr(plhs[0]);
    //////////////////////////////////////////////
    // Run computation
    calculatePointsInds(pointsInds, weight1, weight2, numPoints1, numPoints2);
}
