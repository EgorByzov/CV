#include "mex.h"
#include <math.h>

void makeIrmap(double* irmap, double* traces, int numTraces, double sx, double sy,
        double* x, double* y, int numX, int numY, double dx, double dy) {
    // unpack traces X, Y, Flux
    double* tracesX = traces;
    double* tracesY = tracesX + numTraces;
    double* tracesFlux = tracesY + numTraces;
    // other definitions
    double xcMin, xcMax, ycMin, ycMax;
    int i, j, iiMin, iiMax, jjMin, jjMax, ii, jj, curIrrInd;
    double deltaX, deltaY, multCur;
    // constants
    double mult = 0.5 / 3.14159265358979323846 / sx / sy;
    double sx2m2 = 2 * sx * sx;
    double sy2m2 = 2 * sy * sy;
    
    // let's go!
    for (i = 0; i < numTraces; i++) {
        // borders of ray's area of influence
        xcMin = tracesX[i] - 3*sx;
        xcMax = tracesX[i] + 3*sx;
        ycMin = tracesY[i] - 3*sy;
        ycMax = tracesY[i] + 3*sy;
        
        // compute indices of borders
        iiMin = ceil((xcMin - x[0]) / dx);
        iiMax = floor((xcMax - x[0]) / dx);
        jjMin = ceil((ycMin - y[0]) / dy);
        jjMax = floor((ycMax - y[0]) / dy);
        if (iiMin < 0) iiMin = 0;
        if (jjMin < 0) jjMin = 0;
        if (iiMax >= numX) iiMax = numX - 1;
        if (jjMax >= numY) jjMax = numY - 1;
        
        // update irradiance map
        multCur = mult * tracesFlux[i];
        for (ii = iiMin; ii <= iiMax; ii++)
            for (jj = jjMin; jj <= jjMax; jj++) {
                curIrrInd = ii * numY + jj;
                deltaX = x[ii] - tracesX[i];
                deltaY = y[jj] - tracesY[i];
                irmap[curIrrInd] += multCur *
                        exp( - deltaX * deltaX / sx2m2
                        - deltaY * deltaY / sy2m2 );
            }
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    ////////////////////////////
    // definitions
    size_t mrows, ncols;
    double *traces, *x, *y, *irmap;
    double sx, sy, dx, dy;
    int numTraces, numX, numY;
    
    ////////////////////////////
    // check of input arguments
    // number of arguments
    if (nrhs != 7) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:invalidNumInputs",
                "7 inputs required: traces, sx, sy, x, y, dx, dy");
    } else if( nlhs > 1 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:maxlhs",
                "Too many output arguments");
    }
    // traces
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
            !(mxGetNumberOfDimensions(prhs[0]) == 2) ||
            !(mxGetN(prhs[0]) == 3) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongTraces",
                "traces must be double matrix Nx3");
    }
    // sx, sy, dx, dy
    if ( (mxGetNumberOfDimensions(prhs[1]) != 2) ||
            (mxGetNumberOfDimensions(prhs[2]) != 2) ||
            (mxGetNumberOfDimensions(prhs[5]) != 2) ||
            (mxGetNumberOfDimensions(prhs[6]) != 2) ||
            (mxGetM(prhs[1]) != 1) || (mxGetN(prhs[1]) != 1) ||
            (mxGetM(prhs[2]) != 1) || (mxGetN(prhs[2]) != 1) ||
            (mxGetM(prhs[5]) != 1) || (mxGetN(prhs[5]) != 1) ||
            (mxGetM(prhs[6]) != 1) || (mxGetN(prhs[6]) != 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongSxSyDxDy",
                "sx, sy, dx, dy must be scalars");
    }
    // x, y
    if ( !mxIsDouble(prhs[3])   ||   mxIsComplex(prhs[3])   ||
            !mxIsDouble(prhs[4])   ||   mxIsComplex(prhs[4])   ||
            mxGetNumberOfDimensions(prhs[3]) != 2   ||
            mxGetNumberOfDimensions(prhs[4]) != 2   ||
            mxGetM(prhs[3]) != 1   ||   mxGetM(prhs[4]) != 1 ) {
        mexErrMsgIdAndTxt( "MATLAB:makeIrmapMEX:wrongXY",
                "X and Y must be double rows");
    }
    
    //////////////////////////////////////////////
    // Prepare data for C call
    traces = mxGetPr(prhs[0]);
    numTraces = (int) mxGetM(prhs[0]);
    sx = mxGetScalar(prhs[1]);
    sy = mxGetScalar(prhs[2]);
    x = mxGetPr(prhs[3]);
    y = mxGetPr(prhs[4]);
    numX = (int) mxGetN(prhs[3]);
    numY = (int) mxGetN(prhs[4]);
    dx = mxGetScalar(prhs[5]);
    dy = mxGetScalar(prhs[6]);
    plhs[0] = mxCreateDoubleMatrix(numY, numX, mxREAL);
    irmap = mxGetPr(plhs[0]);
    
    //////////////////////////////////////////////
    // Run computation
    makeIrmap(irmap, traces, numTraces, sx, sy, x, y, numX, numY, dx, dy);
}
