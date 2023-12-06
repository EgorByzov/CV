#include "trees\TriangTree.h"

#ifndef CPPONLY
#include "mex.h"

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){   
    // check of input arguments
    // number of arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:getIntersTriangMEX:invalidNumInputs",
                "2 inputs required: tree, rays");
    } else if( nlhs > 2 ) {
        mexErrMsgIdAndTxt( "MATLAB:getIntersTriangMEX:maxlhs",
                "Too many output arguments");
    }
    
    // tree
    if ( !mxIsUint64(prhs[0]) || mxIsComplex(prhs[0]) ||
            !(mxGetNumberOfDimensions(prhs[0]) == 2) ||
            !(mxGetN(prhs[0]) == 1) || !(mxGetM(prhs[0]) == 1) ) {
        mexErrMsgIdAndTxt( "MATLAB:getIntersTriangMEX:wrongTree",
                "tree must be scalar (pointer to tree, uint64)");
    }
    
    // rays
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            !(mxGetNumberOfDimensions(prhs[1]) == 2) ||
            !(mxGetN(prhs[1]) == 7) ) {
        mexErrMsgIdAndTxt( "MATLAB:getIntersTriangMEX:wrongRays",
                "rays must be double matrix Nx7");
    }

    unsigned long long pointer = (unsigned long long) mxGetPr(prhs[0])[0];
    TriangTree *tree = (TriangTree *) pointer;
    double *ray = mxGetPr(prhs[1]);
    int numRays = mxGetM(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(numRays, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(numRays, 1, mxREAL);
    
    double *t = mxGetPr(plhs[0]);
    double *ind = mxGetPr(plhs[1]);
    
    tree->getInter(ray, numRays, t, ind);
}
#endif