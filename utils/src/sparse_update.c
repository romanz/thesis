#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    double *pr, *qr; 
    mwSize pn, qn, i; 
    
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input arguments required."); 
    }
    if (nlhs > 0) {
        mexErrMsgTxt("Too many output arguments."); 
    } 
    
    pn = mxGetJc(prhs[0])[mxGetN(prhs[0])];
    qn = mxGetNumberOfElements(prhs[1]);
    
    if (pn != qn && qn > 1) {
        mexErrMsgTxt("Incorrect number of elements to update."); 
    }
    
    pr = mxGetPr(prhs[0]);
    qr = mxGetPr(prhs[1]);
    
    if (qn > 1) {
        for (i = 0; i < pn; ++i)
            pr[i] = qr[i];
    } else {
        for (i = 0; i < pn; ++i)
            pr[i] = qr[0];
    }    
    
    return;
    
}
