#include <math.h>
#include "mex.h"

/* Sparse matrix non-zeroes update.
 *
 * This utility can be useful when you need to update a sparse matrix 
 * values without changing its sparsity pattern.
 * 
 * Since MATLAB uses CSC (Compressed Sparse Column) format for sparse 
 * matries, the non-zeroes are stored in an array, which can be safely 
 * overwritten, while preserving the sparsity pattern.
 *
 * NOTE: this function supports only real sparse matrices.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    double *pr, *qr; 
    mwSize pn, qn, i; 
    
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input arguments required."); 
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Number of effective non-zeroes in 1st argument */
    pn = mxGetJc(prhs[0])[mxGetN(prhs[0])];
    /* Number of elements to copy */
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
