#include "mex.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mwSize nzmax; 
    mwSize columns;
    mwIndex i, nnz;
    mxArray *A;
    double *pr;
        
    /* Check for proper number of input and output arguments */    
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    } 
    if(nlhs > 0){
        mexErrMsgTxt("Too many output arguments.");
    }
    A = prhs[0];
    if (!mxIsSparse(A))  {
    	mexErrMsgTxt("Input argument must be a sparse array.");
    } 
    nzmax = mxGetNzmax(A);
    columns = mxGetN(A);
    
    /* NOTE: nnz is the actual number of nonzeros and is stored as the
       last element of the jc array where the size of the jc array is the
       number of columns + 1 */
    nnz = *(mxGetJc(A) + columns);    
    pr = mxGetPr(A);
    for (i = 0; i < nnz; ++i) {
        pr[i] = pr[i] + 1;
        mexPrintf(" %.1f", pr[i]);
    }
}

