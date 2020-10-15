/*
 *
 * L:   Desired parameter lenght
 * N:   Number of nodes
 * K:   Number of iterations
 *
 *
 * A:   N x N
 * C:   N x N
 * w_0  L x N
 * u    L x N x K
 * d    K x N
 * mu   N x 1
 * w*   L x 1
 */


#include "diffusion.h"



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mwSize N,L,K;
    double *A, *C, *w_0, *u, *d, *mu, *w_star, *e, *w_k;

    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:nrhs","Seven inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:nlhs","Two output required.");
    }

    /* check for proper dimensions */

    /* Getting Dimensions */

    N=mxGetM(prhs[0]); // getting agents number
    L=mxGetM(prhs[2]); // geting the data length
    K=mxGetM(prhs[4]); // getting iteration number

    /* make sure A is a square (N x N) matrix */
    if( mxGetM(prhs[0])!=mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","A must be a square matrix.");
    }

    /* make sure A contains double type entries */
    if( !mxIsDouble(prhs[0]) ||
     mxIsComplex(prhs[0])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","A entries must be type double.");
    }

    /* make sure C is a square (N x N) matrix and has the same size as the matrix A */
    if( (mxGetM(prhs[1])!=mxGetN(prhs[1])) || (mxGetM(prhs[1])!=N) ) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","C must be square with the same size as A.");
    }

    /* make sure C contains double type entries */
    if( !mxIsDouble(prhs[1]) ||
     mxIsComplex(prhs[1])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","C entries must be type double.");
    }

    /* make sure w_0 has the correct dimensions */
    if(mxGetN(prhs[2])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble"," w_0 must as much columns as A and C .");
    }

    /* make sure w_0 contains double type entries */
    if( !mxIsDouble(prhs[2]) ||
     mxIsComplex(prhs[2])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","w_0 entries must be type double.");
    }

    /* make sure u has the correct dimensions */
    if(mxGetM(prhs[3])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","u must have as much rows as as w_o");
    }

    if((mxGetNumberOfElements(prhs[3])/(L*K))!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","u must have as much columns as A ");
    }

    if((mxGetNumberOfElements(prhs[3])/(L*N))!=K) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","u must be (L x N X K) where L: is w_star length, N is the number of nodes, and K the number of ieration ");
    }

    /* make sure u contains double type entries */
    if( !mxIsDouble(prhs[3]) ||
     mxIsComplex(prhs[3])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","u entries must be type double.");
    }

    /* make sure d has the correct dimensions */
    if(mxGetN(prhs[4])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","d must have as much columns as A ");
    }

    /* make sure d contains double type entries */
    if( !mxIsDouble(prhs[4]) ||
     mxIsComplex(prhs[4])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","d entries must be type double.");
    }

    /* make sure mu has the correct dimensions */
    if(mxGetM(prhs[5])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","mu must have as much rows as A ");
    }

    /* make sure mu contains double type entries */
    if( !mxIsDouble(prhs[5]) ||
     mxIsComplex(prhs[5])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","mu entries must be type double.");
    }


    /* make sure w_star has the correct dimensions */
    if(mxGetM(prhs[6])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","w_star length mus be the same as w_0 ");
    }

    /* make sure w_star contains double type entries */
    if( !mxIsDouble(prhs[6]) ||
     mxIsComplex(prhs[6])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","w_star entries must be type double.");
    }



//     arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);




    //mexPrintf("N= %d \n K=%d \n L=%d \n",N,K,L);

    /*input pointerq */
    A       = mxGetPr(prhs[0]);
    C       = mxGetPr(prhs[1]);
    w_0     = mxGetPr(prhs[2]);
    u       = mxGetPr(prhs[3]);
    d       = mxGetPr(prhs[4]);
    mu      = mxGetPr(prhs[5]);
    w_star  = mxGetPr(prhs[6]);

    /* output pointers */
    plhs[0] = mxCreateDoubleMatrix(K,N,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(L,N,mxREAL);


    e = mxGetPr(plhs[0]);
    w_k = mxGetPr(plhs[1]);



    /* call the computational routine */
    ATC2(N, L, K, w_star, w_0, A, C, mu, d, u, e, w_k);

}
