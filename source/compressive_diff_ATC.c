/*
 *
 * L:   Desired parameter lenght
 * N:   Number of nodes
 * K:   Number of iterations
 *
 *
 * A:       N x N
 * w_0      L x N
 * u        L x N x K
 * d        K x N
 * mu       N x 1
 * H        L x N x K
 * w        L x 1
 * eta      N x 1
 * delta:   N x 1          initial value for the confidence parameter (it is calculated using equation 78)
 * mu_cvx:  1 x 1          confidence parameter calculation step size
 */


#include "diffusion.h"

#include "mex.h"




/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mwSize N,L,K;
    double *A, *w_0, *u, *d, *mu, *w_star, *e, *w_k ,*H, *eta, *delta, mu_cvx;
    
    /* check for proper number of arguments */
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:nrhs","Nine inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:nlhs","Two output required.");
    }
    
    /* check for proper dimensions */
    
    /* Getting Dimensions */
    
    N=mxGetM(prhs[0]); // getting agents number
    L=mxGetM(prhs[1]); // geting the data length
    K=mxGetM(prhs[3]); // getting iteration number
    
    /* make sure A is a square (N x N) matrix */
    if( mxGetM(prhs[0])!=mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","A must be a square matrix.");
    }
    
    /* make sure A contains double type entries */
    if( !mxIsDouble(prhs[0]) ||
            mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","A entries must be type double.");
    }
    
    
    /* make sure w_0 has the correct dimensions */
    if(mxGetN(prhs[1])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble"," w_0 must as much columns as A and C .");
    }
    
    /* make sure w_0 contains double type entries */
    if( !mxIsDouble(prhs[1]) ||
            mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_0 entries must be type double.");
    }
    
    /* make sure u has the correct dimensions */
    if(mxGetM(prhs[2])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must have as much rows as as w_o");
    }
    
    if((mxGetNumberOfElements(prhs[2])/(L*K))!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must have as much columns as A ");
    }
    
    if((mxGetNumberOfElements(prhs[2])/(L*N))!=K) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must be (L x N X K) where L: is w_star length, N is the number of nodes, and K the number of ieration ");
    }
    
    /* make sure u contains double type entries */
    if( !mxIsDouble(prhs[2]) ||
            mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u entries must be type double.");
    }
    
    /* make sure d has the correct dimensions */
    if(mxGetN(prhs[3])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","d must have as much columns as A ");
    }
    
    /* make sure d contains double type entries */
    if( !mxIsDouble(prhs[3]) ||
            mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","d entries must be type double.");
    }
    /* make sure mu has the correct dimensions */
    if(mxGetM(prhs[4])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","mu must have as much rows as A ");
    }
    
    /* make sure mu contains double type entries */
    if( !mxIsDouble(prhs[4]) ||
            mxIsComplex(prhs[4])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","mu entries must be type double.");
    }
    
    /* make sure H has the correct dimensions */
    if(mxGetM(prhs[5])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble"," H must be the same size as the matrix A ");
    }
    
    if((mxGetNumberOfElements(prhs[5])/(L*K))!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H must be the same size as the matrix A ");
    }
    
    if((mxGetNumberOfElements(prhs[5])/(N*L))!=K) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H must be (N x N X K) where N is the number of nodes, and K the number of ieration ");
    }
    
    /* make sure H contains double type entries */
    if( !mxIsDouble(prhs[5]) ||
            mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H entries must be type double.");
    }
    
    
    /* make sure w_star has the correct dimensions */
    if(mxGetM(prhs[6])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_star length mus be the same as w_0 ");
    }
    
    /* make sure w_star contains double type entries */
    if( !mxIsDouble(prhs[6]) ||
            mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_star entries must be type double.");
    }
    
    /* make sure eta has the correct dimensions */
    if(mxGetM(prhs[7])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","eta length must be equal to N ");
    }
    
    /* make sure eta contains double type entries */
    if( !mxIsDouble(prhs[7]) ||
            mxIsComplex(prhs[7])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","eta entries must be type double.");
    }
    
    /* make sure eta has the correct dimensions */
    if(mxGetM(prhs[8])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","delta length must be equal to N ");
    }
    
    /* make sure delta contains double type entries */
    if( !mxIsDouble(prhs[8]) ||
            mxIsComplex(prhs[8])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","delta entries must be type double.");
    }
    
    /*check if mu_cvx is a scalar */
    if( !mxIsDouble(prhs[9]) ||
            mxIsComplex(prhs[9]) ||
            mxGetNumberOfElements(prhs[9]) != 1 ) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notScalar", "mu_cvx must be a scalar.");
    }
    
    
    
    
    
    
    
    /*input pointers */
    
    A               = mxGetPr(prhs[0]);
    w_0             = mxGetPr(prhs[1]);
    u               = mxGetPr(prhs[2]);
    d               = mxGetPr(prhs[3]);
    mu              = mxGetPr(prhs[4]);
    H               = mxGetPr(prhs[5]);
    w_star          = mxGetPr(prhs[6]);
    eta             = mxGetPr(prhs[7]);
    delta           = mxGetPr(prhs[8]);
    mu_cvx          = mxGetScalar(prhs[9]);
    
    
    
    //mexPrintf("N= %d \n K=%d \n L=%d  \n",N,K,L);
    
    
    /* output pointers */
    plhs[0] = mxCreateDoubleMatrix(K,N,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(L,N,mxREAL);
    
    
    e = mxGetPr(plhs[0]);
    w_k = mxGetPr(plhs[1]);
    
    
    
    /* call the computational routine */
    compressive_diff_ATC(N, L,  K, w_star, w_0, A, mu, d, u, H, eta, delta, mu_cvx, e, w_k);
    
}
