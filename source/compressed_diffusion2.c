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
 * H    L x N x K
 * M    1 x 1
 * w*   L x 1
 */


#include "diffusion.h"

#include "mex.h"




/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mwSize N,L,K,M;
    double *A, *C, *w_0, *u, *d, *mu, *w_star, *e, *w_k ,*H;
    
    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:nrhs","Nine inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:nlhs","Two output required.");
    }
    
    /* check for proper dimensions */
    
    /* Getting Dimensions */
    
    N=mxGetM(prhs[0]); // getting agents number
    L=mxGetM(prhs[2]); // geting the data length
    K=mxGetM(prhs[4]); // getting iteration number
    
    /* make sure A is a square (N x N) matrix */
    if( mxGetM(prhs[0])!=mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","A must be a square matrix.");
    }
    
    /* make sure A contains double type entries */
    if( !mxIsDouble(prhs[0]) || 
     mxIsComplex(prhs[0])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","A entries must be type double.");
    }
    
    /* make sure C is a square (N x N) matrix and has the same size as the matrix A */
    if( (mxGetM(prhs[1])!=mxGetN(prhs[1])) || (mxGetM(prhs[1])!=N) ) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_mex:notDouble","C must be square with the same size as A.");
    }
    
    /* make sure C contains double type entries */
    if( !mxIsDouble(prhs[1]) || 
     mxIsComplex(prhs[1])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","C entries must be type double.");
    }
    /* make sure w_0 has the correct dimensions */
    if(mxGetN(prhs[2])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble"," w_0 must as much columns as A and C .");
    }
    
    /* make sure w_0 contains double type entries */
    if( !mxIsDouble(prhs[2]) || 
     mxIsComplex(prhs[2])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_0 entries must be type double.");
    }
    /* make sure u has the correct dimensions */
    if(mxGetM(prhs[3])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must have as much rows as as w_o");
    }
    
    if((mxGetNumberOfElements(prhs[3])/(L*K))!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must have as much columns as A ");
    }
    
    if((mxGetNumberOfElements(prhs[3])/(L*N))!=K) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u must be (L x N X K) where L: is w_star length, N is the number of nodes, and K the number of ieration ");
    }
    
    /* make sure u contains double type entries */
    if( !mxIsDouble(prhs[3]) || 
     mxIsComplex(prhs[3])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","u entries must be type double.");
    }
    
    /* make sure d has the correct dimensions */
    if(mxGetN(prhs[4])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","d must have as much columns as A ");
    }
    
    /* make sure d contains double type entries */
    if( !mxIsDouble(prhs[4]) || 
     mxIsComplex(prhs[4])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","d entries must be type double.");
    }
    /* make sure mu has the correct dimensions */
    if(mxGetM(prhs[5])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","mu must have as much rows as A ");
    }
    
    /* make sure mu contains double type entries */
    if( !mxIsDouble(prhs[5]) || 
     mxIsComplex(prhs[5])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","mu entries must be type double.");
    }
    
    /* make sure H has the correct dimensions */
    if(mxGetM(prhs[6])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble"," H must have as much rows as as w_o ");
    }
    
    if((mxGetNumberOfElements(prhs[6])/(L*K))!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H must have as much columns as A ");
    }
    
    if((mxGetNumberOfElements(prhs[6])/(L*N))!=K) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H must be (L x N X K) where L: is w_star length, N is the number of nodes, and K the number of ieration ");
    }
    
    /* make sure H contains double type entries */
    if( !mxIsDouble(prhs[6]) || 
     mxIsComplex(prhs[6])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","H entries must be type double.");
    }
    /*check if M is a scalar */
    if( !mxIsDouble(prhs[7]) || 
     mxIsComplex(prhs[7]) ||
     mxGetNumberOfElements(prhs[7]) != 1 ) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notScalar", "M must be a scalar.");
    }
    
    /* make sure w_star has the correct dimensions */
    if(mxGetM(prhs[8])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_star length mus be the same as w_0 ");
    }
    
    /* make sure w_star contains double type entries */
    if( !mxIsDouble(prhs[8]) || 
     mxIsComplex(prhs[8])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:compressed_diffusion:notDouble","w_star entries must be type double.");
    }
      

    
    
    
    
    
    
    /*input pointers */
    A       = mxGetPr(prhs[0]);
    C       = mxGetPr(prhs[1]);
    w_0     = mxGetPr(prhs[2]);
    u       = mxGetPr(prhs[3]);
    d       = mxGetPr(prhs[4]);
    mu      = mxGetPr(prhs[5]);
    H       = mxGetPr(prhs[6]);
    M       = mxGetScalar(prhs[7]);
    w_star  = mxGetPr(prhs[8]);
    
    //mexPrintf("N= %d \n K=%d \n L=%d \n M=%d \n",N,K,L,M);
    
    /* output pointers */
    plhs[0] = mxCreateDoubleMatrix(K,N,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(L,N,mxREAL);
    
    
    e = mxGetPr(plhs[0]);
    w_k = mxGetPr(plhs[1]);
    
    
    
    /* call the computational routine */
    compressed_diffusion2(N, L, K,M, w_star, w_0, A, C, mu, d, u, H, e, w_k);
    
}
