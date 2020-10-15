/* 
 *
 * L:   Desired parameter lenght
 * N:   Number of nodes
 * K:   Number of iterations
 *
 *
 * A:   N x N
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

    const char **fnames;       /* pointers to field names */
    
    
    
    
    int        ifield, nfields;
    mwSize     NStructElems;
    
    
    mwSize N,L;
    int K;
    
    double *A, *w_0, *mu, *w_star, *e, *err, *w_k;
    mxArray *options;
    
    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:nrhs","Six inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:nlhs","Three outputs required.");
    }
    
    /* check for proper dimensions */
    
    /* Getting Dimensions */
     
    N=mxGetM(prhs[0]); // getting agents number
    L=mxGetM(prhs[1]); // geting the data length
    
    
    /* make sure A is a square (N x N) matrix */
    if( mxGetM(prhs[0])!=mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","A must be a square matrix.");
    }
    
    /* make sure A contains double type entries */
    if( !mxIsDouble(prhs[0]) || 
     mxIsComplex(prhs[0])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","A entries must be type double.");
    }
    


    /* make sure w_0 has the correct dimensions */
    if(mxGetN(prhs[1])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble"," w_0 must as much columns as A and C .");
    }
    if(mxGetM(prhs[1])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble"," w_0 must as much lines data model .");
    }
    
    /* make sure w_0 contains double type entries */
    if( !mxIsDouble(prhs[1]) || 
     mxIsComplex(prhs[1])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","w_0 entries must be type double.");
    }
    
    
    
    /* make sure mu has the correct dimensions */
    if(mxGetM(prhs[2])!=N) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","mu must have as much rows as A ");
    }
    
    /* make sure mu contains double type entries */
    if( !mxIsDouble(prhs[2]) || 
     mxIsComplex(prhs[2])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","mu entries must be type double.");
    }
    
    
    /* make sure mu contains double type entries */
    if( !mxIsDouble(prhs[3]) || 
     mxIsComplex(prhs[3])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","mu entries must be type double.");
    }
    
    
    
    /* make sure w_star has the correct dimensions */
    if(mxGetM(prhs[4])!=L) {
        mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","w_star length must be the same as w_0 ");
    }
    

    /* make sure w_star contains double type entries */
    if( !mxIsDouble(prhs[4]) || 
     mxIsComplex(prhs[4])) {
    mexErrMsgIdAndTxt("DiffusionToolbox:ATC_RMT:notDouble","w_star entries must be type double.");
    }
    
    if(!mxIsStruct(prhs[5]))
        mexErrMsgIdAndTxt( "MATLAB:DiffusionToolbox:ATC_RMT:inputNotStruct",
                "options must be a structure.");
  
    
//     /* get input arguments */
//     nfields = mxGetNumberOfFields(prhs[5]);
//     NStructElems = mxGetNumberOfElements(prhs[5]);
//     /* allocate memory  for storing classIDflags */
//     
//     
//     
//     /* allocate memory  for storing pointers */
//     fnames = mxCalloc(nfields, sizeof(*fnames));
//     /* get field name pointers */
//     for (ifield=0; ifield< nfields; ifield++)
//     {
//         fnames[ifield] = mxGetFieldNameByNumber(prhs[5],ifield);
//     }
//     
//     /* create a 1x1 struct matrix for output  */
//     options = mxCreateStructMatrix(1, 1, nfields, fnames);
//     mxFree(fnames);
//     
//     for (ifield=0; ifield< nfields; ifield++)
//     {
//         mxSetFieldByNumber(options, 0, ifield, mxGetFieldByNumber(prhs[5],0,ifield));
//     }
//     
    
   
    
    //mexPrintf("N= %d \n K=%d \n L=%d \n",N,K,L);
    
    /*input pointerq */
    A       = mxGetPr(prhs[0]);
    w_0     = mxGetPr(prhs[1]);
    mu      = mxGetPr(prhs[2]);
    K       = mxGetScalar(prhs[3]);
    w_star  = mxGetPr(prhs[4]);
    
            
    /* output pointers */
    plhs[0] = mxCreateDoubleMatrix(K,N,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(L,N,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(K,N,mxREAL);
    
    
    e       = mxGetPr(plhs[0]);
    w_k     = mxGetPr(plhs[1]);
    err     = mxGetPr(plhs[2]);
    
    
    /* call the computational routine */
    ATC_RMT(N, L, K, w_0, A, mu, e, err, w_k, w_star, prhs[5]);
    
}






