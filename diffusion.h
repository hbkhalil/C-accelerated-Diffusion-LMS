//
//  diffusion.h
//
//
//  Created by Ibrahim Harrane on 11/05/2016.
//
//

#ifndef diffusion_h
#define diffusion_h

/*
 *
 * L:   Desired parameter lenght
 * N:   Number of nodes
 * K:   Number of iterations
 *
 *
 * A:           N x N
 * C:           N x N
 * w_0          L x N
 * u            L x N x K
 * d            K x N
 * mu           N x 1
 * H            L x N x K
 * H_grad       L x N x K
 * M            1 x 1
 * M_grad       1 x 1
 * w*   L x 1
 */


#include "mex.h"
#include "math.h"
#include <string.h>
//#include <omp.h>

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif



double dotProduct(double *x, double *y, mwSize n)
{
    mwSize i;
    double z=0;
    
    for (i=0; i<n; i++)
    {
        z += (*(x+i)) * (*(y+i));
    }
    return z;
}

void arrayProduct(double x, double *y, double *z, mwSize n)
{
    mwSize i;
    /* multiply each element y by x */
    for (i=0; i<n; i++) {
        *(z+i) = x * (*(y+i));
    }
}

void arraySum(double *x, double *y, double *z, mwSize n)
{
    mwSize i;
    /* add each element y to x */
    for (i=0; i<n; i++) {
        *(z+i) = *(x+i) + (*(y+i));
    }
}

void ATC_RMT_theo(mwSize N, mwSize L, int K, double *tw_0, double *A, double *mu,  double *B, double *B_noise, mxArray *options, double *MSD, double *err)
{
    
    
    mwSize i,j,k,m,l;
    
    
    
    
    double *tw_t, *tw, *Met,*Met_add,*Met_mult;
    
    
    
    
    
    
    /* Memory allocation */
    tw_t=mxMalloc(L*N*sizeof(double));
    tw=mxMalloc(L*N*sizeof(double));
    Met=mxMalloc(L*sizeof(double));
    Met_add=mxMalloc(L*sizeof(double));
    Met_mult=mxMalloc(L*sizeof(double));
    //MSD_temp=mxMalloc(L*sizeof(double));
    
    memcpy(tw,tw_0,L*N*sizeof(double));
    
    memset (Met,0,L*sizeof(double));
    memset (MSD,0,K*N*sizeof(double));
    
    for(m=0; m<L ;m++)
    {
        Met_add[m]=1;
        Met_mult[m]=1; 
    }
    
    
//             for(m=0; m<L ;m++)
//         {
//             mexPrintf("Met[%d]= %f \n",m,Met_add[m]);
//         }


    
    for(k=0;k<K;k++)
    {
        for(m=0; m<L ;m++)
        {
            Met[m]=(1-mu[1]*B[m+L*1])*(1-mu[1]*B[m+L*1]);
            Met_mult[m]*=Met[m];
            Met_add[m] += Met_mult[m];
            //Met_add[m + j*L + l*(L*N) + i*(L*N*N) ] += pow(Met[m ],k+1);
        }
        
        //mexPrintf("Met OK \n \n");
//         for(m=0; m<L ;m++)
//         {
//             mexPrintf("Met[%d]= %f \n",m,Met_add[m]);
//         }
        memcpy(tw_t,tw,L*N*sizeof(double));
        
        
        
        for (i=0;i<N;i++)
        {

            err[k+K*i]=0;
            for (m=0;m<L;m++)
            {
                err[k+K*i]+= fabs(tw_t[m+i*L])/L;
            }
                        
            //mexPrintf("err[%d]=%f \n \n",k,err[k+K*i]);
            
            //memset (tw+i*L,0,L*sizeof(double));
            
            for(m=0; m<L ;m++)
            {
                tw[m+i*L] = 0 ;
            }
            
            for (j=0;j<N;j++)
            {
                //err_mean(:,j)= err_mean(:,j) + A(k,j)*(err_mean_t(:,k).*(ones(L,1)-mu(k)*lambda_mean'))
                for(m=0; m<L ;m++)
                    if (IsNonZero(A[j+N*i]))
                    {
                        {
                        tw[m+i*L] += A[j+N*i] *( tw_t[m+L*j]*(1-mu[j]*B[m+L*j]) ) ;
                        }
                    }
                
            }
            

            
            //mexPrintf("afer err[%d]=%f \n \n",k,err[k+K*i]);
            
            //mexPrintf("tw OK \n \n");
            

            
            //mexPrintf("err OK \n \n");
            
            for (j=0;j<N ; j++)
                
            {
                for (l=0; l<N; l++)
                {
                    if (IsNonZero(A[j+N*i]) && IsNonZero(A[l+N*i]) )
                    {
                        for(m=0; m<L ;m++)
                        {
                            Met[m]=(1-mu[j]*B[m+L*j])*(1-mu[l]*B[m+L*l]);
                            MSD[k+K*i] +=  A[j+N*i]*A[l+N*i]*tw_t[m+L*j]*tw_t[m+L*l]*Met[m]+(1.0/N)*A[j+N*i]*A[l+N*i]*mu[j]*mu[l]*B_noise[m+L*i]*Met_add[m];
                        }
                    }
                    
                }
            }
            

            //mexPrintf("MSD OK \n \n");

            
            
            
            mexPrintf("iter: %d \t node %d \n",k,i);
            
            
        }
        
        
        
        //mexPrintf("iter: %d \n",k);
        
        
    }
    

    
    mxFree(tw_t);
    mxFree(tw);
    mxFree(Met);
    mxFree(Met_add);
    mxFree(Met_mult);

}

void ATC_RMT(mwSize N, mwSize L, int K, double *w_0, double *A, double *mu, double *e, double *err, double *w_k, double *w_star, mxArray *options)
{
    
    
    mwSize i,j,k,m,iter_batch;
    mwIndex index;
    
    
    
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i,*Lambda,*noise;
    mxArray *output_get_data[3], *input_get_data[1];
    
    
    
    
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*N*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    //Lambda=mxMalloc(L*sizeof(double));
    //noise=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    output_get_data[2]=options;
    iter_batch=mxGetScalar(mxGetField(options,0,"iter_batch"));
    //mexPrintf("iter_batch=%d \n",iter_batch);
    
    for(k=0;k<K;k++)
    {
        
        
        
        
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
            
            err[k+K*i]=0;
            for (m=0;m<L;m++)
            {
                err[k+K*i]+= fabs(w_star[m] - w_k[L*i+m])/L;
            }
            //mexPrintf("tw[%d]=%f \n \n",k,err[k+K*i]);
        }
        

        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        memset (s,0,L*N*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        
        if ((k % iter_batch)==0)
        {
            mexCallMATLAB( 3 , output_get_data, 1, &output_get_data[2] ,  "get_data_c");
            Lambda=mxGetPr(output_get_data[0]);
            noise=mxGetPr(output_get_data[1]);
        }
        
        //mexPrintf("mexcallMATLAB ok \n");
        for (i=0;i<N;i++)
        {
            //if (IsNonZero(C[j+N*i]))
            //{
            
            //phi=mu(i)*(uk*uk')*(w_star -W(:,i))
            //#pragma omp parallel for
            for(m=0; m<L ;m++)
            {
                phi[m+L*i] = w_k[L*i+m] + mu[i]*Lambda[m+L*(k % iter_batch)]*(w_star[m] - w_k[L*i+m])+(1.0/N)*mu[i]*noise[m+L*(k % iter_batch)];
            }
            
            //}
            
            
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        //#pragma omp parallel for
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                for(m=0;m<L;m++)
                {
                    w_k[L*i+m]+= A[j+N*i]*phi[m+L*j];
                    //w_k[L*i+m]=Lambda[m];
                    //mexPrintf("W_K[%d,%d] = %f \n ",i,j,w_k[L*i+m]);
                }
                }
            }
        }
        
        
        
        
        mexPrintf("iter %d \n",k);
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
    //mxFree(Lambda);
    //mxFree(noise);
}


/* ATC diffusion */
void ATC2(mwSize N, mwSize L, mwSize K, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *e, double *w_k)
{
    
    mwSize i,j,k,m;
    double *phi,*s,*u_j,*phi_j,p,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    
    //w_k=w_0
    for(m=0;m<L;m++)
    {
        for(j=0;j<N;j++)
        {
            *(w_k+m+j*L)= *(w_0+m+j*L);
        }
    }
    
    
    
    for(k=0;k<K;k++)
    {
        
        
        
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            *(e+i*K+k)=0;
            // e(k,i)=norm(w_star-w_k(:,i),2)^2;
            for(m=0;m<L;m++)
            {
                *(e+k+i*K)+=(*(w_star+m)-(*(w_k+L*i+m)))*(*(w_star+m)-(*(w_k+L*i+m)));
            }
        }
        
        /* intermediate estmeate Phi calculation */
        for (i=0;i<N;i++)
        {
            //s=0
            //w_i=w(:,i)
            //phi(:,i)=w_i
            for(m=0;m<L;m++)
            {
                *(s+m)=0;
                *(w_i+m)=*(w_k+L*i+m);
                *(phi+L*i+m)=*(w_i+m);
            }
            
            for(j=0;j<N;j++)
            {
                //u_j=u(:,j,k)
                for(m=0; m<L ;m++)
                {
                    *(u_j+m)=*(u+m+j*L+k*(N*L));
                }
                //p=u_j'*w_i
                p=dotProduct(u_j,w_i,L);
                //p=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                p=(*(mu+i))*(*(C+j+N*i))*(*(d+k+K*j)-p);
                //res=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)
                arrayProduct(p,u_j, res ,L);
                
                //s=s+mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)
                for(m=0;m<L;m++)
                {
                    *(s+m)+= *(res+m);
                }
            }
            
            //phi(:,i)=w_i+s
            for(m=0;m<L;m++)
            {
                *(phi+m+L*i)+=*(s+m);
            }
        }
        
        for(i=0;i<N;i++)
        {
            //w(:,i)=0
            for(m=0;m<L;m++)
            {
                *(w_k+m+i*L)=0;
            }
            
            
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                for(m=0;m<L;m++)
                {
                    *(w_k+L*i+m)+= (*(A+j+N*i))*(*(phi+m+L*j));
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}


/* ATC diffusion */
void ATC(mwSize N, mwSize L, mwSize K, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *e, double *w_k)
{
    
    mwSize i,j,k,m;
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*N*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        memset (s,0,L*N*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        for (i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    
                    //p=u_j'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_k[L*i+m];
                    }
                    
                    //p=En(k,i)*mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=mu[i]*C[j+N*i]*(d[k+K*j]-p);
                    
                    //phi(:,i)+=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)=s+mu(i)*C(j,i)*u_j*p2
                    for(m=0;m<L;m++)
                    {
                        phi[m+L*i] += u[m+j*L+k*(N*L)]*p2;
                    }
                }
            }
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    for(m=0;m<L;m++)
                    {
                        w_k[L*i+m]+= A[j+N*i]*phi[m+L*j];
                    }
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}



/* compressed diffusion2 */

void compressed_diffusion2(mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx;
    double *phi,*s,*u_j,*phi_j,p,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    
    //w_k=w_0
    for(m=0;m<L;m++)
    {
        for(j=0;j<N;j++)
        {
            *(w_k+m+j*L)= *(w_0+m+j*L);
        }
    }
    
    
    
    for(k=0;k<K;k++)
    {
        
        
        
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            *(e+i*K+k)=0;
            // e(k,i)=norm(w_star-w_k(:,i),2)^2;
            for(m=0;m<L;m++)
            {
                *(e+k+i*K)+=(*(w_star+m)-(*(w_k+L*i+m)))*(*(w_star+m)-(*(w_k+L*i+m)));
            }
        }
        
        /* intermediate estmeate Phi calculation */
        for (i=0;i<N;i++)
        {
            //s=0
            //w_i=w(:,i)
            //phi(:,i)=w_i
            for(m=0;m<L;m++)
            {
                *(s+m)=0;
                *(w_i+m)=*(w_k+L*i+m);
                *(phi+L*i+m)=*(w_i+m);
            }
            
            for(j=0;j<N;j++)
            {
                
                //w_i(idx(M+1:end))=w(M+1:end,j)
                for(m=M; m<L ;m++)
                {
                    idx=*(H+m+L*i+k*(N*L))-1;
                    *(w_i+idx)=*(w_k+L*j+idx);
                }
                
                //u_j=u(:,j,k)
                for(m=0; m<L ;m++)
                {
                    *(u_j+m)=*(u+m+j*L+k*(N*L));
                }
                //p=u_j'*w_i
                p=dotProduct(u_j,w_i,L);
                //p=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                p=(*(mu+i))*(*(C+j+N*i))*(*(d+k+K*j)-p);
                //res=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)
                arrayProduct(p,u_j, res ,L);
                
                //s=s+mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)
                for(m=0;m<L;m++)
                {
                    *(s+m)+= *(res+m);
                }
            }
            
            //phi(:,i)=w_i+s
            for(m=0;m<L;m++)
            {
                *(phi+m+L*i)+=*(s+m);
            }
        }
        
        for(i=0;i<N;i++)
        {
            //w(:,i)=0
            for(m=0;m<L;m++)
            {
                *(w_k+m+i*L)=0;
            }
            
            
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                for(m=0;m<L;m++)
                {
                    *(w_k+L*i+m)+= (*(A+j+N*i))*(*(phi+m+L*j));
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}

/* ATC diffusion with energy consideration */
void ATC_En(mwSize N, mwSize L, mwSize K, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u,double *En, double *e, double *w_k)
{
    
    mwSize i,j,k,m;
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*N*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        memset (s,0,L*N*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        for (i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    //p=u_j'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_k[L*i+m];
                    }
                    
                    //p=En(k,i)*mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=En[k+K*i]*En[k+K*j]*mu[i]*C[j+N*i]*(d[k+K*j]-p);
                    
                    //phi(:,i)+=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)=s+mu(i)*C(j,i)*u_j*p2
                    for(m=0;m<L;m++)
                    {
                        phi[m+L*i] += u[m+j*L+k*(N*L)]*p2;
                    }
                }
            }
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    for(m=0;m<L;m++)
                    {
                        w_k[L*i+m]+= A[j+N*i]*phi[m+L*j];
                    }
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}


/*compressed diffusion */
void compressed_diffusion (mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx;
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i, *phi_t, *w_temp;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    phi_t=mxMalloc(L*sizeof(double));
    w_temp=mxMalloc(L*N*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        
        memset (s,0,L*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        
        
        for (i=0;i<N;i++)
        {
            //w_i=w(:,i)
            //memcpy(w_i, w_k+L*i*sizeof(double),L*sizeof(double));
            
            for(m=0;m<L;m++)
            {
                w_i[m]=w_k[L*i+m];
            }
            
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    //w_k(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m<L ;m++)
                    {
                        idx=H[m+L*i+k*(N*L)]-1;
                        w_i[idx]=w_k[L*j+idx];
                    }
                    
                    //p=u_j'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_i[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=mu[i]*C[j+N*i]*(d[k+K*j]-p);
                    
                    //phi(:,i)+=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)=s+mu(i)*C(j,i)*u_j*p2
                    for(m=0;m<L;m++)
                    {
                        phi[m+L*i] += u[m+j*L+k*(N*L)]*p2;
                    }
                }
            }
            
        }
        
        
        memcpy(w_temp, w_k,L*N*sizeof(double));
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            
            
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    for(m=0;m<L;m++)
                    {
                        phi_t[m]=w_temp[L*j+m];
                    }
                    //w_t(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m< L ;m++)
                    {
                        idx=H[m+L*j+k*(N*L)]-1;
                        phi_t[idx]=phi[L*i+idx];
                    }
                    
                    if (i!=j)
                    {
                        for(m=0;m<L;m++)
                        {
                            w_k[L*i+m]+= A[j+N*i]*phi_t[m];
                        }
                    }
                    else
                    {
                        for(m=0;m<L;m++)
                        {
                            w_k[L*i+m]+= A[j+N*i]*phi[L*i+m];
                        }
                    }
                }
                
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
    mxFree(phi_t);
    mxFree(w_temp);
}

/*compressed diffusion with energy consideration */
void compressed_diffusion_En (mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *En, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx;
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        
        memset (s,0,L*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        
        
        for (i=0;i<N;i++)
        {
            //w_i=w(:,i)
            //memcpy(w_i, w_k+L*i*sizeof(double),L*sizeof(double));
            
            for(m=0;m<L;m++)
            {
                w_i[m]=w_k[L*i+m];
            }
            
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    //w_k(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m<L ;m++)
                    {
                        idx=H[m+L*i+k*(N*L)]-1;
                        w_i[idx]=w_k[L*j+idx];
                    }
                    
                    //p=u_j'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_i[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=En[k+K*i]*En[k+K*j]*mu[i]*C[j+N*i]*(d[k+K*j]-p);
                    
                    //phi(:,i)+=mu(i)*C(j,i)*u_j*(d(k,j))-u_j'*w_i)=s+mu(i)*C(j,i)*u_j*p2
                    for(m=0;m<L;m++)
                    {
                        phi[m+L*i] += u[m+j*L+k*(N*L)]*p2;
                    }
                }
            }
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    for(m=0;m<L;m++)
                    {
                        w_k[L*i+m]+= A[j+N*i]*phi[m+L*j];
                    }
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
    
}

/* doubly compressed Diffusion2 */
void doubly_compressed_diffusion2(mwSize N, mwSize L, mwSize K, mwSize M, mwSize M_grad, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *H_grad, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx,idx_grad;
    double *grad, *grad_i, *s, *u_j, *phi_j, p, *res, *w_i;
    
    /* Memory allocation */
    
    grad=mxMalloc(L*N*sizeof(double));
    grad_i=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    
    //w_k=w_0
    for(m=0;m<L;m++)
    {
        for(j=0;j<N;j++)
        {
            *(w_k+m+j*L)= *(w_0+m+j*L);
        }
    }
    
    
    
    for(k=0;k<K;k++)
    {
        
        
        
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            *(e+i*K+k)=0;
            // e(k,i)=norm(w_star-w_k(:,i),2)^2;
            for(m=0;m<L;m++)
            {
                *(e+k+i*K)+=(*(w_star+m)-(*(w_k+L*i+m)))*(*(w_star+m)-(*(w_k+L*i+m)));
            }
        }
        
        
        for (i=0;i<N;i++)
        {
            //s=0
            //w_i=w(:,i)
            for(m=0;m<L;m++)
            {
                *(s+m)=0;
                *(w_i+m)=*(w_k+L*i+m);
            }
            
            for(j=0;j<N;j++)
            {
                
                //w_i(idx(M+1:end))=w(M+1:end,j)
                for(m=M; m<L ;m++)
                {
                    idx=*(H+m+L*i+k*(N*L))-1;
                    *(w_i+idx)=*(w_k+L*j+idx);
                    *(grad+m+j*L)=0;
                }
                
                //u_j=u(:,j,k)
                for(m=0; m<L ;m++)
                {
                    *(u_j+m)=*(u+m+j*L+k*(N*L));
                }
                //p=u_j'*w_i
                p=dotProduct(u_j,w_i,L);
                //p=C(j,i)*(d(k,j))-u_j'*w_i)
                p=(*(C+j+N*i))*(*(d+k+K*j)-p);
                //res=C(j,i)*u_j*(d(k,j))-u_j'*w_i)
                arrayProduct(p, u_j, res, L);
                
                //grad(:,j)=res(:)
                for(m=0;m<L;m++)
                {
                    *(grad+m+j*L)= *(res+m);
                }
            }
            
            //grad_i=grad(:,i)
            for(m=0;m<L;m++)
            {
                *(grad_i+m)=*(grad+m+L*i);
            }
            
            for(j=0;j<N;j++)
            {
                
                
                //grad_i(idx(M+1:end))=w(M+1:end,j)
                for(m=M_grad; m<L ;m++)
                {
                    idx_grad=*(H_grad+m+L*i+k*(N*L))-1;
                    *(grad_i+idx_grad)=*(grad+idx_grad+L*j);
                }
                
                //s=s+mu(i)*C(j,i)*u_j*grad_p
                for(m=0;m<L;m++)
                {
                    *(s+m)+= (*(mu+i))*(*(grad_i+m));
                }
                
                
            }
            //w_k(:,i)=w_k+s
            for(m=0;m<L;m++)
            {
                *(w_k+m+L*i)+=*(s+m);
            }
            
            
        }
        
        
    }
    mxFree(grad);
    mxFree(grad_i);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}

/* doubly compressed Diffusion */
void doubly_compressed_diffusion(mwSize N, mwSize L, mwSize K, mwSize M, mwSize M_grad, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *H_grad, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx,idx_grad;
    double *grad, *grad_i, *phi, *phi_i, *phi_t,*s,*u_j,*phi_j,p,p2,*res,*w_i,*w_t,*grad_l,*grad_e, *w_temp;
    
    /* Memory allocation */
    grad=mxMalloc(L*N*sizeof(double));
    grad_e=mxMalloc(L*sizeof(double));
    grad_l=mxMalloc(L*sizeof(double));
    grad_i=mxMalloc(L*sizeof(double));
    phi=mxMalloc(L*N*sizeof(double));
    phi_i=mxMalloc(L*sizeof(double));
    phi_t=mxMalloc(L*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    w_t=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    w_temp=mxMalloc(N*L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        
        //phi=w_k
        
        
        
        
        memcpy(phi,w_k,L*N*sizeof(double));
        
        for (i=0;i<N;i++)
        {
            
            //w_i=w(:,i)
            //memcpy(w_i, w_k+L*i*sizeof(double),L*sizeof(double));
            for(m=0;m<L;m++)
            {
                w_i[m]=w_k[L*i+m];
                w_t[m]=w_k[L*i+m];
            }
            
            //s=0
            memset(s,0,L*sizeof(double));
            //grad=0;
            memset(grad,0,L*N*sizeof(double));
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    //w_k(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m<L ;m++)
                    {
                        idx=H[m+L*i+k*(N*L)]-1;
                        w_t[idx]=w_k[L*j+idx];
                    }
                    
                    //p=u_j'*w_t
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_t[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=C[j+N*i]*(d[k+K*j]-p);
                    
                    //grad_e=C(j,i)*(u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p));
                    for(m=0;m<L;m++)
                    {
                        grad_e[m] = u[m+j*L+k*(N*L)]*p2;
                    }
                    
                    //p=u_i'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+i*L+k*(N*L)]*w_i[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=C[j+N*i]*(d[k+K*i]-p);
                    
                    //grad_l=C(j,i)*(u(k,:,i)'*(d(k,j) - u(k,:,i)*W_p));
                    for(m=0;m<L;m++)
                    {
                        grad_l[m] = u[m+i*L+k*(N*L)]*p2;
                    }
                    
                    for(m=0;m<L;m++)
                    {
                        grad_i[m]=grad_e[m];
                    }
                    //if (IsNonZero(C[j+N*i]))
                    
                    //grad_i(idx(M+1:end))=grad(M+1:end,j)
                    for(m=M_grad; m<L ;m++)
                    {
                        idx_grad=H_grad[m+L*j+k*(N*L)]-1;
                        grad_i[idx_grad]=grad_l[idx_grad];
                    }
                    //s=s+mu(i)*C(j,i)*u_j*grad_p
                    for(m=0;m<L;m++)
                    {
                        s[m]+= mu[i]*grad_i[m];
                    }
                }
            }
            
            //w_k(:,i)=w_k+s
            for(m=0;m<L;m++)
            {
                phi[m+L*i]+=s[m];
            }
        }
        
        memcpy(w_temp, w_k,L*N*sizeof(double));
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            
            
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    for(m=0;m<L;m++)
                    {
                        phi_t[m]=w_temp[L*j+m];
                    }
                    //w_t(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m< L ;m++)
                    {
                        idx=H[m+L*j+k*(N*L)]-1;
                        phi_t[idx]=phi[L*i+idx];
                    }
                    
                    if (i!=j)
                    {
                        for(m=0;m<L;m++)
                        {
                            w_k[L*i+m]+= A[j+N*i]*phi_t[m];
                        }
                    }
                    else
                    {
                        for(m=0;m<L;m++)
                        {
                            w_k[L*i+m]+= A[j+N*i]*phi[L*i+m];
                        }
                    }
                }
                
            }
        }
        
        
        
        
        
    }
    
    
    
    
    mxFree(grad);
    mxFree(grad_e);
    mxFree(grad_l);
    mxFree(grad_i);
    mxFree(phi);
    mxFree(phi_i);
    mxFree(phi_t);
    mxFree(s);
    mxFree(u_j);
    mxFree(w_i);
    mxFree(w_t);
    mxFree(res);
    mxFree(w_temp);
}



/* doubly compressed Diffusion with energy consideration */
void doubly_compressed_diffusion_En(mwSize N, mwSize L, mwSize K, mwSize M, mwSize M_grad, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *H_grad, double *En, double *e, double *w_k)
{
    
    mwSize i,j,k,m,idx,idx_grad;
    double *grad, *grad_i, *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i,*w_t,*grad_l,*grad_e;
    
    /* Memory allocation */
    grad=mxMalloc(L*N*sizeof(double));
    grad_e=mxMalloc(L*sizeof(double));
    grad_l=mxMalloc(L*sizeof(double));
    grad_i=mxMalloc(L*sizeof(double));
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    w_t=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        
        //phi=w_k
        
        
        
        memcpy(phi,w_k,L*N*sizeof(double));
        
        
        for (i=0;i<N;i++)
        {
            //w_i=w(:,i)
            //w_t=w(:,i)
            //memcpy(w_i, w_k+L*i*sizeof(double),L*sizeof(double));
            for(m=0;m<L;m++)
            {
                w_i[m]=w_k[L*i+m];
                w_t[m]=w_k[L*i+m];
            }
            
            //s=0
            memset(s,0,L*sizeof(double));
            //grad=0;
            memset(grad,0,L*N*sizeof(double));
            for(j=0;j<N;j++)
            {
                
                if (IsNonZero(C[j+N*i]))
                {
                    //w_t(idx(M+1:end),i)=w(idx(M+1:end),j)
                    for(m=M; m<L ;m++)
                    {
                        idx=H[m+L*i+k*(N*L)]-1;
                        w_t[idx]=w_k[L*j+idx];
                    }
                    
                    //p=u_j'*w_t
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_t[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=C[j+N*i]*(d[k+K*j]-p);
                    
                    //grad_e=C(j,i)*(u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p));
                    for(m=0;m<L;m++)
                    {
                        grad_e[m] = u[m+j*L+k*(N*L)]*p2;
                    }
                    
                    //p=u_i'*w_i
                    p=0;
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+i*L+k*(N*L)]*w_i[m];
                    }
                    
                    //p2=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
                    p2=En[k+K*j]*C[j+N*i]*(d[k+K*i]-p);
                    
                    //grad_l=C(j,i)*(u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p));
                    for(m=0;m<L;m++)
                    {
                        grad_l[m] = u[m+i*L+k*(N*L)]*p2;
                    }
                    
                    for(m=0;m<L;m++)
                    {
                        grad_i[m]=grad_e[m];
                    }
                    //if (IsNonZero(C[j+N*i]))
                    //{
                    //grad_i(idx(M+1:end))=grad(M+1:end,j)
                    for(m=M_grad; m<L ;m++)
                    {
                        idx_grad=H_grad[m+L*j+k*(N*L)]-1;
                        grad_i[idx_grad]=grad_l[idx_grad];
                    }
                    //s=s+mu(i)*C(j,i)*u_j*grad_p
                    for(m=0;m<L;m++)
                    {
                        s[m]+= En[k+K*i]* mu[i]*grad_i[m];
                    }
                }
            }
            
            
            
            
            
            
            
            //w_k(:,i)=w_k+s
            for(m=0;m<L;m++)
            {
                w_k[m+L*i]+=s[m];
            }
            
            
        }
        
        
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}


/* Partial model diffusion 2014*/
void ATC_partial_model(mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *mu, double *d, double *u ,double *H, double *e, double *w_k)
{
    
    mwSize i,j,k,m, idx;
    double *phi,*phi_i,*s,*u_j,p,p2,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    phi_i=mxMalloc(L*sizeof(double));
    s=mxMalloc(L*N*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        memset (s,0,L*N*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        for (i=0;i<N;i++)
        {
            
            //p=u_j'*w_i
            p=0;
            for(m=0; m<L ;m++)
            {
                p +=u[m+i*L+k*(N*L)]*w_k[L*i+m];
            }
            
            //p=mu(i)*(d(k,j))-u_j'*w_i)
            p2=mu[i]*(d[k+K*i]-p);
            
            //phi(:,i)+=mu(i)*u_i*(d(k,i))-u_i'*w_i)=s+mu(i)*u_i*p2
            for(m=0;m<L;m++)
            {
                phi[m+L*i] += u[m+i*L+k*(N*L)]*p2;
            }
            
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        
        for(i=0;i<N;i++)
        {
            
            
            //w(:,i)=w(:,i)+ a(j,i)*(H*phi(:,j)+(I-H))*phi(:,i)
            for(j=0;j<N;j++)
            {
                if (IsNonZero(A[j+N*i]))
                {
                    
                    //phi_i=phi(:,i)
                    
                    for(m=0;m<L;m++)
                    {
                        phi_i[m]=phi[L*i+m];
                    }
                    //phi_i(idx)=w(idx,j)
                    for(m=0; m<M ;m++)
                    {
                        idx=H[m+L*i+k*(N*L)]-1;
                        phi_i[idx]=phi[L*j+idx];
                    }
                    for(m=0;m<L;m++)
                    {
                        w_k[L*i+m]+= A[j+N*i]*phi_i[m];
                    }
                }
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(phi_i);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}

/* Reduced Communication ATC Diffusion Arablouei 2015 */
/* The matrix H here is a N x N x K, it containes scrambled indeces vector for each node, we use it to select a subset of nodes to communicate with. we choose the M first indices.*/
void ATC_RCD(mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *mu, double *d, double *u ,double *H, double *e, double *w_k)
{
    
    mwSize i,j,k,m;
    double *phi,*s,*u_j,*phi_j,p,p2,*res,*w_i;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    s=mxMalloc(L*N*sizeof(double));
    u_j=mxMalloc(L*sizeof(double));
    w_i=mxMalloc(L*sizeof(double));
    res=mxMalloc(L*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        /* intermediate estmeate Phi calculation */
        //s=0
        //phi=w_k
        memset (s,0,L*N*sizeof(double));
        memcpy(phi,w_k,L*N*sizeof(double));
        for (i=0;i<N;i++)
        {
            //p=u_j'*w_i
            p=0;
            for(m=0; m<L ;m++)
            {
                p +=u[m+i*L+k*(N*L)]*w_k[L*i+m];
            }
            
            //p=mu(i)*(d(k,j))-u_j'*w_i)
            p2=mu[i]*(d[k+K*i]-p);
            
            //phi(:,i)+=mu(i)*u_j*(d(k,j))-u_j'*w_i)
            for(m=0;m<L;m++)
            {
                phi[m+L*i] += u[m+i*L+k*(N*L)]*p2;
            }
            
            
        }
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
                
                for(m=0;m<L;m++)
                {
                    w_k[L*i+m]+= A[j+N*i]*(H[j+N*i+k*(N*N)]*phi[m+L*j]+(1-H[j+N*i+k*(N*N)])*phi[m+L*i]);
                }
                
            }
        }
        
        
        
    }
    mxFree(phi);
    mxFree(s);
    mxFree(w_i);
    mxFree(u_j);
    mxFree(res);
}

/* Compressive diffusion ATC sayin 2014 */
/* The confidence parameter is calculated using the equation 78 in the same paper */
/* The matrix H here contain projection vectors for K iteration */
void compressive_diff_ATC(mwSize N, mwSize L, mwSize K, double *w_star, double *w_0, double *A, double *mu, double *d, double *u,double *H, double *eta, double* delta, double mu_cvx,  double *e, double *w_k)
{
    
    mwSize i,j,k,m,n;
    double *phi, p, p2, *err, *Gamma, *phi_t, *epsilon ,*alpha;
    
    /* Memory allocation */
    phi=mxMalloc(L*N*sizeof(double));
    phi_t=mxMalloc(L*N*sizeof(double));
    
    Gamma=mxMalloc(L*N*sizeof(double));
    err=mxMalloc(N*sizeof(double));
    epsilon=mxMalloc(N*sizeof(double));
    alpha=mxMalloc(N*sizeof(double));
    
    //w_k=w_0
    memcpy(w_k,w_0,L*N*sizeof(double));
    //alpa=zeros(N,1);
    memset (alpha,0,N*sizeof(double));
    //phi_t=zeros(L,N);
    memset (phi_t,0,L*N*sizeof(double));
    //Gamma=zeros(L,N);
    memset (Gamma,0,L*N*sizeof(double));
    
    
    for(k=0;k<K;k++)
    {
        /* Error calculation*/
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            e[k+i*K]=0;
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                e[k+i*K]+=(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
        }
        
        
        
        /* intermediate estmeate Phi calculation */
        
        //phi=w_k
        memcpy(phi,w_k,L*N*sizeof(double));
        //epsilon=zeros(N,1);
        memset (epsilon,0,N*sizeof(double));
        
        for (i=0;i<N;i++)
        {
            //p=u_j'*w_i
            p=0;
            for(m=0; m<L ;m++)
            {
                p +=u[m+i*L+k*(N*L)]*w_k[L*i+m];
            }
            
            // err(i)=(d(k,i) - u(k,:,i)*W(:,i))
            err[i]=(d[k+K*i]-p);
            
            //p2=mu(i)*(d(k,j))-u_j'*w_i)
            p2=mu[i]*err[i];
            
            //phi(:,i)+=mu(i)*u_j*(d(k,j))-u_j'*w_i)
            for(m=0;m<L;m++)
            {
                phi[m+L*i] += u[m+i*L+k*(N*L)]*p2;
            }
            
            //epsilon(i)=H(:,i,k)'*(Phi(:,i)-Gamma(:,i));
            
            for(m=0;m<L;m++)
            {
                epsilon[i] += H[m+L*i+k*(N*L)]*(phi[m+L*i]- Gamma[m+L*i]);
            }
            
            //Gamma(:,i)=Gamma(:,i)+eta(i)*H(:,i,k)*epsilon(i);
            for(m=0;m<L;m++)
            {
                Gamma[m+i*L] += eta[i]*H[m+i*L+k*(L*N)]*epsilon[i];
            }
            
            
            // p=u(k,:,i)*(Phi(:,i)-Phi_t(:,i))
            p=0;
            for(m=0; m<L ;m++)
            {
                p +=u[m+i*L+k*(N*L)]*(phi[m+L*i]-phi_t[m+L*i]);
            }
            
            
            // alpha(i)=alpha(i)-mu_cvx*err(i)*u(k,:,i)*(Phi(:,i)-Phi_t(:,i))*delta(i)*(1-delta(i));
            alpha[i]-=mu_cvx * err[i]*p*delta[i]*(1-delta[i]);
            
            // delta(i)=1/(1+exp(-alpha(i)));
            delta[i]=1/(1+exp(-alpha[i]));
            
            
            
        }
        
        
        //w_k=0
        memset (w_k,0,L*N*sizeof(double));
        for(i=0;i<N;i++)
        {
            //Phi_t(:,i)=A(i,i)*Phi(:,i);
            for(m=0; m<L ;m++)
            {
                phi_t[m+i*L] =A[i+N*i]*phi[m+i*L];
            }
            
            for(j=0;j<N;j++)
            {
                if(i != j)
                {
                    if (IsNonZero(A[j+N*i]))
                    {
                        
                        //Phi_t(:,i)=Phi_t(:,i)+A(j,i)*Gamma(:,j);
                        for(m=0; m<L ;m++)
                        {
                            phi_t[m+i*L] += A[j+N*i]*Gamma[m+j*L];
                        }
                    }
                }
            }
            //W(:,i)=(1-delta(i))*Phi_t(:,i)+delta(i)*Phi(:,i);
            for(m=0;m<L;m++)
            {
                w_k[L*i+m]= (1-delta[i])*phi_t[m+L*i]+delta[i]*phi[m+L*i];
            }
        }
        
        
        
        
    }
    mxFree(phi);
    mxFree(phi_t);
    mxFree(Gamma);
    mxFree(err);
    mxFree(epsilon);
    mxFree(alpha);
    
}

#endif /* diffusion_h */
