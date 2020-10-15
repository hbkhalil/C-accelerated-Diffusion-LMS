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
#include <omp.h>


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
        //#pragma omp parallel for private(m)
        for(i=0;i<N;i++)
        {
            // e(k,i)=0
            double temp=0;
            //*(e+i*K+k)=0;
            // e(k,i)=norm(w_star-w_k(:,i),2)^2;
            //#pragma omp parallel for reduction(+:temp)
            for(m=0;m<L;m++)
            {
                //#pragma omp critical
                //*(e+k+i*K)+=(*(w_star+m)-(*(w_k+L*i+m)))*(*(w_star+m)-(*(w_k+L*i+m)));
                temp=temp+(w_star[m]-w_k[L*i+m])*(w_star[m]-w_k[L*i+m]);
            }
            //#pragma omp critical
            e[k+i*K]=temp;
        }

        /* intermediate estmeate Phi calculation */
        //s=0
        memset (s,0,L*N*sizeof(double));
        //memset (p,0,N*sizeof(double));
        //#pragma omp parallel for  private(j,p,p2,m)
        memcpy(phi,w_k,L*N*sizeof(double));
        //#pragma omp parallel for  private(j,p,p2,m)
        for (i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                //p=u_j'*w_i
                if (IsNonZero(C[j+N*i]))
                {
                    p=0;
                    //#pragma omp parallel for reduction(+:p)
                    for(m=0; m<L ;m++)
                    {
                        p +=u[m+j*L+k*(N*L)]*w_k[L*i+m];
                    }

                    //p=mu(i)*C(j,i)*(d(k,j))-u_j'*w_i)
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
        //#pragma omp parallel for
        for(i=0;i<N;i++)
        {
            //w(:,i)=w(:,i)+A(j,i)*phi(:,j)
            for(j=0;j<N;j++)
            {
            if (IsNonZero(A[j+N*i]))
                {
                    //#pragma omp parallel for shared(A,phi,w_k) private(m)
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



/* compressed diffusion */
void compressed_diffusion(mwSize N, mwSize L, mwSize K, mwSize M, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *e, double *w_k)
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

/* doubly compressed Diffusion */
void doubly_compressed_diffusion(mwSize N, mwSize L, mwSize K, mwSize M, mwSize M_grad, double *w_star, double *w_0, double *A, double *C, double *mu, double *d, double *u, double *H, double *H_grad, double *e, double *w_k)
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

#endif /* diffusion_h */
