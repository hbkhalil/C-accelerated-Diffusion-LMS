clear all
close all
addpath(genpath('.'));
addpath(genpath('../toolbox'));

%%
rng(356);


L=40;
N=80;
M=3;
N_iter=1e4;
N_exp=100;
step_size=3e-2;
mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=sigma_u*ones(N,1);
Sigma_v=sigma_v*ones(N,1);
a=zeros(N,1);

r=0.075;
w_star=randn(L,1);
w_0=0*randn(L,N);


%%
options.M=L;
options.N=N;
options.N_iter=N_iter;
options.parameter=r;
options.corr=0;
options.a=a;
options.type='radial';
options.wo=w_star;
options.rule='metropolis';
options.sigma_u=sigma_u;
options.sigma_v=sigma_v;

%%
[A,L1,Ac,xy] = get_network( options);
[ Coef ] = get_coef( A,Sigma_u,Sigma_v,mu,options );
C=Coef;
%A=eye(N);
A=C;

%%

[ d,u,Sigma_u,Sigma_v,w_star] = get_data( options );
%d=d';
%u=u(:,:,:,1);
%d=d(:,:,1);
u=permute(u,[2 3 1]);


%%
%mex compressed_diffusion.c

H=zeros(L,N,3);
for i=1:N
    for j=1:N_iter
        H(:,i,j)=randperm(L);
    end
end

%%
tic
%[e,w_k]=ATC_mex2(A,C,w_0,u,d,mu,w_star);
[e,w_k]=ATC_partial_model(A,w_0,u,d,mu,H,M,w_star);
toc
MSD=mean(e,2);



%%
 
u=permute(u,[3 1 2]);

err_w=zeros(N_iter,N);
Psi=zeros(L,N);

tic

W=w_0;

for k=1:N_iter

        for i=1:N
            %Error estimation for diffusion
            err_w(k,i)=norm(w_star-W(:,i),2)^2;
            Psi(:,i)=W(:,i)+mu(i)*u(k,:,i)'*(d(k,i) - u(k,:,i)*W(:,i));

        end



        for i=1:N

            idsel=H(:,i,k);
            
            
            
            Psi_l=Psi(:,i);
            W(:,i)=zeros(L,1);
            
            for j=1:N
                Psi_p=Psi(:,j);
                Psi_p(idsel(M+1:end))=Psi_l(idsel(M+1:end));
                W(:,i)=W(:,i)+A(j,i)*Psi_p;
            end
        end

        %W=reshape(A2_i*Psi(:),[M N]);
end

toc

MSD_1=mean(err_w,2);
%%
file_name='Diff_Partial_model_perf';
%Theo=load(['data/' file_name '.mat']);
%plot(10*log10([MSD Theo.epsilon Theo.MSD_t*ones(Theo.N_iter,1)]))
plot(10*log10([MSD MSD_1]))
legend('mexc','matlab')