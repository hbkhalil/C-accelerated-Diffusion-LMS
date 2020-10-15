clear all
close all
addpath(genpath('.'));
addpath(genpath('../toolbox'));

%%
rng(356);


L=10;
N=20;
N_iter=1e4;
N_exp=100;
step_size=1e-2;
mu=(1.6/L)*step_size*ones(N,1);
eta=0.15*ones(N,1);
mu_cvx=10;
%mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=sigma_u*ones(N,1);
Sigma_v=sigma_v*ones(N,1);
a=zeros(N,1);
%a=[0.1 0.2 0.3];
r=0.075;
w_star=randn(L,1);
%w_1=zeros(L,N);
w_0=0*randn(L,N);
%file_name=['Diff_Partial_model_grad_' num2str(M)];
%file_name='Diff_Partial_model_grad_perf';

delta=zeros(N,1);

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

H=zeros(L,N,N_iter);
for i=1:N
    for j=1:N_iter
        H(:,i,j)=randn(L,1);
    end
end

 %
tic
%[e,w_k]=ATC_mex2(A,C,w_0,u,d,mu,w_star);
[e,w_k]=compressive_diff_ATC(A,w_0,u,d,mu,H,w_star,eta,delta,mu_cvx);
toc
MSD=mean(e,2);



%%

u=permute(u,[3 1 2]);

err_w=zeros(N_iter,N);
Phi=zeros(L,N);
Gamma=zeros(L,N);
epsilon=zeros(N,1);

alpha=zeros(N,1);
e=zeros(N,1);
Phi_t=zeros(L,N);

tic

W=w_0;
% 
% for k=1:N_iter
%     
%     for i=1:N
%         Error estimation for diffusion
%         err_w(k,i)=norm(w_star-W(:,i),2)^2;
%         e(i)=(d(k,i) - u(k,:,i)*W(:,i));
%         Phi(:,i)=W(:,i)+mu(i)*u(k,:,i)'*e(i);
%         epsilon(i)=H(:,i,k)'*(Phi(:,i)-Gamma(:,i));
%         
%         
%         
%         
%         Gamma(:,i)=Gamma(:,i)+eta(i)*H(:,i,k)*epsilon(i);              
%         alpha(i)=alpha(i)-mu_cvx*e(i)*u(k,:,i)*(Phi(:,i)-Phi_t(:,i))*delta(i)*(1-delta(i));
%         delta(i)=1/(1+exp(-alpha(i)));
%     end
%     
%     
%     W=zeros(L,N);
%     for i=1:N        
%         Phi_t(:,i)=A(i,i)*Phi(:,i);
%         for j=1:N            
%             if (i ~=j)
%                 Phi_t(:,i)=Phi_t(:,i)+A(j,i)*Gamma(:,j);
%             end
%         end
%         W(:,i)=(1-delta(i))*Phi_t(:,i)+delta(i)*Phi(:,i);
%     end
%     
%     W=reshape(A2_i*Phi(:),[M N]);
% end
% 
% toc
% 
% MSD_1=mean(err_w,2);
%%
file_name='Diff_Partial_model_perf';
%Theo=load(['data/' file_name '.mat']);
%plot(10*log10([MSD Theo.epsilon Theo.MSD_t*ones(Theo.N_iter,1)]))
plot(10*log10([MSD]))
%legend('mexc','matlab')