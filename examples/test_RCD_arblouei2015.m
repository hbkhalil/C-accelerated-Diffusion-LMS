clear all
close all
addpath(genpath('.'));
addpath(genpath('../../toolbox'));
addpath(genpath('../toolbox_mex'));
%%
rng(356);


L=40;
N=80;
M=1;
N_iter=1e4;
N_exp=100;
step_size=3e-2;
mu=step_size*ones(N,1);
%mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=sigma_u*ones(N,1);
Sigma_v=sigma_v*ones(N,1);
a=zeros(N,1);
%a=[0.1 0.2 0.3];
r=5;
w_star=randn(L,1);
%w_1=zeros(L,N);
w_0=repmat(randn(L,1),1,N);
%file_name=['Diff_Partial_model_grad_' num2str(M)];
%file_name='Diff_Partial_model_grad_perf';

%%
options.M=L;
options.N=N;
options.N_iter=N_iter;
options.parameter=r;
options.corr=0;
options.type='radial';
options.wo=w_star;
options.rule='metropolis';
options.sigma_u=sigma_u;
options.Sigma_u=Sigma_u;
options.sigma_v=sigma_v;
options.one_iter=0;
options.rtype='provided';

%%
[X,L1,Ac,xy] = get_network( options);
print_network( xy,X,options )
[ Coef ] = get_coef( X,Sigma_u,Sigma_v,mu,options );
A=Coef;
%A=eye(N);

C=eye(N);

%%

[ d,u,Sigma_u,Sigma_v,w_star] = get_data( options );
%d=d';
%u=u(:,:,:,1);
%d=d(:,:,1);
u=permute(u,[2 3 1]);


%%
X_t=X-eye(N);
d_k=sum(X_t,2);

H=zeros(N,N,N_iter);
for i=1:N_iter
    H(:,:,i)=X_t;
    for j=1:N
        if(d_k(j)>M)
            idx=find(X_t(j,:)==1);
            id_sel=randperm(length(idx));
            H(j,idx(id_sel(M+1:end)),i)=0;
        end
    end
    H(:,:,i)=H(:,:,i)+eye(N);
end

%%


[e_1,w_k_1]=ATC_RCD(A,w_0,u,d,0.08*mu,H,M,w_star);
MSD=mean(e_1,2);


%%
[e,w_k]=ATC_mex(eye(N),A,w_0,u,d,mu,w_star);


MSD_2=mean(e,2);


%%
% 
% u=permute(u,[3 1 2]);
% 
% err_w=zeros(N_iter,N);
% Psi=zeros(L,N);
% 
% tic
% 
% W=w_0;
% 
% for k=1:N_iter
%     %A_t=H(:,:,k).*A;
%     for i=1:N
%         %Error estimation for diffusion
%         err_w(k,i)=norm(w_star-W(:,i),2)^2;
%         Psi(:,i)=W(:,i)+mu(i)*u(k,:,i)'*(d(k,i) - u(k,:,i)*W(:,i));
%     end
%     
%     
%     W=zeros(L,N);
%     for i=1:N
%         
%         for j=1:N            
%             W(:,i)=W(:,i)+A(j,i)*(H(j,i,k)*Psi(:,j)+(1-H(j,i,k))*Psi(:,i));
%         end
%     end
%     
%     %W=reshape(A2_i*Psi(:),[M N]);
% end
% 
% toc
% 
% MSD_1=mean(err_w,2);
%%
%file_name='Diff_Partial_model_perf';
%Theo=load(['data/' file_name '.mat']);
%plot(10*log10([MSD Theo.epsilon Theo.MSD_t*ones(Theo.N_iter,1)]))
plot(10*log10([MSD MSD_2 ]))
legend('mexc','ATC')