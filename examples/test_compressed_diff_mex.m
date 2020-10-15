clear all
close all
addpath(genpath('.'));
addpath(genpath('/Users/harrane/Dropbox/Diffusion_strategies/toolbox'));

%%
rng(356);


L=5;
N=5;
M=1;
N_iter=2e4;
N_exp=50;
step_size=1e-3;
mu=step_size*ones(N,1);
%mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=sigma_u*ones(N,1);
Sigma_v=sigma_v*ones(N,1);
a=zeros(N,1);
%a=[0.1 0.2 0.3];
r=10;
w_star=randn(L,1);
%w_1=zeros(L,N);
w_0=randn(L,N);
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
options.Sigma_v=Sigma_v;
options.sigma_v=sigma_v;
options.one_iter=0;
options.rtype='provided';
options.ntype='different';
options.sigma_net=50;
options.minAllowableDistance=5;

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
        H(:,i,j)=randperm(L);
    end
end
tic
[e,~]=ATC_mex(eye(N),C,w_0,u,d,mu,w_star);
%[e,w_k]=compressed_diffusion2(A,C,w_0,u,d,mu,H,M,w_star);
toc
MSD_ATC=mean(e,2);

tic
[e,w_k]=compressed_diffusion(eye(N),C,w_0,u,d,mu,H,M,w_star);
toc
MSD=mean(e,2);

%%

% u=permute(u,[3 1 2]);
% 
% err_w=zeros(N_iter,N);
% Psi=zeros(L,N);
% 
% tic
% 
% W=w_0;
% 
% 
% for k=1:N_iter
%     
%     
%     
%     for i=1:N
%         
%         
%         
%         %Error estimation for diffusion
%         err_w(k,i)=norm(w_star-W(:,i),2)^2;
%         idsel=H(:,i,k);
%         
%         
%         J=zeros(L);
%         for m=1:M
%             J(idsel(m,1),idsel(m,1))=1;
%         end
%         
%         s=zeros(L,1);
%         
%         
%         
%         for j=1:N
%             
%             %idsel=H(:,i,k);
% 
%             
%             W_e=J * W(:,i);
%             W_e(idsel(M+1:end))=0;
% 
%             W_l=( eye(L) - J ) * W(:,j);
%             W_l(idsel(1:M))=0;
% 
%             W_p=W_e+W_l;
%             
%             
%             
%             %W_p(idsel(M+1:end))=W_l(idsel(M+1:end));
%             %W_p=J*W(:,i)+(eye(L)-J)*W(:,j);
%             s=s+mu(i)*C(j,i)*u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p);
%         end
%         
%         
%         Psi(:,i)=W(:,i)+s;
%         
%     end
%     
%     %W=Psi;
%     
%     
%     W_temp=W;
%     
%     W=zeros(L,N);
%     
%     for i=1:N
%         %W(:,i)=Psi(:,i);
% 
%         for j=1:N
%             
%             idsel=H(:,j,k);
%             J=zeros(L);
%             
%             for m=1:M
%                 J(idsel(m,1),idsel(m,1))=1;
%             end
%             
%             phi_e=J*W_temp(:,j);
%             phi_i= (eye(L)-J)*Psi(:,i);
%             phi_p=phi_i+phi_e;
%             
%             if (i~=j)
%                 W(:,i)=W(:,i)+A(j,i)*phi_p;
%             else
%                 W(:,i)=W(:,i)+A(j,i)*Psi(:,i);
%             end
%         end
%     end
%     
%     %W=reshape(A2_i*Psi(:),[M N]);
% end
% 
% toc
% 
 %MSD_1=mean(err_w,2);

 %% Mean square stability
% 
% 
% fprintf('Theoritical model \n');
% wt_1=repmat(w_star,1,N)- w_0;
% C_i=kron(C,eye(L));
% R_u=[];
% R=[];
% S=[];
% R_v=[];
% M_mu=[];
% 
% 
% 
% Ru=zeros(L,L,N);
% for k=1:N
%     for i=1:L
%         for j=i:L
%             Ru(i,j,k)=(a(k)^abs(j-i));
%             Ru(j,i,k)=Ru(i,j,k);
%         end
%     end
%     R_u=blkdiag(R_u,Ru(:,:,k));
% end
% 
% 
% 
% for i=1:N
%     R_ut=0;
%     for j=1:N
%         R_ut=R_ut+C(j,i)*Ru(:,:,j);
%     end
%     M_mu=blkdiag(M_mu,mu(i)*eye(L));
%     R=blkdiag(R,R_ut);
%     R_v=blkdiag(R_v,Sigma_v(i));
%     S=blkdiag(S,Sigma_v(i)*Ru(:,:,i));
% end
% 
% 
% M_diag=0;
% H_p=eye(L*N);
% 
% P1_blk=0;
% P1_diag=0;
% 
% P2_blk=0;
% P2_diag=0;
% 
% P3_blk=0;
% P3_diag=0;
% 
% P4_blk=0;
% P4_diag=0;
% 
% for i=1:N*L
%     G=diag(H_p(:,i));
%     
%     P1_diag=P1_diag+kron(G*R*M_mu,G*R*M_mu);
%     
%     
% end
% 
% 
% for i=1:N*L
%     G=diag(H_p(:,i));
%     
% 
%     for j=i:L:L*N
%         Gj=diag(H_p(:,j));
%         P2_diag=P2_diag+kron(G*R_u*C_i*M_mu,Gj*R*M_mu)...`
%             +kron(Gj*R_u*C_i*M_mu,G*R*M_mu);
%         P3_diag=P3_diag+kron(G*R*M_mu,Gj*R_u*C_i*M_mu)...
%             + kron(Gj*R*M_mu,G*R_u*C_i*M_mu);
%         P4_diag=P4_diag+kron(G*R_u*C_i*M_mu,Gj*R_u*C_i*M_mu)...
%             + kron(Gj*R_u*C_i*M_mu,G*R_u*C_i*M_mu);
%     end
%     P2_diag=P2_diag-kron(G*R_u*C_i*M_mu,G*R*M_mu);
%     P3_diag=P3_diag-kron(G*R*M_mu,G*R_u*C_i*M_mu);
%     P4_diag=P4_diag-kron(G*R_u*C_i*M_mu,G*R_u*C_i*M_mu);
% end
% 
% 
% M_blk=0;
% for i=1:N
%     Gt=zeros(N);
%     Gt(i,i)=1;
%     G=kron(Gt,eye(L));
%     
%     P1_blk=P1_blk+kron(G*R*M_mu,G*R*M_mu);
%     P2_blk=P2_blk+kron(G*R_u*C_i*M_mu,G*R*M_mu);
% end
% 
% 
% coef1=(M/L)*((M-1)/(L-1)-(M/L));
% coef2=(M/L)*(1-(M-1)/(L-1));
% coef3=(M/L)^2;
% 
% c1=(M/L)*(1-(M-1)/(L-1));
% c2=(M/L)*(M-1)/(L-1);
% 
% P1= coef1*P1_blk + coef2*P1_diag + coef3*kron(R*M_mu,R*M_mu);
% P2= (M/L - c2) * kron(R_u*C_i*M_mu, R*M_mu) - c1 * P2_diag ;
% P3= (M/L - c2) * kron(R*M_mu, R_u*C_i*M_mu) - c1 * P3_diag ;
% P4= (1 - 2* M/L +c2) * kron(R_u*C_i*M_mu,R_u*C_i*M_mu) + c1 * P4_diag ;
% 
% I_NM=eye(L*N);
% x=M/L;
% 
% 
% F=eye((N*L)^2)-x*kron(R*M_mu,eye(N*L))-(1-x)*kron(R_u*C_i*M_mu,eye(N*L))...
%     -x*kron(I_NM,R*M_mu)...
%     -(1-x)*kron(I_NM,R_u*C_i*M_mu)...
%     +P1+P2+P3+P4;
%     
% 
% 
% F=sparse(F);
% 
% 
% G_i=M_mu*C_i';
% vec_y=reshape(G_i*S'*G_i',1,[]);
% 
% 
% %sig=R_u;
% sig=eye(N*L);
% epsilon=zeros(N_iter,1);
% 
% 
% epsilon_t=(1/N)*wt_1(:)'*reshape(F*sig(:),N*L,N*L)*wt_1(:)+(1/N)*vec_y*sig(:);
% epsilon(1)=epsilon_t;
% 
% Fi=F;
% 
% 
% for i=2:N_iter
%     
%     fprintf('iter n %d  \n',i);
%     
%     Met=(eye((N*L)^2)-F)*Fi*sig(:);
%     Met=reshape(Met,N*L,N*L);
%     
%     epsilon_t=epsilon_t+(1/N)*vec_y*Fi*sig(:)-(1/N)*wt_1(:)'*Met*wt_1(:);
%     epsilon(i)=epsilon_t;
%     
%     Fi=Fi*F;
% end
% 
% MSD_t=(1/N)*vec_y*(eye(size(F))-F)^(-1)*reshape(eye(N*L),[],1);
% %EMSE=(1/N)*vec_y*(eye(size(F))-F)^(-1)*R(:);
%%
%file_name='Diff_Partial_model_perf';
%Theo=load(['data/' file_name '.mat']);
%plot(10*log10([MSD Theo.epsilon Theo.MSD_t*ones(Theo.N_iter,1)]))
plot(10*log10([MSD_ATC MSD MSD]))
legend('ATC','mexc','matlab')