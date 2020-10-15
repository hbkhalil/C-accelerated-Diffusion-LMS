clear all
close all
addpath(genpath('.'));
addpath(genpath('../../toolbox'));
addpath(genpath('../toolbox_mex'));
%%
rng(356);



L=50;
N=80;
M=30;

M_grad=1;
N_iter=1e4;
N_exp=10;
step_size=3e-2;
mu=step_size*ones(N,1);
%mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=ones(N,1);
Sigma_v=sigma_v*rand(N,1);
a=zeros(N,1);
%a=[0.1 0.2 0.3];
r=10;
w_star=randn(L,1);
%w_1=zeros(L,N);
w_0=repmat(randn(L,1),1,N);
%file_name=['Diff_Partial_model_grad_' num2str(M)];
file_name='Diff_Partial_model_grad_perf';

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
%options.ntype='different';
options.sigma_net=50;
options.minAllowableDistance=5;
%%
[X,L1,Ac,xy] = get_network( options);
figure(1)
print_network( xy,X,options )
[ Coef ] = get_coef( X,Sigma_u,Sigma_v,mu,options );
C=Coef;
%A=eye(N);
A=C;

%%

[ d,u,Sigma_u,Sigma_v,w_star] = get_data( options );
%d=d';

u=u(:,:,:,1);
d=d(:,:,1);
u=permute(u,[2 3 1 4]);

H=zeros(L,N,N_iter);
H_grad=zeros(L,N,N_iter);
for i=1:N
    for j=1:N_iter
        H(:,i,j)=randperm(L);
        H_grad(:,i,j)=randperm(L);
    end
end

%%
%mex doubly_compressed_diffusion.c


% tic
% [e,w_k]=doubly_compressed_diffusion(eye(N),C,w_0,u,d,mu,H,H_grad,M,M_grad,w_star);
% toc
% MSD_2=mean(e,2);


%[e,w_k]=doubly_compressed_diffusion(eye(N),C,w_0,u,d,0.035*mu,H,H_grad,M,M_grad,w_star);

%MSD=mean(e,2);
 %%
 [e,w_k]=ATC_mex(eye(N),C,w_0,u,d,0.18*mu,w_star);
 MSD_2=mean(e,2);
%%
 [e,w_k]=compressed_diffusion(eye(N),C,w_0,u,d,0.16*mu,H,M,w_star);
 MSD_cmp=mean(e,2);
%%
figure(2)
plot(10*log10([ MSD_2  MSD_cmp]))
legend('ATC','CMP' )
 

%%

%plot(10*log10([MSD MSD_2]))
%%

% u=permute(u,[3 1 2]);
% grad=zeros(L,N);
% err_w=zeros(N_iter,N);
% W=w_0;

% 
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
%         %             for j=1:K
%         %                 H(idsel(j,1),idsel(j,1))=1;
%         %             end
% 
%         e(k,i)=norm((d(k,i)-u(k,:,i)*W(:,i)),2)^2;
%         err_w(k,i)=norm(w_star-W(:,i),2)^2;
% 
% 
%         s=zeros(L,1);
%         for j=1:N
%             W_e=W(:,i);
%             W_e(idsel(M+1:end))=0;
% 
%             W_l=W(:,j);
%             W_l(idsel(1:M))=0;
% 
%             W_p=W_e+W_l;
%             %W_p=W(:,i);
% 
%             J=zeros(L);
%             idsel=H_grad(:,j,k);
%             
%             for m=1:M_grad
%                 J(idsel(m,1),idsel(m,1))=1;
%             end
% 
%             grad_e=C(j,i)*(u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p));
% 
%             grad_e(idsel(M_grad+1:end))=0;
%             grad_l=C(j,i) *(u(k,:,i)'*(d(k,i) - u(k,:,i)*W_p));
% 
%             grad_l(idsel(1:M_grad))=0;
% 
%             
%             grad_p=grad_e+grad_l;
%             s=s+mu(j)*grad_p;
% 
%         end
% 
%         
% 
% 
% 
%         W(:,i)=W(:,i)+s;
%     end
% end
% 
% 
% MSD_1=mean(err_w,2);
% 
% err_w=zeros(N_iter,N);
% Psi=zeros(L,N);
% W=w_0;
%

phi=zeros(L,N);
for k=1:N_iter
    
    
    
    for i=1:N
        
        
        
        %Error estimation for diffusion
        err_w(k,i)=norm(w_star-W(:,i),2)^2;
        
        
        
        idsel_w=H(:,i,k);
        J_w=zeros(L);
        for m=1:M
            J_w(idsel_w(m,1),idsel_w(m,1))=1;
        end
        
        s=zeros(L,1);
        
        
        for j=1:N
            
            J=zeros(L);
            idsel=H_grad(:,j,k);
            
            

            
            
            for m=1:M_grad
                J(idsel(m,1),idsel(m,1))=1;
                %J_w(idsel_w(m,1),idsel_w(m,1))=1;
            end
            
            %W_e=W(:,i);
            %W_e=J_w * W(:,i);
            %W_e(idsel_w(M+1:end))=0;

            %W_l=W(:,j);
            %W_l=( eye(L) - J_w ) * W(:,j);
            %W_l(idsel_w(1:M))=0;

            %W_p=W_e+W_l;
            
            W_p=J_w * W(:,i)+( eye(L) - J_w ) * W(:,j);
            
            s=s+mu(j)*C(j,i)*J*(u(k,:,j)'*(d(k,j) - u(k,:,j)*W_p))...
                + mu(j)*C(j,i) * (eye(L) - J )*(u(k,:,i)'*(d(k,i) - u(k,:,i)*W(:,i)));
            
        end
        

        
        phi(:,i)=W(:,i)+s;
        %W(:,i)=phi(:,i);
        
    end
    
    W_temp=phi;
    
    W=zeros(L,N);
    
    for i=1:N

        for j=1:N
            
            idsel_w=H(:,j,k);
            J_w=zeros(L);
            
            for m=1:M
                J_w(idsel_w(m,1),idsel_w(m,1))=1;
            end
            
            phi_e=J_w*W_temp(:,j);
            phi_i= (eye(L)-J_w)*phi(:,i);
            phi_p=phi_i+phi_e;
            
            if (i~=j)
                W(:,i)=W(:,i)+A(j,i)*phi_p;
            else
                W(:,i)=W(:,i)+A(j,i)*phi(:,i);
            end
        end
    end
    
    
    

    
    %W=reshape(A2_i*Psi(:),[M N]);
end
 MSD_1=mean(err_w,2);

%%
% plot(10*log10([ MSD MSD_cmp MSD_2]))
% legend('DCMP','CMP','ATC' )
% %%
% figure(2)
% plot(10*log10([MSD MSD_1 MSD_2]))
% legend('mexc','matlab','Diff' )