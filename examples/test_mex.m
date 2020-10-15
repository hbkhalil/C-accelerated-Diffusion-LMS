clear all
close all
addpath(genpath('.'));
addpath(genpath('/Users/harrane/Dropbox/Diffusion_strategies/toolbox'));

%%
rng(360);


L=400;
N=700;

N_iter=5e2;
N_exp=50 ;
step_size=1e-3;
mu=step_size*ones(N,1);
%mu=step_size*ones(N,1);
sigma_u=1e0;
sigma_v=1e-3;
Sigma_u=sigma_u*ones(N,1);
Sigma_v=sigma_v*ones(N,1);
a=zeros(N,1);
%a=[0.1 0.2 0.3];
r=10000;
w_star=100*randn(L,1);
%w_1=zeros(L,N);
w_0=repmat(randn(L,1),1,N);
%file_name=['Diff_Partial_model_grad_' num2str(M)];
path='/Users/harrane/Documents/RMT_for_Diffusion_LMS'; 
file_name=sprintf('%s/N_%d_L_%d_unif',path,N,L); 

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
options.one_iter=1;
options.parameter=5;
options.iter_batch=1e2;
%%
[A,L1,Ac,xy] = get_network( options);
[ Coef ] = get_coef( A,Sigma_u,Sigma_v,mu,options );
C=Coef;
%C=eye(N);
%A=1/N*ones(N);
%A=eye(N);
A=Coef;
%%

% [ d,u,Sigma_u,Sigma_v,w_star,v] = get_data( options );
% %d=d';
% u=permute(u,[2 3 1]);
% tic
%  [ Lambda,noise,options ] = get_data_c(options);
% toc



%%
%%plot(10*log10([MSD_mean]))
%legend('matlab','mexc (new)')

%% Theo
x=randn(L,N); % Normal distribution
c=L/N;
n=50;

s=std(x(:));


% spectral matrix
r= x*x'/N;

% Boundaries 
a=(s^2)*(1-sqrt(c))^2;
b=(s^2)*(1+sqrt(c))^2;

mu=1/(2*(a+b))*ones(N,1);
%mu=0.1*ones(N,1);
% boundaries for R*R'
a_s=a^2;
b_s=b^2;

%% eigenvalues of R
[Q,l]=eig(r);
l=diag(l);

[f,lambda]=hist(l,linspace(a,b,n));
%lambda_1=ones(1,n)-mu*lambda;

% Normalization
f=f/sum(f);

% Theoretical pdf
ft=@(lambda,a,b,c) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
F_R=ft(lambda,a,b,c);

% Processing numerical pdf
F_R=F_R/sum(F_R);
F_R(isnan(F_R))=0;

%% eigenvalues of R*R'

[Q_s,l_s]=eig(r*r);
l_s=diag(l_s);


[f_s,lambda_s]=hist(l_s,linspace(a_s,b_s,n));


% Normalization
f_s=f_s/sum(f_s);

% Theoretical pdf
ft=@(lambda,a,b,c) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
F_R=ft(lambda,a,b,c);

% Processing numerical pdf
F_R=F_R/sum(F_R);
F_R(isnan(F_R))=0;



Fs=ft(sqrt(lambda_s),a,b,c)*2./(sqrt(lambda_s));

% Processing numerical pdf
Fs=Fs/sum(Fs);
Fs(isnan(Fs))=0;
%%
figure(1);
clf
h=bar(lambda,f);
set(h,'FaceColor',[.75 .75 .8]);
set(h,'LineWidth',0.25);
xlabel('Eigenvalues \lambda of (I_N-MR_i) for \mu=\mu_{opt}');
ylabel(' Probability Density Function f(\lambda) ');
title(' Marchenko-Pastur distribution ');
%axis([(1-mu*b) (1-mu*a) 0 1.1*max(f)])

%xlim([0 1.5])

hold on;
plot(lambda,F_R,'g','LineWidth',2);
hold off;

%%
lambda_mean = sort(rand_gen(F_R,lambda,L),'descend');
lambda_mean_sq =sort(rand_gen(Fs,lambda_s,L),'descend');

R_t= diag(lambda_mean) ;
R_s= diag(lambda_mean_sq) ;

Metric= eye(L)- 2*mu(1)*R_t + mu(1)^2*R_s;
Metric_v=Metric(:);

F=Metric;
F=sparse(F);

epsilon=zeros(N_iter,1);

wt_1=(w_star- w_0(:,1));

epsilon_t=wt_1(:)'*F*wt_1(:);%+reshape(Y',[],1)'*sig(:);
epsilon(1)=epsilon_t;

%%
%tic
e=zeros(N_iter,N);
for k=1:N_exp
[ e_t ,w_k, err] = ATC_RMT( A,w_0,mu,N_iter,w_star,options );
e=e+e_t;
end
e=e/N_exp;
MSD=mean(e,2);
save('MSD_sim','MSD');
%%

fprintf('\n \n theo \n')
B=repmat(lambda_mean',1,N);
B_noise=sigma_v*B;
%load('MSD_mean_1.mat');
%load('err_1.mat')
[ MSD_mean_1 ,err_1 ] = ATC_RMT_theo( A,repmat(wt_1,1,N),mu,N_iter,B,B_noise,options);
%%%
%load('sim_200_50.mat')
save(file_name,'MSD_plot','err_mean','err')
%%
clf
figure(1)
set(0, 'defaultTextInterpreter', 'latex')
fs=20; 
h=plot(10*log10([ mean(MSD_mean_1,2) MSD  MSD_mean_1(end)*ones(N_iter,1) ]));
l=legend('Theoretical MSD','Simulated MSD','Steady-state MSD');
set(l,'Interpreter', 'latex','FontSize',0.8*fs) 
title('Mean square deviation','FontSize',1.2*fs)
xlabel('iterations i')
ylabel('MSD ( dB )')
set(gca,'FontSize',0.8*fs)
% set(gca,'Ytick',-10:10:90);
% set(gca,'Xtick',0:20:200);
set(h(1:end),'Linewidth',0.8);
set(h(1),'LineStyle','-.');
set(h(2),'LineStyle','-');
set(h(3),'LineStyle','--');
set(gca,'Ytick',-60:10:80);
set(gca,'Xtick',0:100:N_iter);
grid
print('MSD_50_200','-depsc')


%%
figure(2)
err_plot=[mean(err,2) mean(err_1,2)  ];
plot(err_plot )
legend('sim', 'theo')

title('convergence rate')
xlabel('iteration i')
ylabel('error')
axis([0 2e2 -10 90])












%%

% Met=reshape(eye(L),[],1);
% Metric= (eye(L)- mu(1)*R_t)*(eye(L)- mu(1)*R_t);
% Metric_v=Metric(:);
% 
% err_mean=repmat(Q'*wt_1,1,N);
% MSD_mean=zeros(N_iter,N);
% %MSD_mean(1,:)=err_mean(:,1)'*err_mean(:,1);
% err=zeros(N_iter,N) ; 
% 
% 
% 
% 
% 
% for i=1: N_iter
%     
%    % fprintf('iter  n %d \n',i);
%     Met=Met+Metric_v(:).^(i-1);
%     err_mean_t=err_mean;
%     
%     for j=1:N 
%         err_mean(:,j)=zeros(L,1);
%         for k=1:N
%             err_mean(:,j)= err_mean(:,j) + A(k,j)*(err_mean_t(:,k).*(ones(L,1)-mu(k)*lambda_mean'));%- mu*1/N_k*sum(repmat(v',[L,1]).*x,2)';
%         end
%         err(i,j)=mean(abs(err_mean(:,j)));
%         
%         
%         for k=1:N
%             for n=1:N
%             MSD_mean(i,j)= MSD_mean(i,j) +A(k,j)*A(n,j)*reshape(err_mean_t(:,k)*err_mean_t(:,n)', 1, []  )*Metric_v(:)  + A(k,j)*A(n,j)*mu(j)^2/N *sigma_v* reshape(diag(lambda_mean),1,[])*Met;   
%             end
%         end
%     end
%    
% 
%     
% end

%%



%%

% %mex ATC_mex.c
% tic
% [e2,w_k2]=ATC_mex2(A,C,w_0,u,d,mu,w_star);
% toc
% MSD_2=mean(e2,2);

%mex ATC_mex.c
% tic
% [e,w_k]=ATC_mex(A,C,w_0,u,d,mu,w_star);
% toc
% MSD=mean(e,2);

%mex ATC_mex.c
% tic
% [e_omp,w_k_omp]=ATC_mex_omp(A,C,w_0,u,d,mu,w_star);
% toc
% MSD_omp=mean(e_omp,2);

%plot(10*log10([MSD MSD_2 MSD_omp]))
%legend('mexc (new)','mexc (old)','mex (omp)')
%

% u=permute(u,[3 1 2]);
% 
% tic
%  
% err_w=zeros(N_iter,N);
% Psi=zeros(L,N);
% W=w_0;
% options.seed=1;
% 
% for k=1:N_iter
%         
%         
%         if ( mod ((k-1),options.iter_batch) ==0 )
%         [ Lambda,noise,options ] = get_data_c( options );
%         end
%         
%         for i=1:N
%             
%             
%             
%             %Error estimation for diffusion
%             err_w(k,i)=norm(w_star-W(:,i),2)^2;
%             
% 
%             
%             s=zeros(L,1);
% 
%             
% %             for j=1:N
% %                 %s=s+mu(i)*C(j,i)*u(k,:,j)'*(d(k,j) -
% %                 u(k,:,j)*W(:,i));
% %             end
%                 
%             
%             %Psi(:,i)=W(:,i)+s; 
% 
% %            uk=reshape(u(k,:,:),L,N).*repmat(sqrt(C(:,i)'),[L 1]);
%           %  uk_n= reshape(u(k,:,:),L,N).*repmat(C(:,i)',[L 1]);
%             %Psi(:,i)= W(:,i) + mu(i)*(uk*uk')*(w_star -W(:,i)) + mu(i)*sum(repmat(v(k,:),[L,1]).*uk_n,2);
%             Psi(:,i)= W(:,i) + mu(i)*1/N*Lambda(:,mod ((k-1),options.iter_batch)+1).*(w_star -W(:,i)) + 1/N*mu(i)*noise(:,mod ((k-1),options.iter_batch)+1);
%         end
%         
%         
%         
%         for i=1:N
%             
%             W(:,i)=zeros(L,1);
%             
%             for j=1:N
%                 W(:,i)=W(:,i)+A(j,i)*Psi(:,j);
%             end
%         end
%         
%         %W=reshape(A2_i*Psi(:),[M N]);  
%         fprintf('iter n %d \n',k)
% end
% 
% toc
% 
% MSD_1=mean(err_w,2);

