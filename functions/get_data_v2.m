function [ d,u,Sigma_u,Sigma_v,wo ] = get_data_v2( options )
% This function generates noisy data for N node
%       d=u*wo+v
%   inputs
%      options.N:       node number (default=10)
%      options.M:       regression vector dimention (default=2)
%      options.N_iter:  iteration number (default=1e3)
%      options.N_exp:   exprement number (default=50)
%      options.rtype    regression vector type: (default='uniform')
%                           uniform:     R_u=sigma_u*I for all nodes
%                           different:   R_u_i=sigma_u(i)*I
%      options.sigma_u:     regression variance (default=0.1)
%      options.ntype:   noise type: (default='uniform')
%                           uniform:    same noise variance for all nodes
%                           uniform_random: sigma_v is a random value
%                           diffrent:   diffrent variance for each node
%      options.sigma_v:     noise variance (default=0.1)
%      options.wo:      (default=randn(options.M,1))
%      options.one_iter     Generate data for one iteration (default=0)
%      options.corr         Generate correlated regression vector (default=0)
%      options.a            Correlation parameter (default=0.3)
%      options.sigma_n      Correlation noise variance (default=0.1)


options=initoptions('get_data',options);
N=options.N;
M=options.M;
N_iter=options.N_iter;
sigma_u=options.sigma_u;
sigma_v=options.sigma_v;
wo=options.wo;
corr=options.corr;
a=options.a;
%sigma_n=options.sigma_n;


switch options.rtype
    case 'uniform'
        Sigma_u=sigma_u*ones(N,1);
        
    case 'different'
        Sigma_u=sigma_u*rand(N,1);
end


switch options.ntype
    case 'uniform'
        Sigma_v=sigma_v*ones(N,1);
        
    case 'diffrent'
        Sigma_v=sigma_v*rand(N,1);
end




% u=zeros(N_iter,M,N,N_exp);
% d=zeros(N_iter,N,N_exp);
% R_uk=zeros(M,M,N);
% R=zeros(M,M,N);

if (corr==0)
    if (options.one_iter==0)
        
        u=randn(N_iter,M,N).*repmat(permute(sqrt(Sigma_u),[2,3,1]),N_iter,M);
        v=randn(N_iter,N).*repmat(sqrt(Sigma_v)',N_iter,1);
        d=reshape(permute(u,[1 3 2]),N_iter*N,M)*wo+v(:);
        d=reshape(d,N_iter,N);
    else
        u=randn(N,M).*repmat(sqrt(Sigma_u),1,M);
        v=randn(N,1).*(sqrt(Sigma_v));
        d=u*wo+v;
    end
else
    u=zeros(N_iter,N,M);
    
    Sigma_n=ones(N,1)-a.^2.*Sigma_u;
    %Sigma_n=options.sigma_n*ones(N,1);
    
    u(:,:,1)=randn(N_iter,N).*repmat(sqrt(Sigma_u'),N_iter,1);
    for i=2:M
        for j=1:N
        u(:,j,i)=a(j)*u(:,j,i-1)+ randn(N_iter,1).*repmat(sqrt(Sigma_n(j)'),N_iter,1);
        end
    end
    v=randn(N_iter,N).*repmat(sqrt(Sigma_v)',N_iter,1);
    d=reshape(u,N_iter*N,M)*wo+v(:);
    d=reshape(d,N_iter,N);
    u=permute(u,[1 3 2]);
    
    
    % for l=1:N_exp
    %     for k=1:N_iter
    %         for i=1:N
    %             u(k,:,i,l)=randn(1,M)*sqrt(Sigma_u(i));
    %             d(k,i,l)=u(k,:,i,l)*wo+randn*sqrt(Sigma_v(i));
    %             R(:,:,i)=R(:,:,i)+u(k,:,i,l)'*u(k,:,i,l)/N_iter;
    %         end
    %     end
    %     R_uk=R_uk+R/N_exp;
    % end
    
    
end

