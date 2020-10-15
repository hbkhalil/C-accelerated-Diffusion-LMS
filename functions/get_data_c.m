function [ Lambda, noise, options ] = get_data_c(options)
% This function generates noisy data for N node for mex c RMT
%       d=u*wo+v
%   inputs
%      options.N:           node number (default=700)
%      options.M:           regression vector dimention (default=400)
%      options.N_iter:      iteration number (default=1e3)
%      options.N_exp:       exprement number (default=50)
%      options.rtype        regression vector type: (default='uniform')
%                           uniform:     R_u=sigma_u*I for all nodes
%                           different:   R_u_i=sigma_u(i)*I
%      options.sigma_u:     regression variance (default=1)
%      options.ntype:       noise type: (default='uniform')
%                           uniform:    same noise variance for all nodes
%                           uniform_random: sigma_v is a random value
%                           diffrent:   diffrent variance for each node
%      options.sigma_v:     noise variance (default=0.1)
%      options.wo:          (default=1:options.M)
%      options.one_iter     * 1: Generate data for one iteration (default=1)
%                           * 0: Generate data for all iterations
%                           *-1: Generate data for all iteration and
%                           experiments
%      options.corr         Generate correlated regression vector (default=0)
%      options.a            Correlation parameter (default=0.3)
%      options.sigma_n      Correlation noise variance (default=0.1)
%      options.seed         seed (default=1)
%      options.iter_batch   (default=1)
%      options.det          (default=0)




options=initoptions('get_data_c',options);
N=options.N;
M=options.M;
%N_iter=options.N_iter;
N_iter=options.iter_batch;
N_exp=options.N_exp;
sigma_u=options.sigma_u;
sigma_v=options.sigma_v;
wo=options.wo;
corr=options.corr;
a=options.a;
%C=options.c;
%sigma_n=options.sigma_n;

%rng(options.seed)


switch options.rtype
    case 'uniform'
        Sigma_u=sigma_u*ones(N,1);
        
    case 'different'
        Sigma_u=sigma_u*rand(N,1);
        
    case 'provided'
        Sigma_u=options.Sigma_u;
        
end


switch options.ntype
    case 'uniform'
        Sigma_v=sigma_v*ones(N,1);    
           
    case 'different'
        Sigma_v=sigma_v*rand(N,1);
end




% u=zeros(N_iter,M,N,N_exp);
% d=zeros(N_iter,N,N_exp);
% R_uk=zeros(M,M,N);
% R=zeros(M,M,N);


if N_iter>1

    
    
%     v=randn(N,N_iter).*repmat(Sigma_v,1,N_iter);
%     v=reshape(v,[N,1,N_iter]);
%     u=randn(N,M,N_iter).*repmat(sqrt(Sigma_u),1,M,N_iter);
%     noise=sum(repmat(v,[1,M,1]).*u,1);
%     noise=reshape(noise,M,N_iter);
%     Lambda=zeros(M,N_iter);
%     
%     R=[];
%     for i=1:N_iter
%         R=blkdiag(R,u(:,:,i)'*u(:,:,i));
%     end
%     
%     [Q, Lambda_t]=eig(R,'vector');
%     
%     for i=1:N_iter
%         idx=(Q(i,:)~=0);
%         Qt=Q(:,idx');
%         Qt(Qt==0)=[];
%         Qt=reshape(Qt,[M M]);
%         noise(:,i)=Qt'*noise(:,i);
%         Lambda(:,i)=Lambda_t(idx);      
%     end
    
    




    
    v=randn(N,N_iter).*repmat(sqrt(Sigma_v),1,N_iter);
    v=reshape(v,[N,1,N_iter]);
    u=randn(N,M,N_iter).*repmat(sqrt(Sigma_u),1,M,N_iter);
    noise=sum(repmat(v,[1,M,1]).*u,1);
    noise=reshape(noise,M,N_iter);
    Lambda=zeros(M,N_iter);
      
    for i=1:N_iter
        U=1/N*u(:,:,i)'*u(:,:,i);
        [Q, Lambda_t]=eig(U,'vector');
        %Lambda(:,i)=sort(Lambda_t,'descend');
        Lambda(:,i)=Lambda_t;
        noise(:,i)=Q'*noise(:,i);
    end
    
    
else

if (corr==0)
    switch options.one_iter
        
        
        case -1
            
            u=randn(N_iter,M,N,N_exp).*permute(repmat(permute(sqrt(Sigma_u),[2,3,4,1]),N_iter,M,N_exp),[1,2,4,3]);
            v=randn(N_iter,N,N_exp).*permute(repmat(permute(sqrt(Sigma_v),[2,3,1]),N_iter,N_exp,1),[1,3,2]);
            d=reshape(permute(u,[2 1 3 4]),M,N_iter*N*N_exp)'*wo+v(:);
            d=reshape(d,[N_iter,N, N_exp]);
            
        case 0
            
            u=randn(N_iter,M,N).*repmat(permute(sqrt(Sigma_u),[2,3,1]),N_iter,M);
            v=randn(N_iter,N).*repmat(sqrt(Sigma_v)',N_iter,1);
            d=reshape(permute(u,[1 3 2]),N_iter*N,M)*wo+v(:);
            d=reshape(d,N_iter,N);       
            
        case 1
            u=randn(N,M).*repmat(sqrt(Sigma_u),1,M);
            v=randn(N,1).*(sqrt(Sigma_v));
            d=u*wo+v;
            R= 1/N*u'*u;
            [Q, Lambda_t]=eig(R,'vector');
            Lambda=sort(Lambda_t,'descend');
            noise=Q'*sum(repmat(v',[M,1]).*u',2);
    end
else
    u=zeros(N_iter,N,M);
    
    Sigma_n=ones(N,1)-a^2*sigma_u*ones(N,1);
    %Sigma_n=options.sigma_n*ones(N,1);
    
    u(:,:,1)=randn(N_iter,N).*repmat(sqrt(Sigma_u'),N_iter,1);
    for i=2:M
        u(:,:,i)=a*u(:,:,i-1)+ randn(N_iter,N).*repmat(sqrt(Sigma_n'),N_iter,1);
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

end
options.seed= options.seed +1;
end

