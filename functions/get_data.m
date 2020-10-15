function [ d,u,Sigma_u,Sigma_v,wo,v ] = get_data( options )
% This function generates noisy data for N node
%       d=u*wo+v
%   inputs
%      options.N:               node number (default=10)
%      options.M:               regression vector dimention (default=2)
%      options.N_iter:          iteration number (default=1e3)
%      options.N_exp:           exprement number (default=50)
%      options.rtype            regression vector type: (default='uniform')
%                               uniform:     R_u=sigma_u*I for all nodes
%                               different:   R_u_i=sigma_u(i)*I
%      options.sigma_u:         regression variance (default=0.1)
%      options.ntype:           noise type: (default='uniform')
%                               uniform:    same noise variance for all nodes
%                               uniform_random: sigma_v is a random value
%                               diffrent:   diffrent variance for each node
%      options.sigma_v:         noise variance (default=0.1)
%      options.wo:              (default=randn(options.M,1))
%      options.one_iter         * 1: Generate data for one iteration (default=0)
%                               * 0: Generate data for all iterations
%                               *-1: Generate data for all iteration and
%                               experiments
%      options.corr             Generate correlated regression vector (default=0)
%      options.a                Correlation parameter (default=0.3)
%      options.sigma_n          Correlation noise variance (default=0.1)
%      options.Nb_tasts         Task number in a multitask setting (defailt=1)
%      options.task_jump        swaping tasks (defailt=0)
%      options.wo1              (default=randn(options.M,options.N))
%      options.Nb_tasks         (default=1)
%      options.N_agent_per_task (default=ceil(options.N/options.Nb_tasks)*ones(options.Nb_tasks,1))
%      options.task_jump        (default=0)

options=initoptions('get_data',options);
N=options.N;
M=options.M;
N_iter=options.N_iter;
N_exp=options.N_exp;
sigma_u=options.sigma_u;
sigma_v=options.sigma_v;
wo=options.wo;
corr=options.corr;
a=options.a;
Nb_tasks=options.Nb_tasks;
N_agent_per_task=options.N_agent_per_task;
%sigma_n=options.sigma_n;


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
if Nb_tasks==1
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
    
else
    
    if (options.task_jump==0)
        u=randn(N_iter,M,N).*repmat(permute(sqrt(Sigma_u),[2,3,1]),N_iter,M);
        v=randn(N_iter,N).*repmat(sqrt(Sigma_v)',N_iter,1);
        
        d=[];
        
        %N_task=ceil(N/Nb_tasks);
        N_task=N_agent_per_task;
        
            u_temp=u(:,:,1:N_task(1));
            v_temp=v(:,1:N_task(1));
            
            d_temp=reshape(permute(u_temp,[1 3 2]),N_iter*N_task(1),M)*wo(:,1)+v_temp(:);
            d_temp=reshape(d_temp,N_iter,N_task(1));
            d=[d d_temp];
        
        for i=2:Nb_tasks
            u_temp=u(:,:,N_task(i-1)+1:N_task(i-1)+N_task(i));
            v_temp=v(:,N_task(i-1)+1:N_task(i-1)+N_task(i));
            
            d_temp=reshape(permute(u_temp,[1 3 2]),N_iter*N_task(i),M)*wo(:,i)+v_temp(:);
            d_temp=reshape(d_temp,N_iter,N_task(i));
            d=[d d_temp];
        end
        
        d=d(:,1:N);
    else
        
        wo1=options.wo1;
        u=randn(N_iter,M,N).*repmat(permute(sqrt(Sigma_u),[2,3,1]),N_iter,M);
        v=randn(N_iter,N).*repmat(sqrt(Sigma_v)',N_iter,1);
        
        d=[];
        
        %N_task=ceil(N/Nb_tasks);
        N_task=N_agent_per_task;
        
        N_per_task=N_task(i)/2;
        
        for i=1:Nb_tasks
            idx=1:Nb_tasks;
            idx(idx==i)=[];
            
            u_temp=u(1:N_iter/2,:,(i-1)*N_task(i)+1:i*N_task(i));
            v_temp=v(1:N_iter/2,(i-1)*N_task(i)+1:i*N_task(i));
            
            d_temp=reshape(permute(u_temp,[1 3 2]),N_iter/2*N_task(i),M)*wo(:,i)+v_temp(:);
            d_temp=reshape(d_temp,N_iter/2,N_task);
            
            
            u_temp_1=u(N_iter/2+1:end,:,(i-1)*N_task(i)+1:i*N_task(i));
            v_temp_1=v(N_iter/2+1:end,(i-1)*N_task(i)+1:i*N_task(i));

            d_temp_2=[];
            for j=1:(Nb_tasks-1)

                d_temp_1=reshape(permute(u_temp_1(:,:,(j-1)*N_per_task+1:j*N_per_task),[1 3 2]),N_iter/2*N_per_task,M)*wo1(:,idx(j))+reshape(v_temp_1(:,(j-1)*N_per_task+1:j*N_per_task),[],1);
                d_temp_1=reshape(d_temp_1,N_iter/2,N_per_task);
                d_temp_2=[d_temp_2 d_temp_1];
            end
            
            
            d_temp_3=[d_temp ; d_temp_2];
            
            d=[d d_temp_3];
            
        end
        
        d=d(:,1:N);
    end
    
end



