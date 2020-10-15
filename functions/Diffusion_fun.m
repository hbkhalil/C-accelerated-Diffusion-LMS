function [ res ] = Diffusion_fun( options )
%       Diffusion Function
%   returns teoretical and Simulated MSD, Estimation error, and model's
%   convergence
%
%       Input:
%   options.A:          w weigths matrix
%   options.C:          gradients' weights matrix
%   options.type        Diffision Type (ATC or CTA)
%   options.Tr:         Gradient transforming (0 or 1)
%   options.M_w:        Transformation matrix rank
%   options.M:          w's dimension
%   options.N_iter:     iteration number
%   options.N_exp:      experimentation number
%   options.mu:         step size vector
%   options.Sim:        Simulation on or off (0 or 1)
%   options.Teo:        Teoretical calculation on or off (0 or 1)
%   options.sigma_g:    Transformation matrix variance
%   options.d:          Input data
%   options.u:          regression vectors
%   options.W_1:        W Initialazion
N=size(A,1);

%% Simulation

Psi= zeros(M,N);                % psi column vectors for all agents
e=zeros(N_iter,N);              % error column vectors for all agents in the ATC network
err=zeros(N_iter,N);
MSD=zeros(N_iter,N);
err_w=zeros(N_iter,N);
W_tot=zeros(N_iter,M);
W_1=options_W_1
wt_1=repmat(wo,1,N)-w_1;
d=options.d;
u=options.u
N_iter=options.N_iter
W=W_1;

switch options.type
    case 'ATC'
        if (options.Sim==1)
            
            for k=1:N_iter
                for i=1:N
                    %Error estimation for diffusion
                    e(k,i)=norm((d(k,i)-u(k,:,i)*W(:,i)),2)^2;
                    err_w(k,i)=norm(wo-W(:,i),2)^2;
                    
                    
                    
                    s=zeros(M,1);
                    s_tg=zeros(M,1);
                    s_W_M=zeros(M,1);
                    
                    for j=1:N
                        s=s+mu(j)*C(j,i)*u(k,:,j)'*(d(k,j) - u(k,:,j)*W(:,i));
                    end
                    
                    
                    Psi(:,i)=W(:,i)+s;
                end
                
                
                W=reshape(A2_i*Psi(:),[M N]);
                
                
                W_t(k,:)=mean(W,2)';
                
            end
        end
end


res.W=W;
res.W_t=W_t;
res.e=e;
res.err_w=err_w;



end

