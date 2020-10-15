function [ Res ] = ATC_theo( options )

A2=options.A2;
C=options.C;
N=options.N;
L=options.M;
mu=options.mu;
Sigma_v=options.Sigma_v;
Sigma_u=options.Sigma_u;
a=options.a;
w_1=options.w_1;
wo=options.wo;
N_iter=options.N_iter;





R=[];
R_u=[];
S=[];
R_v=[];
M_mu=[];



Ru=zeros(L,L,N);
switch options.corr
    case 1
        
        for k=1:N
            for i=1:L
                for j=i:L
                    Ru(i,j,k)=(a(k)^abs(j-i))*Sigma_u(k);
                    Ru(j,i,k)=Ru(i,j,k);
                end
            end
            R_u=blkdiag(R_u,Ru(:,:,k));
            
        end
        
    case 0
        
        for k=1:N
            Ru(:,:,k)=Sigma_u(k)*eye(L);
            R_u=blkdiag(R_u,Ru(:,:,k));
        end
end



for i=1:N
    R_ut=0;
    for j=1:N
        R_ut=R_ut+C(j,i)*Ru(:,:,j);
    end
    M_mu=blkdiag(M_mu,mu(i)*eye(L));
    R=blkdiag(R,R_ut);
    R_v=blkdiag(R_v,Sigma_v(i));
    S=blkdiag(S,Sigma_v(i)*Ru(:,:,i));
end








A2_i=kron(A2,eye(L));
C_i=kron(C,eye(L));

B=A2_i'*(eye(N*L) - M_mu*R);
F=kron(B',B');
F=sparse(F);
G=A2_i*M_mu*C_i';
Y=G*S*G';



sig=eye(N*L);
epsilon=zeros(N_iter,1);

wt_1=repmat(wo,1,N)-w_1;

epsilon_t=(1/N)*wt_1(:)'*reshape(F*sig(:),N*L,N*L)*wt_1(:)+(1/N)*reshape(Y',[],1)'*sig(:);
epsilon(1)=epsilon_t;



Fi=F;


for i=2:N_iter
    % Theoritical model for diffusion
    fprintf('iter= %d  \n',i);
    Met=(eye((N*L)^2)-F)*Fi*sig(:);
    Met=reshape(Met,N*L,N*L);
    
    epsilon_t=epsilon_t+(1/N)*reshape(Y',(N*L)^2,1)'*Fi*sig(:)-(1/N)*wt_1(:)'*Met*wt_1(:);
    epsilon(i)=epsilon_t;
    
    Fi=Fi*F;
    
end

MSD_t=(1/N)*reshape(Y',[],1)'*(eye(size(F))-F)^(-1)*reshape(eye(N*L),[],1);
EMSE=(1/N)*reshape(Y',[],1)'*(eye(size(F))-F)^(-1)*R_u(:);


%% Convergence in mean



E_w=zeros(N_iter,L);
E_w_t=w_1(:);
E_w(1,:)=mean(w_1,2);




for i=2:N_iter
    E_w_t=(eye(N*L)-A2_i*(eye(L*N)-M_mu*R_u)*A2_i)*repmat(wo,N,1)+A2_i*(eye(L*N)-M_mu*R_u)*A2_i*E_w_t;
    E_w_m=mean(reshape(E_w_t,L,N),2);
    E_w(i,:)=E_w_m';
end

Net_cost=(1:N_iter)*sum(sum(C))*2*L';

Res.MSD_t=MSD_t;
Res.EMSE=EMSE;
Res.epsilon=epsilon;
Res.Net_cost=Net_cost;


end




