function [ Coef ] = get_coef( A,uv,nv,mu,options )
% This function generates coefficients
%
%       inputs:
%   A:  Topology matrix
%   nv:  noise variance vector for relative variance rule
%   uv:  regression variance vector for relative variance rule
%   nu:  Adaptation step vector (N,1)
%   options.rule:       metropolis, relative variance, uniform (default='uniform')



options=initoptions('get_data',options);
N=max(size(A));

%% Neighbors Number
n=zeros(N,1);

for i=1:N
    n(i)=sum(A(i,:));
end

switch options.rule
    case 'uniform'
        Coef=A;
        for i = 1:N
            Coef(i,:) = Coef(i,:)/sum(Coef(i,:));
        end
        
    case 'metropolis'
        Coef=zeros(N);
        
        for i=1:N
            for j=1:N
                Coef(i,j)=A(i,j)/max(n(i),n(j));
            end
            
            Coef(i,i)=1+Coef(i,i)-sum(Coef(i,:));
        end
        
    case 'rel_var'
        
        gamma=zeros(N,1);
        Coef=zeros(N);
        for i=1:N
            gamma(i)=(mu(i))^2*nv(i)*sum(uv);
            gamma(i)=1./gamma(i);
        end
        
        for i=1:N
            Coef(i,:)=(A(i,:).*gamma')/(A(i,:)*gamma);
        end
end
        
        
        
        
end

