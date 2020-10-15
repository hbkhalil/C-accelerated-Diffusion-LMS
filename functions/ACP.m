function [V,u,Vp] = ACP( X,n )
%       ACP
%   Input
%           X:  data matrixe
%           n:  vector number
%           
%    Output
%           V=X*u
%           u:  eigen vectors
%           Vp: eigen values

%%

Sigma=cov(X);
[u Vp]=eig(Sigma); 
V=X*u(:,end:-1:end-(n-1));


end