function [ e ,w_k, err ] = ATC_RMT( A,w_0,mu,K,w_star,options )


% Please cite the paper referenced below if you use this code.


%%
%           Adapt Then Combine diffusion (ATC)
%          
%               Inputs:
%   L:   Desired parameter lenght
%   N:   Number of nodes
%   K:   Number of iterations
%  
%               
%   A:      N x N          combination matrix
%   C:      N x N          gradient combination matrix
%   w_0     L x N          initial estimate
%   mu      N x 1          step size
%
%
%
%               Outputs:
%
%   e:      K x N          MSD for each node
%   w_k:    L x N          estimated parameter
%   u       L x N x K      regression vectors 
%   d       K x N          refrence signal        
%   w_star  L x 1          desired parameter
end

