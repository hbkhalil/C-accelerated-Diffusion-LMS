function [ e ,w_k ] = compressive_diff_ATC(A,w_0,u,d,mu,H,w_star,eta,delta,mu_cvx)
%
% Please cite the paper refrenced below if you use this code.
% @article{harranedoubly,
%   title={Doubly compressed diffusion LMS over adaptive networks},
%   author={Harrane, Ibrahim El Khalil and Flamary, R{\'e}mi and Richard, C{\'e}dric}
%
%%
%
%
%           Reduced Communication Diffusion Arablouei 2015
%          
%               Inputs:
%   L:   Desired parameter lenght
%   N:   Number of nodes
%   K:   Number of iterations
%   
%  
%               
%   A:      N x N          combination matrix
%   w_0:    L x N          initial estimate
%   u:      L x N x K      regression vectors 
%   d:      K x N          refrence signal        
%   mu:     N x 1          step size
%   H:      L x N x K      A 3D matrix containing projection vectors for N
%                           nodes at every time instance k
%   w_star: L x 1          desired parameter
%   eta:    N x 1          construction step size
%   delta:  N x 1          initial value for the confidence parameter (it is calculated using equation 78)
%   mu_cvx: 1 x 1          confidence parameter calculation step size
%
%
%               Outputs:
%
%   e:      K x N           MSD for each node
%   w_k:    L x N           estimated parameter

end


