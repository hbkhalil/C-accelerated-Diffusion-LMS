function [ e ,w_k ] = ATC_En( A,C,w_0,u,d,mu,w_star,En )


% Please cite the paper referenced below if you use this code.
% @article{harranedoubly,
%   title={Doubly compressed diffusion LMS over adaptive networks},
%   author={Harrane, Ibrahim El Khalil and Flamary, R{\'e}mi and Richard, C{\'e}dric}

%%
%           Adapt Then Combine diffusion (ATC)
%
%   L:   Desired parameter lenght
%   N:   Number of nodes
%   K:   Number of iterations
%
%               Inputs:
%
%  
%               
%   A:      N x N          combination matrix
%   C:      N x N          gradient combination matrix
%   w_0     L x N          initial estimate
%   u       L x N x K      regression vectors 
%   d       K x N          refrence signal        
%   mu      N x 1          step size
%   w_star  L x 1          desired parameter
%   En:     K x N          This parameter allows to enable an agent  
%
%
%
%               Outputs:
%
%   e:      K x N           MSD for each node
%   w_k:    L x N           estimated parameter


end

