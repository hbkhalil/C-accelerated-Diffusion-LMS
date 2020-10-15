function [ e ,w_k ] = ATC_RCD(A,w_0,u,d,mu,H,M,w_star)
%
% Please cite the paper refrenced below if you use this code.
% @article{harranedoubly,
%   title={Doubly compressed diffusion LMS over adaptive networks},
%   author={Harrane, Ibrahim El Khalil and Flamary, R{\'e}mi and Richard, C{\'e}dric}
%
%%
%
%           Reduced Communication Diffusion Arablouei 2015
%          
%               Inputs:
%   L:   Desired parameter lenght
%   N:   Number of nodes
%   K:   Number of iterations
%   M:   Number of estimates' shared entries
%  
%               
%   A:      N x N          combination matrix
%   w_0:    L x N          initial estimate
%   u:      L x N x K      regression vectors 
%   d:      K x N          refrence signal        
%   mu:     N x 1          step size
%   H:      N x N x K      scrambled nodes indeces going from 1 to |N_k|:where |N_k| is the number of node k neighbors
%                          We use it to select a subset of nodes to
%                          communicate with. we choose the M first indices.
%   w_star: L x 1          desired parameter
%
%
%
%               Outputs:
%
%   e:      K x N           MSD for each node
%   w_k:    L x N           estimated parameter

end


