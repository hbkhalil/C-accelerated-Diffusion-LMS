function [ e ,w_k ] = doubly_compressed_diffusion(A,C,w_0,u,d,mu,H,H_grad,M,M_grad,w_star)
%
%
% Please cite the paper referenced below if you use this code.
% @article{harranedoubly,
%   title={Doubly compressed diffusion LMS over adaptive networks},
%   author={Harrane, Ibrahim El Khalil and Flamary, R{\'e}mi and Richard, C{\'e}dric}
%
%%
%
%
%           Compressed Adapt Then Combine diffusion (ATC)
%          
%                   Inputs:
%
%   L:          Desired parameter lenght
%   N:          Number of nodes
%   K:          Number of iterations
%   M:          Number of estimates' shared entries
%   M_grad:     Number of gradients' shared entries
%  
%               
%   A:          N x N          combination matrix
%   C:          N x N          gradient combination matrix
%   w_0:        L x N          initial estimate
%   u:          L x N x K      regression vectors 
%   d:          K x N          refrence signal        
%   mu:         N x 1          step size
%   H:          L x N x K      scrambled index vectors for estimates' entries selection
%   H_grad:     L x N x K      scrambled index vectors for the gradients' entries selection
%   w_star:     L x 1          desired parameter
%
%
%
%               Outputs:
%
%   e:      K x N           MSD for each node
%   w_k:    L x N           estimated parameter

end





