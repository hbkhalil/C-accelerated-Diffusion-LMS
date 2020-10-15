function [ e ,w_k ] = compressed_diffusion2(A,C,w_0,u,d,mu,H,M,w_star)
%
%           Compressed Adapt Then Combine diffusion (ATC)
%          
%               Inputs:
%   L:   Desired parameter lenght
%   N:   Number of nodes
%   K:   Number of iterations
%   M:   Number of estimates' shared entries
%  
%               
%   A:      N x N          combination matrix
%   C:      N x N          gradient combination matrix
%   w_0:    L x N          initial estimate
%   u:      L x N x K      regression vectors 
%   d:      K x N          refrence signal        
%   mu:     N x 1          step size
%   H:      L x N x K      scrambled index vectors for entry selection
%   w_star: L x 1          desired parameter
%
%
%
%               Outputs:
%
%   e:      K x N           MSD for each node
%   w_k:    L x N           estimated parameter

end




