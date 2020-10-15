%% Compressed Diffusion LMS
clear all
close all
addpath(genpath('.'));
addpath(genpath('../toolbox'));


%%

x=randn(10,1)';
y=rand(10,1)';


[estimate,nbias,sigma,descriptor]=information(y,y);
estimate