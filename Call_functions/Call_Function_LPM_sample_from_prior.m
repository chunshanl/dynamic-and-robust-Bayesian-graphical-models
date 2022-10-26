function [results]...
    =Call_Function_LPM_sample_from_prior(burnin, nmc, S, v0, v1, pii, p)

%% Sample precision matrices, adjacency matrices and Phi from the linked spike-and-slab prior
%% in order to compute prior edge inclusion probabilities

%% Outputs
% results: structure with the following fields:
%   method: 'LPM_sample_from_prior'
%   running time: in seconds
%   Saved values across MCMC iterations:
%     C_save : p x p x S x nmc array of precision matrices across iterations
%     adj_save : p x p x S x nmc array of adjacency matrices across iterations
%     Sig_save: p x p x S x nmc array of covariance matrices across iterations
%     ppi_edges : p x p x S array of posterior probability of each edge
%     Phi_prop_save: S x S x nmc array of the proposed state similarity matrix Phi across iterations
%     Phi_save: S x S x nmc array of the state similarity matrix Phi across iterations
%     ar_Phi: acceptance rate for Phi
%   Initial values:
%     C:  p x p x S array of initial value of precision matrices
%     adj:  p x p x S array of initial value of  adjacency matrices 
%     Sig: p x p x S array of initial value of covariance matrices
%     Phi: S x S array of initial value of the state similarity matrix
%   Hyperparameter values:
%     burnin, nmc, S, v0, v1, pii, p
%     as specified in the input

%% Inputs
% burnin : number of burnin iterations
% nmc : number of Monte Carlo samples to save
% S : number of hidden states
% v0 : v0^2 in the continuous spike and slab prior  
% v1 : v1 in the continuous spike and slab prior  
% pii : pi in the continuous spike and slab prior  
% p: number of variables

%%
Sig = zeros(p,p,S);
C = Sig;
for s = 1:S 
    Sig(:,:,s) = eye(p);
	C(:,:,s) = eye(p);
end
adj = abs(C) > 1e-5;
C = C .* adj;

%% Initial value for Phi
Phi = eye(S);


%% Run the algorithm
tic
[C_save, Sig_save, adj_save, Phi_save, ar_Phi,...
    log_mh_ratio_Phi_save, Phi_prop_save] = ...
    MCMC_Algorithm_LPM_sample_from_prior(burnin, nmc, ...
    S,...   % number of hidden states 
    v0, v1, pii,...   % parameters for the continuous spike-and-slab prior
    Sig, C, adj, Phi, ...   % initialization of the linked graphical model
    p);  % number of variables
running_time=toc;

%% 
% Edge PPIs for each graph
ppi_edges = mean(adj_save, 4);

% Get 95% credible intervals for omega (precision matrix)
CI_omega_lower = quantile(C_save, 0.025, 4);
CI_omega_upper = quantile(C_save, 0.975, 4);

% Create struct containing the results
w = whos;
for a = 1:length(w) 
    results.(w(a).name) = eval(w(a).name); 
end

results.method = 'LPM_sample_from_prior';

clearvars -except results
