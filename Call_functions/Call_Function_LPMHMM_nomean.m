function [results, states_save]...
    = Call_Function_LPMHMM_nomean(data, burnin, nmc, S, v0, v1, pii, disp_result)

%% Outputs
% states_save:  T x nmc array of hidden states across iterations

% results: structure with the following fields:
%   method: "LPMHMM": dynamic Gaussian graphical model
%           "CMTHMM": dynamic classical-t graphical model
%           "TDTHMM": (truncated) dynamic Dirichlet-t graphical model 
%   running time: in seconds
%   Saved values across MCMC iterations:
%     ppi_HMM : S x T array of posterior probability of each state across
%     time
%     C_save : p x p x S x nmc array of precision matrices across iterations
%     adj_save : p x p x S x nmc array of adjacency matrices across iterations
%     ppi_edges : p x p x S array of posterior probability of each edge
%     Phi_prop_save: S x S x nmc array of the proposed state similarity matrix Phi across iterations
%     Phi_save: S x S x nmc array of the state similarity matrix Phi across iterations
%     ar_Phi: acceptance rate for Phi
%   Initial values:
%     states: 1 x T array of initial hidden states
%     C:  p x p x S array of initial value of precision matrices
%     adj:  p x p x S array of initial value of  adjacency matrices 
%     Sig: p x p x S array of initial value of covariance matrices
%     Phi: S x S array of initial value of the state similarity matrix
%   Hyperparameter values:
%     burnin, nmc, S, v0, v1, pii, disp_result
%     as specified in the input
%     p: number of variables; T: number of time points

%% Inputs
% data : a structure containing the follow elements 
%   data.TC : T x p array of time courses
%   Other fields are not required in this function
% burnin : number of burnin iterations
% nmc : number of Monte Carlo samples to save
% S : number of hidden states
% v0 : v0^2 in the continuous spike and slab prior  
% h : v1 = h*v0 in the continuous spike and slab prior  
% pii : pi in the continuous spike and slab prior  
% disp_result : true/false, whether to display simulation process during iterations


%%
Y = data.TC;
% p is number of variables; T is the number of time points 
[T, p] = size(Y);


%% Initial hidden states
states=repelem(1:S, floor(T/S)); % equally divide into S parts
states=[states, S*ones(1,S+2)];
states=states(1:T);


%% Initalize covariance and precision matrices 
Sig = zeros(p,p,S); 
C = Sig;
for s = 1:S 
    Sig(:,:,s) = robustcov(Y(states == s,:)); % use Minimum Covariance Determinant method to initialize covariance
	C(:,:,s) = inv(Sig(:,:,s));
end
adj = abs(C) > 1e-5;  % initial adjacency matrices as defined by intial precision matrices 
C = C .* adj;


%% Initial value for Phi
Phi = eye(S);


%% Run the algorithm
tic
[C_save, adj_save, ppi_HMM, states_save, Phi_save, ar_Phi,...
    log_mh_ratio_Phi_save, Phi_prop_save]=...
    MCMC_Algorithm_LPMHMM_nomean(burnin, nmc, Y, ...
    S,...    % number of hidden states 
    v0, v1, pii,...   % parameters for the continuous spike-and-slab prior
    Sig, C, adj, Phi,...   % initialization of the linked graphical model
    states,...   % initialization of hidden states
    disp_result); 
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

results.Y=[];
results.data=[];
results.states_save=[];
results.method = 'LPMHMM';

clearvars -except results states_save
