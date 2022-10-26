function [results, states_save]...
    =Call_Function_TDTHMM_nomean_application_givenstate(data, burnin, nmc, S, v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result, states)

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
%     tau_t_save: p x T x nmc array of tau_ti across iterations
%     alpha_save: 1 x nmc array of alpha in Dirichlet-process across
%     iterations
%     w_kt_save: K_sb x T x nmc array of stick-breaking weights across
%     iterations
%   Initial values:
%     states: 1 x T array of initial hidden states
%     C:  p x p x S array of initial value of precision matrices
%     adj:  p x p x S array of initial value of  adjacency matrices 
%     Sig: p x p x S array of initial value of covariance matrices
%     Phi: S x S array of initial value of the state similarity matrix
%     tau_t: p x T array of initial value of tau_ti
%     alpha: initial value of alpha in the Dirichlet-process
%     w_kt: K_sb x T array of initial value of stick-breaking weights
%     Z_it: p x T array of initial value of the clusters for tau_ti
%     eta_kt: K_sb x T array of initial value of tau_ti in each cluster
%     n_kt: K_sb x T array of initial value of number of tau_ti in each cluster
%   Hyperparameter values:
%     burnin, nmc, S, v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result
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
% nu : nu in the prior of tau 
% a_alpha, b_alpha: parameters in the prior of alpha in the Dirichlet-t model
% K_sb: level of truncation in the Dirichlet-t model
% disp_result : true/false, whether to display simulation process during iterations
% states: T x 1 states given to the model and fixed during MCMC

%%
Y = data.TC;
% p is number of variables; T is the number of time points 
[T, p] = size(Y);

% Mean standardize and subtract out linear and quadratic trends
std_input = mad(Y,1);
mean_input = median(Y);
Y = Y - repmat(mean_input,size(Y,1),1);
Y = Y./repmat(std_input,T,1);


%% Initial hidden states
% states=repelem(1:S, floor(T/S)); % equally divide into S parts
% states=[states, S*ones(1,S+2)];
% states=states(1:T);


%% Initalize covariance and precision matrices 
Sig = zeros(p,p,S);
C = Sig;
for s = 1:S 
    if sum(states == s)>100
        Sig(:,:,s) = robustcov(Y(states == s,:));  % use Minimum Covariance Determinant method to initialize covariance
	    C(:,:,s) = inv(Sig(:,:,s));
    else 
        Sig(:,:,s) = eye(p);
        C(:,:,s) = eye(p);
    end
end
adj = abs(C) > 1e-5;  % initial adjacency matrices as defined by intial precision matrices 
C = C .* adj;


%% Initial value for Phi
Phi = eye(S);

%% Initial tau_t and alpha
alpha = 1;
w_kt = 1/K_sb* ones(K_sb, T); % equal weight for all clusters
tau_t=ones(p,T); % tau_ti are all ones
eta_kt =  gamrnd(nu, 1/nu, [K_sb, T]); % randomly generate eta
Z_it = nan(p,T);  % latent cluster indicator Z
n_kt = zeros(K_sb, T); % number of variables in each cluster
for t = 1:T
    temp = unique(tau_t(:,t));
    eta_kt(1:length(temp), t) = temp;
    for k = 1:length(temp)
        Z_it( tau_t(:, t)== eta_kt(k, t), t) = k;
        n_kt(k, t) = sum(tau_t(:, t)== eta_kt(k, t));
    end
end


%% Run the algorithm
tic
[C_save, adj_save, ppi_HMM, states_save, Phi_save, ar_Phi,...
     log_mh_ratio_Phi_save, Phi_prop_save, tau_t_save, w_kt_save, alpha_save ]=...
     MCMC_Algorithm_TDTHMM_nomean_givenstate(burnin, nmc, Y, ...
     S,...      % number of hidden states
     v0, v1, pii,...        % parameters for the continuous spike-and-slab prior
     nu, a_alpha, b_alpha, K_sb, ...        % parameters for the t-distribution
     Sig, C, adj, Phi,...  % initialization of the linked graphical model
     states,...     % initialization of hidden states
     tau_t, eta_kt, alpha, w_kt,...  % initialization of the t-distribution
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
results.n_kt = [];
results.method = 'TDTHMM';

clearvars -except results states_save
