%%% Simulation study for p = 10
%%%  Add the paths containing the relevant files
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions
%%% Create folders to save results
mkdir('Synthetic_data_appendix') % save synthetic data from the simulation study
mkdir('Simulation_results_appendix') % save MCMC results 

%%%%%%%%%%%%%%%%%
%% Generate data and plot data, p = 10
%%%%%%%%%%%%%%%%%
Generate_and_plot_data_appendix
% data are saved in the Synthetic_data_appendix folder

%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MCMC: an example
%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Select one of three types of data in the simulation study
data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;
% cmthmm: classical-t data
% 5spikeshmm: slightly contaminated data
% 10spikeshmm: highly contaminated data
data_method = data_method_list(3);
% Load 
rep = 10; % select a replication from 25 replications
data = load(strcat('Synthetic_data_appendix/Synthetic_data_', data_method, '_rep_', num2str(rep)));
data = data.data;

%%% Select one of the three proposed models
model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"]; 
% LPMHMM: dynamic Gaussian graphical model
% CMTHMM: dynamic classical-t graphical model
% TDTMHH: (truncated) dynamic Dirichlet-t graphical model
model_method = model_method_list(3);

%%% Set up parameters
S = 3; % number of hidden states
v0 = 0.02^2; % v0^2 in the continuous spike and slab prior  
h = 50^2; 
v1 = h * v0; % v1^2 in the continuous spike and slab prior  
pii = 3/(size(data.TC,2)-1); % pi in the continuous spike and slab prior  
nu = 3; % nu in the prior of tau
a_alpha = 1; b_alpha = 1;  %  parameters in the prior of alpha in the Dirichlet-t model
K_sb = 7; % level of truncation in the Dirichlet-t model
if model_method == "TDTHMM" 
    burnin=700; % number if burn-ins
    nmc=1000;  % number of number of Monte Carlo samples to save
else
    burnin = 2000;
    nmc = 8000;
end
disp_result = true; % whether to display simulation process during iterations

%%% Run MCMC
results=[]; states_save = [];
rng(12345 + rep, 'twister')
% Examples of how to run the three proposed models
% Run Dynamic Gaussian Graphical Model
% [results, states_save] = Call_Function_LPMHMM_nomean(data, burnin, nmc, S, v0, v1, pii, disp_result);
% Run Dynamic Classical-t Graphical Model
% [results, states_save] = Call_Function_CMTHMM_nomean(data, burnin, nmc, S, v0, v1, pii, nu, disp_result);
% Run Dynamic Dirhchlet-t Graphical Model
% [results, states_save]= Call_Function_TDTHMM_nomean(data, burnin, nmc, S, v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result);

% Example of running Call_functions with model_method as an input
fname = strcat('Call_Function_', model_method, '_nomean');
        if model_method == "TDTHMM"
            [results, states_save] = feval(fname,...
                data, burnin, nmc, S, v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result);
        elseif model_method == "CMTHMM"
            [results, states_save] = feval(fname,...
                 data, burnin, nmc, S, v0, v1, pii, nu, disp_result);          
        else
            [results, states_save] = feval(fname,...
                 data, burnin, nmc, S, v0, v1, pii, disp_result);          
        end

% Save results
save(strcat('Simulation_results_appendix/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)), 'results')
save(strcat('Simulation_results_appendix/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)), 'states_save')


%%%%%%%%%%%%%%%%
%% Analyze MCMC results and compute model performance
%%%%%%%%%%%%%%%%
%% Load data and results
% Load data
rep = 1; 
data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;
data_method = data_method_list(3);
data = load(strcat('Synthetic_data_appendix/Synthetic_data_', data_method, '_rep_', num2str(rep)));
data = data.data;
% Load results
model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"]; 
model_method = model_method_list(3);
load(strcat('Simulation_results_appendix/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)))
load(strcat('Simulation_results_appendix/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)))
results.states_save = states_save;

%% Plot MCMC results
Plot_simulation_results(data, results)
% plots are saved in the Simulation_analysis folder

%% Compute performance measures
plot_disp = true;
rate_disp = true;
[performance_state_match_state_avg, performance_state_match_by_state,...
performance_graph_state_avg, performance_graph_by_state, performance_graph_roc_state_avg] = ...
  Compute_simulation_performance( data, results,  plot_disp, rate_disp);
% performance of states match: true positive rate (tpr), false positive rate (fpr), mcc; 
% performance of graph estimation: tpr, fpr, mcc, auc, frobinous loss (fl),
% fl considering tau (fl_tau)

%% MCMC diagonosis
plot_disp = true;
[geweke_num_edges, geweke_state, prob_of_last_cluster ] = ...
      Analyze_MCMC_results(data, results, plot_disp);
display(geweke_num_edges)  % geweke test of number of edges in each state
display(geweke_state)  % geweke test of number of observations in each state
display(prob_of_last_cluster)  % probability of tau_ti in the Kth cluster


%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate performances of four models from 25 replicates
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three kinds of data: classical-t, slightly contaminated, highly
% contaminated
% Fit four kinds of model: the dynamic Gaussian graphical model, the dynamic classical-t
% graphical model, the Dinamic Dirichlet-t graphical model and the joint
% graphical lasso
% 25 replicates for each data and model combination

%% Run MCMC 
% Three kinds of Bayesian dynamic models, parallel computing:
open('Helper_functions/Simulation_25rep_appendix.m')
% Joint graphical lasso model, R code, parallel computing:
open('Helper_functions/Simulation_JGL_25rep_appendix.R')

%% Analyze results of 25 replications
open('Helper_functions/Simulation_analyze_results_25rep_appendix.m')


    
