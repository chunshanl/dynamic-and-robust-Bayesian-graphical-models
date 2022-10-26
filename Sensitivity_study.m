%% Sensitivity study
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%%%%%%%%%%%%%%%%%%
%% Run MCMC on a grid of parameter settings, parallel computing
%%%%%%%%%%%%%%%%%%%
open('Sensitivity_study_parallel')
% Sensitivity analysis of three proposed dynamic Bayesian graphical models  
% Use slightly contaminated data
% Run MCMC on a grid of parameter settings
% parallel computing on 13 cores


%%%%%%%%%%%%%%%%%%%
%% Load 
%%%%%%%%%%%%%%%%%%%
mkdir('Synthetic_data') % save synthetic data from the simulation study
Generate_and_plot_data
data = load('Synthetic_data/Synthetic_data_5spikeshmm_rep_1');
data = data.data;

p = 20;
pii_list = [2:5] / (p-1);
nu_list = [3, 6, 9];
K_sb_list = [4, 7, 10];

%% Sensitivity of the dynamic Gaussian graphical model - not shown in the paper
%%% pi = 2/(p-1),3/(p-1),4/(p-1),5/(p-1)
%%% Load MCMC results and compute performance rates
% model_method = "LPMHMM";
% perf_gaussian_state = zeros(1, length(pii_list));
% perf_gaussian_graph = zeros(2, length(pii_list));
% for rep = 1:length(pii_list)
%     % load results
%     load(strcat('Sensitivity_results/sens_results_',model_method,...
%                 '_pii', num2str(rep)), 'results')
%     load(strcat('Sensitivity_results/sens_results_',model_method,...
%         '_pii', num2str(rep), '_states'))
%     results.states_save = states_save;
%     states_save = [];
%     % get mcc for state estimation and tpr and fpr for graph estimation
%     plot_disp = false;
%     rate_disp = false;
%     [performance_state_match_state_avg, ~,...
%         performance_graph_state_avg, ~, ~] = ...
%         Compute_simulation_performance( data, results,  plot_disp, rate_disp);
%     perf_gaussian_state(rep) = performance_state_match_state_avg(3);
%     perf_gaussian_graph(1, rep) = performance_graph_state_avg(1);
%     perf_gaussian_graph(2, rep) = performance_graph_state_avg(2);   
% end
% save('Sensitivity_results/perf_gaussian_state','perf_gaussian_state');
% save('Sensitivity_results/perf_gaussian_graph','perf_gaussian_graph');

%%% Display performance rates
load('Sensitivity_results/perf_gaussian_state')
load('Sensitivity_results/perf_gaussian_graph')
% mcc for state estimation under different pi = 2/(p-1),3/(p-1),4/(p-1),5/(p-1)
disp(perf_gaussian_state);
% tpr and fpr for graph estimation under different pi = 2/(p-1),3/(p-1),4/(p-1),5/(p-1)
disp(perf_gaussian_graph);


%% Sensitivity of the dynamic classical-t graphical model - generate Table 3
%%% pi = 2/(p-1),3/(p-1),4/(p-1),5/(p-1)
%%% nu = 3, 6, 9
%%% Load MCMC results and compute performance rates
% model_method = "CMTHMM";
% perf_classical_state_tpr = zeros(length(nu_list), length(pii_list));
% perf_classical_state_fpr = perf_classical_state_tpr;
% perf_classical_state_mcc = perf_classical_state_tpr;
% perf_classical_graph_tpr = zeros(length(nu_list), length(pii_list));
% perf_classical_graph_fpr = perf_classical_graph_tpr ;
% perf_classical_graph_mcc = perf_classical_graph_tpr;
% for ii = 1:length(nu_list)
%     nu = nu_list(ii);
%     for rep = 1:length(pii_list)
%         % load results
%         load(strcat('Sensitivity_results/sens_results_',model_method, '_nu', num2str(nu) ,...
%             '_pii', num2str(rep)))
%         load(strcat('Sensitivity_results/sens_results_',model_method,'_nu', num2str(nu) ,...
%             '_pii', num2str(rep), '_states'))
%         results.states_save = states_save;
%         states_save = [];
%         % get tpr, fpr and mcc of state estimation and graph estimation
%         plot_disp = false;
%         rate_disp = false;
%         [performance_state_match_state_avg, ~,...
%             performance_graph_state_avg, ~, ~] = ...
%             Compute_simulation_performance( data, results,  plot_disp, rate_disp);
%         perf_classical_state_tpr(ii, rep) = performance_state_match_state_avg(1);
%         perf_classical_state_fpr(ii, rep) = performance_state_match_state_avg(2);
%         perf_classical_state_mcc(ii, rep) = performance_state_match_state_avg(3);
%         perf_classical_graph_tpr(ii, rep) = performance_graph_state_avg(1);
%         perf_classical_graph_fpr(ii, rep) = performance_graph_state_avg(2);
%         perf_classical_graph_mcc(ii, rep) = performance_graph_state_avg(3);
%     end
% end
% save('Sensitivity_results/perf_classical_state_tpr', 'perf_classical_state_tpr')
% save('Sensitivity_results/perf_classical_state_fpr', 'perf_classical_state_fpr')
% save('Sensitivity_results/perf_classical_state_mcc', 'perf_classical_state_mcc')
% save('Sensitivity_results/perf_classical_graph_tpr', 'perf_classical_graph_tpr')
% save('Sensitivity_results/perf_classical_graph_fpr', 'perf_classical_graph_fpr')
% save('Sensitivity_results/perf_classical_graph_mcc', 'perf_classical_graph_mcc')

%%% Display performance rates
load('Sensitivity_results/perf_classical_state_tpr')
load('Sensitivity_results/perf_classical_state_fpr')
load('Sensitivity_results/perf_classical_state_mcc')
load('Sensitivity_results/perf_classical_graph_tpr')
load('Sensitivity_results/perf_classical_graph_fpr')
load('Sensitivity_results/perf_classical_graph_mcc')
% Performance of graph estimation - Table  3
a = perf_classical_graph_tpr;
a = a(:)';
b = perf_classical_graph_fpr;
b = b(:)';
c = perf_classical_graph_mcc;
c = c(:)';
input.data = [a; b; c];
input.dataForma = {'%.2f'};
latext = latexTable(input);
% Performance of state estimation, not shown in the paper
% rows: nu = 3, 6, 9; columns: pi = 2/(p-1),...,5/(p-1)
disp(perf_classical_state_tpr)
disp(perf_classical_state_fpr)
disp(perf_classical_state_mcc)


%% Sensitivity of the dynamic Dirichlet-t graphical model - Table 3
%%% Load MCMC results and compute performance rates
%%% pi = 2/(p-1),...,5/(p-1)
%%% nu = 3, 6, 9
%%% K_sb = 4, 7, 10
% model_method = "TDTHMM";
% perf_dirichlet_state_tpr = zeros(length(nu_list), length(pii_list),  length(K_sb_list));
% perf_dirichlet_state_fpr = perf_dirichlet_state_tpr;
% perf_dirichlet_state_mcc = perf_dirichlet_state_tpr;
% perf_dirichlet_graph_tpr = zeros(length(nu_list), length(pii_list), length(K_sb_list));
% perf_dirichlet_graph_fpr = perf_dirichlet_graph_tpr ;
% perf_dirichlet_graph_mcc = perf_dirichlet_graph_tpr ;
% for jj = 1:length(K_sb_list)
%     K_sb = K_sb_list(jj);
%     for ii = 1:length(nu_list)
%         nu = nu_list(ii);
%         for rep = 1:length(pii_list)
%             % load results
%             load(strcat('Sensitivity_results/sens_results_',model_method, '_nu', num2str(nu) ,...
%                 '_K', num2str(K_sb),'_pii', num2str(rep)))
%             load(strcat('Sensitivity_results/sens_results_',model_method,'_nu', num2str(nu) ,...
%                 '_K', num2str(K_sb),'_pii', num2str(rep), '_states'))
%             results.states_save = states_save;
%             states_save = [];
%             % get mcc of graph estimation
%             plot_disp = false;
%             rate_disp = false;
%             [performance_state_match_state_avg, ~,...
%                 performance_graph_state_avg, ~, ~] = ...
%                 Compute_simulation_performance( data, results,  plot_disp, rate_disp);
%             perf_dirichlet_state_tpr(ii, rep, jj) = performance_state_match_state_avg(1);
%             perf_dirichlet_state_fpr(ii, rep, jj) = performance_state_match_state_avg(2);
%             perf_dirichlet_state_mcc(ii, rep, jj) = performance_state_match_state_avg(3);
%             perf_dirichlet_graph_tpr(ii, rep, jj) = performance_graph_state_avg(1);
%             perf_dirichlet_graph_fpr(ii, rep, jj) = performance_graph_state_avg(2);
%             perf_dirichlet_graph_mcc(ii, rep, jj) = performance_graph_state_avg(3);
%         end
%     end
% end
% save('Sensitivity_results/perf_dirichlet_state_tpr', 'perf_dirichlet_state_tpr')
% save('Sensitivity_results/perf_dirichlet_state_fpr', 'perf_dirichlet_state_fpr')
% save('Sensitivity_results/perf_dirichlet_state_mcc', 'perf_dirichlet_state_mcc')
% save('Sensitivity_results/perf_dirichlet_graph_tpr', 'perf_dirichlet_graph_tpr')
% save('Sensitivity_results/perf_dirichlet_graph_fpr', 'perf_dirichlet_graph_fpr')
% save('Sensitivity_results/perf_dirichlet_graph_mcc', 'perf_dirichlet_graph_mcc')

load('Sensitivity_results/perf_dirichlet_state_tpr')
load('Sensitivity_results/perf_dirichlet_state_fpr')
load('Sensitivity_results/perf_dirichlet_state_mcc')
load('Sensitivity_results/perf_dirichlet_graph_tpr')
load('Sensitivity_results/perf_dirichlet_graph_fpr')
load('Sensitivity_results/perf_dirichlet_graph_mcc')
jj = 2; % for each K_sb
% performance of graph estimation - Table 3
a = perf_dirichlet_graph_tpr(:,:,jj);
a = a(:)';
b = perf_dirichlet_graph_fpr(:,:,jj);
b = b(:)';
c = perf_dirichlet_graph_mcc(:,:,jj);
c = c(:)';
input.data = [a; b; c];
input.dataForma = {'%.2f'};
latext = latexTable(input);
% performance of state estimation - not shown in the paper
% rows: nu = 3, 6, 9; columns: pi = 2/(p-1),...,5/(p-1)
disp(perf_dirichlet_state_tpr(:,:,jj))
disp(perf_dirichlet_state_fpr(:,:,jj))
disp(perf_dirichlet_state_mcc(:,:,jj))


%%%%%%%%%%%%%%%%%%%
%% Compute prior edge inclusion probabilities for different v0 and pi
%%%%%%%%%%%%%%%%%%%
open('Helper_functions/Prior_edge_inclusion_probability_study')

