function [performance_state_match_state_avg, performance_state_match_by_state,...
    performance_graph_state_avg, performance_graph_by_state, performance_graph_roc_state_avg] = ...
      Compute_simulation_performance( data, results,  plot_disp, rate_disp)
%% Compute performance of graph estimation and state estimation

%% Outputs:
% performance_state_match_state_avg: performance of state match averaged
% across states
% performance_state_match_by_state: performance of state match in each state
% performance_graph_state_avg: performance of graph estimation averaged
% across states
% performance_graph_by_state: performance of graph estimation in each state
% performance_graph_roc_state_avg: ROC of graph estimation averaged across
% states

% performance of states match: true positive rate (tpr), false positive rate (fpr), mcc; 
% performance of graph estimation: tpr, fpr, mcc, auc, frobinous loss (fl),
% frobinous loss considering tau (fl_tau)

%% Inputs:
% data: a structure of simulated data and true parameters generated by
% Generate_and_plot_data.m
% results: a structure of MCMC results, output of the call functions
% plot_disp: true or false, whether to plot the ROC curves
% rate_disp: true or false, whether to display the performance rates

%% Load 
TC=data.TC;
S_true=data.S_true;
states_true=data.states_true;
S=results.S;
[T,p]=size(TC);

%% Get the most possible posterior states 
mm=results.ppi_HMM';
states_est=zeros(1, T);
for i=1:T
    [~,c]= max(mm(i,:));
    states_est(i)=c;
end

%% Performace of state match -- by state 
% Match posterior states to true states
[post_order, states_est_sort] = match_state(states_true, states_est);
% compute performance rates by state
performance_state_match_by_state = zeros(S_true, 3);
for s = 1:S_true
    temp_est = states_est_sort == s;
    temp_true = states_true ==s;
    tp = sum(temp_est& temp_true);
    tn = sum((~temp_est)& (~temp_true));
    fp = sum((temp_est)& (~temp_true));
    fn = sum((~temp_est)& (temp_true));
    tpr = tp / (tp + fn);
    fpr = fp / (fp + tn);
    mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
    performance_state_match_by_state(s,:) = [tpr, fpr, mcc];
end
if rate_disp
    display(performance_state_match_by_state)
end

%% Performace of state match -- averaged aross states
% Compute performance of state estimation averaged across states
[tpr, fpr, mcc] = get_perf_state(states_true, states_est_sort);
performance_state_match_state_avg = [tpr, fpr, mcc];
if rate_disp
    display(performance_state_match_state_avg)
end

%% Performance of graph estimation -- by state 
%%% Compute true precision 
C_true = data.C_true;
adj_true = abs(C_true)>0;
%%% Compute posterior precision
ppi = mean(results.adj_save, 4);
C_est = mean(results.C_save, 4);
if ~isfield(results, 'std_input')
    results.std_input = ones(p, 1);
end
for s = 1:S
    C_temp = diag(1./ results.std_input) * C_est(:,:,s) * diag(1 ./ results.std_input);
    C_temp( ppi(:,:,s) < 0.5 ) = 0;
    C_est(:,:,s) = C_temp;
end
%%% Consider tau_t
if isfield(data,'tau_t_true')
    tau_t_true = data.tau_t_true;
end
if isfield(results,'tau_t_save')
    tau_t_est = mean(results.tau_t_save, length(size(results.tau_t_save)));
    tau_t_est = tau_t_est';
end
%%%
C_tau_true = zeros(p,p,S_true);
C_est_sort = zeros(p,p,S_true);
C_tau_est_sort = zeros(p,p,S_true);
ppi_sort = zeros(p,p,S_true);

performance_graph_by_state = [];

for s_true_temp = 1:S_true
    %%% Match states
    s_post_temp = post_order(s_true_temp);
    C_est_sort(:,:,s_true_temp) = C_est(:,:,s_post_temp);
    ppi_sort(:,:,s_true_temp) = ppi(:,:,s_post_temp);
    
    %%% Performance without considering tau
    [tpr_cur, fpr_cur, mcc_cur, auc_cur, ~, fl_cur] = ...
     get_perf_graph(...
                adj_true(:,:,s_true_temp), ...
                ppi(:,:, s_post_temp), ...
                C_true(:,:, s_true_temp), ...
                C_est(:,:, s_post_temp));

    %%% Performance Considering tau
    if strcmp(data.method, 'dmthmm')
        temp = mean(sqrt(tau_t_true(states_true == s_true_temp, :))) ;
        C_tau_true_temp = diag(temp)* C_true(:, :, s_true_temp) * diag(temp);
    elseif strcmp(data.method, 'cmthmm' )
        C_tau_true_temp = mean(tau_t_true(states_true == s_true_temp)) * C_true(:,:, s_true_temp);
    else
        C_tau_true_temp = C_true(:, :, s_true_temp);
    end
    C_tau_true(:,:,s_true_temp) = C_tau_true_temp;
    
    if strcmp(results.method, 'TDTHMM')
        temp = mean(sqrt(tau_t_est(states_est == s_post_temp,:))) ;
        C_tau_est_temp = diag(temp)* C_est(:, :, s_post_temp) * diag(temp);
    elseif strcmp(results.method, 'CMTHMM')
        C_tau_est_temp = mean(tau_t_est(states_est == s_post_temp)) *...
            C_est(:,:, s_post_temp);
    elseif strcmp(results.method, 'LPMHMM')
        C_tau_est_temp = C_est(:,:, s_post_temp);
    end
    C_tau_est_sort(:,:,s_true_temp) = C_tau_est_temp;
    
    [~, ~, ~, ~, ~, fl_tau_cur] = ...
     get_perf_graph(...
                adj_true(:,:,s_true_temp), ...
                ppi(:,:, s_post_temp), ...
                C_tau_true_temp, ...
                C_tau_est_temp);
   
    % Save
    performance_graph_by_state = [performance_graph_by_state;...
       [tpr_cur, fpr_cur, mcc_cur, auc_cur, fl_cur, fl_tau_cur]];
end

if rate_disp 
    display(performance_graph_by_state)
end

%% Performance of graph estimation - averaged across states
%%% Performance without considering tau
[tpr_cur, fpr_cur, mcc_cur, auc_cur, roc_curve, fl_cur] = ...
 get_perf_graph(...
            adj_true, ...
            ppi_sort, ...
            C_true, ...
            C_est_sort);
%%% Performance Considering tau
[~, ~, ~, ~, ~, fl_tau_cur] = ...
 get_perf_graph(...
            adj_true, ...
            ppi_sort, ...
            C_tau_true, ...
            C_tau_est_sort);
   
% Save
performance_graph_state_avg = [tpr_cur, fpr_cur, mcc_cur, auc_cur, fl_cur, fl_tau_cur];
performance_graph_roc_state_avg = roc_curve;

if rate_disp 
    display(performance_graph_state_avg)
end
if plot_disp
    figure
    plot(performance_graph_roc_state_avg(:,1), performance_graph_roc_state_avg(:,2))
    title('ROC')
    saveas(gcf, 'Simulation_analysis/sim_post_roc.png')
end

