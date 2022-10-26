function [geweke_num_edges, geweke_state, prob_of_last_cluster ] = ...
      Analyze_MCMC_results( data, results, plot_disp)
%% Output:
% geweke_num_edges: S x p array of p-values, geweke test for number of
% edges in each state
% geweek_state:  S x p array of p-values, geweke test for number of
% observations in each state
% prob_of_last_cluster: averaged probability of tau_ti in the last cluster in
% the dynamic Dirichlet-t graphical model

%% Input:
% data : a structure containing the follow elements 
%   data.TC : T x p array of time courses
%   Other fields are not required in this function
% results: a structure of MCMC results, output of the call functions
% plot_disp: true or false, whether to display convergence plots

%% Load
[T, p] = size(data.TC);

%% Geweke test for convergence:  number of edges in each state across iterations
num_edges = (squeeze(sum(sum(results.adj_save)))-p)/2;
csvwrite('Helper_functions/geweke_input.csv',num_edges)

if plot_disp
    figure()
    hold on
    for s = 1:results.S
        subplot(results.S,1,s)  
        plot(num_edges(s,:), 'color', 'black')
        if s < results.S
            set(gca,'Xticklabel',[])
        end
        set(gca,'FontName', 'Times New Roman','Xticklabel',[], 'box', 'off')
        title(['State ', num2str(s)], 'FontWeight','Normal', 'FontSize', 15)       
    end
    hold off
    set(gcf,'Position',[10,10,1200,800])
    % saveas(gcf,'Application_results/futures_post_nedges.png')
end

system('R CMD BATCH Helper_functions/geweke_test.R');
geweke_num_edges = csvread('Helper_functions/geweke_pvalues.csv');


%% Geweke test for convergence:  number of observations in each state across iterations
S = results.S;
nmc = results.nmc;
num_obs = zeros(S, nmc);
for s = 1:S
    num_obs(s,:) = sum(results.states_save == s);
end
csvwrite('Helper_functions/geweke_input.csv',num_obs)

if plot_disp
    figure()
    plot(num_obs')
    title('Number of observations in each state across iterations')
end

system('R CMD BATCH Helper_functions/geweke_test.R');
geweke_state = csvread('Helper_functions/geweke_pvalues.csv');


%% Probability that the last cluster in the truncated Dirichlet-t is empty
prob_of_last_cluster = nan;
if strcmp(results.method, 'TDTHMM')
    probs = zeros(results.K_sb, T);
    for t = 1:T 
        for ii = 1:results.nmc
            temp = size(unique(results.tau_t_save(:, t, ii)),1);
            probs(1:temp, t) = probs(1:temp, t) + 1;
        end
    end

    probs = probs / results.nmc;
    prob_of_last_cluster = mean(probs(end, :));
end

