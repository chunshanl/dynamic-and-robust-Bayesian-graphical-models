function [geweke_num_edges, geweke_state, waic, prob_of_last_cluster ] = ...
      Analyze_MCMC_results_waic( data, results, plot_disp)
%% Output:
% geweke_num_edges: S x p array of p-values, geweke test for number of
% edges in each state
% geweek_state:  S x p array of p-values, geweke test for number of
% observations in each state
% waic
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

%% Geweke test for convergence: number of edges in each state across iterations
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
    saveas(gcf,'Application_results/post_nedges.png')
end

system('R CMD BATCH Helper_functions/geweke_test.R');
geweke_num_edges = csvread('Helper_functions/geweke_pvalues.csv');
% disp(geweke_num_edges)

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
    for s = 1:S
        subplot(S, 1, s)
        plot(num_obs(s,:))
    end
    title('Number of observations in each state across iterations')
    saveas(gcf,'Application_results/post_nobs.png')
end

system('R CMD BATCH Helper_functions/geweke_test.R');
geweke_state = csvread('Helper_functions/geweke_pvalues.csv');

% disp(geweke_state)

%% Compute WAIC
TC_scaled = data.TC;
if isfield(results, 'mean_input')
    TC_scaled = TC_scaled - results.mean_input;
end
if isfield(results, 'std_input')
    TC_scaled = TC_scaled./results.std_input;
end
log_likelihoods = zeros(results.nmc, T);
if strcmp(results.method, 'TDTHMM')
    for t = 1:T
        for ii = 1:results.nmc
            C_temp = results.C_save(:,:,results.states_save(t, ii), ii);
            %adj_temp = ~results.adj_save(:,:,results.states_save(t, ii),ii);
            %C_temp(adj_temp) = 0;
            C_temp = diag(sqrt(results.tau_t_save(:,t,ii))) * C_temp * diag(sqrt(results.tau_t_save(:,t,ii)));
            temp = -p/2 * log(2*pi) + 1/2 * log(det(C_temp)) - ...
                    1/2 * TC_scaled(t,:) * C_temp * TC_scaled(t,:)';
            log_likelihoods(ii, t) = temp;
%             if ~isreal(temp)
%                 log_likelihoods(ii, t) = 0;
%             end
        end
    end
elseif  strcmp(results.method, 'CMTHMM')
    for t = 1:T
        for ii = 1:results.nmc
            C_temp = results.C_save(:,:,results.states_save(t, ii), ii);
            % adj_temp = ~results.adj_save(:,:,results.states_save(t, ii),ii);
            % C_temp(adj_temp) = 0;
            C_temp = results.tau_t_save(t,ii)* C_temp;
            temp = -p/2 * log(2*pi) + 1/2 * log(det(C_temp)) - ...
                    1/2 * TC_scaled(t,:) * C_temp * TC_scaled(t,:)';
            log_likelihoods(ii, t) = temp;
%             if ~isreal(temp)
%                 log_likelihoods(ii, t) = 0;
%             end
        end
    end
end

lppd = sum(log(mean(exp(log_likelihoods))));
pwaic2 = sum(var(log_likelihoods));
waic = - 2 * lppd + 2 * pwaic2;
% disp(waic)

%% Probability that the last cluster in the truncated Dirichlet-t is not empty
prob_of_last_cluster = nan;
if strcmp(results.method, 'TDTHMM')
    probs = zeros(results.K_sb, T);
    for t = 1:T 
        for ii = (results.nmc-1000):results.nmc
            temp = size(unique(results.tau_t_save(:, t, ii)),1);
            probs(1:temp, t) = probs(1:temp, t) + 1;
        end
    end
    probs = probs / results.nmc;
    prob_of_last_cluster = mean(probs(end, :));
end

%% Trace plots
if plot_disp
    
    % trace plots of elements in the precision matrix
    for s_temp=1:results.S
        figure()
        for i = 1:5 % 1:(p-1)
            for j= i + 1 %(i+1):p
                plot(squeeze(results.C_save(i,j,s_temp,:)))
                hold on
            end
        end
        title(['Precision matrix of state' string(s_temp) 'across iter'])
        hold off
    end
       
end

%% Geweke test for convergence: elements in Phi
temp = [];
for i = 1:(results.S-1)
    for j = (i+1):results.S
        temp = [temp, squeeze(results.Phi_save(i,j,:))];
    end
end
temp = temp';
csvwrite('Helper_functions/geweke_input.csv', temp)

% Trace plots of elements in Phi
ntemp = size(temp,1);
figure()
for i = 1:ntemp
    subplot(ceil(ntemp/2),2,i)  
    plot(1:results.nmc, temp(i,:), 'color', 'black')
    ylim([-1,1])
    if i < ntemp - 1
        set(gca,'Xticklabel',[])
    end
end
set(gcf,'Position',[10,10,1200,800])
saveas(gcf,'Application_results/trace_phi.png')
    
system('R CMD BATCH Helper_functions/geweke_test.R');
geweke_Phi = csvread('Helper_functions/geweke_pvalues.csv');


