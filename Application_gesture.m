addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%
%mkdir('Application_results')

%%%%%%%%%%%%%%%%%%
%% Get data
%%%%%%%%%%%%%%%%%%
%% Load and plot gesture data - Figure 5
load('Application_data/a1_va3.mat')
data.TC = TC;
TC = [];

figure
set(gcf,'Position',[10,10,1200,250])
plot(data.TC)
xlim([0, length(data.TC)])
%saveas(gcf,'Application_results/gesture_data.png')
exportgraphics(gcf,'Application_results/gesture_data.jpeg','Resolution',300)


%% Load and plot the phases provided by the original data set
load('Application_data/a1_va3_phase.mat')

figure();
set(gcf,'Position',[10,10,1200,250])
plot(phase, '-o','MarkerSize',1, 'color', 'black');
ylim([0.5, 5+0.5])
xlim([0, length(data.TC)])
title('Phases', 'FontWeight','Normal', 'FontSize', 15)
%saveas(gcf,'Application_results/gesture_true_phase.png')
exportgraphics(gcf,'Application_results/gesture_true_phase.jpeg','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%
%% Run MCMC for the Bayesian models
%%%%%%%%%%%%%%%%%%%%%%
S = 2;
burnin = 8000; 
nmc = 8000;
h = 50^2; 
v0 = 0.02^2; 
v1 = h * v0;
pii = 3/(size(data.TC,2)-1);  
K_sb = 6; 
nu=3;
a_alpha = 1;
b_alpha = 1;

disp_result = true;
results =[]; states_save = [];
rng(12453, 'Threefry')
% Dynamic Dirichlet-t graphical model
[results,states_save] = Call_Function_TDTHMM_nomean_application(data, burnin, nmc, S , v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result);
% Linked Dirichlet-t graphical model without HMM, given phases from the data
%[results,states_save] = Call_Function_TDTHMM_nomean_application_givenstate(data, burnin, nmc, S , v0, v1, pii, nu, a_alpha, b_alpha, K_sb, disp_result, ...
%    phase);
% Linked Gaussian graphical model without HMM, given phases from the
% dynamic Dirichlet graphical model
%[results, states_save] = Call_Function_LPMHMM_nomean_application_givenstate(data, burnin, nmc, S, v0, v1, pii, disp_result, ...
%    states_est);
% Dynamic Gaussian graphical model
% [results, states_save] = Call_Function_LPMHMM_nomean_application(data, burnin, nmc, S, v0, v1, pii, disp_result);
%% Save MCMC results
tau_t_save = results.tau_t_save; results.tau_t_save = [];

save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states'],'results', '-v7.3')
save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_tau'],'tau_t_save', '-v7.3')
save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_states'],'states_save', '-v7.3')
%save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate'],'results', '-v7.3')
%save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate_tau'],'tau_t_save', '-v7.3')
%save(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate_states'],'states_save', '-v7.3')

%% Load MCMC results
S = 2;
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states'])
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_tau'])
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_states'])
%load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate'])
%load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate_tau'])
%load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_givenstate_states'])

results.states_save = states_save;
results.tau_t_save = tau_t_save;
states_save = []; tau_t_save = [];

%% Plot MCMC results - Figure 5-8
Plot_application_results_gesture(data, results)

%% MCMC diaganosis
plot_disp = false;
[geweke_num_edges, geweke_state, waic, prob_of_last_cluster ] = ...
       Analyze_MCMC_results_waic(data, results, plot_disp);
display(geweke_num_edges)  % geweke test of number of edges in each state
display(geweke_state)  % geweke test of number of observations in each state
display(waic)  % WAIC
display(prob_of_last_cluster)  % probability of tau_ti in the Kth cluster

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MCMC of 16000 iterations under different numder of hidden states
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel computing
open('Helper_functions/Application_gesture_varyingS.m')

%% Plot the posterior probability of each state across time under each S
% S_all = 2:4;
% for S = S_all
%     results = [];
%     load(['Application_results/gesture_result_', num2str(S),'states'])
%     ppi_HMM = results.ppi_HMM;
%     save(['Application_results/gesture_ppiHMM_', num2str(S),'states'], 'ppi_HMM')
% end
data.X_t = 1:size(data.TC,1);
S_all = 2:5;
for S = S_all
    load(['Application_results/gesture_ppiHMM_', num2str(S),'states'])
    figure()
    for i = 1:S
        subplot(S,1,i)  
        plot(data.X_t,ppi_HMM(i,:),'-o','MarkerSize', 1, 'color', 'black');
        ylim([0,1])
        xlim([min(data.X_t), max(data.X_t)])
        title(['State ', num2str(i)])
    end
    set(gcf,'Position',[10,10,1200,800])
    saveas(gcf,['Application_results/gesture_sensitivity_', num2str(S),'states.png'])
end

%% Choose the number of hidden states based on WAIC
%%% Compute WAIC for different number of hidden states
% plot_disp = true;
% S_all = 2:5;
% waic_all = zeros(1,length(S_all));
% for ii = 1:length(S_all)
%     S = S_all(ii);
%     load(['Application_results/futures_result_', num2str(S),'states'])
%     load(['Application_results/futures_result_', num2str(S),'states_tau'])
%     load(['Application_results/futures_result_', num2str(S),'states_states'])
%     results.states_save = states_save;
%     results.tau_t_save = tau_t_save;
%     states_save = [];
%     tau_t_save = [];
% 
%     [geweke_num_edges, geweke_state, waic, prob_of_last_cluster ] = ...
%        Analyze_MCMC_results_new( data, results, plot_disp);    
%     waic_all(ii) = waic;
%     display(waic)  % WAIC
% end
% save('Application_results/waic_varyingS.mat', 'waic_all')

%%%  WAIC for S = 2, 3, 4, 5
load('Application_results/waic_varyingS.mat')
disp(waic_all)


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MCMC of 16000 iterations for different pi
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel computing
open('Helper_functions/Application_gesture_varyingPi.m')
%% Plot the posterior state probability under each pi
%%% The probabilities are provided in the Application_results folder
pii_list = [2,4,5];
S=2;
% data.X_t = 1:size(data.TC,1);
% for pii = pii_list
%     results = [];
%     load(['Application_results/gesture_result_pii', num2str(pii)])
%     ppi_HMM = results.ppi_HMM;
%     save(['Application_results/gesture_ppiHMM_pii', num2str(pii)], 'ppi_HMM')
% end
for pii = pii_list
    load(['Application_results/gesture_ppiHMM_pii', num2str(pii)])
    figure()
    for i = 1:S
        subplot(S,1,i)  
        plot(data.X_t, ppi_HMM(i,:),'-o','MarkerSize', 1, 'color', 'black');
        ylim([0,1])
        xlim([min(data.X_t), max(data.X_t)])
        title(['State ', num2str(i)])
    end
    set(gcf,'Position',[10,10,1200,800])
    saveas(gcf,['Application_results/gesture_sensitivity_pii', num2str(pii),'.png'])
end
%% Compute posterior edge inclusion probabilities under each pi
% pii_list = [2, 4, 5];
% probs_all = [];
% thresh = 0.5;
% p = size(data.TC, 2);
% for pii = pii_list
%     results = [];
%     load(['Application_results/gesture_result_pii', num2str(pii)])
%     ppi = mean(results.adj_save, 4);
%     edge_include_all = ppi > thresh;
%     % edge inclusion probabilities
%     temp = 1 - ( squeeze(sum(sum(~edge_include_all))) /2 ) / (p*(p-1)/2);
%     probs_all = [probs_all; temp'];
% end
% save('Application_results/gesture_piigrph_allpii', 'probs_all')

%%%
load('Application_results/gesture_piigrph_allpii')
% rows: pi = 2/(p-1), 4/(p-1), 5/(p-1)
% columns: 2 states 
disp(probs_all)

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Fused graphical lasso model
%%%%%%%%%%%%%%%%%%%%%%%%%
% Load states from the Dirichlet-t model
S = 2;
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states'])
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_tau'])
load(['Application_results/gesture_result_TDTHMM_', num2str(S),'states_states'])
results.states_save = states_save;
results.tau_t_save = tau_t_save;
states_save = []; tau_t_save = [];
T = size(data.TC,1);
mm=results.ppi_HMM';
states_est=zeros(1, T);
for i=1:T
    [~,c]= max(mm(i,:));
    states_est(i)=c;
end
states_dirichlet = 3 - states_est;  % re-order
csvwrite('states_dirichlet.csv', states_dirichlet')
%% Run fused graphical lasso
open('Helper_functions/Application_JGL.R')
%% Load results and plot
p = 32;
C_est = csvread('Application_results/gesture_result_FGL_C.csv');
adj_est = csvread('Application_results/gesture_result_FGL_adj.csv');
CC_all = zeros(p,p,S);
edge_include_all = zeros(p, p, S);
for s = 1:S
    dims = ((s-1)*p + 1):(s*p);
    CC_all(:,:,s) = C_est(dims,:);
    edge_include_all(:,:,s) = adj_est(dims,:);
end
%% Plot parcial correlations
Pcc_all = zeros(p, p, S);
for s =1:S
    CC = CC_all(:,:,s);
    CC(~edge_include_all(:, :, s))=0;  
    temp=CC;
    a=sqrt(diag(temp));  
    temp=temp./repmat(a,1,p);
    temp=-temp./repmat(a,1,p)';    
    temp= (temp+temp')/2;
    Pcc=temp+eye(p);
    Pcc_all(:,:,s) = Pcc;
end

% plot graphs
for s =1:S
   Pcc = Pcc_all(:,:,s);
   % plot graph
   temp=Pcc;
   temp=temp-diag(diag(temp));
   G=graph(temp);
   figure
   h=plot(G,'Layout','circle');
   h.NodeFontWeight = 'bold';
   h.NodeFontSize = 18;
     
   % change line width according to the absolute value of partial correlations
   m=max(abs(Pcc(:)));  % rescale edge weights so the maximum edge width in plot is 7 
   h.LineWidth=abs(G.Edges.Weight)/m*7;
  
   % change edge style if the edge is unique to this state
   unique_edges = zeros(p, p);
   for i = 1:(p-1)
       for j = 1:p
           if edge_include_all(i, j, s) == 1 && sum(edge_include_all(i,j,:))==1
               unique_edges(i, j) = 1;
           end
       end
   end
   temp = edge_include_all(:,:,s);
   temp(upperLogic(p))= 0;
   temp = temp - eye(p);
   unique_edges = unique_edges(logical(temp));
   a=cell(length(G.Edges.Weight),1);
   a(:)={'-'};
   a(logical(unique_edges))={'--'};
   h.LineStyle=a;
   
   % show edge width on graph
%    h.EdgeLabel=round(G.Edges.Weight,2);
%    h.EdgeFontWeight = 'bold';
%    h.EdgeFontSize = 8;
   
   % change edge color if weight is negative
   a= repmat([0,0.477,0.741],[length(G.Edges.Weight),1]);
   a(G.Edges.Weight<0,:)= repmat( [0.77,0,0.1], [sum(G.Edges.Weight<0),1]);
   h.EdgeColor=a;
   
   title(['State ', num2str(s)], 'FontSize', 19)
   set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
   set(gcf,'Position',[10,10,800,600])
   %saveas(gcf, ['Application_results/gesture_JGL_graph_state_', num2str(s), '.png'])
   exportgraphics(gcf, ['Application_results/gesture_JGL_graph_state_', num2str(s), '.jpeg'], Resolution=300)

end

