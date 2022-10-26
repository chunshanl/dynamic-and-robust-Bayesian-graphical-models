%% Generate 3 kinds of data, each with 25 replications and plot data

%% Generate data 
data=[];
p = 10; % number of variables
S_true = 3; % number of hidden states
T = S_true * 250;  % length of the time series
nblock = 25;  % divide the time series equally into nblock number of blocks
n_perturb_delete = [5,5,5];  % number of edges to delete from s-1 to s
n_perturb_add = [5,5,5];  % number of edges to add from s-1 to s
% n_perturb_delete = [5,5];  % number of edges to delete from s-1 to s
% n_perturb_add = [5,5];  % number of edges to add from s-1 to s

% classical-t data
Generate_data_cmthmm_rep_appendix
% slightly contaminated data
nspikes = 5;  % number of spikes of each signal
Generate_data_spikeshmm_rep_appendix
% highly contaminated data
nspikes = 10;  % number of spikes of each signal
Generate_data_spikeshmm_rep_appendix
clear
% data are saved in the Synthetic_data_appendix folder

%% Plot data
% data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;  % three kinds of synthetic data in the simulation study
% for ind = 1:length(data_method_list)
%     data_method = data_method_list (ind);
%     load(strcat('Synthetic_data_appendix/Synthetic_data_', data_method, '_rep_1'))
%     Plot_synthetic_data(data)  % function to plot synthetic data
%     % plots of the data are saved in the Simulation_analysis folder
% end

