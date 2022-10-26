%% Generate 3 kinds of data, each with 25 replications and plot data

%% Generate data 
data=[];
% classical-t data
Generate_data_cmthmm_rep
% slightly contaminated data
nspikes = 5;  % number of spikes of each signal
Generate_data_spikeshmm_rep
% highly contaminated data
nspikes = 10;  % number of spikes of each signal
Generate_data_spikeshmm_rep
clear
% data are saved in the Synthetic_data folder

%% Plot data
data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;  % three kinds of synthetic data in the simulation study
for ind = 1:length(data_method_list)
    data_method = data_method_list (ind);
    load(strcat('Synthetic_data/Synthetic_data_', data_method, '_rep_1'))
    Plot_synthetic_data(data)  % function to plot synthetic data
end
% plots of the data are saved in the Simulation_analysis folder
