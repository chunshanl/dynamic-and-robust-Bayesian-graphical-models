%% Compare model performances
% Three kinds of data: classical-t, slightly contaminated, highly
% contaminated
% Three models: the dynamic Gaussian graphical model, the dynamic classical-t
% graphical model, the dynamic Dirichlet-t graphical model
% (the joint graphical lasso is in R)
% 25 replicates for each data and model combination
% Parallel computing on 9 cores

% set work directory to home
mkdir('Synthetic_data')  % save synthetic data from the simulation study
mkdir('Simulation_results') % MCMC results are saved in the Simulation_results folder
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%%%%%%%%%%%%%%%%
%% Generate data
%%%%%%%%%%%%%%%%%
Generate_and_plot_data
clear

%% Set up parallel computing
ncore = 9;
parpool(ncore)
% load functions 
poolobj = gcp;
addAttachedFiles(poolobj,{'Call_Function_TDTHMM_nomean.m',...
    'MCMC_Algorithm_TDTHMM_nomean.m',...
    'Call_Function_CMTHMM_nomean.m',...
    'MCMC_Algorithm_CMTHMM_nomean.m',...
    'Call_Function_LPMHMM_nomean.m',...
    'MCMC_Algorithm_LPMHMM_nomean.m',...
})

%%
parfor coreind =1:ncore

    data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;
    model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"];
    [data_vec, model_vec] = meshgrid(data_method_list, model_method_list);
    
    data_method = data_vec(coreind);
    model_method = model_vec(coreind);
    disp(['data:', data_method])
    disp(['model:', model_method])
    

    %% Fit model
    for rep = 1:25

        fprintf('current_rep = %d\n', rep);

        %% Load data
        data = [];
        data = load(strcat('Synthetic_data/Synthetic_data_', data_method, '_rep_', num2str(rep)));
        data = data.data;

        %% MCMC
        S = 3;
        h = 50^2;
        v0 = 0.02^2;  
        v1 = h * v0;
        pii = 3/(size(data.TC,2)-1);     
        nu=3;
        a_alpha = 1;
        b_alpha = 1;
        K_sb = 7;
        if model_method == "TDTHMM"
            burnin=500;
            nmc=1000;
        else
            burnin = 2000;
            nmc = 8000;
        end
              
        results=[]; states_save = [];
        disp_result = false;
        rng(12345 + rep, 'twister')
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
        
        %% Save results
        results.states_save = [];
        parsave(results, states_save, rep, data_method, model_method)

    end
end

%%
function parsave(results, states_save, rep, data_method, model_method)
    save(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)), 'results')
    save(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)), 'states_save')
end