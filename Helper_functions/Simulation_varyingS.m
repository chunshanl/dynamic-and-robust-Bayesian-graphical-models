%% Simulation study, select the number of hidden states S
% 25 randomly generated highly contaminated data
% The dynamic Dirichlet-t graphical model
% (the joint graphical lasso is in R)
% 25 replicates for each data and model combination
% Parallel computing on 9 cores

% set work directory to home
mkdir('Simulation_results_varyingS') % MCMC results are saved in the Simulation_results_varyingS folder
addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%%%%%%%%%%%%%%%%
%% Generate data
%%%%%%%%%%%%%%%%%
mkdir('Synthetic_data')
Generate_and_plot_data
clear

%% Set up parallel computing
ncore = 12;
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

    data_method = "10spikeshmm";
    model_method = "TDTHMM";
    
    if coreind < 13
       temp = [coreind, coreind + 13];
    else
        temp = coreind
    end

    %% Fit model
    for rep = temp
        fprintf('current_rep = %d\n', rep);

        %% Load data
        data = [];
        data = load(strcat('Synthetic_data/Synthetic_data_', data_method, '_rep_', num2str(rep)));
        data = data.data;

        %% MCMC
        S = 2;  % number of hidden states
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
        parsave(results, states_save, rep, data_method, model_method, S)

    end
end

%%
function parsave(results, states_save, rep, data_method, model_method, S)
    % save(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)), 'results')
    % save(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)), 'states_save')
    save(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method, '_', num2str(S),'states_rep_', num2str(rep)), 'results')
    save(strcat('Simulation_results_varyingS/sim_results_',data_method,'_',model_method, '_', num2str(S), 'states_states_rep_', num2str(rep)), 'states_save')
end

