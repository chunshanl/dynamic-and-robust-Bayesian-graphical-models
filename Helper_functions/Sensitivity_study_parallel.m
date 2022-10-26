%% Sensitivity study
% Sensitivity analysis of three proposed dynamic Bayesian graphical models  
% Use slightly contaminated data
% Run MCMC on a grid of parameter settings
% parallel computing on 13 cores

mkdir('Sensitivity_results')

addpath Call_functions
addpath MCMC_algorithms
addpath Helper_functions

%%%%%%%%%%%%%%%%%
%% Generate data
%%%%%%%%%%%%%%%%%
% mkdir('Synthetic_data')
% nspikes = 5;  % number of spikes of each signal
% Generate_data_spikeshmm_rep
% clear

%% Set up parameter grid
model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"];

p = 20;
pii_list = 2:5 / (p-1);
nu_list = [3, 6, 9];
K_sb_list = [4, 7, 10];
%
nu_vec_2 = nu_list';
K_sb_vec_2 = nan( length(nu_list),1);
model_method_vec_2 = repelem("CMTHMM", length(nu_vec_2),1);
[nu_vec_3, K_sb_vec_3] = ...
    ndgrid(nu_list,  K_sb_list);
nu_vec_3 = nu_vec_3(:);
K_sb_vec_3 = K_sb_vec_3(:);
model_method_vec_3 = repelem("TDTHMM", length(nu_vec_3),1);
%
model_method_vec = ["LPMHMM"; model_method_vec_2; model_method_vec_3];
nu_vec = [nan; nu_vec_2; nu_vec_3];
K_sb_vec = [nan; K_sb_vec_2; K_sb_vec_3];

ncore = length(nu_vec);

%% Set up parallel computing
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
% for  coreind =1:ncore
    
    p = 20;
    pii_list = [2:5] / (p-1);
    nu_list = [3, 6, 9];
    K_sb_list = [4, 7, 10];
    %
    nu_vec_2 = nu_list';
    K_sb_vec_2 = nan( length(nu_list),1);
    model_method_vec_2 = repelem("CMTHMM", length(nu_vec_2),1);
    [nu_vec_3, K_sb_vec_3] = ...
        ndgrid(nu_list,  K_sb_list);
    nu_vec_3 = nu_vec_3(:);
    K_sb_vec_3 = K_sb_vec_3(:);
    model_method_vec_3 = repelem("TDTHMM", length(nu_vec_3),1);
    %
    model_method_vec = ["LPMHMM"; model_method_vec_2; model_method_vec_3];
    nu_vec = [nan; nu_vec_2; nu_vec_3];
    K_sb_vec = [nan; K_sb_vec_2; K_sb_vec_3];
    
    model_method = model_method_vec(coreind);
    nu = nu_vec(coreind);
    K_sb = K_sb_vec(coreind);
    disp(['model:', model_method])
    disp(['nu:', num2str(nu)])
    disp(['K:', num2str(K_sb)])
    

    %% MCMC
    for rep = 1 : length(pii_list)

        fprintf('current_rep = %d\n', rep);
        pii = pii_list(rep);

        %% Load data
        data = load('Synthetic_data/Synthetic_data_5spikeshmm_rep_1');
        data = data.data;

        %% MCMC
        S = 3;
        h = 50^2;
        v0 = 0.02^2;  
        v1 = h * v0;
        a_alpha = 1;
        b_alpha = 1;
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
        parsave(results, states_save, rep, nu, K_sb,  model_method)

    end
end

%%
function parsave(results, states_save, rep, nu, K_sb, model_method)
    if strcmp(model_method, "LPMHMM")
        save(strcat('Sensitivity_results/sens_results_',model_method,...
            '_pii', num2str(rep)), 'results')
        save(strcat('Sensitivity_results/sens_results_',model_method,...
            '_pii', num2str(rep), '_states'), 'states_save')
    elseif strcmp(model_method, "CMTHMM")
        save(strcat('Sensitivity_results/sens_results_',model_method, '_nu', num2str(nu) ,...
            '_pii', num2str(rep)), 'results')
        save(strcat('Sensitivity_results/sens_results_',model_method,'_nu', num2str(nu) ,...
            '_pii', num2str(rep), '_states'), 'states_save')
    else
        save(strcat('Sensitivity_results/sens_results_',model_method, '_nu', num2str(nu) ,...
            '_K', num2str(K_sb),...
            '_pii', num2str(rep)), 'results')
        save(strcat('Sensitivity_results/sens_results_',model_method,'_nu', num2str(nu) ,...
            '_K', num2str(K_sb),...
            '_pii', num2str(rep), '_states'), 'states_save')
    end
end